use thiserror::Error;

use std::fs::File;
use std::io::BufReader;
use std::path::Path;

use std::io::prelude::*;
use std::vec::Vec;

#[derive(Debug, Default)]
pub struct FastQRead {
    pub letters: Vec<u8>,
    pub qualities: Vec<u8>,
    pub title: String,
    pub sub_title: String,
}

impl FastQRead {
    pub fn new(title: &str) -> Self {
        let mut read = Self::default();

        read.title = String::from(title);
        read
    }

    pub fn reverse_complement_nucleotides(&mut self) {
        let letters: Vec<u8> = self
            .letters
            .iter()
            .rev()
            .map(|n| match n {
                0 => 0,
                1 => 2,
                2 => 1,
                3 => 4,
                4 => 3,
                _ => panic!("Invalid nuceotide {} found!", n),
            })
            .collect();

        self.letters = letters;
    }
}

#[derive(Error, Debug)]
pub enum FastQFileError {
    #[error("IO error while reading fastq file")]
    IO {
        #[from]
        source: std::io::Error,
    },
    #[error("Did not find title (titles should start with '@')")]
    NoTitleLine,
    #[error("Did not find expected line starting with '+'")]
    NoDescriptionLine,
    #[error("The quality sequence has unexpected characters")]
    InvalidQualityLetter,
    #[error("The nucleotide sequence and the quality sequence are different lengths")]
    MismatchedSequenceLength,
    #[error("Found nucleotide {c} that is not |ATCGNatcgn|")]
    InvalidNucleotideLetter { c: char },
    #[error("EOF caused Incomplete record")]
    IncompleteRecord,
    #[error("Found FASTA style title (title started with a '>'). Expected FASTQ files.")]
    FastATitleLine,
    #[error("Missing read from pair, file truncated")]
    MissingPairedRead,
}

fn nuc_string_to_vec(letters: &str) -> Result<Vec<u8>, FastQFileError> {
    let mut ret: Vec<u8> = Vec::with_capacity(letters.len());
    for n in letters.chars() {
        match n {
            'n' | 'N' => ret.push(0),
            'a' | 'A' => ret.push(1),
            't' | 'T' => ret.push(2),
            'c' | 'C' => ret.push(3),
            'g' | 'G' => ret.push(4),
            _ => return Err(FastQFileError::InvalidNucleotideLetter { c: n }),
        }
    }

    Ok(ret)
}

/// Fastq file things
pub trait FastQFileReader {
    fn read_next(&mut self, buf: &mut FastQRead) -> Result<bool, FastQFileError>;
}

#[derive(Debug)]
pub struct FastQFile<R: Read> {
    file_obj: BufReader<R>,
}

impl<R: Read> FastQFile<R> {
    pub fn new(file_obj: BufReader<R>) -> Self {
        FastQFile { file_obj: file_obj }
    }
}

impl<R: Read> FastQFileReader for FastQFile<R> {
    fn read_next(&mut self, buf: &mut FastQRead) -> Result<bool, FastQFileError> {
        let mut title = String::new();
        let mut nucleotides = String::new();
        let mut sub_title = String::new();
        let mut quality_letters = String::new();

        if self.file_obj.read_line(&mut title)? == 0 {
            return Ok(false);
        }

        if !title.starts_with("@") {
            if title.starts_with(">") {
                return Err(FastQFileError::FastATitleLine);
            }
            return Err(FastQFileError::NoTitleLine);
        }
        title = title[1..].trim_end().to_string();

        if self.file_obj.read_line(&mut nucleotides)? == 0 {
            return Err(FastQFileError::IncompleteRecord);
        }

        nucleotides = nucleotides.trim_end().to_string();
        let letters: Vec<u8> = nuc_string_to_vec(&nucleotides)?;
        if self.file_obj.read_line(&mut sub_title)? == 0 {
            return Err(FastQFileError::IncompleteRecord);
        }

        if !sub_title.starts_with("+") {
            return Err(FastQFileError::NoDescriptionLine);
        }
        sub_title = sub_title[1..].trim_end().to_string();

        if self.file_obj.read_line(&mut quality_letters)? == 0 {
            return Err(FastQFileError::IncompleteRecord);
        }

        quality_letters = quality_letters.trim_end().to_string();
        if quality_letters.bytes().any(|c| !c.is_ascii_graphic()) {
            return Err(FastQFileError::InvalidQualityLetter);
        }

        if nucleotides.len() != quality_letters.len() {
            return Err(FastQFileError::MismatchedSequenceLength);
        }
        let qualities: Vec<u8> = quality_letters.bytes().map(|v| v - 32).collect();

        *buf = FastQRead {
            letters: letters,
            qualities: qualities,
            title: title,
            sub_title: sub_title,
        };
        return Ok(true);
    }
}

impl FastQFile<File> {
    pub fn open<P: AsRef<Path>>(path: &P) -> Result<Self, std::io::Error> {
        let file = File::open(path)?;
        Ok(FastQFile::new(BufReader::new(file)))
    }
}

impl FastQFile<std::io::Stdin> {
    pub fn from_stdin() -> Self {
        FastQFile::new(BufReader::new(std::io::stdin()))
    }
}

pub trait PairedFastQReader {
    fn read_next(
        &mut self,
        buf_r1: &mut FastQRead,
        buf_r2: &mut FastQRead,
    ) -> Result<bool, FastQFileError>;
}

pub struct FastQPairedFiles {
    r1_reader: Box<dyn FastQFileReader>,
    r2_reader: Box<dyn FastQFileReader>,
    reverse_complement_r2_nucleotides: bool,
}

impl FastQPairedFiles {
    pub fn new(
        stream_r1: Box<dyn FastQFileReader>,
        stream_r2: Box<dyn FastQFileReader>,
        reverse_complement_r2_nucleotides: bool,
    ) -> Self {
        FastQPairedFiles {
            r1_reader: stream_r1,
            r2_reader: stream_r2,
            reverse_complement_r2_nucleotides: reverse_complement_r2_nucleotides,
        }
    }
}

impl PairedFastQReader for FastQPairedFiles {
    fn read_next(
        &mut self,
        buf_r1: &mut FastQRead,
        buf_r2: &mut FastQRead,
    ) -> Result<bool, FastQFileError> {
        if !self.r1_reader.read_next(buf_r1)? {
            return Ok(false);
        }

        if !self.r2_reader.read_next(buf_r2)? {
            return Err(FastQFileError::MissingPairedRead);
        }

        if self.reverse_complement_r2_nucleotides {
            buf_r2.reverse_complement_nucleotides();
        }

        Ok(true)
    }
}

pub struct FastQInterleavedFile {
    reader: Box<dyn FastQFileReader>,
    reverse_complement_r2_nucleotides: bool,
}

impl FastQInterleavedFile {
    pub fn new(stream: Box<dyn FastQFileReader>, reverse_complement_r2_nucleotides: bool) -> Self {
        FastQInterleavedFile {
            reader: stream,
            reverse_complement_r2_nucleotides: reverse_complement_r2_nucleotides,
        }
    }
}

impl PairedFastQReader for FastQInterleavedFile {
    fn read_next(
        &mut self,
        buf_r1: &mut FastQRead,
        buf_r2: &mut FastQRead,
    ) -> Result<bool, FastQFileError> {
        if !self.reader.read_next(buf_r1)? {
            return Ok(false);
        }

        if !self.reader.read_next(buf_r2)? {
            return Err(FastQFileError::MissingPairedRead);
        }

        if self.reverse_complement_r2_nucleotides {
            buf_r2.reverse_complement_nucleotides();
        }

        Ok(true)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_reverse_complement_nucleotides() {
        let mut read = FastQRead::default();
        read.letters = [
            2, 2, 1, 1, 2, 2, 4, 4, 2, 1, 1, 1, 2, 1, 1, 1, 2, 3, 2, 3, 3, 2, 1, 1, 2, 1, 4, 3, 2,
            2, 1, 4, 1, 2, 0, 2, 2, 1, 3, 3, 2, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 1, 4, 2, 2, 2,
            3, 2, 2, 4, 1, 4, 1, 2, 2, 2, 4, 2, 2, 4, 4, 4, 4, 4, 1, 4, 1, 3, 1, 2, 2, 2, 2, 2, 4,
            2, 4, 1, 2, 2, 4, 3, 3, 2, 2, 4, 1, 2,
        ]
        .to_vec();
        read.reverse_complement_nucleotides();

        assert_eq!(
            read.letters,
            [
                1, 2, 3, 1, 1, 4, 4, 3, 1, 1, 2, 3, 1, 3, 1, 1, 1, 1, 1, 2, 4, 2, 3, 2, 3, 3, 3, 3,
                3, 1, 1, 3, 1, 1, 1, 2, 3, 2, 3, 1, 1, 4, 1, 1, 1, 3, 2, 1, 0, 0, 0, 0, 0, 0, 0, 0,
                0, 0, 1, 1, 4, 4, 2, 1, 1, 0, 1, 2, 3, 2, 1, 1, 4, 3, 2, 1, 2, 2, 1, 4, 4, 1, 4, 1,
                2, 2, 2, 1, 2, 2, 2, 1, 3, 3, 1, 1, 2, 2, 1, 1
            ]
        );
    }

    const FASTQ_RECORD: &str = concat!(
        "@HWI-EAS209_0006_FC706VJ:5:58:5894:21141#ATCACG/1\n",
        "TTAATTGGTAAATAAATCTCCTAATAGCTTAGATNTTACCTTNNNNNNNNNNTAGTTTCTTGAGATTTGTTGGGGGAGACATTTTTGTGATTGCCTTGAT\n",
        "+HWI-EAS209_0006_FC706VJ:5:58:5894:21141#ATCACG/1\n",
        "efcfffffcfeefffcffffffddf`feed]`]_B__^__[YBBBBBBBBBBRTT\\]][]dddd`ddd^dddadd^BBBBBBBBBBBBBBBBBBBBBBBB\n"
    );

    #[test]
    fn test_correct_read() -> Result<(), FastQFileError> {
        let mut reader = FastQFile::new(BufReader::new(FASTQ_RECORD.as_bytes()));
        let mut seq = FastQRead::default();

        assert_eq!(true, reader.read_next(&mut seq)?);

        assert_eq!(
            seq.title,
            "HWI-EAS209_0006_FC706VJ:5:58:5894:21141#ATCACG/1"
        );
        assert_eq!(
            seq.sub_title,
            "HWI-EAS209_0006_FC706VJ:5:58:5894:21141#ATCACG/1"
        );
        assert_eq!(
            seq.letters,
            //T  T  A  A  T  T  G  G  T  A  A  A  T  A  A  A  T  C  T  C  C  T  A  A  T  A  G  C
            [
                2, 2, 1, 1, 2, 2, 4, 4, 2, 1, 1, 1, 2, 1, 1, 1, 2, 3, 2, 3, 3, 2, 1, 1, 2, 1, 4, 3,
                //T  T  A  G  A  T  N  T  T  A  C  C  T  T  N  N  N  N  N  N  N  N  N  N  T  A  G  T
                2, 2, 1, 4, 1, 2, 0, 2, 2, 1, 3, 3, 2, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 1, 4, 2,
                //T  T  C  T  T  G  A  G  A  T  T  T  G  T  T  G  G  G  G  G  A  G  A  C  A  T  T  T
                2, 2, 3, 2, 2, 4, 1, 4, 1, 2, 2, 2, 4, 2, 2, 4, 4, 4, 4, 4, 1, 4, 1, 3, 1, 2, 2, 2,
                //T  T  G  T  G  A  T  T  G  C  C  T  T  G  A  T
                2, 2, 4, 2, 4, 1, 2, 2, 4, 3, 3, 2, 2, 4, 1, 2
            ]
        );
        assert_eq!(
            seq.qualities,
            //e   f   c   f   f   f   f   f   c   f   e   e   f   f   f   c   f   f   f   f   f
            [
                69, 70, 67, 70, 70, 70, 70, 70, 67, 70, 69, 69, 70, 70, 70, 67, 70, 70, 70, 70, 70,
                //f   d   d   f   `   f   e   e   d   ]   `   ]   _   B   _   _   ^   _   _   [   Y
                70, 68, 68, 70, 64, 70, 69, 69, 68, 61, 64, 61, 63, 34, 63, 63, 62, 63, 63, 59, 57,
                //B   B   B   B   B   B   B   B   B   B   R   T   T   \   ]   ]   [   ]   d   d   d
                34, 34, 34, 34, 34, 34, 34, 34, 34, 34, 50, 52, 52, 60, 61, 61, 59, 61, 68, 68, 68,
                //d   `   d   d   d   ^   d   d   d   a   d   d   ^   B   B   B   B   B   B   B   B
                68, 64, 68, 68, 68, 62, 68, 68, 68, 65, 68, 68, 62, 34, 34, 34, 34, 34, 34, 34, 34,
                //B   B   B   B   B   B   B   B   B   B   B   B   B   B   B   B
                34, 34, 34, 34, 34, 34, 34, 34, 34, 34, 34, 34, 34, 34, 34, 34
            ]
        );

        assert_eq!(false, reader.read_next(&mut seq)?);

        //let mut reader_rc = FastQFile::new(BufReader::new(FASTQ_RECORD.as_bytes()));
        //let mut seq_rc = FastQRead::default();

        //assert_eq!(true, reader_rc.read_next(&mut seq_rc)?);
        //seq_rc.reverse_complement_nucleotides();
        //assert_eq!(seq_rc.title, seq.title);
        //assert_eq!(seq_rc.sub_title, seq.sub_title);
        //assert_eq!(seq_rc.qualities, seq.qualities);
        //assert_eq!(seq_rc.letters, seq.letters);

        Ok(())
    }

    const FASTQ_RECORD_INVALID_SEQUENCE_LENGTH: &str = concat!(
        "@HWI-EAS209_0006_FC706VJ:5:58:5894:21141#ATCACG/1\n",
        "TTAATTGGTAAATAAATCTCCTAATAGCTTAGATNTTACCTTNNNNNNNNNNTAGTTTCTTGAGATTTGTTGGGGGAGACATTTTTGTGATTGCCTTGA\n",
        "+HWI-EAS209_0006_FC706VJ:5:58:5894:21141#ATCACG/1\n",
        "efcfffffcfeefffcffffffddf`feed]`]_B__^__[YBBBBBBBBBBRTT\\]][]dddd`ddd^dddadd^BBBBBBBBBBBBBBBBBBBBBBBB\n"
    );

    #[test]
    fn test_invalid_sequence_length() {
        let mut reader = FastQFile::new(BufReader::new(
            FASTQ_RECORD_INVALID_SEQUENCE_LENGTH.as_bytes(),
        ));
        let mut seq = FastQRead::default();

        let actual = reader.read_next(&mut seq);
        assert!(actual.is_err());
        assert!(matches!(
            actual.unwrap_err(),
            FastQFileError::MismatchedSequenceLength
        ));
    }

    const FASTQ_RECORD_INVALID_QUALITY_LENGTH: &str = concat!(
        "@HWI-EAS209_0006_FC706VJ:5:58:5894:21141#ATCACG/1\n",
        "TTAATTGGTAAATAAATCTCCTAATAGCTTAGATNTTACCTTNNNNNNNNNNTAGTTTCTTGAGATTTGTTGGGGGAGACATTTTTGTGATTGCCTTGAT\n",
        "+HWI-EAS209_0006_FC706VJ:5:58:5894:21141#ATCACG/1\n",
        "efcfffffcfeefffcffffffddf`feed]`]_B__^__[YBBBBBBBBBBRTT\\]][]dddd`ddd^dddadd^BBBBBBBBBBBBBBBBBBBBBBB\n"
    );

    #[test]
    fn test_invalid_quality_length() {
        let mut reader = FastQFile::new(BufReader::new(
            FASTQ_RECORD_INVALID_QUALITY_LENGTH.as_bytes(),
        ));
        let mut seq = FastQRead::default();

        let actual = reader.read_next(&mut seq);
        assert!(actual.is_err());
        assert!(matches!(
            actual.unwrap_err(),
            FastQFileError::MismatchedSequenceLength
        ));
    }

    const FASTQ_RECORD_NO_DESCRIPTION: &str = concat!(
        "@HWI-EAS209_0006_FC706VJ:5:58:5894:21141#ATCACG/1\n",
        "TTAATTGGTAAATAAATCTCCTAATAGCTTAGATNTTACCTTNNNNNNNNNNTAGTTTCTTGAGATTTGTTGGGGGAGACATTTTTGTGATTGCCTTGAT\n",
        "HWI-EAS209_0006_FC706VJ:5:58:5894:21141#ATCACG/1\n",
        "efcfffffcfeefffcffffffddf`feed]`]_B__^__[YBBBBBBBBBBRTT\\]][]dddd`ddd^dddadd^BBBBBBBBBBBBBBBBBBBBBBB\n"
    );

    #[test]
    fn test_no_description_line() {
        let mut reader = FastQFile::new(BufReader::new(FASTQ_RECORD_NO_DESCRIPTION.as_bytes()));
        let mut seq = FastQRead::default();

        let actual = reader.read_next(&mut seq);
        assert!(actual.is_err());
        assert!(matches!(
            actual.unwrap_err(),
            FastQFileError::NoDescriptionLine
        ));
    }

    const FASTQ_RECORD_NO_TITLE: &str = concat!(
        "HWI-EAS209_0006_FC706VJ:5:58:5894:21141#ATCACG/1\n",
        "TTAATTGGTAAATAAATCTCCTAATAGCTTAGATNTTACCTTNNNNNNNNNNTAGTTTCTTGAGATTTGTTGGGGGAGACATTTTTGTGATTGCCTTGAT\n",
        "+HWI-EAS209_0006_FC706VJ:5:58:5894:21141#ATCACG/1\n",
        "efcfffffcfeefffcffffffddf`feed]`]_B__^__[YBBBBBBBBBBRTT\\]][]dddd`ddd^dddadd^BBBBBBBBBBBBBBBBBBBBBBB\n"
    );

    #[test]
    fn test_no_title_line() {
        let mut reader = FastQFile::new(BufReader::new(FASTQ_RECORD_NO_TITLE.as_bytes()));
        let mut seq = FastQRead::default();

        let actual = reader.read_next(&mut seq);
        assert!(actual.is_err());
        assert!(matches!(actual.unwrap_err(), FastQFileError::NoTitleLine));
    }

    const FASTQ_RECORD_INVALID_NUCLEOTIDE: &str = concat!(
        "@HWI-EAS209_0006_FC706VJ:5:58:5894:21141#ATCACG/1\n",
        "TTAATTGGTAAATAAATCTCCTAATAGCTTAGATNTTACCTTNNNNNNzNNTAGTTTCTTGAGATTTGTTGGGGGAGACATTTTTGTGATTGCCTTGAT\n",
        "+HWI-EAS209_0006_FC706VJ:5:58:5894:21141#ATCACG/1\n",
        "efcfffffcfeefffcffffffddf`feed]`]_B__^__[YBBBBBBBBBBRTT\\]][]dddd`ddd^dddadd^BBBBBBBBBBBBBBBBBBBBBBB\n"
    );

    #[test]
    fn test_invalid_nucleotide() {
        let mut reader = FastQFile::new(BufReader::new(FASTQ_RECORD_INVALID_NUCLEOTIDE.as_bytes()));
        let mut seq = FastQRead::default();

        let actual = reader.read_next(&mut seq);
        assert!(actual.is_err());
        assert!(matches!(
            actual.unwrap_err(),
            FastQFileError::InvalidNucleotideLetter { c: 'z' }
        ));
    }

    const FASTQ_RECORD_INVALID_QUALITY: &str = concat!(
        "@HWI-EAS209_0006_FC706VJ:5:58:5894:21141#ATCACG/1\n",
        "TTAATTGGTAAATAAATCTCCTAATAGCTTAGATNTTACCTTNNNNNNNNNTAGTTTCTTGAGATTTGTTGGGGGAGACATTTTTGTGATTGCCTTGAT\n",
        "+HWI-EAS209_0006_FC706VJ:5:58:5894:21141#ATCACG/1\n",
        "efcfffffcfeefffcffffffddf`feed]`]_B__^__[YBBBBBB BBBRTT\\]][]dddd`ddd^dddadd^BBBBBBBBBBBBBBBBBBBBBBB\n"
    );

    #[test]
    fn test_invalid_quality() {
        let mut reader = FastQFile::new(BufReader::new(FASTQ_RECORD_INVALID_QUALITY.as_bytes()));
        let mut seq = FastQRead::default();

        let actual = reader.read_next(&mut seq);
        assert!(actual.is_err());
        assert!(matches!(
            actual.unwrap_err(),
            FastQFileError::InvalidQualityLetter
        ));
    }

    const FASTQ_RECORD_TRUNCATED_QUALITY: &str = concat!(
        "@HWI-EAS209_0006_FC706VJ:5:58:5894:21141#ATCACG/1\n",
        "TTAATTGGTAAATAAATCTCCTAATAGCTTAGATNTTACCTTNNNNNNNNNTAGTTTCTTGAGATTTGTTGGGGGAGACATTTTTGTGATTGCCTTGAT\n",
        "+HWI-EAS209_0006_FC706VJ:5:58:5894:21141#ATCACG/1\n",
    );

    const FASTQ_RECORD_TRUNCATED_DESCRIPTION: &str = concat!(
        "@HWI-EAS209_0006_FC706VJ:5:58:5894:21141#ATCACG/1\n",
        "TTAATTGGTAAATAAATCTCCTAATAGCTTAGATNTTACCTTNNNNNNNNNTAGTTTCTTGAGATTTGTTGGGGGAGACATTTTTGTGATTGCCTTGAT\n",
    );

    const FASTQ_RECORD_TRUNCATED_NUCLEOTIDES: &str =
        concat!("@HWI-EAS209_0006_FC706VJ:5:58:5894:21141#ATCACG/1\n",);

    #[test]
    fn test_truncated() {
        let mut reader = FastQFile::new(BufReader::new(FASTQ_RECORD_TRUNCATED_QUALITY.as_bytes()));
        let mut seq = FastQRead::default();

        let actual = reader.read_next(&mut seq);
        assert!(actual.is_err());
        assert!(matches!(
            actual.unwrap_err(),
            FastQFileError::IncompleteRecord
        ));

        let mut reader = FastQFile::new(BufReader::new(
            FASTQ_RECORD_TRUNCATED_DESCRIPTION.as_bytes(),
        ));
        let mut seq = FastQRead::default();

        let actual = reader.read_next(&mut seq);
        assert!(actual.is_err());
        assert!(matches!(
            actual.unwrap_err(),
            FastQFileError::IncompleteRecord
        ));

        let mut reader = FastQFile::new(BufReader::new(
            FASTQ_RECORD_TRUNCATED_NUCLEOTIDES.as_bytes(),
        ));
        let mut seq = FastQRead::default();

        let actual = reader.read_next(&mut seq);
        assert!(actual.is_err());
        assert!(matches!(
            actual.unwrap_err(),
            FastQFileError::IncompleteRecord
        ));
    }

    const FASTQ_RECORD_FASTA_RECORD: &str = concat!(
        "@HWI-EAS209_0006_FC706VJ:5:58:5894:21141#ATCACG/1\n",
        "TTAATTGGTAAATAAATCTCCTAATAGCTTAGATNTTACCTTNNNNNNNNNNTAGTTTCTTGAGATTTGTTGGGGGAGACATTTTTGTGATTGCCTTGAT\n",
        "+HWI-EAS209_0006_FC706VJ:5:58:5894:21141#ATCACG/1\n",
        "efcfffffcfeefffcffffffddf`feed]`]_B__^__[YBBBBBBBBBBRTT\\]][]dddd`ddd^dddadd^BBBBBBBBBBBBBBBBBBBBBBBB\n",
        ">HWI-EAS209_0006_FC706VJ:5:58:5894:21141#ATCACG/1\n",
        "TTAATTGGTAAATAAATCTCCTAATAGCTTAGATNTTACCTTNNNNNNNNNNTAGTTTCTTGAGATTTGTTGGGGGAGACATTTTTGTGATTGCCTTGAT\n",
    );

    #[test]
    fn test_fasta_record() -> Result<(), FastQFileError> {
        let mut reader = FastQFile::new(BufReader::new(FASTQ_RECORD_FASTA_RECORD.as_bytes()));
        let mut seq = FastQRead::default();

        assert_eq!(true, reader.read_next(&mut seq)?);
        assert_eq!(
            seq.title,
            "HWI-EAS209_0006_FC706VJ:5:58:5894:21141#ATCACG/1"
        );

        let actual = reader.read_next(&mut seq);
        assert!(actual.is_err());
        assert!(matches!(
            actual.unwrap_err(),
            FastQFileError::FastATitleLine
        ));

        Ok(())
    }
}
