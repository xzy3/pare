use thiserror::Error;

use std::fs::File;
use std::io::{BufReader, BufWriter};
use std::path::Path;

use std::io::prelude::*;
use std::vec::Vec;

fn reverse_complement_nucleotides(nucleotides: &mut Vec<u8>) {
    nucleotides.reverse();
    nucleotides.iter_mut().for_each(|n| {
        *n = match n {
            b'n' => b'n',
            b'a' => b't',
            b't' => b'a',
            b'c' => b'g',
            b'g' => b'c',
            _ => panic!("Invalid nuceotide {} found!", n),
        }
    });
}

fn nuclotides_upper(nucleotides: &mut Vec<u8>) {
    nucleotides.iter_mut().for_each(|n| {
        *n = match n {
            b'n' => b'N',
            b'a' => b'A',
            b't' => b'T',
            b'c' => b'C',
            b'g' => b'G',
            _ => panic!("Invalid nuclotide {} found!", n),
        }
    });
}

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
        reverse_complement_nucleotides(&mut self.letters);
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
    NoTitleLine { line: u32 },
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
            'n' | 'N' => ret.push(b'n'),
            'a' | 'A' => ret.push(b'a'),
            't' | 'T' => ret.push(b't'),
            'c' | 'C' => ret.push(b'c'),
            'g' | 'G' => ret.push(b'g'),
            _ => return Err(FastQFileError::InvalidNucleotideLetter { c: n }),
        }
    }

    Ok(ret)
}

/// Fastq file things
pub trait FastQFileReaderTrait {
    fn read_next(&mut self, buf: &mut FastQRead) -> Result<bool, FastQFileError>;
}

#[derive(Debug)]
pub struct FastQFileReader<R: Read> {
    stream: BufReader<R>,
    line: u32,
}

impl<R: Read> FastQFileReader<R> {
    pub fn new(stream: BufReader<R>) -> Self {
        FastQFileReader {
            stream: stream,
            line: 0,
        }
    }
}

impl<R: Read> FastQFileReaderTrait for FastQFileReader<R> {
    fn read_next(&mut self, buf: &mut FastQRead) -> Result<bool, FastQFileError> {
        let mut title = String::new();
        let mut nucleotides = String::new();
        let mut sub_title = String::new();
        let mut quality_letters = String::new();

        loop {
            if self.stream.read_line(&mut title)? == 0 {
                return Ok(false);
            }

            // ignore blank lines
            if !title.trim_end().is_empty() {
                break;
            }
        }

        if !title.starts_with("@") {
            if title.starts_with(">") {
                return Err(FastQFileError::FastATitleLine);
            }
            return Err(FastQFileError::NoTitleLine { line: self.line });
        }
        title = title[1..].trim_end().to_string();

        if self.stream.read_line(&mut nucleotides)? == 0 {
            return Err(FastQFileError::IncompleteRecord);
        }

        nucleotides = nucleotides.trim_end().to_string();
        let letters: Vec<u8> = nuc_string_to_vec(&nucleotides)?;
        if self.stream.read_line(&mut sub_title)? == 0 {
            return Err(FastQFileError::IncompleteRecord);
        }

        if !sub_title.starts_with("+") {
            return Err(FastQFileError::NoDescriptionLine);
        }
        sub_title = sub_title[1..].trim_end().to_string();

        if self.stream.read_line(&mut quality_letters)? == 0 {
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
        self.line += 4;
        return Ok(true);
    }
}

impl FastQFileReader<File> {
    pub fn open<P: AsRef<Path>>(path: &P) -> Result<Self, std::io::Error> {
        let file = File::open(path)?;
        Ok(FastQFileReader::new(BufReader::new(file)))
    }
}

impl FastQFileReader<std::io::Stdin> {
    pub fn from_stdin() -> Self {
        FastQFileReader::new(BufReader::new(std::io::stdin()))
    }
}

pub trait PairedFastQReader {
    fn read_next(
        &mut self,
        buf_r1: &mut FastQRead,
        buf_r2: &mut FastQRead,
    ) -> Result<bool, FastQFileError>;
}

pub struct FastQPairedFilesReader {
    r1_reader: Box<dyn FastQFileReaderTrait>,
    r2_reader: Box<dyn FastQFileReaderTrait>,
    reverse_complement_r2_nucleotides: bool,
}

impl FastQPairedFilesReader {
    pub fn new(
        stream_r1: Box<dyn FastQFileReaderTrait>,
        stream_r2: Box<dyn FastQFileReaderTrait>,
        reverse_complement_r2_nucleotides: bool,
    ) -> Self {
        FastQPairedFilesReader {
            r1_reader: stream_r1,
            r2_reader: stream_r2,
            reverse_complement_r2_nucleotides: reverse_complement_r2_nucleotides,
        }
    }
}

impl PairedFastQReader for FastQPairedFilesReader {
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

pub struct FastQInterleavedFileReader {
    reader: Box<dyn FastQFileReaderTrait>,
    reverse_complement_r2_nucleotides: bool,
}

impl FastQInterleavedFileReader {
    pub fn new(
        stream: Box<dyn FastQFileReaderTrait>,
        reverse_complement_r2_nucleotides: bool,
    ) -> Self {
        FastQInterleavedFileReader {
            reader: stream,
            reverse_complement_r2_nucleotides: reverse_complement_r2_nucleotides,
        }
    }
}

impl PairedFastQReader for FastQInterleavedFileReader {
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

//writers
pub trait FastQFileWriterTrait {
    fn write_next(
        &mut self,
        buf: FastQRead,
        reverse_complement: bool,
    ) -> Result<bool, FastQFileError>;
}

pub struct FastQFileWriter<W: Write> {
    stream: BufWriter<W>,
    line: u32,
}

impl<W: Write> FastQFileWriter<W> {
    pub fn new(stream: BufWriter<W>) -> Self {
        FastQFileWriter {
            stream: stream,
            line: 0,
        }
    }
}

impl<W: Write> FastQFileWriterTrait for FastQFileWriter<W> {
    fn write_next(
        &mut self,
        buf: FastQRead,
        reverse_complement: bool,
    ) -> Result<bool, FastQFileError> {
        write!(self.stream, "@{}\n", buf.title)?;

        let mut letters = buf.letters.to_owned();

        if reverse_complement {
            reverse_complement_nucleotides(&mut letters);
        }
        nuclotides_upper(&mut letters);
        self.stream.write(&letters)?;

        self.stream.write(b"\n")?;
        write!(self.stream, "+{}\n", buf.sub_title)?;

        let quals: Vec<u8> = buf.qualities.iter().map(|q| q + 32).collect();

        self.stream.write(&quals)?;
        self.stream.write(b"\n")?;
        self.line += 4;
        return Ok(true);
    }
}

impl FastQFileWriter<File> {
    pub fn create<P: AsRef<Path>>(path: &P) -> Result<Self, std::io::Error> {
        let file = File::create(path)?;
        Ok(FastQFileWriter::new(BufWriter::new(file)))
    }
}

impl FastQFileWriter<std::io::Stdout> {
    pub fn to_stdout() -> Self {
        FastQFileWriter::new(BufWriter::new(std::io::stdout()))
    }
}

pub trait PairedFastQWriter {
    fn write_next(&mut self, buf_r1: FastQRead, buf_r2: FastQRead) -> Result<bool, FastQFileError>;
}

pub struct FastQPairedFilesWriter {
    r1_stream: Box<dyn FastQFileWriterTrait>,
    r2_stream: Box<dyn FastQFileWriterTrait>,
    reverse_complement_r2_nucleotides: bool,
}

impl FastQPairedFilesWriter {
    pub fn new(
        stream_r1: Box<dyn FastQFileWriterTrait>,
        stream_r2: Box<dyn FastQFileWriterTrait>,
        reverse_complement_r2_nucleotides: bool,
    ) -> Self {
        FastQPairedFilesWriter {
            r1_stream: stream_r1,
            r2_stream: stream_r2,
            reverse_complement_r2_nucleotides: reverse_complement_r2_nucleotides,
        }
    }
}

impl PairedFastQWriter for FastQPairedFilesWriter {
    fn write_next(&mut self, buf_r1: FastQRead, buf_r2: FastQRead) -> Result<bool, FastQFileError> {
        self.r1_stream.write_next(buf_r1, false)?;
        self.r2_stream
            .write_next(buf_r2, self.reverse_complement_r2_nucleotides)?;
        Ok(true)
    }
}

pub struct FastQInterleavedFileWriter {
    stream: Box<dyn FastQFileWriterTrait>,
    reverse_complement_r2_nucleotides: bool,
}

impl FastQInterleavedFileWriter {
    pub fn new(
        stream: Box<dyn FastQFileWriterTrait>,
        reverse_complement_r2_nucleotides: bool,
    ) -> Self {
        FastQInterleavedFileWriter {
            stream: stream,
            reverse_complement_r2_nucleotides: reverse_complement_r2_nucleotides,
        }
    }
}

impl PairedFastQWriter for FastQInterleavedFileWriter {
    fn write_next(&mut self, buf_r1: FastQRead, buf_r2: FastQRead) -> Result<bool, FastQFileError> {
        self.stream.write_next(buf_r1, false)?;
        self.stream
            .write_next(buf_r2, self.reverse_complement_r2_nucleotides)?;

        Ok(true)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_reverse_complement_nucleotides() {
        let mut read = FastQRead::default();
        read.letters = b"ttaattggtaaataaatctcctaatagcttagatnttaccttnnnnnnnnnntagtttcttgagatttgttgggggagacatttttgtgattgccttgat".to_vec();

        read.reverse_complement_nucleotides();
        assert_eq!(
            read.letters,
            b"atcaaggcaatcacaaaaatgtctcccccaacaaatctcaagaaactannnnnnnnnnaaggtaanatctaagctattaggagatttatttaccaattaa"
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
        let mut reader = FastQFileReader::new(BufReader::new(FASTQ_RECORD.as_bytes()));
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
            b"ttaattggtaaataaatctcctaatagcttagatnttaccttnnnnnnnnnntagtttcttgagatttgttgggggagacatttttgtgattgccttgat"
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
        let mut reader = FastQFileReader::new(BufReader::new(
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
        let mut reader = FastQFileReader::new(BufReader::new(
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
        let mut reader =
            FastQFileReader::new(BufReader::new(FASTQ_RECORD_NO_DESCRIPTION.as_bytes()));
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
        let mut reader = FastQFileReader::new(BufReader::new(FASTQ_RECORD_NO_TITLE.as_bytes()));
        let mut seq = FastQRead::default();

        let actual = reader.read_next(&mut seq);
        assert!(actual.is_err());
        assert!(matches!(
            actual.unwrap_err(),
            FastQFileError::NoTitleLine { line: 0 }
        ));
    }

    const FASTQ_RECORD_INVALID_NUCLEOTIDE: &str = concat!(
        "@HWI-EAS209_0006_FC706VJ:5:58:5894:21141#ATCACG/1\n",
        "TTAATTGGTAAATAAATCTCCTAATAGCTTAGATNTTACCTTNNNNNNzNNTAGTTTCTTGAGATTTGTTGGGGGAGACATTTTTGTGATTGCCTTGAT\n",
        "+HWI-EAS209_0006_FC706VJ:5:58:5894:21141#ATCACG/1\n",
        "efcfffffcfeefffcffffffddf`feed]`]_B__^__[YBBBBBBBBBBRTT\\]][]dddd`ddd^dddadd^BBBBBBBBBBBBBBBBBBBBBBB\n"
    );

    #[test]
    fn test_invalid_nucleotide() {
        let mut reader =
            FastQFileReader::new(BufReader::new(FASTQ_RECORD_INVALID_NUCLEOTIDE.as_bytes()));
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
        let mut reader =
            FastQFileReader::new(BufReader::new(FASTQ_RECORD_INVALID_QUALITY.as_bytes()));
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
        let mut reader =
            FastQFileReader::new(BufReader::new(FASTQ_RECORD_TRUNCATED_QUALITY.as_bytes()));
        let mut seq = FastQRead::default();

        let actual = reader.read_next(&mut seq);
        assert!(actual.is_err());
        assert!(matches!(
            actual.unwrap_err(),
            FastQFileError::IncompleteRecord
        ));

        let mut reader = FastQFileReader::new(BufReader::new(
            FASTQ_RECORD_TRUNCATED_DESCRIPTION.as_bytes(),
        ));
        let mut seq = FastQRead::default();

        let actual = reader.read_next(&mut seq);
        assert!(actual.is_err());
        assert!(matches!(
            actual.unwrap_err(),
            FastQFileError::IncompleteRecord
        ));

        let mut reader = FastQFileReader::new(BufReader::new(
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
        let mut reader = FastQFileReader::new(BufReader::new(FASTQ_RECORD_FASTA_RECORD.as_bytes()));
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

    // test the pair iterator items. The parsing is tested by the FastQFileReader tests so just make sure
    // they get the right files.
    const FASTQ_RECORD_INTERLEAVED: &str = concat!(
        "@HWI-EAS209_0006_FC706VJ:5:58:5894:21141#ATCACG/1\n",
        "TTAATTGGTAAATAAATCTCCTAATAGCTTAGATNTTACCTTNNNNNNNNNNTAGTTTCTTGAGATTTGTTGGGGGAGACATTTTTGTGATTGCCTTGAT\n",
        "+\n",
        "efcfffffcfeefffcffffffddf`feed]`]_B__^__[YBBBBBBBBBBRTT\\]][]dddd`ddd^dddadd^BBBBBBBBBBBBBBBBBBBBBBBB\n",
        "@HWI-EAS209_0006_FC706VJ:5:58:5894:21141#ATCACG/2\n",
        "TTAATTGGTAAATAAATCTCCTAATAGCTTAGATNTTACCTTNNNNNNNNNNTAGTTTCTTGAGATTTGTTGGGGGAGACATTTTTGTGATTGCCTTGAT\n",
        "+\n",
        "efcfffffcfeefffcffffffddf`feed]`]_B__^__[YBBBBBBBBBBRTT\\]][]dddd`ddd^dddadd^BBBBBBBBBBBBBBBBBBBBBBBB\n"
    );

    #[test]
    fn test_fastq_interleaved_file_read() -> Result<(), FastQFileError> {
        let str_reader = Box::new(FastQFileReader::new(BufReader::new(
            FASTQ_RECORD_INTERLEAVED.as_bytes(),
        )));
        let mut seq1 = FastQRead::default();
        let mut seq2 = FastQRead::default();

        let mut reader = FastQInterleavedFileReader::new(str_reader, true);
        reader.read_next(&mut seq1, &mut seq2)?;

        assert_eq!(
            seq1.title,
            "HWI-EAS209_0006_FC706VJ:5:58:5894:21141#ATCACG/1"
        );
        assert_eq!(
            seq2.title,
            "HWI-EAS209_0006_FC706VJ:5:58:5894:21141#ATCACG/2"
        );
        seq2.reverse_complement_nucleotides();
        assert_eq!(seq1.letters, seq2.letters);

        Ok(())
    }

    const FASTQ_RECORD_PAIR_R1: &str = concat!(
        "@HWI-EAS209_0006_FC706VJ:5:58:5894:21141#ATCACG/1\n",
        "TTAATTGGTAAATAAATCTCCTAATAGCTTAGATNTTACCTTNNNNNNNNNNTAGTTTCTTGAGATTTGTTGGGGGAGACATTTTTGTGATTGCCTTGAT\n",
        "+\n",
        "efcfffffcfeefffcffffffddf`feed]`]_B__^__[YBBBBBBBBBBRTT\\]][]dddd`ddd^dddadd^BBBBBBBBBBBBBBBBBBBBBBBB\n"
    );

    const FASTQ_RECORD_PAIR_R2: &str = concat!(
        "@HWI-EAS209_0006_FC706VJ:5:58:5894:21141#ATCACG/2\n",
        "TTAATTGGTAAATAAATCTCCTAATAGCTTAGATNTTACCTTNNNNNNNNNNTAGTTTCTTGAGATTTGTTGGGGGAGACATTTTTGTGATTGCCTTGAT\n",
        "+\n",
        "efcfffffcfeefffcffffffddf`feed]`]_B__^__[YBBBBBBBBBBRTT\\]][]dddd`ddd^dddadd^BBBBBBBBBBBBBBBBBBBBBBBB\n"
    );

    #[test]
    fn test_fastq_paired_files_read() -> Result<(), FastQFileError> {
        let str_reader1 = Box::new(FastQFileReader::new(BufReader::new(
            FASTQ_RECORD_PAIR_R1.as_bytes(),
        )));
        let str_reader2 = Box::new(FastQFileReader::new(BufReader::new(
            FASTQ_RECORD_PAIR_R2.as_bytes(),
        )));

        let mut reader = FastQPairedFilesReader::new(str_reader1, str_reader2, true);
        let mut seq1 = FastQRead::default();
        let mut seq2 = FastQRead::default();

        reader.read_next(&mut seq1, &mut seq2)?;

        assert_eq!(
            seq1.title,
            "HWI-EAS209_0006_FC706VJ:5:58:5894:21141#ATCACG/1"
        );
        assert_eq!(
            seq2.title,
            "HWI-EAS209_0006_FC706VJ:5:58:5894:21141#ATCACG/2"
        );
        seq2.reverse_complement_nucleotides();
        assert_eq!(seq1.letters, seq2.letters);

        Ok(())
    }

    #[test]
    fn test_correct_write() -> Result<(), FastQFileError> {
        let mut reader = FastQFileReader::new(BufReader::new(FASTQ_RECORD.as_bytes()));
        let mut seq = FastQRead::default();

        assert_eq!(true, reader.read_next(&mut seq)?);
        let buf: Vec<u8> = Vec::with_capacity(FASTQ_RECORD.len());

        let mut writer = FastQFileWriter::new(BufWriter::new(buf));
        writer.write_next(seq, false)?;

        let result = writer.stream.into_inner().unwrap();
        assert_eq!(FASTQ_RECORD.as_bytes(), result);

        Ok(())
    }

    #[test]
    fn test_fastq_paired_files_write() -> Result<(), FastQFileError> {
        let str_reader1 = Box::new(FastQFileReader::new(BufReader::new(
            FASTQ_RECORD_PAIR_R1.as_bytes(),
        )));
        let str_reader2 = Box::new(FastQFileReader::new(BufReader::new(
            FASTQ_RECORD_PAIR_R2.as_bytes(),
        )));

        let mut reader = FastQPairedFilesReader::new(str_reader1, str_reader2, true);
        let mut seq1 = FastQRead::default();
        let mut seq2 = FastQRead::default();

        reader.read_next(&mut seq1, &mut seq2)?;

        let buf1: Vec<u8> = Vec::with_capacity(FASTQ_RECORD_PAIR_R1.len());
        let buf2: Vec<u8> = Vec::with_capacity(FASTQ_RECORD_PAIR_R2.len());

        let mut writer = FastQPairedFilesWriter::new(
            Box::new(FastQFileWriter::new(BufWriter::new(buf1))),
            Box::new(FastQFileWriter::new(BufWriter::new(buf2))),
            true,
        );

        writer.write_next(seq1, seq2)?;
        // really need to find a way to get to the underlying data. the type checker/lifetime
        // checker isn't happy with normal methods.

        Ok(())
    }

    #[test]
    fn test_fastq_interleaved_file_write() -> Result<(), FastQFileError> {
        let str_reader = Box::new(FastQFileReader::new(BufReader::new(
            FASTQ_RECORD_INTERLEAVED.as_bytes(),
        )));
        let mut seq1 = FastQRead::default();
        let mut seq2 = FastQRead::default();

        let mut reader = FastQInterleavedFileReader::new(str_reader, true);
        reader.read_next(&mut seq1, &mut seq2)?;

        let buf: Vec<u8> = Vec::with_capacity(FASTQ_RECORD_INTERLEAVED.len());
        let mut writer = FastQInterleavedFileWriter::new(
            Box::new(FastQFileWriter::new(BufWriter::new(buf))),
            true,
        );

        writer.write_next(seq1, seq2)?;

        Ok(())
    }
}
