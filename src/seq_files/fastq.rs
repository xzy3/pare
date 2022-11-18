use thiserror::Error;

use std::io::BufReader;

use std::io::prelude::*;
use std::vec::Vec;

#[derive(Debug, Default)]
pub struct FastQRead {
    letters: Vec<u8>,
    qualities: Vec<u8>,
    title: String,
    sub_title: String,
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
    MismatchedSequnceLength,
    #[error("Found nucleotide that is not |ATCGNatcgn|")]
    InvalidNucleotideLetter,
    #[error("EOF caused Incomplete record")]
    IncompleteRecord,
}

/// Fastq file things
pub struct FastQFile<R: BufRead> {
    file_obj: R,
}

impl<R> FastQFile<R>
        where R: BufRead {

    pub fn read_next(&mut self, buf: &mut FastQRead) -> Result<bool, FastQFileError> {
        let mut title = String::new();
        let mut nucleotides = String::new();
        let mut sub_title = String::new();
        let mut quality_letters = String::new();

        if self.file_obj.read_line(&mut title)? == 0 {
            return Ok(false);
        }

        if !title.starts_with("@") {
            return Err(FastQFileError::NoTitleLine);
        }
        title = title[1..].trim_end().to_string();

        if self.file_obj.read_line(&mut nucleotides)? == 0 {
            return Err(FastQFileError::IncompleteRecord);
        }

        nucleotides = nucleotides.trim_end().to_string();
        let letters: Result<Vec<u8>, FastQFileError> = nucleotides.trim_end().chars().map(|n| {
            match n {
                'n' | 'N' => Ok(0),
                'a' | 'A' => Ok(1),
                't' | 'T' => Ok(2),
                'c' | 'C' => Ok(3),
                'g' | 'G' => Ok(4),
                _ => return Err(FastQFileError::InvalidNucleotideLetter),
            }
        }).collect();
        if letters.is_err() {
            return Err(letters.unwrap_err());
        }

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
            return Err(FastQFileError::MismatchedSequnceLength);
        }
        let qualities: Vec<u8> = quality_letters.bytes().map(|v| v-32).collect();

        *buf = FastQRead {
            letters: letters.unwrap(),
            qualities: qualities,
            title: title,
            sub_title: sub_title,
        };
        return Ok(true);
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    const FASTQ_RECORD: &str = concat!(
        "@HWI-EAS209_0006_FC706VJ:5:58:5894:21141#ATCACG/1\n",
        "TTAATTGGTAAATAAATCTCCTAATAGCTTAGATNTTACCTTNNNNNNNNNNTAGTTTCTTGAGATTTGTTGGGGGAGACATTTTTGTGATTGCCTTGAT\n",
        "+HWI-EAS209_0006_FC706VJ:5:58:5894:21141#ATCACG/1\n",
        "efcfffffcfeefffcffffffddf`feed]`]_Ba_^__[YBBBBBBBBBBRTT\\]][]dddd`ddd^dddadd^BBBBBBBBBBBBBBBBBBBBBBBB\n"
    );


    #[test]
    fn test_correct_read() -> Result<(), FastQFileError> {
        let mut reader = FastQFile { file_obj: BufReader::new(FASTQ_RECORD.as_bytes()) };
        let mut seq = FastQRead::default();

        assert_eq!(true, reader.read_next(&mut seq)?);

        assert_eq!(seq.title, "HWI-EAS209_0006_FC706VJ:5:58:5894:21141#ATCACG/1");
        assert_eq!(seq.sub_title, "HWI-EAS209_0006_FC706VJ:5:58:5894:21141#ATCACG/1");
        assert_eq!(seq.letters,
           //T  T  A  A  T  T  G  G  T  A  A  A  T  A  A  A  T  C  T  C  C  T  A  A  T  A  G  C
            [2, 2, 1, 1, 2, 2, 4, 4, 2, 1, 1, 1, 2, 1, 1, 1, 2, 3, 2, 3, 3, 2, 1, 1, 2, 1, 4, 3,
           //T  T  A  G  A  T  N  T  T  A  C  C  T  T  N  N  N  N  N  N  N  N  N  N  T  A  G  T
             2, 2, 1, 4, 1, 2, 0, 2, 2, 1, 3, 3, 2, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 1, 4, 2,
           //T  T  C  T  T  G  A  G  A  T  T  T  G  T  T  G  G  G  G  G  A  G  A  C  A  T  T  T
             2, 2, 3, 2, 2, 4, 1, 4, 1, 2, 2, 2, 4, 2, 2, 4, 4, 4, 4, 4, 1, 4, 1, 3, 1, 2, 2, 2,
           //T  T  G  T  G  A  T  T  G  C  C  T  T  G  A  T
             2, 2, 4, 2, 4, 1, 2, 2, 4, 3, 3, 2, 2, 4, 1, 2]);
        assert_eq!(seq.qualities,
            //e   f   c   f   f   f   f   f   c   f   e   e   f   f   f   c   f   f   f   f   f
            [69, 70, 67, 70, 70, 70, 70, 70, 67, 70, 69, 69, 70, 70, 70, 67, 70, 70, 70, 70, 70,
            //f   d   d   f   `   f   e   e   d   ]   `   ]   _   B   a   _   ^   _   _   [   Y
             70, 68, 68, 70, 64, 70, 69, 69, 68, 61, 64, 61, 63, 34, 65, 63, 62, 63, 63, 59, 57,
            //B   B   B   B   B   B   B   B   B   B   R   T   T   \   ]   ]   [   ]   d   d   d
             34, 34, 34, 34, 34, 34, 34, 34, 34, 34, 50, 52, 52, 60, 61, 61, 59, 61, 68, 68, 68,
            //d   `   d   d   d   ^   d   d   d   a   d   d   ^   B   B   B   B   B   B   B   B
             68, 64, 68, 68, 68, 62, 68, 68, 68, 65, 68, 68, 62, 34, 34, 34, 34, 34, 34, 34, 34,
            //B   B   B   B   B   B   B   B   B   B   B   B   B   B   B   B
             34, 34, 34, 34, 34, 34, 34, 34, 34, 34, 34, 34, 34, 34, 34, 34]
        );


        assert_eq!(false, reader.read_next(&mut seq)?);

        Ok(())
    }

    const FASTQ_RECORD_INVALID_SEQUENCE_LENGTH: &str = concat!(
        "@HWI-EAS209_0006_FC706VJ:5:58:5894:21141#ATCACG/1\n",
        "TTAATTGGTAAATAAATCTCCTAATAGCTTAGATNTTACCTTNNNNNNNNNNTAGTTTCTTGAGATTTGTTGGGGGAGACATTTTTGTGATTGCCTTGA\n",
        "+HWI-EAS209_0006_FC706VJ:5:58:5894:21141#ATCACG/1\n",
        "efcfffffcfeefffcffffffddf`feed]`]_Ba_^__[YBBBBBBBBBBRTT\\]][]dddd`ddd^dddadd^BBBBBBBBBBBBBBBBBBBBBBBB\n"
    );

    #[test]
    fn test_invalid_sequence_length() {
        let mut reader = FastQFile {
            file_obj: BufReader::new(FASTQ_RECORD_INVALID_SEQUENCE_LENGTH.as_bytes())
        };
        let mut seq = FastQRead::default();

        let actual = reader.read_next(&mut seq);
        assert!(actual.is_err());
        assert!(matches!(actual.unwrap_err(), FastQFileError::MismatchedSequnceLength));
    }

    const FASTQ_RECORD_INVALID_QUALITY_LENGTH: &str = concat!(
        "@HWI-EAS209_0006_FC706VJ:5:58:5894:21141#ATCACG/1\n",
        "TTAATTGGTAAATAAATCTCCTAATAGCTTAGATNTTACCTTNNNNNNNNNNTAGTTTCTTGAGATTTGTTGGGGGAGACATTTTTGTGATTGCCTTGAT\n",
        "+HWI-EAS209_0006_FC706VJ:5:58:5894:21141#ATCACG/1\n",
        "efcfffffcfeefffcffffffddf`feed]`]_Ba_^__[YBBBBBBBBBBRTT\\]][]dddd`ddd^dddadd^BBBBBBBBBBBBBBBBBBBBBBB\n"
    );

    #[test]
    fn test_invalid_quality_length() {
        let mut reader = FastQFile {
            file_obj: BufReader::new(FASTQ_RECORD_INVALID_QUALITY_LENGTH.as_bytes())
        };
        let mut seq = FastQRead::default();

        let actual = reader.read_next(&mut seq);
        assert!(actual.is_err());
        assert!(matches!(actual.unwrap_err(), FastQFileError::MismatchedSequnceLength));
    }

    const FASTQ_RECORD_NO_DESCRIPTION: &str = concat!(
        "@HWI-EAS209_0006_FC706VJ:5:58:5894:21141#ATCACG/1\n",
        "TTAATTGGTAAATAAATCTCCTAATAGCTTAGATNTTACCTTNNNNNNNNNNTAGTTTCTTGAGATTTGTTGGGGGAGACATTTTTGTGATTGCCTTGAT\n",
        "HWI-EAS209_0006_FC706VJ:5:58:5894:21141#ATCACG/1\n",
        "efcfffffcfeefffcffffffddf`feed]`]_Ba_^__[YBBBBBBBBBBRTT\\]][]dddd`ddd^dddadd^BBBBBBBBBBBBBBBBBBBBBBB\n"
    );

    #[test]
    fn test_no_description_line() {
        let mut reader = FastQFile {
            file_obj: BufReader::new(FASTQ_RECORD_NO_DESCRIPTION.as_bytes())
        };
        let mut seq = FastQRead::default();

        let actual = reader.read_next(&mut seq);
        assert!(actual.is_err());
        assert!(matches!(actual.unwrap_err(), FastQFileError::NoDescriptionLine));
    }

    const FASTQ_RECORD_NO_TITLE: &str = concat!(
        "HWI-EAS209_0006_FC706VJ:5:58:5894:21141#ATCACG/1\n",
        "TTAATTGGTAAATAAATCTCCTAATAGCTTAGATNTTACCTTNNNNNNNNNNTAGTTTCTTGAGATTTGTTGGGGGAGACATTTTTGTGATTGCCTTGAT\n",
        "+HWI-EAS209_0006_FC706VJ:5:58:5894:21141#ATCACG/1\n",
        "efcfffffcfeefffcffffffddf`feed]`]_Ba_^__[YBBBBBBBBBBRTT\\]][]dddd`ddd^dddadd^BBBBBBBBBBBBBBBBBBBBBBB\n"
    );

    #[test]
    fn test_no_title_line() {
        let mut reader = FastQFile {
            file_obj: BufReader::new(FASTQ_RECORD_NO_TITLE.as_bytes())
        };
        let mut seq = FastQRead::default();

        let actual = reader.read_next(&mut seq);
        assert!(actual.is_err());
        assert!(matches!(actual.unwrap_err(), FastQFileError::NoTitleLine));
    }

    const FASTQ_RECORD_INVALID_NUCLEOTIDE: &str = concat!(
        "@HWI-EAS209_0006_FC706VJ:5:58:5894:21141#ATCACG/1\n",
        "TTAATTGGTAAATAAATCTCCTAATAGCTTAGATNTTACCTTNNNNNNzNNTAGTTTCTTGAGATTTGTTGGGGGAGACATTTTTGTGATTGCCTTGAT\n",
        "+HWI-EAS209_0006_FC706VJ:5:58:5894:21141#ATCACG/1\n",
        "efcfffffcfeefffcffffffddf`feed]`]_Ba_^__[YBBBBBBBBBBRTT\\]][]dddd`ddd^dddadd^BBBBBBBBBBBBBBBBBBBBBBB\n"
    );

    #[test]
    fn test_invalid_nucleotide() {
        let mut reader = FastQFile {
            file_obj: BufReader::new(FASTQ_RECORD_INVALID_NUCLEOTIDE.as_bytes())
        };
        let mut seq = FastQRead::default();

        let actual = reader.read_next(&mut seq);
        assert!(actual.is_err());
        assert!(matches!(actual.unwrap_err(), FastQFileError::InvalidNucleotideLetter));
    }

    const FASTQ_RECORD_INVALID_QUALITY: &str = concat!(
        "@HWI-EAS209_0006_FC706VJ:5:58:5894:21141#ATCACG/1\n",
        "TTAATTGGTAAATAAATCTCCTAATAGCTTAGATNTTACCTTNNNNNNNNNTAGTTTCTTGAGATTTGTTGGGGGAGACATTTTTGTGATTGCCTTGAT\n",
        "+HWI-EAS209_0006_FC706VJ:5:58:5894:21141#ATCACG/1\n",
        "efcfffffcfeefffcffffffddf`feed]`]_Ba_^__[YBBBBBB BBBRTT\\]][]dddd`ddd^dddadd^BBBBBBBBBBBBBBBBBBBBBBB\n"
    );

    #[test]
    fn test_invalid_quality() {
        let mut reader = FastQFile {
            file_obj: BufReader::new(FASTQ_RECORD_INVALID_QUALITY.as_bytes())
        };
        let mut seq = FastQRead::default();

        let actual = reader.read_next(&mut seq);
        assert!(actual.is_err());
        assert!(matches!(actual.unwrap_err(), FastQFileError::InvalidQualityLetter));
    }

    const FASTQ_RECORD_TRUNCATED: &str = concat!( "@HWI-EAS209_0006_FC706VJ:5:58:5894:21141#ATCACG/1\n",
        "TTAATTGGTAAATAAATCTCCTAATAGCTTAGATNTTACCTTNNNNNNNNNTAGTTTCTTGAGATTTGTTGGGGGAGACATTTTTGTGATTGCCTTGAT\n",
        "+HWI-EAS209_0006_FC706VJ:5:58:5894:21141#ATCACG/1\n",
    );

    #[test]
    fn test_truncated() {
        let mut reader = FastQFile {
            file_obj: BufReader::new(FASTQ_RECORD_TRUNCATED.as_bytes())
        };
        let mut seq = FastQRead::default();

        let actual = reader.read_next(&mut seq);
        assert!(actual.is_err());
        assert!(matches!(actual.unwrap_err(), FastQFileError::IncompleteRecord));
    }
}
