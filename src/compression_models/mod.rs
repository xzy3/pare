pub mod lzma_multi_stream;
pub mod lzma_single_file;

use std::string::FromUtf8Error;

use thiserror::Error;

use crate::seq_files::fastq::{FastQFileError, PairedFastQReader, PairedFastQWriter};

#[derive(Error, Debug)]
pub enum CompressionModelError {
    #[error("IO error while reading fastq file")]
    IO {
        #[from]
        source: std::io::Error,
    },
    #[error("bad value found file likely corrupted")]
    EncodingError {
        #[from]
        source: FromUtf8Error,
    },

    #[error("Invalid fastq file")]
    FastQError {
        #[from]
        source: FastQFileError,
    },

    #[error("EOF caused Incomplete record")]
    IncompleteRecord,
    #[error("Could not find the expected version string")]
    MissingVersion,
}

pub trait EncoderModel {
    fn compress(
        &mut self,
        reader: &mut Box<dyn PairedFastQReader>,
    ) -> Result<(), CompressionModelError>;
}

pub trait DecoderModel {
    fn decompress(
        &mut self,
        writer: &mut Box<dyn PairedFastQWriter>,
    ) -> Result<(), CompressionModelError>;
}
