pub mod lzma_multi_stream;
pub mod lzma_single_file;

use std::io::{Cursor, Read, Seek, SeekFrom, Write};
use std::string::FromUtf8Error;

use serde_json;
use tar::{Builder, Header};
use tempfile::SpooledTempFile;
use thiserror::Error;
use xz2::read::XzDecoder;
use xz2::write::XzEncoder;

use crate::seq_files::fastq::{FastQFileError, PairedFastQReader, PairedFastQWriter};

type Result<T> = std::result::Result<T, CompressionModelError>;

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

    #[error("metadata error in parsing JSON")]
    MetadataError {
        #[from]
        source: serde_json::Error,
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

pub struct PareArchiveEncoder<W: Write> {
    arch: Builder<W>,
}

impl<W: Write> PareArchiveEncoder<W> {
    pub fn new(sink: W) -> Self {
        PareArchiveEncoder {
            arch: Builder::new(sink),
        }
    }

    pub fn write_metadata(&mut self, metadata: serde_json::Value) -> Result<()> {
        let mut buffer: Vec<u8> = Vec::new();
        serde_json::to_writer(&mut buffer, &metadata)?;

        self.write_stream(&mut Cursor::new(buffer), "metadata")?;
        Ok(())
    }

    pub fn write_xz_spool(&mut self, spool: XzEncoder<SpooledTempFile>, path: &str) -> Result<()> {
        let mut finished_spool = spool.finish()?;
        self.write_stream(&mut finished_spool, &path)?;

        Ok(())
    }

    pub fn write_stream<T: Read + Seek>(&mut self, source: &mut T, path: &str) -> Result<()> {
        let mut header = Header::new_gnu();
        header.set_size(source.seek(SeekFrom::End(0))?);
        header.set_path(path)?;
        header.set_mode(0600);
        header.set_cksum();

        source.rewind()?;
        self.arch.append(&header, source)?;
        Ok(())
    }

    pub fn finish(&mut self) -> Result<()> {
        self.arch.finish()?;
        Ok(())
    }
}

pub trait EncoderModel {
    fn compress(&mut self, reader: &mut Box<dyn PairedFastQReader>) -> Result<()>;
}

pub trait DecoderModel {
    fn decompress(&mut self, writer: &mut Box<dyn PairedFastQWriter>) -> Result<()>;
}
