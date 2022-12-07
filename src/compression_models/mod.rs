pub mod lzma_multi_stream;
pub mod lzma_single_file;

use std::fs;
use std::fs::File;
use std::io::{Cursor, Read, Seek, SeekFrom, Write};
use std::string::FromUtf8Error;

use bson::Document;
use bson::{de, document, ser};
use tar::{Archive, Builder, Header};
use tempfile::{tempdir, SpooledTempFile, TempDir};
use thiserror::Error;
use xz2::read::XzDecoder;
use xz2::write::XzEncoder;

use crate::seq_files::fastq::{FastQFileError, PairedFastQReader, PairedFastQWriter};

type Result<T> = std::result::Result<T, CompressionModelError>;

#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord)]
enum CompressionModel {
    LZMASingle,
    LZMAMulti,
}

impl CompressionModel {
    fn as_str(&self) -> &'static str {
        match self {
            CompressionModel::LZMASingle => "lzma_single_stream",
            CompressionModel::LZMAMulti => "lzma_multi_stream",
        }
    }
}

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

    #[error("metadata error in parsing BSON")]
    MetadataError {
        #[from]
        source: de::Error,
    },

    #[error("metadata error in serializing BSON")]
    MetadataSerError {
        #[from]
        source: ser::Error,
    },

    #[error("error parsing metadata BSON")]
    MetadataAccessError {
        #[from]
        source: document::ValueAccessError,
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
    #[error("The wrong model was used to open the file")]
    OpenedWithWrongModel,
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

    pub fn write_metadata(&mut self, metadata: Document) -> Result<()> {
        let mut buffer: Vec<u8> = Vec::new();
        metadata.to_writer(&mut buffer)?;

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
        header.set_mode(0o600);
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

pub struct PareArchiveDecoder<R: Read> {
    arc: Archive<R>,
    tmpdir: TempDir,
}

impl<R: Read> PareArchiveDecoder<R> {
    pub fn new(source: R) -> Result<Self> {
        let dir = tempdir()?;
        let mut arc = Archive::<R>::new(source);

        arc.unpack(dir.path())?;
        Ok(PareArchiveDecoder {
            arc: arc,
            tmpdir: dir,
        })
    }

    pub fn get_stream(&mut self, path: &str) -> Result<File> {
        let file_path = self.tmpdir.path().join(path);
        Ok(File::open(file_path)?)
    }

    pub fn get_xz_stream(&mut self, path: &str) -> Result<XzDecoder<File>> {
        Ok(XzDecoder::new(self.get_stream(path)?))
    }

    pub fn get_metadata(&mut self) -> Result<Document> {
        //TODO: handle not finding the metadata file
        let cont = fs::read(self.tmpdir.path().join("metadata"))?;
        Ok(Document::from_reader(&mut Cursor::new(cont))?)
    }
}

pub trait EncoderModel {
    fn compress(&mut self, reader: &mut Box<dyn PairedFastQReader>) -> Result<()>;
}

pub trait DecoderModel {
    fn decompress(&mut self, writer: &mut Box<dyn PairedFastQWriter>) -> Result<()>;
}
