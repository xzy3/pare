use thiserror::Error;

use std::fs::File;
use std::io;
use std::io::{BufRead, BufReader, BufWriter, Read, Write};
use std::path::Path;
use std::string::FromUtf8Error;

use xz2::read::XzDecoder;
use xz2::write::XzEncoder;

use crate::seq_files::fastq::{FastQFileError, FastQRead, PairedFastQReader, PairedFastQWriter};

#[derive(Error, Debug)]
pub enum XZSingleFileError {
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
}

pub struct XZSingleFileWriter<W: Write> {
    encoder: BufWriter<XzEncoder<W>>,
}

impl<W: Write> XZSingleFileWriter<W> {
    pub fn new(sink: W) -> Self {
        XZSingleFileWriter {
            encoder: BufWriter::new(XzEncoder::<W>::new(sink, 9)),
        }
    }

    fn write_string(&mut self, record: &String) -> io::Result<()> {
        self.write_u8(record.as_bytes())?;
        Ok(())
    }

    fn write_u8(&mut self, record: &[u8]) -> io::Result<()> {
        //let len = record.len() as u16;
        //self.encoder.write(&len.to_be_bytes())?;
        self.encoder.write(record)?;
        self.encoder.write(b"\xFF")?;

        Ok(())
    }

    pub fn compress(
        &mut self,
        reader: &mut Box<dyn PairedFastQReader>,
    ) -> Result<(), XZSingleFileError> {
        let mut r1 = FastQRead::default();
        let mut r2 = FastQRead::default();

        loop {
            if !reader.read_next(&mut r1, &mut r2)? {
                break;
            }

            self.write_string(&r1.title)?;
            self.write_string(&r2.title)?;
            self.write_u8(&r1.letters)?;
            self.write_u8(&r2.letters)?;
            self.write_string(&r1.sub_title)?;
            self.write_string(&r2.sub_title)?;
            self.write_u8(&r1.qualities)?;
            self.write_u8(&r2.qualities)?;
        }
        Ok(())
    }
}

impl XZSingleFileWriter<std::io::Stdout> {
    pub fn to_stdout() -> Self {
        XZSingleFileWriter::new(std::io::stdout())
    }
}

impl XZSingleFileWriter<File> {
    pub fn create<P: AsRef<Path>>(path: &P) -> Result<Self, std::io::Error> {
        let file = File::create(path)?;
        Ok(XZSingleFileWriter::new(file))
    }
}

// readers
pub struct XZSingleFileReader<R: Read> {
    decoder: BufReader<XzDecoder<R>>,
}

impl<R: Read> XZSingleFileReader<R> {
    pub fn new(source: R) -> Self {
        XZSingleFileReader {
            decoder: BufReader::new(XzDecoder::<R>::new(source)),
        }
    }

    fn read_string(&mut self, record: &mut String) -> Result<bool, XZSingleFileError> {
        let mut buffer = vec![];

        let ret = self.read_u8(&mut buffer)?;
        record.clear();
        *record = String::from_utf8(buffer)?;
        Ok(ret)
    }

    fn read_u8(&mut self, record: &mut Vec<u8>) -> Result<bool, XZSingleFileError> {
        record.clear();
        if self.decoder.read_until(b'\xFF', record)? == 0 {
            return Ok(false);
        }

        match record.pop() {
            Some(b'\xFF') => {}
            _ => {
                return Err(XZSingleFileError::IncompleteRecord);
            }
        }

        Ok(true)
    }

    fn read_next(
        &mut self,
        r1: &mut FastQRead,
        r2: &mut FastQRead,
    ) -> Result<bool, XZSingleFileError> {
        if !self.read_string(&mut r1.title)? {
            return Ok(false);
        }

        if !self.read_string(&mut r2.title)? {
            return Err(XZSingleFileError::IncompleteRecord);
        }

        if !self.read_u8(&mut r1.letters)? {
            return Err(XZSingleFileError::IncompleteRecord);
        }

        if !self.read_u8(&mut r2.letters)? {
            return Err(XZSingleFileError::IncompleteRecord);
        }

        if !self.read_string(&mut r1.sub_title)? {
            return Err(XZSingleFileError::IncompleteRecord);
        }

        if !self.read_string(&mut r2.sub_title)? {
            return Err(XZSingleFileError::IncompleteRecord);
        }

        if !self.read_u8(&mut r1.qualities)? {
            return Err(XZSingleFileError::IncompleteRecord);
        }

        if !self.read_u8(&mut r2.qualities)? {
            return Err(XZSingleFileError::IncompleteRecord);
        }

        Ok(true)
    }

    pub fn decompress(
        &mut self,
        writer: &mut Box<dyn PairedFastQWriter>,
    ) -> Result<(), XZSingleFileError> {
        let mut r1 = FastQRead::default();
        let mut r2 = FastQRead::default();

        loop {
            if !self.read_next(&mut r1, &mut r2)? {
                break;
            }
            writer.write_next(&r1, &r2)?;
        }
        Ok(())
    }
}

impl XZSingleFileReader<std::io::Stdin> {
    pub fn from_stdin() -> Self {
        XZSingleFileReader::new(std::io::stdin())
    }
}

impl XZSingleFileReader<File> {
    pub fn open<P: AsRef<Path>>(path: &P) -> Result<Self, std::io::Error> {
        let file = File::open(path)?;
        Ok(XZSingleFileReader::new(file))
    }
}
