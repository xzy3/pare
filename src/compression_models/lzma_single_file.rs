use thiserror::Error;

use std::fs::File;
use std::io;
use std::io::{BufWriter, Write};
use std::path::Path;

use xz2::write::XzEncoder;

use crate::seq_files::fastq::{FastQFileError, FastQRead, PairedFastQReader};

#[derive(Error, Debug)]
pub enum XZSingleFileError {
    #[error("IO error while reading fastq file")]
    IO {
        #[from]
        source: std::io::Error,
    },

    #[error("Invalid fastq file")]
    FastQError {
        #[from]
        source: FastQFileError,
    },
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

    pub fn write_string(&mut self, record: &String) -> io::Result<()> {
        self.write_u8(record.as_bytes())?;
        Ok(())
    }

    pub fn write_u8(&mut self, record: &[u8]) -> io::Result<()> {
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
