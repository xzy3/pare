use std::fs::File;
use std::io::{BufRead, BufReader, Read, Write};
use std::path::Path;

use serde_json::json;
use xz2::read::XzDecoder;
use xz2::write::XzEncoder;

use crate::compression_models::*;
use crate::seq_files::fastq::{FastQRead, PairedFastQReader, PairedFastQWriter};

const FILE_VERSION: &'static [u8] = &*b"PARE lzma_single_file v1\xFF";

pub struct XZSingleFileWriter<W: Write> {
    sink: PareArchiveEncoder<W>,
}

impl<W: Write> XZSingleFileWriter<W> {
    pub fn new(sink: W) -> Self {
        XZSingleFileWriter {
            sink: PareArchiveEncoder::<W>::new(sink),
        }
    }

    fn write_string(
        &mut self,
        spool: &mut XzEncoder<SpooledTempFile>,
        record: &String,
    ) -> Result<()> {
        self.write_u8(spool, record.as_bytes())?;
        Ok(())
    }

    fn write_u8(&mut self, spool: &mut XzEncoder<SpooledTempFile>, record: &[u8]) -> Result<()> {
        spool.write(record)?;
        spool.write(b"\xFF")?;

        Ok(())
    }
}

impl<W: Write> EncoderModel for XZSingleFileWriter<W> {
    fn compress(&mut self, reader: &mut Box<dyn PairedFastQReader>) -> Result<()> {
        let mut r1 = FastQRead::default();
        let mut r2 = FastQRead::default();

        let mut spool = XzEncoder::new(SpooledTempFile::new(4096), 9);

        self.sink.write_metadata(json!({
            "model": "lzma_single_stream",
            "version": 1,
        }))?;

        loop {
            if !reader.read_next(&mut r1, &mut r2)? {
                break;
            }

            self.write_string(&mut spool, &r1.title)?;
            self.write_string(&mut spool, &r2.title)?;
            self.write_u8(&mut spool, &r1.letters)?;
            self.write_u8(&mut spool, &r2.letters)?;
            spool.write(&r1.qualities)?;
            spool.write(&r2.qualities)?;
        }
        self.sink.write_stream(&mut spool.finish()?, "data")?;
        Ok(())
    }
}

impl XZSingleFileWriter<std::io::Stdout> {
    pub fn to_stdout() -> Self {
        XZSingleFileWriter::new(std::io::stdout())
    }
}

impl XZSingleFileWriter<File> {
    pub fn create<P: AsRef<Path>>(path: &P) -> Result<Self> {
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

    fn check_magic(&mut self) -> Result<()> {
        let mut buffer = vec![];
        if self.decoder.read_until(b'\xFF', &mut buffer)? == 0 {
            return Err(CompressionModelError::MissingVersion);
        }

        if buffer != FILE_VERSION {
            return Err(CompressionModelError::MissingVersion);
        }

        Ok(())
    }

    fn read_string(&mut self, record: &mut String) -> Result<bool> {
        let mut buffer = vec![];

        let ret = self.read_u8(&mut buffer)?;
        record.clear();
        *record = String::from_utf8(buffer)?;
        Ok(ret)
    }

    fn read_u8(&mut self, record: &mut Vec<u8>) -> Result<bool> {
        record.clear();
        if self.decoder.read_until(b'\xFF', record)? == 0 {
            return Ok(false);
        }

        match record.pop() {
            Some(b'\xFF') => {}
            _ => {
                return Err(CompressionModelError::IncompleteRecord);
            }
        }

        Ok(true)
    }

    fn read_next(&mut self, r1: &mut FastQRead, r2: &mut FastQRead) -> Result<bool> {
        if !self.read_string(&mut r1.title)? {
            return Ok(false);
        }

        if !self.read_string(&mut r2.title)? {
            return Err(CompressionModelError::IncompleteRecord);
        }

        if !self.read_u8(&mut r1.letters)? {
            return Err(CompressionModelError::IncompleteRecord);
        }

        if !self.read_u8(&mut r2.letters)? {
            return Err(CompressionModelError::IncompleteRecord);
        }

        r1.qualities.clear();
        r1.qualities.resize(r1.letters.len(), 0);

        match self.decoder.read_exact(&mut r1.qualities[..]) {
            Ok(()) => {}
            _ => {
                return Err(CompressionModelError::IncompleteRecord);
            }
        }

        r2.qualities.clear();
        r2.qualities.resize(r2.letters.len(), 0);
        match self.decoder.read_exact(&mut r2.qualities[..]) {
            Ok(()) => {}
            _ => {
                return Err(CompressionModelError::IncompleteRecord);
            }
        }

        Ok(true)
    }
}

impl<R: Read> DecoderModel for XZSingleFileReader<R> {
    fn decompress(&mut self, writer: &mut Box<dyn PairedFastQWriter>) -> Result<()> {
        let mut r1 = FastQRead::default();
        let mut r2 = FastQRead::default();

        self.check_magic()?;

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
    pub fn open<P: AsRef<Path>>(path: &P) -> Result<Self> {
        let file = File::open(path)?;
        Ok(XZSingleFileReader::new(file))
    }
}
