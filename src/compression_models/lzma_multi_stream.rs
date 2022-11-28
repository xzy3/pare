use std::fs::File;
use std::io::{BufRead, BufReader, Cursor, Read, Seek, Write};
use std::path::Path;

use serde_json::json;
use tar::{Builder, Header};
use tempfile::SpooledTempFile;
use xz2::read::XzDecoder;
use xz2::write::XzEncoder;

use crate::compression_models::*;
use crate::seq_files::fastq::{FastQRead, PairedFastQReader, PairedFastQWriter};

const FILE_VERSION: &'static [u8] = &*b"PARE lzma_multi_stream v1\xFF";

pub struct XZMultiStreamWriter<W: Write> {
    sink: PareArchiveEncoder<W>,
}

impl<W: Write> XZMultiStreamWriter<W> {
    pub fn new(sink: W) -> Self {
        XZMultiStreamWriter {
            sink: PareArchiveEncoder::<W>::new(sink),
        }
    }
}

impl<W: Write> EncoderModel for XZMultiStreamWriter<W> {
    fn compress(&mut self, reader: &mut Box<dyn PairedFastQReader>) -> Result<()> {
        let mut r1 = FastQRead::default();
        let mut r2 = FastQRead::default();

        let mut title_spool = XzEncoder::new(SpooledTempFile::new(4096), 9);
        let mut nucleotides_spool = XzEncoder::new(SpooledTempFile::new(4096), 9);
        let mut qualities_spool = XzEncoder::new(SpooledTempFile::new(4096), 9);

        self.sink.write_metadata(json!({
            "model": "lzma_multi_stream",
            "version": 1,
        }))?;

        loop {
            if !reader.read_next(&mut r1, &mut r2)? {
                break;
            }

            write!(title_spool, "{}\n", r1.title)?;
            write!(title_spool, "{}\n", r2.title)?;

            nucleotides_spool.write(&r1.letters)?;
            nucleotides_spool.write(b"\n")?;

            nucleotides_spool.write(&r2.letters)?;
            nucleotides_spool.write(b"\n")?;

            qualities_spool.write(&r1.qualities)?;
            qualities_spool.write(&r2.qualities)?;
        }

        self.sink.write_xz_spool(title_spool, "titles")?;
        self.sink.write_xz_spool(nucleotides_spool, "nucleotides")?;
        self.sink.write_xz_spool(qualities_spool, "qualities")?;
        self.sink.finish()?;

        Ok(())
    }
}

impl XZMultiStreamWriter<std::io::Stdout> {
    pub fn to_stdout() -> Self {
        XZMultiStreamWriter::new(std::io::stdout())
    }
}

impl XZMultiStreamWriter<File> {
    pub fn create<P: AsRef<Path>>(path: &P) -> Result<Self> {
        let file = File::create(path)?;
        Ok(XZMultiStreamWriter::new(file))
    }
}

// readers
pub struct XZMultiStreamReader<R: Read> {
    decoder: BufReader<XzDecoder<R>>,
}

impl<R: Read> XZMultiStreamReader<R> {
    pub fn new(source: R) -> Self {
        XZMultiStreamReader {
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

impl<R: Read> DecoderModel for XZMultiStreamReader<R> {
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

impl XZMultiStreamReader<std::io::Stdin> {
    pub fn from_stdin() -> Self {
        XZMultiStreamReader::new(std::io::stdin())
    }
}

impl XZMultiStreamReader<File> {
    pub fn open<P: AsRef<Path>>(path: &P) -> Result<Self> {
        let file = File::open(path)?;
        Ok(XZMultiStreamReader::new(file))
    }
}
