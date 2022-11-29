use std::fs::File;
use std::io::{BufRead, BufReader, Cursor, Read, Seek, Write};
use std::path::Path;

use serde_json::json;
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
    arc: PareArchiveDecoder<R>,
}

type XzBufReader = BufReader<XzDecoder<File>>;

impl<R: Read> XZMultiStreamReader<R> {
    pub fn new(source: R) -> Result<Self> {
        let arc = PareArchiveDecoder::<R>::new(source)?;

        Ok(XZMultiStreamReader { arc: arc })
    }

    fn check_magic(&mut self) -> Result<()> {
        let metadata = self.arc.get_metadata()?;
        if metadata["model"] != "lzma_multi_stream" || metadata["version"] != 1 {
            return Err(CompressionModelError::OpenedWithWrongModel);
        }

        Ok(())
    }

    fn read_line(&mut self, source: &mut XzBufReader, record: &mut String) -> Result<bool> {
        record.clear();
        if source.read_line(record)? == 0 {
            return Ok(false);
        }

        record.pop();

        Ok(true)
    }

    fn read_u8(&mut self, source: &mut XzBufReader, record: &mut Vec<u8>) -> Result<bool> {
        record.clear();
        if source.read_until(b'\n', record)? == 0 {
            return Ok(false);
        }

        match record.pop() {
            Some(b'\n') => {}
            _ => {
                return Err(CompressionModelError::IncompleteRecord);
            }
        }

        Ok(true)
    }

    fn read_exact(
        &mut self,
        source: &mut XzDecoder<File>,
        l: usize,
        record: &mut Vec<u8>,
    ) -> Result<bool> {
        record.clear();
        record.resize(l, 0);

        source.read_exact(&mut record[..])?;
        Ok(true)
    }
}

impl<R: Read> DecoderModel for XZMultiStreamReader<R> {
    fn decompress(&mut self, writer: &mut Box<dyn PairedFastQWriter>) -> Result<()> {
        let mut r1 = FastQRead::default();
        let mut r2 = FastQRead::default();

        self.check_magic()?;

        let mut title_stream = BufReader::new(self.arc.get_xz_stream("titles")?);
        let mut nuc_stream = BufReader::new(self.arc.get_xz_stream("nucleotides")?);
        let mut qual_stream = self.arc.get_xz_stream("qualities")?;

        loop {
            if !self.read_line(&mut title_stream, &mut r1.title)? {
                break;
            }

            if !self.read_line(&mut title_stream, &mut r2.title)? {
                return Err(CompressionModelError::IncompleteRecord);
            }

            if !self.read_u8(&mut nuc_stream, &mut r1.letters)? {
                return Err(CompressionModelError::IncompleteRecord);
            }

            if !self.read_u8(&mut nuc_stream, &mut r2.letters)? {
                return Err(CompressionModelError::IncompleteRecord);
            }

            self.read_exact(&mut qual_stream, r1.letters.len(), &mut r1.qualities)?;
            self.read_exact(&mut qual_stream, r2.letters.len(), &mut r2.qualities)?;

            writer.write_next(&r1, &r2)?;
        }
        Ok(())
    }
}

impl XZMultiStreamReader<std::io::Stdin> {
    pub fn from_stdin() -> Result<Self> {
        XZMultiStreamReader::new(std::io::stdin())
    }
}

impl XZMultiStreamReader<File> {
    pub fn open<P: AsRef<Path>>(path: &P) -> Result<Self> {
        let file = File::open(path)?;
        XZMultiStreamReader::new(file)
    }
}
