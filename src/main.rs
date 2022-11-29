mod compression_models;
mod seq_files;

use std::ffi::OsString;

use clap::{Parser, Subcommand, ValueEnum};

use crate::compression_models::lzma_multi_stream::*;
use crate::compression_models::lzma_single_file::*;
use crate::compression_models::*;

use crate::seq_files::fastq::*;

#[derive(Debug, Parser)]
#[command(name = "pare")]
#[command(about = "Compression algorithm for fastq files", long_about = None, version, author)]
struct Cli {
    #[command(subcommand)]
    command: Commands,
}

#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord, ValueEnum)]
enum Model {
    LZMA,
    LZMAMulti,
}

#[derive(Debug, Subcommand)]
enum Commands {
    /// compress fastq file(s)
    #[command(arg_required_else_help = true)]
    Compress {
        #[arg(short, long, default_value = "-")]
        output: Option<OsString>,
        #[arg(required = true, num_args(1..3))]
        files: Vec<OsString>,
        #[arg(short, long, action, help = "Don't reverse complement R2")]
        reverse_r2: bool,
        #[arg(short, long, value_enum, help = "which model to use")]
        model: Option<Model>,
    },
    #[command(arg_required_else_help = true)]
    Decompress {
        #[arg(default_value = "-")]
        file: OsString,
        #[arg(default_value = "-", num_args(1..3))]
        outputs: Vec<OsString>,
        #[arg(short, long, action, help = "Don't reverse complement R2")]
        reverse_r2: bool,
        #[arg(short, long, value_enum, help = "which model to use")]
        model: Option<Model>,
    },
    #[command()]
    Cite {},
}

fn compress(
    files: &Vec<OsString>,
    output: Option<OsString>,
    reverse_r2: bool,
    model: Option<Model>,
) -> Result<(), CompressionModelError> {
    let mut sequence_reader: Box<dyn PairedFastQReader>;
    match files.len() {
        1 => {
            //println!("interleaved {:?} {:?} {}", files[0], output, reverse_r2);
            let in_file: Box<dyn FastQFileReaderTrait> = match files[0].to_str() {
                Some("-") => Box::new(FastQFileReader::from_stdin()),
                _ => Box::new(FastQFileReader::open(&files[0])?),
            };

            sequence_reader = Box::new(FastQInterleavedFileReader::new(in_file, reverse_r2));
        }
        2 => {
            //eprintln!("paired files {:?} {:?} {}", files, output, reverse_r2);
            let in_file_r1: Box<dyn FastQFileReaderTrait> = match files[0].to_str() {
                Some("-") => Box::new(FastQFileReader::from_stdin()),
                _ => Box::new(FastQFileReader::open(&files[0])?),
            };

            let in_file_r2: Box<dyn FastQFileReaderTrait> = match files[1].to_str() {
                Some("-") => Box::new(FastQFileReader::from_stdin()),
                _ => Box::new(FastQFileReader::open(&files[1])?),
            };

            sequence_reader = Box::new(FastQPairedFilesReader::new(
                in_file_r1, in_file_r2, reverse_r2,
            ));
        }
        _ => panic!("Too many input files! programming error."),
    }

    let mut writer: Box<dyn EncoderModel>;
    match (output.as_ref().and_then(|o| o.to_str()), model) {
        (Some("-") | None, Some(Model::LZMA) | None) => {
            writer = Box::new(XZSingleFileWriter::to_stdout());
        }
        (_, Some(Model::LZMA) | None) => {
            writer = Box::new(XZSingleFileWriter::create(
                &output.expect("Programming error! output should be Some"),
            )?);
        }
        (Some("-") | None, Some(Model::LZMAMulti)) => {
            writer = Box::new(XZMultiStreamWriter::to_stdout());
        }
        (_, Some(Model::LZMAMulti)) => {
            writer = Box::new(XZMultiStreamWriter::create(
                &output.expect("Programming error! output should be Some"),
            )?);
        }
    }

    writer.compress(&mut sequence_reader)?;
    Ok(())
}

fn decompress(
    file: OsString,
    outputs: Vec<OsString>,
    reverse_r2: bool,
    model: Option<Model>,
) -> Result<(), CompressionModelError> {
    let mut sequence_writer: Box<dyn PairedFastQWriter>;
    match outputs.len() {
        1 => {
            eprintln!("interleaved {:?} {:?} {:?}", file, outputs, reverse_r2);
            let out_file: Box<dyn FastQFileWriterTrait> = match outputs[0].to_str() {
                Some("-") => Box::new(FastQFileWriter::to_stdout()),
                _ => Box::new(FastQFileWriter::create(&outputs[0])?),
            };

            sequence_writer = Box::new(FastQInterleavedFileWriter::new(out_file, reverse_r2));
        }
        2 => {
            eprintln!("paired files {:?} {:?} {:?}", file, outputs, reverse_r2);
            let out_file_r1: Box<dyn FastQFileWriterTrait> = match outputs[0].to_str() {
                Some("-") => Box::new(FastQFileWriter::to_stdout()),
                _ => Box::new(FastQFileWriter::create(&outputs[0])?),
            };

            let out_file_r2: Box<dyn FastQFileWriterTrait> = match outputs[1].to_str() {
                Some("-") => Box::new(FastQFileWriter::to_stdout()),
                _ => Box::new(FastQFileWriter::create(&outputs[1])?),
            };

            sequence_writer = Box::new(FastQPairedFilesWriter::new(
                out_file_r1,
                out_file_r2,
                reverse_r2,
            ));
        }
        _ => panic!("Too many output files! programming error."),
    }

    let mut writer: Box<dyn DecoderModel>;
    match (file.to_str(), model) {
        (Some("-") | None, Some(Model::LZMA) | None) => {
            writer = Box::new(XZSingleFileReader::from_stdin()?);
        }
        (_, Some(Model::LZMA) | None) => {
            writer = Box::new(XZSingleFileReader::open(&file)?);
        }
        (Some("-") | None, Some(Model::LZMAMulti)) => {
            writer = Box::new(XZMultiStreamReader::from_stdin()?);
        }
        (_, Some(Model::LZMAMulti)) => {
            writer = Box::new(XZMultiStreamReader::open(&file)?);
        }
    }
    writer.decompress(&mut sequence_writer)?;

    Ok(())
}

fn main() -> Result<(), CompressionModelError> {
    let args = Cli::parse();

    match args.command {
        Commands::Compress {
            files,
            output,
            reverse_r2,
            model,
        } => compress(&files, output, reverse_r2, model)?,
        Commands::Decompress {
            file,
            outputs,
            reverse_r2,
            model,
        } => decompress(file, outputs, reverse_r2, model)?,
        Commands::Cite {} => {
            println!("print out a citation here");
        }
    }
    Ok(())
}
