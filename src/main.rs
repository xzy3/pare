mod compression_models;
mod seq_files;

use std::ffi::OsString;

use clap::{Parser, Subcommand};

use crate::compression_models::lzma_single_file::*;
use crate::seq_files::fastq::*;

#[derive(Debug, Parser)]
#[command(name = "pare")]
#[command(about = "Compression algorithm for fastq files", long_about = None, version, author)]
struct Cli {
    #[command(subcommand)]
    command: Commands,
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
    },
    #[command(arg_required_else_help = true)]
    Decompress {
        #[arg(short, long, default_value = "-", num_args(1..3))]
        outputs: Vec<OsString>,
        file: OsString,
    },
    #[command()]
    Cite {},
}

fn compress(
    files: &Vec<OsString>,
    output: Option<OsString>,
    reverse_r2: bool,
) -> Result<(), XZSingleFileError> {
    let mut sequence_reader: Box<dyn PairedFastQReader>;
    match files.len() {
        1 => {
            //println!("interleaved {:?} {:?} {}", files[0], output, reverse_r2);
            let in_file: Box<dyn FastQFileReader> = match files[0].to_str() {
                Some("-") => Box::new(FastQFile::from_stdin()),
                _ => Box::new(FastQFile::open(&files[0])?),
            };

            sequence_reader = Box::new(FastQInterleavedFile::new(in_file, reverse_r2));
        }
        2 => {
            //println!("paired files {:?} {:?} {}", files, output, reverse_r2);
            let in_file_r1: Box<dyn FastQFileReader> = match files[0].to_str() {
                Some("-") => Box::new(FastQFile::from_stdin()),
                _ => Box::new(FastQFile::open(&files[0])?),
            };

            let in_file_r2: Box<dyn FastQFileReader> = match files[1].to_str() {
                Some("-") => Box::new(FastQFile::from_stdin()),
                _ => Box::new(FastQFile::open(&files[1])?),
            };

            sequence_reader = Box::new(FastQPairedFiles::new(in_file_r1, in_file_r2, reverse_r2));
        }
        _ => panic!("Too many input files! programming error."),
    }

    match output.as_ref().and_then(|o| o.to_str()) {
        Some("-") | None => {
            let mut writer = XZSingleFileWriter::to_stdout();
            writer.compress(&mut sequence_reader)?;
        }
        _ => {
            let mut writer = XZSingleFileWriter::create(
                &output.expect("Programming error! output should be Some"),
            )?;
            writer.compress(&mut sequence_reader)?;
        }
    }

    Ok(())
}

fn main() -> Result<(), XZSingleFileError> {
    let args = Cli::parse();

    match args.command {
        Commands::Compress {
            files,
            output,
            reverse_r2,
        } => compress(&files, output, reverse_r2)?,
        Commands::Decompress { file, outputs } => match outputs.len() {
            1 => {
                println!("interleaved {:?} {:?}", file, outputs);
            }
            2 => {
                println!("paired files {:?} {:?}", file, outputs[0]);
            }
            _ => panic!("Too many output files! programming error."),
        },
        Commands::Cite {} => {
            println!("print out a citation here");
        }
    }
    Ok(())
}
