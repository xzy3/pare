mod seq_files;

use std::ffi::OsString;

use clap::{Parser, Subcommand};

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

fn main() -> Result<(), FastQFileError> {
    let args = Cli::parse();

    match args.command {
        Commands::Compress { files, output } => match files.len() {
            1 => {
                println!("interleaved {:?} {:?}", files[0], output);
            }
            2 => {
                println!("paired files {:?} {:?}", files, output);
            }
            _ => panic!("Too many input files! programming error."),
        },
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
