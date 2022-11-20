mod seq_files;

use std::path::PathBuf;

use clap::{Parser, Subcommand};

#[derive(Debug, Parser)]
#[command(name = "pare")]
#[command(about = "Compression algorithm for fastq files", long_about = None)]
struct Cli {
    #[command(subcommand)]
    command: Commands,
}

#[derive(Debug, Subcommand)]
enum Commands {
    /// compress fastq file(s)
    #[command(arg_required_else_help = true)]
    Compress {
        #[arg(short, long, default_missing_value = "-")]
        output: Option<PathBuf>,
        #[arg(required = true)]
        files: Vec<PathBuf>,
    },
    #[command(arg_required_else_help = true)]
    Decompress {
        #[arg(short, long, default_missing_value = "-")]
        output: Option<PathBuf>,
        file: PathBuf,
    },
}

fn main() {
    let args = Cli::parse();

    match args.command {
        Commands::Compress { files, output } => {
            println!("files: {:?} {:?}", files, output);
        }
        Commands::Decompress { file, output } => {
            println!("decompress files {:?} {:?}", file, output);
        }
    }
}
