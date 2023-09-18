use std::io::{BufReader, Write};
use std::{fs::File, time::Instant};

use bincode::deserialize_from;

use clap::Subcommand;
use clap::{Args, Parser};
use owo_colors::{OwoColorize, Stream::Stdout, Style};

/// WebGestalt CLI.
/// ORA and GSEA enrichment tool.
/// Created by Bing Zhang Lab.
#[derive(Parser)]
#[command(author, version, about, long_about = None)]
struct CliArgs {
    #[command(subcommand)]
    command: Option<Commands>,
}

#[derive(Subcommand)]
enum Commands {
    /// Benchmark different file formats for gmt.  TODO: Remove later
    Benchmark,
    /// Run provided examples for various types of analyses
    Example(ExampleArgs),
    /// Run GSEA on the provided files
    Gsea(GseaArgs),
}

#[derive(Debug, Args)]
#[command(args_conflicts_with_subcommands = true)]
struct ExampleArgs {
    #[command(subcommand)]
    commands: Option<ExampleOptions>,
}

#[derive(Subcommand, Debug)]
enum ExampleOptions {
    /// Run a GSEA example
    Gsea,
    /// Run an ORA example
    Ora,
}

#[derive(Args)]
struct GseaArgs {
    /// Path to the GMT file of interest
    gmt: Option<String>,
    /// Path to the rank file of interest
    rnk: Option<String>,
}

fn main() {
    let args = CliArgs::parse();
    match &args.command {
        Some(Commands::Benchmark) => {
            benchmark();
        }
        Some(Commands::Example(ex)) => match &ex.commands {
            Some(ExampleOptions::Gsea) => {
                let gene_list = webgestalt_lib::readers::read_rank_file(
                    "webgestalt_lib/data/test.rnk".to_owned(),
                );
                let gmt = webgestalt_lib::readers::read_gmt_file(
                    "webgestalt_lib/data/ktest.gmt".to_owned(),
                );
                let start = Instant::now();
                webgestalt_lib::methods::gsea::gsea(gene_list.unwrap(), gmt.unwrap());
                let duration = start.elapsed();
                println!("New Hash\nTime took: {:?}", duration);
            }
            Some(ExampleOptions::Ora) => {
                println!("ORA example not implemented.")
            }
            _ => todo!("HELLO"),
        },
        Some(Commands::Gsea(gsea_args)) => {
            let style = Style::new().red().bold();
            if gsea_args.gmt.is_none() || gsea_args.rnk.is_none() {
                println!(
                    "{}: DID NOT PROVIDE PATHS FOR GMT AND RANK FILE.",
                    "ERROR".if_supports_color(Stdout, |text| text.style(style))
                );
                return;
            }
            let gene_list =
                webgestalt_lib::readers::read_rank_file(gsea_args.rnk.clone().unwrap()).unwrap();
            let gmt =
                webgestalt_lib::readers::read_gmt_file(gsea_args.gmt.clone().unwrap()).unwrap();
            webgestalt_lib::methods::gsea::gsea(gene_list, gmt);
            println!("Done with GSEA");
        }
        _ => {
            println!("Please select a command. Run --help for options.")
        }
    }
}

fn benchmark() {
    let mut bin_durations: Vec<f64> = Vec::new();
    for _i in 0..1000 {
        let start = Instant::now();
        let mut r = BufReader::new(File::open("test.gmt.wga").unwrap());
        let _x: Vec<webgestalt_lib::readers::utils::Item> = deserialize_from(&mut r).unwrap();
        let duration = start.elapsed();
        bin_durations.push(duration.as_secs_f64())
    }
    let mut gmt_durations: Vec<f64> = Vec::new();
    for _i in 0..1000 {
        let start = Instant::now();
        let _x = webgestalt_lib::readers::read_gmt_file("webgestalt_lib/data/ktest.gmt".to_owned())
            .unwrap();
        let duration = start.elapsed();
        gmt_durations.push(duration.as_secs_f64())
    }
    let gmt_avg: f64 = gmt_durations.iter().sum::<f64>() / gmt_durations.len() as f64;
    let bin_avg: f64 = bin_durations.iter().sum::<f64>() / bin_durations.len() as f64;
    let improvement: f64 = 100.0 * (gmt_avg - bin_avg) / gmt_avg;
    println!(
        " GMT time: {}\tGMT.WGA time: {}\n Improvement: {:.1}%",
        gmt_avg, bin_avg, improvement
    );
    let mut whole_file: Vec<String> = Vec::new();
    whole_file.push("type\ttime".to_string());
    for line in bin_durations {
        whole_file.push(format!("bin\t{:?}", line));
    }
    for line in gmt_durations {
        whole_file.push(format!("gmt\t{:?}", line));
    }
    let mut ftsv = File::create("format_benchmarks.tsv").unwrap();
    writeln!(ftsv, "{}", whole_file.join("\n")).unwrap();
}
