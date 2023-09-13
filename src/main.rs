use std::io::{BufReader, Write};
use std::{fs::File, time::Instant};

use bincode::deserialize_from;

use clap::Parser;
use clap::Subcommand;

/// WebGestalt CLI.
/// ORA and GSEA enrichment tool.
/// Created by Bing Zhang Lab.
#[derive(Parser)]
#[command(author, version, about, long_about = None)]
struct Args {
    // /// Gene List Path
    // #[arg(short, long)]
    // gene_list: String,
    // /// Test argumnet
    // #[arg(short, long)]
    // test: u32,
    #[command(subcommand)]
    command: Option<Commands>,
}

#[derive(Subcommand)]
enum Commands {
    /// Benchmark different file formats for gmt.  TODO: Remove later
    Benchmark,
}

fn main() {
    let args = Args::parse();
    match &args.command {
        Some(Commands::Benchmark) => {
            benchmark();
        }
        None => {
            let gene_list =
                webgestalt_lib::readers::read_rank_file("webgestalt_lib/data/test.rnk".to_owned());
            let gmt =
                webgestalt_lib::readers::read_gmt_file("webgestalt_lib/data/ktest.gmt".to_owned());
            let start = Instant::now();
            webgestalt_lib::methods::gsea::gsea(gene_list.unwrap(), gmt.unwrap());
            let duration = start.elapsed();
            println!("New Hash\nTime took: {:?}", duration);
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
    let improvment: f64 = 100.0 * (gmt_avg - bin_avg) / gmt_avg;
    println!(
        " GMT time: {}\tGMT.WGA time: {}\n Improvment: {:.1}%",
        gmt_avg, bin_avg, improvment
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
