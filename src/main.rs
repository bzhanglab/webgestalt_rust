use std::{time::Instant, fs::File};

use bincode::serialize_into;
use clap::Parser;
use std::io::BufWriter;

/// WebGestalt CLI.
/// ORA and GSEA enrichment tool.
/// Created by Bing Zhang Lab.
#[derive(Parser, Debug)]
#[command(author, version, about, long_about = None)]
struct Args {
    // /// Gene List Path
    // #[arg(short, long)]
    // gene_list: String,

    // /// Test argumnet
    // #[arg(short, long)]
    // test: u32,
}

fn main() {
    let args = Args::parse();
    let gene_list = webgestalt_lib::readers::read_rank_file("test.rnk".to_owned());
    let gmt = webgestalt_lib::readers::read_gmt_file("test.gmt".to_owned());
    // let mut f = BufWriter::new(File::create("test.gmt.wga").unwrap());
    // serialize_into(&mut f, &gmt.unwrap()).unwrap();
    let start = Instant::now();
    webgestalt_lib::methods::gsea::gsea(gene_list.unwrap(), gmt.unwrap());
    let duration = start.elapsed();
    println!("New Hash\nTime took: {:?}", duration);
}
