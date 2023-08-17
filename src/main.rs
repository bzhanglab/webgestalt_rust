use std::time::Instant;

use clap::Parser;

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
    let start = Instant::now();
    webgestalt_lib::methods::gsea::get_gsea(gene_list.unwrap(), gmt.unwrap());
    let duration = start.elapsed();
    println!("New Hash\nTime took: {:?}", duration);
}
