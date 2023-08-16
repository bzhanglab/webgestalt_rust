use clap::Parser;

/// WebGestalt CLI
#[derive(Parser, Debug)]
#[command(author, version, about, long_about = None)]
struct Args {
    /// Gene List Path
    #[arg(short, long)]
    gene_list: String,

    /// Test argumnet
    #[arg(short, long)]
    test: u32,
}

fn main() {
    let args = Args::parse();
    println!("Hello, {:?}", webgestalt_lib::add(1, 2));
}
