use clap::{Args, Parser};
use clap::{Subcommand, ValueEnum};
use owo_colors::{OwoColorize, Stream::Stdout, Style};
use std::io::Write;
use std::{fs::File, time::Instant};
use webgestalt_lib::methods::gsea::GSEAConfig;
use webgestalt_lib::methods::multilist::{combine_gmts, MultiListMethod, NormalizationMethod};
use webgestalt_lib::methods::nta::NTAConfig;
use webgestalt_lib::methods::ora::ORAConfig;
use webgestalt_lib::readers::utils::Item;
use webgestalt_lib::readers::{read_gmt_file, read_rank_file};

/// WebGestalt CLI.
/// ORA and GSEA enrichment tool.
/// Created by John Elizarraras from the Bing Zhang Lab.
#[derive(Parser)]
#[command(author, version, about, long_about = None)]
struct CliArgs {
    #[command(subcommand)]
    command: Option<Commands>,
}

#[derive(Subcommand)]
enum Commands {
    /// Run provided examples for various types of analyses
    Example(ExampleArgs),
    /// Run GSEA on the provided files
    Gsea(GseaArgs),
    /// Run ORA using the provided files
    Ora(ORAArgs),
    /// Run NTA on the provided files
    Nta(NtaArgs),
    /// Combine multiple files into a single file
    Combine(CombineArgs),
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

#[derive(Parser)]
struct NtaArgs {
    /// Path to the file containing the network (tab separated file with two columns: source and target)
    #[arg(short, long)]
    network: String,
    /// Path to the file containing the seeds (one per line)
    #[arg(short, long)]
    seeds: String,
    /// Output path for the results
    #[arg(short, long)]
    output: String,
    /// Probability of random walk resetting
    #[arg(short, long, default_value = "0.5")]
    reset_probability: f64,
    /// Convergence tolerance
    #[arg(short, long, default_value = "0.000001")]
    tolerance: f64,
    /// Number of nodes to prioritize or expand to
    #[arg(short = 'z', long = "size", default_value = "50")]
    neighborhood_size: usize,
    /// Method to use for NTA
    /// Options: prioritize, expand
    #[arg(short, long, default_value = "prioritize")]
    method: NTAMethodClap,
}

#[derive(ValueEnum, Clone)]
enum NTAMethodClap {
    Prioritize,
    Expand,
}

#[derive(Args)]
struct GseaArgs {
    /// Path to the GMT file of interest
    #[arg(short, long)]
    gmt: String,
    /// Path to the rank file of interest
    #[arg(short, long)]
    rnk: String,
    /// Output path for the results
    #[arg(short, long, default_value = "out.json")]
    output: String,
}
#[derive(Parser)]
struct ORAArgs {
    /// Path to the GMT file of interest
    #[arg(short, long)]
    gmt: String,
    /// Path to the file containing the interesting analytes
    #[arg(short, long)]
    interest: String,
    /// Output path for the results
    #[arg(short, long, default_value = "out.json")]
    output: String,
    /// Path the file containing the reference list
    #[arg(short, long)]
    reference: String,
}

#[derive(Args)]
struct CombineArgs {
    #[command(subcommand)]
    combine_type: Option<CombineType>,
}

#[derive(Subcommand)]
enum CombineType {
    Gmt(CombineGmtArgs),
    List(CombineListArgs),
}

#[derive(Args)]
struct CombineGmtArgs {
    out: Option<String>,
    /// Paths to the files to combine
    files: Vec<String>,
}

#[derive(ValueEnum, Clone)]
enum NormMethods {
    MedianRank,
    MedianValue,
    MeanValue,
    None,
}

#[derive(ValueEnum, Clone)]
enum CombinationMethods {
    Max,
    Mean,
}

#[derive(Args)]
struct CombineListArgs {
    combination: Option<CombinationMethods>,
    normalization: Option<NormMethods>,
    out: Option<String>,
    files: Vec<String>,
}

fn prompt_yes_no(question: &str) -> bool {
    loop {
        print!("{} (y/n): ", question);
        std::io::stdout().flush().expect("Could not flush stdout!"); // Ensure the prompt is displayed

        let mut input = String::new();
        std::io::stdin()
            .read_line(&mut input)
            .expect("Could not read line");
        print!("\x1B[2J\x1B[1;1H");
        std::io::stdout().flush().expect("Could not flush stdout!");
        match input.trim().to_lowercase().as_str() {
            "y" => return true,
            "n" => return false,
            _ => println!("Invalid input. Please enter 'y' or 'n'."),
        }
    }
}

fn check_and_overwrite(file_path: &str) {
    // Check if the file exists
    if std::path::Path::new(file_path).exists() {
        // Check if the user wants to overwrite the file
        return if !prompt_yes_no(&format!(
            "File at {} already exists. Do you want to overwrite it?",
            file_path
        )) {
            println!("Stopping analysis.");
            std::process::exit(1);
        };
    }
}

fn main() {
    println!("WebGestalt CLI v{}", env!("CARGO_PKG_VERSION"));
    let args = CliArgs::parse();
    match &args.command {
        Some(Commands::Example(ex)) => match &ex.commands {
            Some(ExampleOptions::Gsea) => {
                let gene_list = webgestalt_lib::readers::read_rank_file(
                    "webgestalt_lib/data/test.rnk".to_owned(),
                );
                let gmt = webgestalt_lib::readers::read_gmt_file(
                    "webgestalt_lib/data/ktest.gmt".to_owned(),
                );
                let start = Instant::now();
                let _res = webgestalt_lib::methods::gsea::gsea(
                    gene_list.unwrap(),
                    gmt.unwrap(),
                    GSEAConfig::default(),
                    None,
                );
                let duration = start.elapsed();
                println!("GSEA\nTime took: {:?}", duration);
            }
            Some(ExampleOptions::Ora) => {
                let (gmt, gene_list, reference) = webgestalt_lib::readers::read_ora_files(
                    "webgestalt_lib/data/ktest.gmt".to_owned(),
                    "webgestalt_lib/data/genelist.txt".to_owned(),
                    "webgestalt_lib/data/reference.txt".to_owned(),
                );
                let gmt_count = gmt.len();
                let start = Instant::now();
                let x: Vec<webgestalt_lib::methods::ora::ORAResult> =
                    webgestalt_lib::methods::ora::get_ora(
                        &gene_list,
                        &reference,
                        gmt,
                        ORAConfig::default(),
                    );
                let mut count = 0;
                for i in x {
                    if i.p < 0.05 && i.fdr < 0.05 {
                        println!("{}: {}, {}, {}", i.set, i.p, i.fdr, i.overlap);
                        count += 1;
                    }
                }
                let duration = start.elapsed();
                println!(
                    "ORA\nTime took: {:?}\nFound {} significant pathways out of {} pathways",
                    duration, count, gmt_count
                );
            }
            _ => {
                println!("Please select a valid example: ora or gsea.");
            }
        },
        Some(Commands::Gsea(gsea_args)) => {
            check_and_overwrite(&gsea_args.output);
            let gene_list = webgestalt_lib::readers::read_rank_file(gsea_args.rnk.clone())
                .unwrap_or_else(|_| {
                    panic!("File {} not found", gsea_args.rnk.clone());
                });
            let gmt = webgestalt_lib::readers::read_gmt_file(gsea_args.gmt.clone()).unwrap_or_else(
                |_| {
                    panic!("File {} not found", gsea_args.gmt.clone());
                },
            );
            let res =
                webgestalt_lib::methods::gsea::gsea(gene_list, gmt, GSEAConfig::default(), None);
            let mut count = 0;
            for i in res {
                if i.p < 0.05 && i.fdr < 0.05 {
                    println!("{}: {}, {}", i.set, i.p, i.fdr);
                    count += 1;
                }
            }
            println!("Done with GSEA: {}", count);
        }
        Some(Commands::Ora(ora_args)) => {
            check_and_overwrite(&ora_args.output);
            let start = Instant::now();
            let (gmt, interest, reference) = webgestalt_lib::readers::read_ora_files(
                ora_args.gmt.clone(),
                ora_args.interest.clone(),
                ora_args.reference.clone(),
            );
            println!("Reading Took {:?}", start.elapsed());
            let start = Instant::now();
            let res = webgestalt_lib::methods::ora::get_ora(
                &interest,
                &reference,
                gmt,
                ORAConfig::default(),
            );
            let output_file =
                File::create(&ora_args.output).expect("Could not create output file!");
            serde_json::to_writer(output_file, &res).expect("Could not create JSON file!");
            println!("Analysis Took {:?}", start.elapsed());
            let mut count = 0;
            for row in res.iter() {
                if row.p < 0.05 && row.fdr < 0.05 {
                    count += 1;
                }
            }
            println!(
                "Found {} significant pathways out of {} pathways",
                count,
                res.len()
            );
        }
        Some(Commands::Nta(nta_args)) => {
            check_and_overwrite(&nta_args.output);
            let network = webgestalt_lib::readers::read_edge_list(nta_args.network.clone());
            let start = Instant::now();
            let nta_method = match nta_args.method {
                NTAMethodClap::Prioritize => {
                    webgestalt_lib::methods::nta::NTAMethod::Prioritize(nta_args.neighborhood_size)
                }
                NTAMethodClap::Expand => {
                    webgestalt_lib::methods::nta::NTAMethod::Expand(nta_args.neighborhood_size)
                }
            };
            let config: NTAConfig = NTAConfig {
                edge_list: network,
                seeds: webgestalt_lib::readers::read_seeds(nta_args.seeds.clone()),
                reset_probability: nta_args.reset_probability,
                tolerance: nta_args.tolerance,
                method: Some(nta_method),
            };
            let res = webgestalt_lib::methods::nta::get_nta(config);
            println!("Analysis Took {:?}", start.elapsed());
            webgestalt_lib::writers::save_nta(nta_args.output.clone(), res).unwrap();
        }
        Some(Commands::Combine(args)) => match &args.combine_type {
            Some(CombineType::Gmt(gmt_args)) => {
                let style = Style::new().blue().bold();
                println!(
                    "{}: READING GMTS",
                    "INFO".if_supports_color(Stdout, |text| text.style(style))
                );
                let mut gmts: Vec<Vec<Item>> = Vec::new();
                let mut tot_length: usize = 0;
                for path in gmt_args.files.clone() {
                    let gmt = read_gmt_file(path).unwrap();
                    tot_length += gmt.len();
                    gmts.push(gmt);
                }
                let combined_gmt = combine_gmts(&gmts);
                println!(
                    "Found {} overlapping sets out of {}",
                    tot_length - combined_gmt.len(),
                    combined_gmt.len()
                );
                println!(
                    "{}: CREATING COMBINED GMT AT {}",
                    "INFO".if_supports_color(Stdout, |text| text.style(style)),
                    gmt_args.out.clone().unwrap()
                );
                let mut file = File::create(gmt_args.out.clone().unwrap()).unwrap();
                for row in combined_gmt {
                    writeln!(file, "{}\t{}\t{}", row.id, row.url, row.parts.join("\t")).unwrap();
                }
            }
            Some(CombineType::List(ora_args)) => {
                let style = Style::new().blue().bold();
                println!(
                    "{}: READING LISTS",
                    "INFO".if_supports_color(Stdout, |text| text.style(style))
                );
                let mut lists = Vec::new();
                for file in ora_args.files.iter() {
                    lists.push(read_rank_file(file.clone()).unwrap());
                }
                let norm_method: NormalizationMethod = match ora_args.normalization {
                    Some(NormMethods::None) => NormalizationMethod::None,
                    Some(NormMethods::MeanValue) => NormalizationMethod::MeanValue,
                    Some(NormMethods::MedianRank) => NormalizationMethod::MedianRank,
                    Some(NormMethods::MedianValue) => NormalizationMethod::MedianValue,
                    None => panic!("No normalization method chosen."),
                };
                let method: MultiListMethod = match ora_args.combination {
                    Some(CombinationMethods::Mean) => MultiListMethod::Mean(norm_method),
                    Some(CombinationMethods::Max) => MultiListMethod::Max(norm_method),
                    None => panic!("No combination method chosen."),
                };
                let mut combined_list =
                    webgestalt_lib::methods::multilist::combine_lists(lists, method);
                combined_list.sort_by(|a, b| b.rank.partial_cmp(&a.rank).unwrap());
                let mut file = File::create(ora_args.out.clone().unwrap()).unwrap();
                println!(
                    "{}: CREATING COMBINED LIST AT {}",
                    "INFO".if_supports_color(Stdout, |text| text.style(style)),
                    ora_args.out.clone().unwrap()
                );
                for row in combined_list {
                    writeln!(file, "{}\t{}", row.analyte, row.rank).unwrap();
                }
            }
            _ => {
                println!("Please select a valid combine type");
            }
        },
        _ => {
            println!("Please select a valid command. Run --help for options.")
        }
    }
}
