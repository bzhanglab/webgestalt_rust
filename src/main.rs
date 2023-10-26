use bincode::deserialize_from;
use clap::{Args, Parser};
use clap::{Subcommand, ValueEnum};
use owo_colors::{OwoColorize, Stream::Stdout, Style};
use std::io::{BufReader, Write};
use std::{fs::File, time::Instant};
use webgestalt_lib::methods::gsea::GSEAConfig;
use webgestalt_lib::methods::multiomics::{combine_gmts, MultiOmicsMethod, NormalizationMethod};
use webgestalt_lib::methods::ora::ORAConfig;
use webgestalt_lib::readers::utils::Item;
use webgestalt_lib::readers::{read_gmt_file, read_rank_file};
use webgestalt_lib::{MalformedError, WebGestaltError};

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
    /// Benchmark different file formats for gmt.  TODO: Remove later
    Benchmark,
    /// Run provided examples for various types of analyses
    Example(ExampleArgs),
    /// Run GSEA on the provided files
    Gsea(GseaArgs),
    /// Run ORA using the provided files
    Ora(ORAArgs),
    /// Run a test
    Test,
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

#[derive(Args)]
struct GseaArgs {
    /// Path to the GMT file of interest
    gmt: Option<String>,
    /// Path to the rank file of interest
    rnk: Option<String>,
}

#[derive(Args)]
struct ORAArgs {
    /// Path to the GMT file of interest
    gmt: Option<String>,
    /// Path to the file containing the interesting analytes
    interest: Option<String>,
    /// Path the file containing the reference list
    reference: Option<String>,
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
                    "webgestalt_lib/data/test.gmt".to_owned(),
                );
                let start = Instant::now();
                webgestalt_lib::methods::gsea::gsea(
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
                    "webgestalt_lib/data/test.gmt".to_owned(),
                    "webgestalt_lib/data/genelist.txt".to_owned(),
                    "webgestalt_lib/data/reference.txt".to_owned(),
                );
                let gmtcount = gmt.len();
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
                    duration, count, gmtcount
                );
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
            let gene_list = webgestalt_lib::readers::read_rank_file(gsea_args.rnk.clone().unwrap())
                .unwrap_or_else(|_| {
                    panic!("File {} not found", gsea_args.rnk.clone().unwrap());
                });
            let gmt = webgestalt_lib::readers::read_gmt_file(gsea_args.gmt.clone().unwrap())
                .unwrap_or_else(|_| {
                    panic!("File {} not found", gsea_args.gmt.clone().unwrap());
                });
            webgestalt_lib::methods::gsea::gsea(gene_list, gmt, GSEAConfig::default(), None);
            println!("Done with GSEA");
        }
        Some(Commands::Ora(ora_args)) => {
            let style = Style::new().red().bold();
            if ora_args.gmt.is_none() || ora_args.interest.is_none() || ora_args.reference.is_none()
            {
                println!(
                    "{}: DID NOT PROVIDE PATHS FOR GMT, INTEREST, AND REFERENCE FILE.",
                    "ERROR".if_supports_color(Stdout, |text| text.style(style))
                );
                return;
            }
            let start = Instant::now();
            let (gmt, interest, reference) = webgestalt_lib::readers::read_ora_files(
                ora_args.gmt.clone().unwrap(),
                ora_args.interest.clone().unwrap(),
                ora_args.reference.clone().unwrap(),
            );
            println!("Reading Took {:?}", start.elapsed());
            let start = Instant::now();
            let res = webgestalt_lib::methods::ora::get_ora(
                &interest,
                &reference,
                gmt,
                ORAConfig::default(),
            );
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
        Some(Commands::Test) => will_err(1).unwrap_or_else(|x| println!("{}", x)),
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
                let method: MultiOmicsMethod = match ora_args.combination {
                    Some(CombinationMethods::Mean) => MultiOmicsMethod::Mean(norm_method),
                    Some(CombinationMethods::Max) => MultiOmicsMethod::Max(norm_method),
                    None => panic!("No combination method chosen."),
                };
                let mut combined_list =
                    webgestalt_lib::methods::multiomics::combine_lists(lists, method);
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
                panic!("Please select a valid combine type");
            }
        },
        _ => {
            panic!("Please select a valid command. Run --help for options.")
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

fn will_err(x: i32) -> Result<(), WebGestaltError> {
    if x == 0 {
        Ok(())
    } else {
        Err(WebGestaltError::MalformedFile(MalformedError {
            path: String::from("ExamplePath.txt"),
            kind: webgestalt_lib::MalformedErrorType::WrongFormat {
                found: String::from("GMT"),
                expected: String::from("rank"),
            },
        }))
    }
}
