pub mod utils;
use crate::methods::gsea::RankListItem;
use ahash::AHashSet;
use std::{
    fs::File,
    io::{prelude::*, BufReader},
};
use utils::Item;

/// Read GMT file from specified path. For format description, see [broadinstitute.org](https://software.broadinstitute.org/cancer/software/gsea/wiki/index.php/Data_formats#GMT:_Gene_Matrix_Transposed_file_format_.28.2A.gmt.29)
///
/// # Parameters
///
/// - `path` - A [`String`] of the path of the GMT to read.
///
/// # Panics
///
/// Panics if there is not file at `path`.
///
/// # Returns
///
/// If result is `Ok`, returns a [`Vec<Item>`] containing the elements of the GMT
pub fn read_gmt_file(path: String) -> Result<Vec<Item>, Box<std::io::Error>> {
    let file = File::open(path)?;
    let mut rdr = csv::ReaderBuilder::new()
        .delimiter(b'\t')
        .flexible(true)
        .has_headers(false)
        .from_reader(file);
    let mut items: Vec<utils::Item> = Vec::new();
    for r in rdr.records() {
        let result = r
            .unwrap()
            .iter()
            .map(|x| x.to_string())
            .collect::<Vec<String>>();
        let id = result.first().unwrap().to_owned();
        let url = result.get(1).unwrap().to_owned();
        let parts = result[2..].to_vec();
        let item = Item { id, url, parts };
        items.push(item);
    }
    Ok(items)
}

pub fn read_rank_file(path: String) -> Result<Vec<RankListItem>, Box<std::io::Error>> {
    let file = File::open(path)?;
    let mut rdr = csv::ReaderBuilder::new()
        .delimiter(b'\t')
        .flexible(true)
        .has_headers(false)
        .from_reader(file);
    let mut items: Vec<RankListItem> = Vec::new();
    for r in rdr.records() {
        let result = r
            .unwrap()
            .iter()
            .map(|x| x.to_string())
            .collect::<Vec<String>>();
        let phenotype = result.first().unwrap().to_owned();
        let rank = result.get(1).unwrap().to_owned().parse::<f64>().unwrap();
        let item = RankListItem {
            analyte: phenotype,
            rank,
        };
        items.push(item);
    }
    Ok(items)
}

pub fn read_single_list(path: String) -> AHashSet<String> {
    let file = File::open(path).expect("no such file");
    let buf = BufReader::new(file);
    let mut h: AHashSet<String> = AHashSet::default();
    let v: Vec<String> = buf
        .lines()
        .map(|l| l.expect("Could not parse line"))
        .collect();
    for i in v {
        h.insert(i);
    }
    h
}

pub fn read_ora_files(
    gmt_path: String,
    interest_path: String,
    ref_path: String,
) -> (Vec<Item>, AHashSet<String>, AHashSet<String>) {
    let file = File::open(gmt_path).unwrap();
    let mut rdr = csv::ReaderBuilder::new()
        .delimiter(b'\t') // TODO: Add option to use different delimiter
        .flexible(true)
        .has_headers(false)
        .from_reader(file);
    let mut items: Vec<utils::Item> = Vec::new();
    let mut annotated_genes: AHashSet<String> = AHashSet::default();
    for r in rdr.records() {
        let result = r
            .unwrap()
            .iter()
            .map(|x| x.to_string())
            .collect::<Vec<String>>();
        let id = result.first().unwrap().to_owned();
        let url = result.get(1).unwrap().to_owned();
        let parts = result[2..].to_vec();
        for analyte in parts.clone().into_iter() {
            annotated_genes.insert(analyte);
        }
        let item = Item { id, url, parts };
        items.push(item);
    }
    let reference_list = read_intersection_list(ref_path, &annotated_genes);
    let analyte_list = read_intersection_list(interest_path, &reference_list);
    (items, analyte_list, reference_list)
}

pub fn read_intersection_list(path: String, ref_list: &AHashSet<String>) -> AHashSet<String> {
    let file = File::open(path).expect("no such file");
    let buf = BufReader::new(file);
    let mut h = AHashSet::default();
    let v: Vec<String> = buf
        .lines()
        .map(|l| l.expect("Could not parse line"))
        .collect();
    for i in v {
        if ref_list.contains(&i) {
            h.insert(i);
        }
    }
    h
}

/// Read edge list from specified path. Separated by whitespace with no support for weights
/// 
/// # Parameters
/// path - A [`String`] of the path of the edge list to read.
/// 
/// # Returns
/// A [`Vec<Vec<String>>`] containing the edge list
pub fn read_edge_list(path: String) -> Vec<Vec<String>> {
    let file = File::open(path).expect("no such file");
    let buf = BufReader::new(file);
    let mut v: Vec<Vec<String>> = Vec::new();
    for line in buf.lines() {
        let l = line.expect("Could not parse line");
        let parts: Vec<String> = l.split_whitespace().map(|s| s.to_string()).collect();
        v.push(parts);
    }
    v
}

pub fn read_seeds(path: String) -> Vec<String> {
    let file = File::open(path).expect("no such file");
    let buf = BufReader::new(file);
    let v: Vec<String> = buf
        .lines()
        .map(|l| l.expect("Could not parse line"))
        .filter(|line| !line.is_empty())
        .collect();
    v
}
