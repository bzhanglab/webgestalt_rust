pub mod utils;
use crate::methods::gsea::RankListItem;
use ahash::AHashSet;
use std::{
    fs::File,
    io::{prelude::*, BufReader},
};
use utils::Item;

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
        let id = result.get(0).unwrap().to_owned();
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
        let phenotype = result.get(0).unwrap().to_owned();
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
    let mut h = AHashSet::default();
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
        .delimiter(b'\t') // TODO: Add option to use different delimeter
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
        let id = result.get(0).unwrap().to_owned();
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
