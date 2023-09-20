pub mod utils;
use crate::methods::gsea::RankListItem;
use rustc_hash::FxHashSet;
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

pub fn read_single_list(path: String) -> FxHashSet<String> {
    let file = File::open(path).expect("no such file");
    let buf = BufReader::new(file);
    let mut h = rustc_hash::FxHashSet::default();
    let v: Vec<String> = buf
        .lines()
        .map(|l| l.expect("Could not parse line"))
        .collect();
    for i in v {
        h.insert(i);
    }
    h
}
