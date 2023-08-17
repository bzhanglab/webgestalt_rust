use std::{sync::{Arc, Mutex}, time::Instant};

use crate::readers::utils::Item;
use rayon::prelude::*;
use rustc_hash::FxHashSet;

pub struct RankListItem {
    pub phenotype: String,
    pub rank: f64,
}

impl RankListItem {
    pub fn to_vecs(input: Vec<RankListItem>) -> (Vec<String>, Vec<f64>) {
        let mut phenotypes: Vec<String> = Vec::new();
        let mut ranks: Vec<f64> = Vec::new();
        for item in input {
            phenotypes.push(item.phenotype);
            ranks.push(item.rank)
        }
        (phenotypes, ranks)
    }
} 

fn gene_set_p(genes: &Vec<String>, ranks: &[f64], item: &Item, p: f64, permutations: usize) {
    let gene_set = FxHashSet::from_iter(item.parts.iter());
    let mut n_r: f64 = 0.0;
    let mut n_miss: f64 = 0.0;

    let inverse_size_dif: f64 = 1.0 / ((genes.len() - item.parts.len()) as f64); // Inverse now,
    let gene_size = genes.len();
    for j in 0..gene_size {
        if gene_set.contains(&genes[j]) {
            n_r += ranks[j].powf(p);
        } else {
            n_miss += 1.0;
        }
    }
    let initial_max_score: f64 = 1.0 - n_miss * inverse_size_dif;
    let inverse_nr = 1.0 / n_r;
    let original_order = 0..(gene_size - 1);
    let has_gene: Vec<bool> = genes.par_iter().map(|x| gene_set.contains(x)).collect();
    let new_ranks: Vec<f64> = ranks
        .par_iter()
        .enumerate()
        .map(|(i, x)| x.powf(p))
        .collect();
    let real_es = enrichment_score(
        &gene_set,
        &has_gene,
        &new_ranks,
        original_order.collect(),
        p,
        inverse_size_dif,
        inverse_nr,
        initial_max_score,
    );
    let perm_es = Arc::new(Mutex::new(Vec::new()));
    (0..permutations).into_par_iter().for_each(|_i| {
        let new_order = fastrand::choose_multiple(0..(gene_size), gene_size);
        perm_es.lock().unwrap().push(enrichment_score(
            &gene_set,
            &has_gene,
            &new_ranks,
            new_order,
            p,
            inverse_size_dif,
            inverse_nr,
            initial_max_score,
        ));
    });
}

fn enrichment_score(
    genes: &Vec<bool>,
    ranks: &[f64],
    order: Vec<usize>,
    inverse_size_dif: f64,
    inverse_nr: f64,
    initial_max: f64,
) -> f64 {
    let mut max_score = initial_max;
    let mut sum_hits = 0.0;
    let mut sum_miss = 0.0;
    for i in 0..(genes.len() - 1) {
        if genes[order[i]] {
            sum_hits += ranks[i];
        }
        else {
            sum_miss += 1.0;
        }
        let es = sum_hits * inverse_nr - sum_miss * inverse_size_dif;
        if es.abs() > max_score.abs() {
            max_score = es;
        }
    }
    max_score
}

pub fn get_gsea(mut gene_list: Vec<RankListItem>, gmt: Vec<Item>) {
    println!("Starting GSEA Calculation.");
    gene_list.sort_by(|a, b| b.rank.partial_cmp(&a.rank).unwrap());
    let (phenotypes, ranks) = RankListItem::to_vecs(gene_list);
    // let mut res: Vec<GSEAResult<f64, Vec<String>>> = Vec::new();
    gmt.par_iter().for_each(|x| {
        gene_set_p(&phenotypes, &ranks, x, 1.0, 1000);
    });
}
