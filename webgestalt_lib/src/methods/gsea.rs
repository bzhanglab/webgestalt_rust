use crate::readers::utils::Item;

use rand::prelude::SliceRandom;
use rand::SeedableRng;
use rayon::prelude::*;
use rustc_hash::FxHashSet;
use std::sync::{Arc, Mutex};

pub struct RankListItem {
    pub analyte: String,
    pub rank: f64,
}

pub struct GSEAResult {
    // TODO: Look at adding enrichment and normalized enrichment score
    pub phenotype: String,
    pub p: f64,
    pub es: f64,
    pub nes: f64,
    overlap: i32,
}

impl GSEAResult {
    pub fn add_fdr(&self, fdr: f64) -> FullGSEAResult {
        FullGSEAResult {
            phenotype: self.phenotype.clone(),
            p: self.p,
            fdr,
            es: self.es,
            nes: self.nes,
        }
    }
}

pub struct FullGSEAResult {
    pub phenotype: String,
    pub p: f64,
    pub fdr: f64,
    pub es: f64,
    pub nes: f64,
}

impl RankListItem {
    pub fn to_vecs(input: Vec<RankListItem>) -> (Vec<String>, Vec<f64>) {
        let mut analytes: Vec<String> = Vec::new();
        let mut ranks: Vec<f64> = Vec::new();
        for item in input {
            analytes.push(item.analyte);
            ranks.push(item.rank)
        }
        (analytes, ranks)
    }
}

fn median(array: &Vec<f64>) -> f64 {
    if array.is_empty() {
        0.0
    } else if (array.len() % 2) == 0 {
        let ind_left = array.len() / 2 - 1;
        let ind_right = array.len() / 2;
        (array[ind_left] + array[ind_right]) as f64 / 2.0
    } else {
        array[array.len() / 2] as f64
    }
}

fn gene_set_p(
    genes: &Vec<String>,
    ranks: &[f64],
    item: &Item,
    p: f64,
    permutations_vec: &Vec<Vec<usize>>,
) -> (GSEAResult, Vec<f64>) {
    let permutations = permutations_vec.len();
    let gene_set = FxHashSet::from_iter(item.parts.iter());
    let mut n_r: f64 = 0.0;
    let inverse_size_dif: f64 = 1.0 / ((genes.len() - gene_set.len()) as f64); // Inverse now,
    let gene_size = genes.len();
    let mut overlap: i32 = 0;
    for j in 0..gene_size {
        if gene_set.contains(&genes[j]) {
            n_r += ranks[j].abs().powf(p);
            overlap += 1;
        }
    }
    if overlap < 20 || overlap > 500 {
        // No GSEA needed.
        (
            GSEAResult {
                phenotype: item.id.clone(),
                p: 1.0,
                nes: 0.0,
                es: 0.0,
                overlap,
            },
            Vec::new(),
        )
    } else {
        let inverse_nr = 1.0 / n_r;
        let original_order = 0..gene_size;
        let has_gene: Vec<bool> = genes.iter().map(|x| gene_set.contains(x)).collect();
        let new_ranks: Vec<f64> = if p != 1.0 {
            ranks.par_iter().map(|x| x.abs().powf(p)).collect()
        } else {
            ranks.par_iter().map(|x| x.abs()).collect()
        };
        let real_es = enrichment_score(
            &has_gene,
            &new_ranks,
            &original_order.collect(),
            inverse_size_dif,
            inverse_nr,
        );
        let perm_es = Arc::new(Mutex::new(Vec::new()));
        (0..permutations).into_par_iter().for_each(|i| {
            perm_es.lock().unwrap().push(enrichment_score(
                &has_gene,
                &new_ranks,
                &permutations_vec[i],
                inverse_size_dif,
                inverse_nr,
            ));
        });
        let es_iter = perm_es.lock().unwrap();
        // es_iter.sort_by(|a, b| a.partial_cmp(b).unwrap());
        let up: Vec<f64> = es_iter
            .par_iter()
            .filter(|&x| *x >= 0_f64)
            .copied()
            .collect();
        let down: Vec<f64> = es_iter
            .par_iter()
            .filter(|&x| *x < 0_f64)
            .copied()
            .collect();
        if up.is_empty() || item.id == "hsa04723" || inverse_nr < 0.01 {
            println!("{:?}, {:?}, {:?}", item.id, inverse_nr, inverse_size_dif);
        }
        let up_len = up.len();
        let down_len = down.len();
        let up_avg: f64 = up.iter().sum::<f64>() / (up_len as f64 + 0.000001) + 0.000001;
        let down_avg: f64 = down.iter().sum::<f64>() / (down_len as f64 + 0.000001) - 0.000001;
        // up.sort_by(|a, b| a.partial_cmp(b).unwrap());
        // down.sort_by(|a, b| a.partial_cmp(b).unwrap());
        // let up_avg = median(&up);
        // let down_avg = median(&down);
        let mut nes_es: Vec<f64> = up.par_iter().map(|x| x / up_avg).collect();
        nes_es.extend(down.par_iter().map(|x| -x / down_avg).collect::<Vec<f64>>());
        let norm_es: f64 = if real_es >= 0_f64 {
            real_es / up_avg
        } else {
            -real_es / down_avg
        };
        let side: Vec<&f64> = if real_es >= 0_f64 {
            es_iter.iter().filter(|x| *x >= &0_f64).collect()
        } else {
            es_iter.iter().filter(|x| *x < &0_f64).collect()
        };
        let tot = side.len();
        let p: f64 = if tot != 0 {
            side.into_iter()
                .filter(|x| x.abs() >= real_es.abs())
                .count() as f64
                / tot as f64
        } else {
            0.0
        };
        (
            GSEAResult {
                phenotype: item.id.clone(),
                p,
                nes: norm_es,
                es: real_es,
                overlap,
            },
            nes_es,
        )
    }
}

fn enrichment_score(
    genes: &Vec<bool>,
    ranks: &[f64],
    order: &Vec<usize>,
    inverse_size_dif: f64,
    inverse_nr: f64,
) -> f64 {
    let mut max_score: f64 = 0.0;
    let mut sum_hits: f64 = 0.0;
    let mut sum_miss: f64 = 0.0;
    for i in 0..genes.len() {
        if genes[order[i]] {
            sum_hits += ranks[i];
        } else {
            sum_miss += 1.0;
        }
        let es = (sum_hits * inverse_nr) - (sum_miss * inverse_size_dif);
        if es.abs() > max_score.abs() {
            max_score = es;
        }
    }
    max_score
}

pub fn gsea(mut gene_list: Vec<RankListItem>, gmt: Vec<Item>) {
    println!("Starting GSEA Calculation.");
    gene_list.sort_by(|a, b| b.rank.partial_cmp(&a.rank).unwrap());
    let (phenotypes, ranks) = RankListItem::to_vecs(gene_list);
    let mut smallrng = rand::rngs::SmallRng::from_entropy();
    let mut permutations: Vec<Vec<usize>> = Vec::new();
    (0..1000).for_each(|_i| {
        let mut new_order = (0..(phenotypes.len())).collect::<Vec<usize>>();
        new_order.shuffle(&mut smallrng);
        permutations.push(new_order);
    });
    let all_nes = Arc::new(Mutex::new(Vec::new()));
    let set_nes = Arc::new(Mutex::new(Vec::new()));
    let all_res = Arc::new(Mutex::new(Vec::new()));
    gmt.par_iter().for_each(|gene_set| {
        let (y, nes_iter) = gene_set_p(&phenotypes, &ranks, gene_set, 1.0, &permutations);
        if y.overlap >= 20 && y.overlap <= 500 {
            all_nes.lock().unwrap().extend(nes_iter);
            set_nes.lock().unwrap().push(y.nes);
            all_res.lock().unwrap().push(y);
        }
    });
    let null_distribution = all_nes.lock().unwrap();
    let observed_distribution = set_nes.lock().unwrap();
    let partial_results: std::sync::MutexGuard<'_, Vec<GSEAResult>> = all_res.lock().unwrap();
    let mut final_gsea: Vec<FullGSEAResult> = Vec::new();
    for i in 0..partial_results.len() {
        let nes = partial_results[i].nes;
        let top_side: Vec<&f64> = if nes > 0_f64 {
            null_distribution
                .par_iter()
                .filter(|&x| x > &0_f64)
                .collect()
        } else {
            null_distribution
                .par_iter()
                .filter(|&x| x < &0_f64)
                .collect()
        };
        let top_len = if top_side.is_empty() {
            0.000001
        } else {
            top_side.len() as f64
        };
        let nes_abs = nes.abs();
        let top_val = top_side.par_iter().filter(|&x| x.abs() >= nes_abs).count() as f64;
        let bottom_side: Vec<&f64> = if nes >= 0_f64 {
            observed_distribution
                .par_iter()
                .filter(|&x| x >= &0_f64)
                .collect()
        } else {
            observed_distribution
                .par_iter()
                .filter(|&x| x <= &0_f64)
                .collect()
        };
        let bottom_len = if bottom_side.is_empty() {
            0.000001
        } else {
            bottom_side.len() as f64
        };
        let bottom_val = bottom_side
            .par_iter()
            .filter(|&x| x.abs() >= nes_abs)
            .count() as f64;
        let fdr: f64 = (top_val / top_len) / (bottom_val / bottom_len);
        if fdr < 0.7 {
            println!(
                "{:?}/{:?}",
                top_val / top_len as f64,
                bottom_val / bottom_len as f64
            );
        }
        final_gsea.push(partial_results[i].add_fdr(fdr));
    }
    let mut sigs: i32 = 0;
    for res in final_gsea {
        if res.p <= 0.05 {
            println!(
                "{:?}: p: {:?}, fdr: {:?}, es: {:?}, nes: {:?}",
                res.phenotype, res.p, res.fdr, res.es, res.nes
            );
            sigs += 1;
        }
    }
    println!("Found {:?} significant pathways.", sigs);
}
