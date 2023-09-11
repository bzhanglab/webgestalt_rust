use crate::readers::utils::Item;

use rand::prelude::SliceRandom;
use rand::SeedableRng;
use rayon::prelude::*;
use rustc_hash::FxHashSet;
use std::sync::atomic::AtomicI32;
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
    let inverse_size_dif: f64 = 1.0 / ((genes.len() - item.parts.len()) as f64); // Inverse now,
    let gene_size = genes.len();
    for j in 0..gene_size {
        if gene_set.contains(&genes[j]) {
            n_r += ranks[j].powf(p).abs();
        }
    }
    if n_r == 0.0 {
        // No GSEA needed.
        (
            GSEAResult {
                phenotype: item.id.clone(),
                p: 1.0,
                nes: 1.0,
                es: 1.0,
            },
            Vec::new(),
        )
    } else {
        let inverse_nr = 1.0 / n_r;
        let original_order = 0..gene_size;
        let has_gene: Vec<bool> = genes.par_iter().map(|x| gene_set.contains(x)).collect();
        let new_ranks: Vec<f64> = if p != 1.0 {
            ranks.par_iter().map(|x| x.powf(p)).collect()
        } else {
            ranks.to_vec()
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
            // let new_order = permutations_vec[i].clone();
            perm_es.lock().unwrap().push(enrichment_score(
                &has_gene,
                &new_ranks,
                &permutations_vec[i],
                inverse_size_dif,
                inverse_nr,
            ));
        });
        let nes_iter = perm_es.lock().unwrap();
        let len: f64 = nes_iter.len() as f64;
        let avg: f64 = nes_iter.iter().sum::<f64>() / len;
        let nes_es: Vec<f64> = nes_iter.iter().map(|x| x / avg.abs()).collect();
        let norm_es: f64 = real_es / avg.abs();
        let side: Vec<&f64> = if real_es >= 0_f64 {
            nes_iter.iter().filter(|x| *x >= &0_f64).collect()
        } else {
            nes_iter.iter().filter(|x| *x <= &0_f64).collect()
        };
        let tot = side.len();
        let p: f64 =
            side.into_iter().filter(|x| x.abs() > real_es.abs()).count() as f64 / tot as f64;
        (
            GSEAResult {
                phenotype: item.id.clone(),
                p,
                nes: norm_es,
                es: real_es,
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
            sum_hits += ranks[i].abs();
        } else {
            sum_miss += 1.0;
        }
        let es = sum_hits * inverse_nr - sum_miss * inverse_size_dif;
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
        let new_order = (0..(phenotypes.len()))
            .collect::<Vec<usize>>()
            .choose_multiple(&mut smallrng, phenotypes.len())
            .copied()
            .collect();
        permutations.push(new_order);
    });
    let all_nes = Arc::new(Mutex::new(Vec::new()));
    let set_nes = Arc::new(Mutex::new(Vec::new()));
    let all_res = Arc::new(Mutex::new(Vec::new()));
    gmt.par_iter().for_each(|x| {
        let (y, nes_iter) = gene_set_p(&phenotypes, &ranks, x, 1.0, &permutations);
        all_nes.lock().unwrap().extend(nes_iter);
        set_nes.lock().unwrap().push(y.nes);
        all_res.lock().unwrap().push(y);
    });
    let all_nes_l = all_nes.lock().unwrap();
    let set_nes_l = set_nes.lock().unwrap();
    let all_res_l = all_res.lock().unwrap();
    let mut final_gsea: Vec<FullGSEAResult> = Vec::new();
    for i in 0..set_nes_l.len() {
        let top_side: Vec<&f64> = if set_nes_l[i] >= 0_f64 {
            all_nes_l.par_iter().filter(|&x| x >= &0_f64).collect()
        } else {
            all_nes_l.par_iter().filter(|&x| x <= &0_f64).collect()
        };
        let top_len = top_side.len();
        let top_val = top_side
            .par_iter()
            .filter(|&x| x.abs() > set_nes_l[i].abs())
            .count() as f64
            / top_len as f64;
        let bottom_side: Vec<&f64> = if set_nes_l[i] >= 0_f64 {
            set_nes_l.par_iter().filter(|&x| x >= &0_f64).collect()
        } else {
            set_nes_l.par_iter().filter(|&x| x <= &0_f64).collect()
        };
        let bottom_len = bottom_side.len();
        let bottom_val = bottom_side
            .par_iter()
            .filter(|&x| x.abs() > set_nes_l[i].abs())
            .count() as f64
            / bottom_len as f64;
        let fdr: f64 = top_val / bottom_val;
        final_gsea.push(all_res_l[i].add_fdr(fdr));
    }
    let mut sigs: i32 = 0;
    for res in final_gsea {
        if res.p < 0.05 {
            println!(
                "{:?}: p: {:?}, fdr: {:?}, es: {:?}, nes: {:?}",
                res.phenotype, res.p, res.fdr, res.es, res.nes
            );
            sigs += 1;
        }
    }
    println!("Found {:?} significant pathways.", sigs)
}
