use crate::readers::utils::Item;

use rand::prelude::SliceRandom;
use rand::rngs::SmallRng;
use rand::{Rng, SeedableRng};
use rayon::prelude::*;
use rustc_hash::FxHashSet;
use statrs::distribution::{ContinuousCDF, Empirical};
use std::sync::atomic::AtomicI32;
use std::{
    sync::{Arc, Mutex},
    time::Instant,
};

pub struct RankListItem {
    pub phenotype: String,
    pub rank: f64,
}

pub struct GSEAResult<T, A, B> {
    // TODO: Look at adding enrichment and normalized enrichment score
    pub phenotype: T,
    pub p: A,
    pub fdr: B,
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

fn gene_set_p(
    genes: &Vec<String>,
    ranks: &[f64],
    item: &Item,
    p: f64,
    permutations: usize,
) -> GSEAResult<String, f64, f64> {
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
    if n_r == 0.0 {
        GSEAResult {
            phenotype: item.id.clone(),
            p: 1.0,
            fdr: 1.0,
        }
    } else {
        let inverse_nr = 1.0 / n_r;
        let original_order = 0..(gene_size - 1);
        let has_gene: Vec<bool> = genes.par_iter().map(|x| gene_set.contains(x)).collect();
        let new_ranks: Vec<f64> = ranks.par_iter().map(|x| x.powf(p)).collect();
        let real_es = enrichment_score(
            &has_gene,
            &new_ranks,
            original_order.collect(),
            inverse_size_dif,
            inverse_nr,
        );
        let perm_es = Arc::new(Mutex::new(Vec::new()));
        (0..permutations).into_par_iter().for_each(|_i| {
            let mut smallrng = rand::rngs::SmallRng::from_entropy();
            let new_order = (0..(gene_size))
                .collect::<Vec<usize>>()
                .choose_multiple(&mut smallrng, gene_size)
                .map(|x| x.clone())
                .collect();
            perm_es.lock().unwrap().push(enrichment_score(
                &has_gene,
                &new_ranks,
                new_order,
                inverse_size_dif,
                inverse_nr,
            ));
        });
        let nes_iter = perm_es.lock().unwrap();
        let len: f64 = nes_iter.len() as f64;
        let avg: f64 = nes_iter.iter().sum::<f64>() / len;
        let nes_es: Vec<f64> = nes_iter.iter().map(|x| x / avg).collect();
        let norm_es = real_es / avg;
        let x = if real_es < 0_f64 {
            let p: f64 = 1.0
                - Empirical::from_vec(
                    nes_iter
                        .par_iter()
                        .filter(|x| x < &&0_f64)
                        .map(|x| x.clone())
                        .collect::<Vec<f64>>(),
                )
                .cdf(real_es);
            let fdr = 0.0;
            GSEAResult {
                phenotype: item.id.clone(),
                p,
                fdr,
            }
        } else {
            let p: f64 = 1.0
                - Empirical::from_vec(
                    nes_iter
                        .par_iter()
                        .filter(|x| x >= &&0_f64)
                        .map(|x| x.clone())
                        .collect::<Vec<f64>>(),
                )
                .cdf(real_es);
            let fdr = 0.0;
            GSEAResult {
                phenotype: item.id.clone(),
                p,
                fdr,
            }
        };
        x
    }
}

fn enrichment_score(
    genes: &Vec<bool>,
    ranks: &[f64],
    order: Vec<usize>,
    inverse_size_dif: f64,
    inverse_nr: f64,
) -> f64 {
    let mut max_score: f64 = 0.0;
    let mut sum_hits = 0.0;
    let mut sum_miss = 0.0;
    for i in 0..(genes.len() - 1) {
        if genes[order[i]] {
            sum_hits += ranks[i];
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

pub fn get_gsea(mut gene_list: Vec<RankListItem>, gmt: Vec<Item>) {
    println!("Starting GSEA Calculation.");
    gene_list.sort_by(|a, b| b.rank.partial_cmp(&a.rank).unwrap());
    let (phenotypes, ranks) = RankListItem::to_vecs(gene_list);
    // let mut res: Vec<GSEAResult<String, f64, f64>> = Vec::new();
    let sigs = AtomicI32::new(0);
    gmt.par_iter().for_each(|x| {
        let y = gene_set_p(&phenotypes, &ranks, x, 1.0, 1000);
        if y.p < 0.05 {
            sigs.fetch_add(1, std::sync::atomic::Ordering::Relaxed);
            println!("{:?}: {:?}, {:?}", x.id, y.p, y.fdr);
        }
    });
    println!("Found {:?} significant pathways.", sigs)
}
