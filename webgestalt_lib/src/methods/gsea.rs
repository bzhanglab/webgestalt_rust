use crate::readers::utils::Item;

use rand::prelude::SliceRandom;
use rand::SeedableRng;
use rayon::prelude::*;
use rustc_hash::FxHashSet;
use std::sync::{Arc, Mutex};

#[derive(Clone)]
pub struct RankListItem {
    pub analyte: String,
    pub rank: f64,
}

pub struct GSEAResult {
    // TODO: Look at adding enrichment and normalized enrichment score
    pub set: String,
    pub p: f64,
    pub es: f64,
    pub nes: f64,
    overlap: i32,
    leading_edge: i32,
}

impl GSEAResult {
    pub fn add_fdr(&self, fdr: f64) -> FullGSEAResult {
        FullGSEAResult {
            set: self.set.clone(),
            p: self.p,
            fdr,
            es: self.es,
            nes: self.nes,
            leading_edge: self.leading_edge,
        }
    }
}

#[derive(Clone)]
pub struct FullGSEAResult {
    pub set: String,
    pub p: f64,
    pub fdr: f64,
    pub es: f64,
    pub nes: f64,
    pub leading_edge: i32,
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

/// Run GSEA for one analyte set.
///
/// Returns a [`GSEAResult`], which does not have FDR.
/// # Parameters
///
/// - `analytes` - Vector containing the names of the analytes
/// - `ranks` - Slice of [`f64`] containing the rank values corresponding to the `analytes` vector
/// - `item` - [`Item`] of the analyte set
/// - `p` - The power to raise the ranks. **Not the statistical p-value**
/// - `permutations_vec` - Vector of Vectors of random permutations of [`usize`] that are indices
/// used to shuffle analytes
///
/// # Panics
///
/// Panics if the `ranks` and `analytes` parameters are not the same length.
pub fn analyte_set_p(
    analytes: &Vec<String>,
    ranks: &[f64],
    item: &Item,
    p: f64,
    permutations_vec: &Vec<Vec<usize>>,
) -> (GSEAResult, Vec<f64>) {
    let permutations = permutations_vec.len();
    let analyte_set = FxHashSet::from_iter(item.parts.iter());
    let mut n_r: f64 = 0.0;
    let inverse_size_dif: f64 = 1.0 / ((analytes.len() - analyte_set.len()) as f64); // Inverse now
    let analyte_count = analytes.len();
    let mut overlap: i32 = 0;
    for j in 0..analyte_count {
        // Calculate N_r
        if analyte_set.contains(&analytes[j]) {
            n_r += ranks[j].abs().powf(p);
            overlap += 1;
        }
    }
    if overlap < 15 || overlap > 500 {
        // TODO: Make numbers parameters (maybe a struct)
        (
            GSEAResult {
                // No GSEA needed
                set: item.id.clone(),
                p: 1.0,
                nes: 0.0,
                es: 0.0,
                overlap,
                leading_edge: 0,
            },
            Vec::new(),
        )
    } else {
        let inverse_nr = 1.0 / n_r; // Invert n_r for the enrichment score
        let original_order = (0..analyte_count).collect::<Vec<usize>>(); // get regular order
        let has_analyte: Vec<bool> = analytes.iter().map(|x| analyte_set.contains(x)).collect(); // create
                                                                                                 // vector of booleans where T means analyte set has analyte at that index
        let new_ranks: Vec<f64> = if p != 1.0 {
            // raise to power p if needed and take abs
            ranks.par_iter().map(|x| x.abs().powf(p)).collect()
        } else {
            ranks.par_iter().map(|x| x.abs()).collect()
        };
        let (real_es, max_hits) = enrichment_score(
            // get normal es and leading edge max_hits
            &has_analyte,
            &new_ranks,
            &original_order,
            inverse_size_dif,
            inverse_nr,
            false,
        );
        let mut es_iter = Vec::new(); // Not parallelized because locking is expensive
        (0..permutations).for_each(|i| {
            // get es for the permutations
            let (p_es, _) = enrichment_score(
                &has_analyte,
                &new_ranks,
                &permutations_vec[i],
                inverse_size_dif,
                inverse_nr,
                true,
            );
            es_iter.push(p_es);
        });
        let side: Vec<&f64> = if real_es >= 0_f64 {
            // get side of distribution for p value
            es_iter.iter().filter(|x| *x >= &0_f64).collect()
        } else {
            es_iter.iter().filter(|x| *x < &0_f64).collect()
        };
        let tot = side.len();
        let p: f64 = if tot != 0 {
            // calculate p value
            side.into_iter()
                .filter(|x| x.abs() >= real_es.abs())
                .count() as f64
                / tot as f64
        } else {
            // no higher values found, so p is '0.0'. Previously < 2.2e-16 on the R version
            0.0
        };
        let up: Vec<f64> = es_iter // get positive (up) ES
            .par_iter()
            .filter(|&x| *x >= 0_f64)
            .copied()
            .collect();
        let down: Vec<f64> = es_iter // down scores
            .par_iter()
            .filter(|&x| *x < 0_f64)
            .copied()
            .collect();
        let up_len = up.len();
        let down_len = down.len();
        let up_avg: f64 = up.iter().sum::<f64>() / (up_len as f64 + 0.000001) + 0.000001; // up
                                                                                          // average
        let down_avg: f64 = down.iter().sum::<f64>() / (down_len as f64 + 0.000001) - 0.000001; // down
                                                                                                // average
        let mut nes_es: Vec<f64> = up.par_iter().map(|x| x / up_avg).collect(); // get all
                                                                                // normalized scores for up
        nes_es.extend(down.par_iter().map(|x| -x / down_avg).collect::<Vec<f64>>()); // extend with
                                                                                     // down scores
        let norm_es: f64 = if real_es >= 0_f64 {
            // get normalized score for the real run
            real_es / up_avg
        } else {
            -real_es / down_avg
        };

        (
            GSEAResult {
                set: item.id.clone(),
                p,
                nes: norm_es,
                es: real_es,
                overlap,
                leading_edge: max_hits,
            },
            nes_es,
        )
    }
}

/// Calculates the enrichment score for the specified list.
fn enrichment_score(
    analytes: &Vec<bool>,
    ranks: &[f64],
    order: &[usize],
    inverse_size_dif: f64,
    inverse_nr: f64,
    is_perm: bool,
) -> (f64, i32) {
    let mut max_score: f64 = 0.0;
    let mut hits: i32 = 0;
    let mut max_hits: i32 = 0;
    let mut sum_hits: f64 = 0.0;
    let mut sum_miss: f64 = 0.0;
    let inv_nr = if is_perm {
        // calculate new N_r if it is a permutation
        let mut temp_n: f64 = 0.0;
        for i in 0..analytes.len() {
            if analytes[order[i]] {
                temp_n += ranks[i];
            }
        }
        1.0 / temp_n
    } else {
        inverse_nr
    };
    for i in 0..analytes.len() {
        // find max ES score
        if analytes[order[i]] {
            // found in gene set
            sum_hits += ranks[i];
            hits += 1;
        } else {
            // not in set
            sum_miss += 1.0;
        }
        let es = (sum_hits * inv_nr) - (sum_miss * inverse_size_dif); // iteration score
        if es.abs() > max_score.abs() {
            // if bigger deviation from zero, store
            max_score = es;
            max_hits = hits;
        }
    }
    (max_score, max_hits)
}

/// Run GSEA andf return a [`Vec<FullGSEAResult`] for all analayte sets.
///
/// # Parameters
///
/// - `analyte_list` - [`Vec<RankListItem>`] of the rank list
/// - `gmt` - [`Vec<Item>`] of gmt file
pub fn gsea(mut analyte_list: Vec<RankListItem>, gmt: Vec<Item>) -> Vec<FullGSEAResult> {
    println!("Starting GSEA Calculation.");
    analyte_list.sort_by(|a, b| b.rank.partial_cmp(&a.rank).unwrap()); // sort list
    let (analytes, ranks) = RankListItem::to_vecs(analyte_list.clone()); // seperate into vectors
    let mut smallrng = rand::rngs::SmallRng::from_entropy();
    let mut permutations: Vec<Vec<usize>> = Vec::new();
    (0..1000).for_each(|_i| {
        // get random permutations that are shared for all analyte sets
        let mut new_order: Vec<usize> = (0..analytes.len()).collect();
        new_order.shuffle(&mut smallrng);
        permutations.push(new_order);
    });
    let all_nes = Arc::new(Mutex::new(Vec::new()));
    let set_nes = Arc::new(Mutex::new(Vec::new()));
    let all_res = Arc::new(Mutex::new(Vec::new()));
    gmt.par_iter().for_each(|analyte_set| {
        // parallelized scoring of all sets
        let (y, nes_iter) = analyte_set_p(&analytes, &ranks, analyte_set, 1.0, &permutations);
        if y.overlap >= 15 && y.overlap <= 500 {
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
        // get all FDR values
        let nes = partial_results[i].nes;
        let top_side: Vec<&f64> = if nes > 0_f64 {
            // positive null distribution
            null_distribution
                .par_iter()
                .filter(|&x| x > &0_f64)
                .collect()
        } else {
            // negative null distribution
            null_distribution
                .par_iter()
                .filter(|&x| x < &0_f64)
                .collect()
        };
        let top_len = if top_side.is_empty() {
            // avoid dividing by 0
            0.000001
        } else {
            top_side.len() as f64
        };
        let nes_abs = nes.abs();
        let top_val = top_side.par_iter().filter(|&x| x.abs() >= nes_abs).count() as f64; // get
                                                                                          // count of scores higher than current NES
        let bottom_side: Vec<&f64> = if nes >= 0_f64 {
            observed_distribution
                .par_iter()
                .filter(|&x| x >= &0_f64)
                .collect()
        } else {
            observed_distribution
                .par_iter()
                .filter(|&x| x < &0_f64)
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
        let fdr: f64 = (top_val * bottom_len) / (bottom_val * top_len); // get FDR value
        final_gsea.push(partial_results[i].add_fdr(fdr));
    }
    let mut sigs: i32 = 0;
    for res in final_gsea.iter() {
        if res.p <= 0.05 {
            // basic reporting  TODO: REMOVE
            println!(
                "{}: p: {:?}, fdr: {:?}, es: {:?}, nes: {:?}",
                res.set, res.p, res.fdr, res.es, res.nes
            );
            sigs += 1;
        }
    }
    println!("Found {:?} significant pathways.", sigs);
    final_gsea
}
