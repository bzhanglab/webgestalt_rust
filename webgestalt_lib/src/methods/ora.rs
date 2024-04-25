use crate::{readers::utils::Item, stat};
use ahash::AHashSet;
use rayon::prelude::*;
use serde::Serialize;
use statrs::distribution::{DiscreteCDF, Hypergeometric};

#[derive(Clone)]
pub struct ORAConfig {
    pub min_overlap: i64,
    pub min_set_size: usize,
    pub max_set_size: usize,
    pub fdr_method: stat::AdjustmentMethod,
}

impl Default for ORAConfig {
    fn default() -> Self {
        ORAConfig {
            min_overlap: 5,
            min_set_size: 5,
            max_set_size: 500,
            fdr_method: stat::AdjustmentMethod::BH,
        }
    }
}

#[derive(Debug, Serialize, Clone)]
pub struct ORAResult {
    pub set: String,
    pub p: f64,
    pub fdr: f64,
    pub overlap: i64,
    pub expected: f64,
    pub enrichment_ratio: f64,
}

#[derive(Debug, Clone)]
struct PartialORAResult {
    set: String,
    p: f64,
    overlap: i64,
    expected: f64,
}

pub fn ora_p(m: i64, j: i64, n: i64, k: i64) -> f64 {
    let result = Hypergeometric::new(m as u64, j as u64, n as u64).unwrap();
    result.sf((k - 1) as u64)
}

/// Get ORA results for the provided interest list and reference list against the GMT file.
/// Requires both the interest list and the reference list to be filtered.
///
/// # Parameters
/// - `interest_list` - A [`AHashSet<String>`] of the interesting analytes
/// - `reference` - A [`AHashSet<String>`] of the reference list
/// - `gmt` - A [`Vec<Item>`] of the gmt file
pub fn get_ora(
    interest_list: &AHashSet<String>,
    reference: &AHashSet<String>,
    gmt: Vec<Item>,
    config: ORAConfig,
) -> Vec<ORAResult> {
    let m: i64 = reference.len() as i64;
    let n: i64 = interest_list.len() as i64;
    let partials: Vec<PartialORAResult> = gmt
        .par_iter()
        .map(|i| {
            let mut j: i64 = 0;
            let mut enriched_parts: AHashSet<String> = AHashSet::default();
            let mut k: i64 = 0;
            for analyte in i.parts.iter() {
                if interest_list.contains(analyte) {
                    k += 1;
                    enriched_parts.insert(analyte.to_owned());
                }
                if reference.contains(analyte) {
                    j += 1;
                }
            }
            let p = if k == 0 { 1.0 } else { ora_p(m, j, n, k) };
            PartialORAResult {
                set: i.id.clone(),
                p,
                overlap: k,
                expected: j as f64 * n as f64 / m as f64,
            }
        })
        .collect();
    let p_vals: Vec<f64> = partials.iter().map(|x| x.p).collect();
    let fdrs: Vec<f64> = stat::adjust(&p_vals, config.fdr_method);
    let mut final_res = Vec::new();
    for (i, row) in partials.clone().into_iter().enumerate() {
        final_res.push(ORAResult {
            set: row.set,
            p: row.p,
            fdr: fdrs[i],
            overlap: row.overlap,
            expected: row.expected,
            enrichment_ratio: row.overlap as f64 / row.expected,
        })
    }
    final_res
}
