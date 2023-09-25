use crate::readers::utils::Item;
use adjustp::{adjust, Procedure};
use rayon::prelude::*;
use rustc_hash::FxHashSet;
use statrs::distribution::{DiscreteCDF, Hypergeometric};
use std::sync::{Arc, Mutex};
#[derive(Debug)]
pub struct ORAResult {
    pub set: String,
    pub p: f64,
    pub fdr: f64,
    pub overlap: i64,
}

#[derive(Debug, Clone)]
pub struct PartialORAResult {
    pub set: String,
    pub p: f64,
    pub overlap: i64,
}

pub fn ora_p(m: i64, j: i64, n: i64, k: i64) -> f64 {
    let result = Hypergeometric::new(m as u64, j as u64, n as u64).unwrap();
    result.sf((k - 1) as u64)
}

/// Get ORA results for the provided interest list and reference list against the GMT file.
/// Requires both the interest list and the reference list to be filtered.
///
/// # Parameters
/// - `interest_list` - A [`FxHashSet<String>`] of the interesting analytes
/// - `reference` - A [`FxHashSet<String>`] of the reference list
/// - `gmt` - A [`Vec<Item>`] of the gmt file
///
/// # Panics
///
/// Panics if the [`Arc`] struggles to lock during parallelization.
///
/// # Errors
///
/// This function will return an error if .
pub fn get_ora(
    interest_list: &FxHashSet<String>,
    reference: &FxHashSet<String>,
    gmt: Vec<Item>,
) -> Vec<ORAResult> {
    let m: i64 = reference.len() as i64;
    let n: i64 = interest_list.len() as i64;
    let res = Arc::new(Mutex::new(Vec::new()));
    gmt.par_iter().for_each(|i| {
        if i.parts.len() >= 5 {
            let mut j: i64 = 0;
            let mut enriched_parts: FxHashSet<String> = FxHashSet::default();
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
            if k >= 5 {
                if i.id == "hsa05221" {
                    println!("{}, {}, {}, {}", m, j, n, k);
                }
                let p = ora_p(m, j, n, k);
                res.lock().unwrap().push(PartialORAResult {
                    set: i.id.clone(),
                    p,
                    overlap: k,
                });
            }
        }
    });

    let partials = res.lock().unwrap();
    let p_vals: Vec<f64> = partials.iter().map(|x| x.p).collect();
    let fdrs: Vec<f64> = adjust(&p_vals, Procedure::BenjaminiHochberg)
        .iter()
        .map(|x| x * 2.0)
        .collect();
    let mut final_res = Vec::new();
    for (i, row) in partials.clone().into_iter().enumerate() {
        final_res.push(ORAResult {
            set: row.set,
            p: row.p,
            fdr: fdrs[i],
            overlap: row.overlap,
        })
    }
    final_res
}
