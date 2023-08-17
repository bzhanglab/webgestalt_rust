use crate::readers::utils::Item;
use rayon::prelude::*;
use rustc_hash::{FxHashMap, FxHashSet};
use statrs::distribution::{Discrete, Hypergeometric};
use std::sync::{Arc, Mutex};

pub fn ora_p(big_n: i64, m: i64, n: i64, k: i64) -> f64 {
    let result = Hypergeometric::new((big_n + n) as u64, n as u64, (m + k) as u64).unwrap();
    let mut p: f64 = 0.0;
    for i in k..=(m + k) {
        let x = result.pmf(i as u64);
        if !x.is_nan() {
            p += x
        }
    }
    p
}

pub fn get_ora(
    gene_list: FxHashMap<String, bool>,
    reference: FxHashMap<String, bool>,
    gmt: Vec<Item>,
) -> Vec<f64> {
    let big_n: i64 = reference.len() as i64;
    let n: i64 = gene_list.len() as i64;
    let res = Arc::new(Mutex::new(Vec::new()));
    gmt.par_iter().for_each(|i| {
        let mut m: i64 = 0;
        let mut enriched_parts: FxHashSet<String> = FxHashSet::default();
        let mut k: i64 = 0;
        for j in i.parts.iter() {
            if gene_list.contains_key(j) {
                k += 1;
                enriched_parts.insert(j.to_owned());
            }
            if reference.contains_key(j) {
                m += 1;
            }
        }
        let p = ora_p(big_n, m, n, k);
        res.lock().unwrap().push(p);
    });
    Arc::try_unwrap(res).unwrap().into_inner().unwrap()
}
