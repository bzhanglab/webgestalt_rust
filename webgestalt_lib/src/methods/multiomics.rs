use ahash::{AHashMap, AHashSet};
use statrs::distribution::{Continuous, ContinuousCDF, Normal};

use super::{
    gsea::{GSEAConfig, GSEAResult, RankListItem},
    ora::{get_ora, ORAConfig, ORAResult},
};
use crate::{methods::gsea::gsea, readers::utils::Item};

pub enum MultiOmicsMethod {
    /// Get the max median ratio of the analyte from any list
    Max(NormalizationMethod),
    /// Get the average median ratio of analyte from all the lists
    Mean(NormalizationMethod),
    /// Run each list separately and calculate a meta-p value
    Meta(MetaAnalysisMethod),
}

pub enum MetaAnalysisMethod {
    Stouffer,
    Fisher,
}

pub enum AnalysisType {
    /// Gene Set Enrichment Analysis
    GSEA,
    /// Over-representation Analysis
    ORA,
}

pub struct GSEAJob {
    pub gmt: Vec<Item>,
    pub rank_list: Vec<RankListItem>,
    pub config: GSEAConfig,
}

pub struct ORAJob {
    pub gmt: Vec<Item>,
    pub interest_list: AHashSet<String>,
    pub reference_list: AHashSet<String>,
    pub config: ORAConfig,
}

#[derive(Copy, Clone)]
pub enum NormalizationMethod {
    MedianRank,
    MedianValue,
    MeanValue,
    None,
}

/// Run a multiomics analysis, using either the max/mean median ratio or a typical meta analysis
/// method
///
/// # Parameters
///
/// - `jobs` - A [`Vec<GSEAJob>`] containing all of the separates 'jobs' or analysis to combine
/// - `method` - A [`MultiOmicsMethod`] enum detailing the analysis method to combine the runs together (meta-analysis, mean median ration, or max median ratio).
///
/// # Returns
///
/// Returns a [`Vec<Vec<FullGSEAResult>>`] containing the results of each analysis. If the method was not meta-analysis, then the outer vector will only have one element.
/// If the method was meta-analysis, then the first element will be the results of the meta-analysis, and the rest of the elements will be the results of each analysis run individually.
pub fn multiomic_gsea(jobs: Vec<GSEAJob>, method: MultiOmicsMethod) -> Vec<Vec<GSEAResult>> {
    if let MultiOmicsMethod::Meta(meta_method) = method {
        let mut phash: AHashMap<String, Vec<f64>> = AHashMap::default();
        let mut results: Vec<Vec<GSEAResult>> = Vec::new();
        for job in jobs {
            let res = gsea(job.rank_list, job.gmt, job.config, None);
            for row in res.iter() {
                let set = row.set.clone();
                phash.entry(set).or_default().push(row.p);
            }
            results.push(res);
        }
        let mut final_result: Vec<GSEAResult> = Vec::new();
        match meta_method {
            MetaAnalysisMethod::Stouffer => {
                let normal = Normal::new(0.0, 1.0).unwrap();
                for set in phash.keys() {
                    final_result.push(GSEAResult {
                        set: set.clone(),
                        p: stouffer_with_normal(&phash[set], &normal),
                        fdr: 0.0,
                        nes: 0.0,
                        es: 0.0,
                        running_sum: Vec::new(),
                        leading_edge: 0,
                    });
                }
            }
            MetaAnalysisMethod::Fisher => {
                for set in phash.keys() {
                    final_result.push(GSEAResult {
                        set: set.clone(),
                        p: fisher(&phash[set]),
                        fdr: 0.0,
                        nes: 0.0,
                        es: 0.0,
                        running_sum: Vec::new(),
                        leading_edge: 0,
                    });
                }
            }
        }
        results.insert(0, final_result);
        results
    } else {
        let lists = jobs.iter().map(|x| x.rank_list.clone()).collect();
        let combined_list = combine_lists(lists, method);
        let gmts = jobs.iter().map(|x| x.gmt.clone()).collect();
        let combined_gmt = combine_gmts(&gmts);
        vec![gsea(
            combined_list,
            combined_gmt,
            jobs.first().unwrap().config.clone(),
            None,
        )]
    }
}

pub fn multiomic_ora(jobs: Vec<ORAJob>, method: MultiOmicsMethod) -> Vec<Vec<ORAResult>> {
    match method {
        MultiOmicsMethod::Meta(meta_method) => {
            let mut phash: AHashMap<String, Vec<f64>> = AHashMap::default();
            let mut results: Vec<Vec<ORAResult>> = Vec::new();
            for job in jobs {
                let res = get_ora(&job.interest_list, &job.reference_list, job.gmt, job.config);
                for row in res.iter() {
                    let set = row.set.clone();
                    phash.entry(set).or_default().push(row.p);
                }
                results.push(res);
            }
            let mut final_result: Vec<ORAResult> = Vec::new();
            match meta_method {
                MetaAnalysisMethod::Stouffer => {
                    let normal = Normal::new(0.0, 1.0).unwrap();
                    for set in phash.keys() {
                        final_result.push(ORAResult {
                            set: set.clone(),
                            p: stouffer_with_normal(&phash[set], &normal),
                            fdr: 0.0,
                            overlap: 0,
                            expected: 0.0,
                            enrichment_ratio: 0.0,
                        });
                    }
                }
                MetaAnalysisMethod::Fisher => {
                    for set in phash.keys() {
                        final_result.push(ORAResult {
                            set: set.clone(),
                            p: fisher(&phash[set]),
                            fdr: 0.0,
                            overlap: 0,
                            expected: 0.0,
                            enrichment_ratio: 0.0,
                        });
                    }
                }
            }
            results.insert(0, final_result);
            results
        }
        _ => {
            panic!("Multi-Omics ORA can only be run with meta-analysis");
        }
    }
}

pub fn combine_lists(
    lists: Vec<Vec<RankListItem>>,
    combination_method: MultiOmicsMethod,
) -> Vec<RankListItem> {
    match combination_method {
        MultiOmicsMethod::Max(normalization_method) => max_combine(lists, normalization_method),
        MultiOmicsMethod::Mean(normalization_method) => mean_combine(lists, normalization_method),
        MultiOmicsMethod::Meta(_) => panic!("Lists can not be combined for meta-analysis"),
    }
}

fn max_combine(
    lists: Vec<Vec<RankListItem>>,
    normalization_method: NormalizationMethod,
) -> Vec<RankListItem> {
    let normalized_lists: Vec<Vec<RankListItem>> = lists
        .into_iter()
        .map(|mut list| normalize(&mut list, normalization_method))
        .collect();
    let mut batches: AHashMap<String, f64> = AHashMap::default();
    for list in normalized_lists {
        for item in list {
            if let Some(val) = batches.get_mut(&item.analyte) {
                if item.rank.abs() > *val {
                    *val = item.rank;
                }
            } else {
                batches.insert(item.analyte, item.rank);
            }
        }
    }
    let mut final_list: Vec<RankListItem> = Vec::new();
    for key in batches.keys() {
        final_list.push(RankListItem {
            analyte: key.clone(),
            rank: batches[key],
        });
    }
    final_list
}

fn mean_combine(
    lists: Vec<Vec<RankListItem>>,
    normalization_method: NormalizationMethod,
) -> Vec<RankListItem> {
    let normalized_lists: Vec<Vec<RankListItem>> = lists
        .into_iter()
        .map(|mut list| normalize(&mut list, normalization_method))
        .collect();
    let mut batches: AHashMap<String, Vec<f64>> = AHashMap::default();
    for list in normalized_lists {
        for item in list {
            if let Some(val) = batches.get_mut(&item.analyte) {
                val.push(item.rank);
            } else {
                batches.insert(item.analyte, vec![item.rank]);
            }
        }
    }
    let mut final_list: Vec<RankListItem> = Vec::new();
    for key in batches.keys() {
        final_list.push(RankListItem {
            analyte: key.clone(),
            rank: batches[key].iter().sum::<f64>() / (batches[key].len() as f64),
        })
    }
    final_list
}

fn normalize(list: &mut Vec<RankListItem>, method: NormalizationMethod) -> Vec<RankListItem> {
    match method {
        NormalizationMethod::None => list.clone(),
        NormalizationMethod::MedianRank => {
            list.sort_by(|a, b| {
                a.rank
                    .partial_cmp(&b.rank)
                    .expect("Invalid float comparison during normalization")
            });
            let median = list.len() as f64 / 2.0;
            let mut final_list: Vec<RankListItem> = Vec::new();
            for (i, item) in list.iter().enumerate() {
                final_list.push(RankListItem {
                    analyte: item.analyte.clone(),
                    rank: (i as f64 - median) / median,
                });
            }
            final_list
        }
        NormalizationMethod::MedianValue => {
            list.sort_by(|a, b| {
                b.rank
                    .partial_cmp(&a.rank)
                    .expect("Invalid float comparison during normalization")
            });
            let min = list.last().unwrap().rank;
            let median = list[list.len() / 2].rank - min;
            let shift = min / median;
            let mut final_list: Vec<RankListItem> = Vec::new();
            for item in list.iter() {
                final_list.push(RankListItem {
                    analyte: item.analyte.clone(),
                    rank: (item.rank - min) / median + shift,
                });
            }
            final_list
        }
        NormalizationMethod::MeanValue => {
            list.sort_by(|a, b| {
                b.rank
                    .partial_cmp(&a.rank)
                    .expect("Invalid float comparison during normalization")
            });
            let min = list.last().unwrap().rank;
            let mean: f64 = list.iter().map(|x| x.rank - min).sum::<f64>() / (list.len() as f64)
                - min / (list.len() as f64);
            let shift = min / mean;
            let mut final_list: Vec<RankListItem> = Vec::new();
            for item in list.iter() {
                final_list.push(RankListItem {
                    analyte: item.analyte.clone(),
                    rank: (item.rank - min) / mean + shift,
                });
            }
            final_list
        }
    }
}

pub fn combine_gmts(gmts: &Vec<Vec<Item>>) -> Vec<Item> {
    let mut combined_parts: AHashMap<String, Vec<String>> = AHashMap::default();
    let mut combined_urls: AHashMap<String, String> = AHashMap::default();
    for gmt in gmts {
        for item in gmt {
            if combined_parts.contains_key(&item.id) {
                combined_parts
                    .get_mut(&item.id)
                    .unwrap()
                    .extend(item.parts.clone());
            } else {
                combined_parts.insert(item.id.clone(), item.parts.clone());
                combined_urls.insert(item.id.clone(), item.url.clone());
            }
        }
    }
    let mut final_gmt: Vec<Item> = Vec::new();
    for (key, parts) in combined_parts {
        final_gmt.push(Item {
            id: key.clone(),
            parts,
            url: combined_urls[&key].clone(),
        })
    }
    final_gmt
}

/// Calculates meta-p values using the Stouffer method ([DOI:10.1037/h0051438](https://doi.org/10.1037/h0051438)) of `vals`
///
/// # Arguments
/// - `val` - `Vec<f64>` of p-values to combine
///
/// # Examples
///
/// ```rust
/// use webgestalt_lib::methods::multiomics::stouffer;
/// let vals: Vec<f64> = vec![0.1, 0.01, 0.11, 0.23];
/// let metap: f64 = stouffer(vals);
/// ```
pub fn stouffer(vals: Vec<f64>) -> f64 {
    let n = Normal::new(0.0, 1.0).unwrap();
    let k = vals.len();
    n.cdf(vals.iter().map(|x| n.inverse_cdf(*x)).sum::<f64>() / f64::sqrt(k as f64))
}

fn stouffer_with_normal(vals: &Vec<f64>, normal: &Normal) -> f64 {
    let k = vals.len();
    normal.cdf(vals.iter().map(|x| normal.inverse_cdf(*x)).sum::<f64>() / f64::sqrt(k as f64))
}

pub fn fisher(vals: &Vec<f64>) -> f64 {
    let k = vals.len();
    let pt = -2.0 * vals.iter().map(|x| x.ln()).sum::<f64>();
    let dist = statrs::distribution::ChiSquared::new(2_f64.powi(k as i32 - 1)).unwrap();
    dist.pdf(pt)
}

/// Calculates meta-p values using the Stouffer weighted method ([10.1214/aoms/1177698861](https://doi.org/10.1214/aoms/1177698861)) of `vals` with weights in `weights`
///
/// # Arguments
/// - `val` - [`Vec<f64>`] of p-values to combine
/// - `weights` - [`Vec<f64>`] of weights corresponding to each p-value
///
/// # Examples
///
/// ```rust
/// use webgestalt_lib::methods::multiomics::stouffer_weighted;
/// let vals: Vec<f64> = vec![0.1, 0.01, 0.11, 0.23];
/// let weights: Vec<f64> = vec![0.1, 0.2, 0.3, 0.4];
/// let metap: f64 = stouffer_weighted(vals, weights);
/// ```
pub fn stouffer_weighted(vals: Vec<f64>, weights: Vec<f64>) -> f64 {
    let n = Normal::new(0.0, 1.0).unwrap();
    n.cdf(
        vals.iter()
            .enumerate()
            .map(|(i, x)| weights[i] * n.inverse_cdf(*x))
            .sum::<f64>()
            / f64::sqrt(weights.iter().map(|x| x * x).sum::<f64>()),
    )
}
