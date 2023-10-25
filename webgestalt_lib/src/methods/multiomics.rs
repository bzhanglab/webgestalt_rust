use ahash::{AHashMap, AHashSet};

use super::{
    gsea::{GSEAConfig, RankListItem},
    ora::ORAConfig,
};
use crate::readers::utils::Item;

pub enum MultiOmicsMethod {
    /// Get the max median ratio of the analyte from any list
    Max,
    /// Get the average median ratio of analyte from all the lists
    Mean,
    /// Run each list separately and calculate a meta-p value
    Meta(MetaAnalysisMethod),
}

pub enum MetaAnalysisMethod {
    Stouffer,
    Fisher,
}

pub enum AnalysisType {
    /// Gene Set Enrichment Analysips
    GSEA,
    /// Over-representation Analysis
    ORA,
}

pub enum AnalysisJob<'a> {
    GSEAJob(GSEAJob<'a>),
    ORAJob(ORAJob<'a>),
}

pub struct GSEAJob<'a> {
    pub gmt: &'a Vec<Item>,
    pub rank_list: &'a Vec<RankListItem>,
    pub config: GSEAConfig,
}

pub struct ORAJob<'a> {
    pub gmt: &'a Vec<Item>,
    pub interest_list: &'a AHashSet<String>,
    pub reference_list: &'a AHashSet<String>,
    pub config: ORAConfig,
}

#[derive(Copy, Clone)]
pub enum NormalizationMethod {
    MedianRank,
    MedianValue,
    MeanValue,
    None,
}

/// Run a multiomics analysis, using iehter the max/mean median ratio or a typical meta analysis
/// method
///
/// # Parameters
///
/// - `jobs` - A [`Vec<AnalysisJob>`] containing all of the seperates 'jobs' or analysis to combine
/// - `analysis_type` - A [`AnalysisType`] enum of the analysis type (GSEA, ORA, or NTA) to run
/// - `method` - A [`MultiOmicsMethod`] enum detailing the analysis method to combine the runs
/// together (meta-analysis, mean median ration, or max median ratio).
pub fn multiomic_analysis(
    _jobs: Vec<AnalysisJob>,
    _analysis_type: AnalysisType,
    method: MultiOmicsMethod,
) {
    if let MultiOmicsMethod::Meta(meta_method) = method {
    } else {
    }
}

pub fn combine_lists(
    lists: Vec<Vec<RankListItem>>,
    combination_method: MultiOmicsMethod,
    normalization_method: NormalizationMethod,
) -> Vec<RankListItem> {
    match combination_method {
        MultiOmicsMethod::Max => max_combine(lists, normalization_method),
        MultiOmicsMethod::Mean => mean_combine(lists, normalization_method),
        MultiOmicsMethod::Meta(_x) => panic!("Lists can not be combine for meta-analysis"),
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
                b.rank
                    .partial_cmp(&a.rank)
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
            let mut final_list: Vec<RankListItem> = Vec::new();
            for item in list.iter() {
                final_list.push(RankListItem {
                    analyte: item.analyte.clone(),
                    rank: (item.rank - min) / median,
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
            let mean: f64 = list.iter().map(|x| x.rank).sum::<f64>() / (list.len() as f64) - min;
            let mut final_list: Vec<RankListItem> = Vec::new();
            for item in list.iter() {
                final_list.push(RankListItem {
                    analyte: item.analyte.clone(),
                    rank: (item.rank - min) / mean,
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
