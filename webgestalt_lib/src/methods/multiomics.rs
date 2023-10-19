use rustc_hash::FxHashSet;

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
    pub interest_list: &'a FxHashSet<String>,
    pub reference_list: &'a FxHashSet<String>,
    pub config: ORAConfig,
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
    jobs: Vec<AnalysisJob>,
    analysis_type: AnalysisType,
    method: MultiOmicsMethod,
) -> () {
    if let MultiOmicsMethod::Meta(meta_method) = method {
    } else {
    }
}
