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
  Fisher
}
