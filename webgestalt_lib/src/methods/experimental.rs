// Experimental methods for WebGestalt

use crate::readers::read_rank_file;

pub fn integrated_gsea(
    gmt_file: &str,
    rank_files: &[String],
) -> Result<(), Box<dyn std::error::Error>> {
    let gmt = crate::readers::read_gmt_file(gmt_file.to_string())
        .expect("Could not sucessfully read GMT file!");
    let mut ranks: Vec<(usize, String, f64)> = Vec::new();
    // For each rank file, create (list_index, analyte, score) tuple
    for (list_index, file_path) in rank_files.iter().enumerate() {
        let rank_file = read_rank_file(file_path.to_string())
            .unwrap_or_else(|_| panic!("Could not read {}", file_path));
        for item in rank_file {
            ranks.push((list_index, item.analyte.to_string(), item.rank));
        }
    }
    // Sort ranks from large to small
    ranks.sort_by(|a, b| b.2.total_cmp(&a.2));
    Ok(())
}
