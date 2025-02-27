// Experimental methods for WebGestalt

pub fn integrated_gsea(
    gmt_file: &str,
    rank_files: &[String],
) -> Result<(), Box<dyn std::error::Error>> {
    let gmt = crate::readers::read_gmt_file(gmt_file.to_string())
        .expect("Could not sucessfully read GMT file!");
    Ok(())
}
