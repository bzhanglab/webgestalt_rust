/// Collection of utilities to save results to a file/folder


use std::fs::File;
use std::io::prelude::*;


pub fn save_nta(
    path: String,
    result: crate::methods::nta::NTAResult,
) -> Result<(), Box<std::io::Error>> {
    let mut file = File::create(path)?;
    let json = serde_json::to_string(&result).unwrap();
    file.write_all(json.as_bytes())?;
    Ok(())
}