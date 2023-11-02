use pretty_assertions::assert_eq;
use statrs::assert_almost_eq;
use webgestalt_lib;

const THRESHOLD: f64 = 0.0001;

#[test]
fn read_gmt() {
    let gmt = webgestalt_lib::readers::read_gmt_file("data/ktest.gmt".to_string()).unwrap();
    assert_eq!(gmt.len(), 330)
}
