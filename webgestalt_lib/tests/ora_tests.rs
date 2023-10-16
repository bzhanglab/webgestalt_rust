use pretty_assertions::assert_eq;
use statrs::assert_almost_eq;
use webgestalt_lib::{self, methods::ora::ORAConfig};
const THRESHOLD: f64 = 0.0001;
#[test]
fn ora() {
    let (gmt, gene_list, reference) = webgestalt_lib::readers::read_ora_files(
        "data/test.gmt".to_owned(),
        "data/genelist.txt".to_owned(),
        "data/reference.txt".to_owned(),
    );
    let gmtcount = gmt.len();
    let x: Vec<webgestalt_lib::methods::ora::ORAResult> =
        webgestalt_lib::methods::ora::get_ora(&gene_list, &reference, gmt, ORAConfig::default());
    let res = x.iter().find(|x| x.set == "GO:2000147").unwrap();
    assert_eq!(gmtcount, 850);
    assert_almost_eq!(res.p, 0.004516370110462129, THRESHOLD);
}
