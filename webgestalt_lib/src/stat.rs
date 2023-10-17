struct Carrier {
    p: f64,
    original_order: usize,
}

pub fn adjust(p_vals: &[f64]) -> Vec<f64> {
    let mut carriers: Vec<Carrier> = p_vals
        .iter()
        .enumerate()
        .map(|(i, p)| Carrier {
            p: *p,
            original_order: i,
        })
        .collect();
    carriers.sort_by(|a, b| a.p.partial_cmp(&b.p).unwrap());
    let m = carriers.len();
    let mut fdr_vals = vec![0.0; m];
    for (i, carrier) in carriers.iter().enumerate() {
        let mut fdr = carrier.p * m as f64 / (i + 1) as f64;
        if fdr > 1.0 {
            fdr = 1.0;
        }
        fdr_vals[carrier.original_order] = fdr;
    }
    let mut prev_fdr = fdr_vals[0];
    for i in 1..fdr_vals.len() {
        if prev_fdr > fdr_vals[i] {
            fdr_vals[i] = prev_fdr;
        } else {
            prev_fdr = fdr_vals[i];
        }
    }
    fdr_vals
}
