struct Carrier {
    p: f64,
    original_order: usize,
}

#[derive(Clone)]
pub enum AdjustmentMethod {
    BH,
    None,
}

pub fn adjust(p_vals: &[f64], method: AdjustmentMethod) -> Vec<f64> {
    match method {
        AdjustmentMethod::BH => benjamini_hochberg(p_vals),
        AdjustmentMethod::None => p_vals.to_vec(),
    }
}

fn benjamini_hochberg(p_vals: &[f64]) -> Vec<f64> {
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
    let mut prev_fdr = 1.0;
    for (i, carrier) in carriers.iter().enumerate().rev() {
        let mut fdr = carrier.p * m as f64 / (i + 1) as f64;
        if fdr > 1.0 {
            fdr = 1.0;
        }
        if fdr > prev_fdr {
            fdr = prev_fdr;
        } else {
            prev_fdr = fdr;
        }
        fdr_vals[carrier.original_order] = fdr;
    }
    fdr_vals
}
