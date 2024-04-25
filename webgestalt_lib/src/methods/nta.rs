use ndarray::{Array2, Axis, Zip};
use serde::Serialize;
use std::ops::Div;

#[derive(Debug, Clone)]
/// A struct representing the options for the NTA algorithm
pub struct NTAConfig {
    /// A vector of vectors of strings representing the edge list of the graph
    pub edge_list: Vec<Vec<String>>,
    /// A vector of strings representing the seeds
    pub seeds: Vec<String>,
    /// A float representing the reset probability during random walk (default: 0.5)
    pub reset_probability: f64,
    /// A float representing the tolerance for probability calculation
    pub tolerance: f64,
    /// The [`NTAMethod`] to use for the analysis
    pub method: Option<NTAMethod>,
}

#[derive(Debug, Clone)]
pub enum NTAMethod {
    Prioritize(usize),
    Expand(usize),
}

impl Default for NTAConfig {
    fn default() -> Self {
        NTAConfig {
            edge_list: vec![],
            seeds: vec![],
            reset_probability: 0.5,
            tolerance: 0.000001,
            method: None,
        }
    }
}

#[derive(Debug, Serialize)]
pub struct NTAResult {
    pub neighborhood: Vec<String>,
    pub scores: Vec<f64>,
    pub candidates: Vec<String>,
}

/// Performs network topology-based analysis using random walk to identify important nodes in a network
///
/// ## Parameters
///
/// - `config`: A [`NTAConfig`] struct containing the parameters for the analysis.
///
/// ## Returns
///
/// Returns a [`NTAResult`] struct containing the results from the analysis. Is [serde](https://serde.rs/) compatible.
pub fn get_nta(config: NTAConfig) -> NTAResult {
    let mut method = config.clone().method;
    if method.is_none() {
        method = Some(NTAMethod::Expand(10));
    }
    let mut nta_res = process_nta(config.clone());
    match method {
        Some(NTAMethod::Prioritize(size)) => {
            let only_seeds = nta_res
                .iter()
                .filter(|(node, _)| config.seeds.contains(node))
                .cloned()
                .collect::<Vec<(String, f64)>>();
            let mut neighborhood: Vec<String> = Vec::new();
            let mut candidates: Vec<String> = Vec::new();
            let mut scores: Vec<f64> = Vec::new();
            for (node, score) in only_seeds.iter().take(size) {
                scores.push(*score);
                neighborhood.push(node.clone());
                if neighborhood.len() < size {
                    candidates.push(node.clone());
                }
            }
            return NTAResult {
                neighborhood,
                scores,
                candidates,
            };
        }
        Some(NTAMethod::Expand(size)) => {
            nta_res = nta_res
                .iter()
                .filter(|(node, _)| !config.seeds.contains(node))
                .cloned()
                .collect::<Vec<(String, f64)>>();
            let mut neighborhood: Vec<String> = Vec::new();
            let mut scores: Vec<f64> = Vec::new();
            for (node, score) in nta_res.iter().take(size) {
                neighborhood.push(node.clone());
                scores.push(*score);
            }
            let candidates: Vec<String> = Vec::new();
            return NTAResult {
                neighborhood,
                scores,
                candidates,
            };
        }
        _ => {
            panic!("Invalid method");
        }
    }
}

/// Uses random walk to calculate the probabilities of each node being walked through
/// Returns [`Vec<String>`] representing the nodes in the neighborhood
///
/// ## Parameters
/// - `config` - A [`NTAOptions`] struct containing the edge list, seeds, neighborhood size, reset probability, and tolerance
///
/// ## Returns
///
/// Returns a [`Vec<(String, f64)>`] where the [`String`] is the original node name, and the following value is the random walk probability (higher is typically better)
pub fn process_nta(config: NTAConfig) -> Vec<(String, f64)> {
    println!("Building Graph");
    let unique_nodes = ahash::AHashSet::from_iter(config.edge_list.iter().flatten().cloned());
    let mut node_map: ahash::AHashMap<String, usize> = ahash::AHashMap::default();
    let mut reverse_map: ahash::AHashMap<usize, String> = ahash::AHashMap::default();
    for (i, node) in unique_nodes.iter().enumerate() {
        node_map.insert(node.clone(), i);
        reverse_map.insert(i, node.clone());
    }
    let mut graph = Array2::<f64>::zeros((unique_nodes.len(), unique_nodes.len()));
    for edge in config.edge_list.iter() {
        let node1 = node_map.get(&edge[0]).unwrap();
        let node2 = node_map.get(&edge[1]).unwrap();
        graph[[*node1, *node2]] = 1.0;
        graph[[*node2, *node1]] = 1.0;
    }
    println!("Calculating NTA");
    let node_indices: Vec<usize> = config
        .seeds
        .iter()
        .map(|seed| *node_map.get(seed).unwrap())
        .collect();
    let walk_res = random_walk_probability(
        &graph,
        &node_indices,
        config.reset_probability,
        config.tolerance,
    );
    let mut walk = walk_res.iter().enumerate().collect::<Vec<(usize, &f64)>>();
    walk.sort_by(|a, b| b.1.partial_cmp(a.1).unwrap());
    walk.iter()
        .map(|(i, p)| (reverse_map.get(&i).unwrap().clone(), **p))
        .collect()
}

/// calculates the probability each node will be walked when starting from the one of the seeds
///
/// ## Parameters
///
/// - `adj_matrix` - A 2d adjacency matrix, where 1 means the node at the row and column indices are connected
/// - `seed_indices` - a [`Vec<usize>`] of the indices of the seeds (starting points)
/// - `r` - a [`f64`] of the reset probability (default in WebGestaltR is 0.5)
/// - `tolerance` - the tolerance/threshold value in [`f64`] (WebGestaltR default is `1e-6`)
///
/// ## Output
///
/// Returns 1d array containing the probability for each node
fn random_walk_probability(
    adj_matrix: &ndarray::Array2<f64>,
    seed_indices: &Vec<usize>,
    r: f64,
    tolerance: f64,
) -> ndarray::Array1<f64> {
    let num_nodes = seed_indices.len() as f64;
    let de = adj_matrix.sum_axis(Axis(0));
    // de to 2d array
    let de = de.insert_axis(Axis(1));
    let temp = adj_matrix.t().div(de);
    let w = temp.t();
    let mut p0 = ndarray::Array1::from_elem(w.shape()[0], 0.0);
    for i in seed_indices {
        p0[*i] = 1.0 / num_nodes;
    }
    let mut pt = p0.clone();
    let mut pt1 = w.dot(&pt) * (1.0 - r) + (r * &p0);
    while Zip::from(&pt1)
        .and(&pt)
        .par_map_collect(|a, b| (a - b).abs())
        .sum()
        > tolerance
    {
        pt = pt1;
        pt1 = w.dot(&pt) * (1.0 - r) + (r * &p0);
    }
    pt1
}
