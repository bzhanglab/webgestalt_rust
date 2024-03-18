use ndarray::{Array2, Axis, Zip};
use std::ops::Div;

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
}

impl Default for NTAConfig {
    fn default() -> Self {
        NTAConfig {
            edge_list: vec![],
            seeds: vec![],
            reset_probability: 0.5,
            tolerance: 0.000001,
        }
    }
}

pub struct NTAResult {
    pub neighborhood: Vec<String>,
    pub scores: Vec<f64>,
    pub candidates: Vec<String>,
}

/// Uses random walk to calculate the neighborhood of a set of nodes
/// Returns [`Vec<String>`]representing the nodes in the neighborhood
///
/// # Parameters
/// - `config` - A [`NTAOptions`] struct containing the edge list, seeds, neighborhood size, reset probability, and tolerance
pub fn nta(config: NTAConfig) -> Vec<(String, f64)> {
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
        config.reset_probability,
    );
    let mut walk = walk_res.iter().enumerate().collect::<Vec<(usize, &f64)>>();
    walk.sort_by(|a, b| b.1.partial_cmp(a.1).unwrap());
    walk.iter()
        .map(|(i, p)| (reverse_map.get(&i).unwrap().clone(), **p))
        .collect()
}

fn random_walk_probability(
    adj_matrix: &ndarray::Array2<f64>,
    node_indices: &Vec<usize>,
    r: f64,
    tolerance: f64,
) -> ndarray::Array1<f64> {
    let num_nodes = node_indices.len() as f64;
    let de = adj_matrix.sum_axis(Axis(0));
    // de to 2d array
    let de = de.insert_axis(Axis(1));
    let temp = adj_matrix.t().div(de);
    let w = temp.t();
    let mut p0 = ndarray::Array1::from_elem(w.shape()[0], 0.0);
    for i in node_indices {
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
