# A Gene Regulatory Network Approach to Cellular Persistence in S. cerevisiae

This repository contains the source code and data for the Master's Thesis titled "A gene regulatory network approach to cellular persistence in *S. cerevisiae*".

## Abstract

Cellular persistence is a phenomenon where a subpopulation of cells transiently tolerates lethal drug doses, posing a major clinical challenge. This non-heritable form of drug tolerance is thought to arise from the stochastic nature of gene expression. We hypothesized that expression noise drives key regulators across critical thresholds, leading cells into a drug-tolerant state. To test this hypothesis, we employed an *in silico* approach. First, persistence-candidate genes were selected from single-cell proteomics data. Next, we analyzed a state-of-the-art *S. cerevisiae* Gene Regulatory Network (GRN), identifying its key architectural features. Finally, the network's dynamics were simulated using a Boolean model under both proliferative and quiescent conditions. The system consistently converged to a single, condition-specific attractor, showing no evidence of stable phenotypic heterogeneity. This result refutes the idea that persistence corresponds to an alternative stable state in this framework. However, we observed complex and extended transient dynamics en route to the quiescent-state attractor. This suggests that persistence may not be a stable attractor but rather a long-lived, noise-driven, non-equilibrium state.

## Project Description

This project is divided into two main analytical parts:

1.  **Topological Network Analysis**: This involves the structural characterization of the *S. cerevisiae* Gene Regulatory Network (GRN) to understand the framework that governs information flow and cellular dynamics. We analyze local properties (e.g., degree distribution, clustering) and global features (e.g., bow-tie architecture, functional communities).

2.  **Dynamical Simulation**: This part implements a threshold Boolean model to simulate the network's behavior under different environmental conditions (rich media vs. starvation). The goal is to investigate whether the model can generate the phenotypic heterogeneity required for persistence by searching for multiple stable states (attractors).

## Repository Structure

The project is organized as follows:

```
.
├── data/
│   ├── S_cerevisiae_GRN.gml
│   └── persistence_candidates.csv
├── notebooks/
│   ├── 01_Topological_Analysis.ipynb
│   └── 02_Dynamical_Simulation.ipynb
├── figures/
│   └── (folder for generated figures)
├── requirements.txt
└── README.md
```

## Installation

To run the analyses, you will need Python 3.8 or higher. Using a virtual environment is recommended.

1.  Clone the repository:
    ```bash
    git clone https://github.com/jaimemartinezcazon/S.-Cerevisae-Gene-Regulatory-Network-as-an-approach-to-cellular-persistence
    cd S.-Cerevisae-Gene-Regulatory-Network-as-an-approach-to-cellular-persistence
    ```

2.  Create and activate a virtual environment:
    ```bash
    python -m venv venv
    source venv/bin/activate  # On Windows, use: venv\Scripts\activate
    ```

3.  Install the required dependencies:
    ```bash
    pip install -r requirements.txt
    ```

## Usage

The analyses are implemented in Jupyter Notebooks for clarity and ease of execution.

### 1. Topological Analysis

Open and run the notebook `notebooks/01_Topological_Analysis.ipynb`. This notebook performs the following steps:
-   Loads the network from the `data/` directory.
-   Analyzes local network properties such as degree distribution, clustering, and assortativity.
-   Analyzes global network properties, including the bow-tie decomposition and community detection.
-   Generates figures summarizing the network's topology.

### 2. Dynamical Simulation

Open and run the notebook `notebooks/02_Dynamical_Simulation.ipynb`. This notebook:
-   Implements the Boolean network model with a majority-rule update scheme and asynchronous updates.
-   Sets up environmental conditions by fixing the states of input transcription factors.
-   Runs thousands of simulations from random initial states to explore the state space.
-   Analyzes the resulting attractors and transient dynamics to identify network behavior.

## Data Sources

-   **Gene Regulatory Network (GRN):** The network is a state-of-the-art reconstruction for *S. cerevisiae* from [Jackson et al. (2020), *eLife*](https://elifesciences.org/articles/51254).
-   **Candidate Gene Proteomics Data:** Candidate genes were selected from single-cell proteomics data available in the [YeastRGB](https://www.nature.com/articles/nchembio.2147) and [LoQate](https://www.jcb.org/content/200/6/839) databases.

## Key Findings

1.  **Single Attractor per Condition:** The deterministic Boolean model converges to a single, globally dominant attractor for each simulated environment. It does not produce the stable phenotypic heterogeneity required to explain persistence as an alternative cell state.
2.  **The Importance of Transients:** The convergence to the quiescent attractor is significantly slower and more heterogeneous than the convergence to the proliferative state. This suggests that persistence could be a long-lived transient state within a rugged dynamical landscape, rather than a fixed attractor.

## Citation

If you use this code or data in your research, please cite the original thesis:

> Martínez, J. (2025). *A gene regulatory network approach to cellular persistence in S. cerevisiae*. Master's Thesis, Physics of Complex Systems and Biophysics, Universitat de Barcelona.

## License

This project is distributed under the MIT License.

## Contact

Jaime Martínez Cazón – `jmartica662@alumnes.ub.edu`
