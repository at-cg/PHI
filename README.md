<div align="center">
  <img src="test/logo/logo_phi.png" alt="PHI Logo" width="200">
</div>

## <div align="center"><span style="color:red;"><b>PHI</b></span> (<span style="color:red;"><b>P</b></span>angenome-based <span style="color:red;"><b>H</b></span>aplotype <span style="color:red;"><b>I</b></span>nference)</div>


## <a name="started"></a>Getting Started

### Prerequisites

Before using PHI, ensure that **Miniforge** is installed: [Miniforge Installation Guide](https://github.com/conda-forge/miniforge).

## <a name="get_phi"></a>Get PHI

```bash
git clone https://github.com/at-cg/PHI
cd PHI
# Install dependencies (Miniforge is required)
./Installdeps
export PATH="$(pwd)/extra/bin:$PATH"
export LD_LIBRARY_PATH="$(pwd)/extra/lib:$LD_LIBRARY_PATH"
make

# test run with IQP (default)
./PHI -t32 -g test/MHC_4.gfa.gz -r test/CHM13_reads.fq.gz -o CHM13.fa

# test run with ILP
./PHI -t32 -q0 -g test/MHC_4.gfa.gz -r test/CHM13_reads.fq.gz -o CHM13.fa

# test run with vcf and reference
./vcf2gfa.py -v test/MHC_4.vcf.gz -r test/MHC-CHM13.0.fa.gz | bgzip > test/MHC_4_vcf.gfa.gz
./PHI -t32 -g test/MHC_4_vcf.gfa.gz -r test/CHM13_reads.fq.gz -o CHM13.fa
```

#### Adding Binary and Library Paths to `.bashrc`
To ensure that the `extra/bin` and `extra/lib` directories are automatically loaded for every terminal session, you can export them to your `~/.bashrc`. This will make sure the required binaries and libraries for `PHI` are available.

```bash
# Add extra/bin and extra/lib to .bashrc
echo 'export PATH="$(pwd)/extra/bin:$PATH"' >> ~/.bashrc
echo 'export LD_LIBRARY_PATH="$(pwd)/extra/lib:$LD_LIBRARY_PATH"' >> ~/.bashrc
source ~/.bashrc
```

## Table of Contents

- [Getting Started](#started)
- [Get PHI](#get_phi)
- [Introduction](#intro)
- [Results](#results)
- [Future work](#future)
- [Publications](#pub)

## <a name="intro"></a>Introduction
PHI reconstructs haploid haplotypes from low-coverage sequencing data (short-reads or long-reads) using a pangenome graph represented as a Directed Acyclic Graph (DAG). Users can provide a pangenome reference as either:
- Graph Format ([GFA v1.1](http://gfa-spec.github.io/GFA-spec/GFA1.html#gfa-11)): A sequence graph-based representation of the pangenome graph.
- Variant Call Format ([VCF](https://samtools.github.io/hts-specs/VCFv4.2.pdf)): A list of multi-sample, multi-allelic phased variants along with a reference genome.

PHI outputs the haplotype sequence associated with the optimal inferred path from the graph in FASTA format.

PHI employs two integer programming formulations, ILP and IQP, and uses the [Gurobi optimizer](https://www.gurobi.com) to solve these formulations optimally. Details of these formulations are described in our [paper](#publications).


## <a name="results"></a>Results
We benchmarked PHI (v1.0) using real Illumina short-reads from five MHC haplotypes: APD, DBB, MANN, QBL, and SSTO from [Houwaart et al. (2022)](https://doi.org/10.1111/tan.15020). These datasets were downsampled to various coverages ranging from 0.1x to 18.2x. Using [Minigraph-Cactus](https://github.com/ComparativeGenomicsToolkit/cactus/tree/master), we constructed a haplotype-resolved pangenome graph comprising 49 complete [MHC sequences](https://doi.org/10.5281/zenodo.6617246). To assess the accuracy of PHI, we evaluated the edit distance between the imputed haplotypes and the ground-truth haplotypes. Additionally, we benchmarked the runtime and memory usage of the ILP and IQP-based formulations. The results show that the IQP method runs faster than the ILP method but requires approximately 1.5x more memory.

<p align="center" id="F1-score">
    <img src="data/edit_distances.jpg" width="700" alt="F1-score"/>
</p>

> Edit distance between ground-truth and imputed MHC haplotypes generated from real Illumina reads at various sequencing coverages (0.1× to full coverage) using different tools (PHI, VG, and PanGenie).

<p align="center" id="F1-score">
    <img src="data/phi_vs_phi_ilp.jpg" width="700" alt="F1-score"/>
</p>

> Performance comparison between ILP and IQP, illustrating runtime and memory usage across different coverage levels and haplotypes.

In our experiments, we used the [Gurobi Academic WLS License](https://www.gurobi.com/academia/academic-program-and-licenses/). The scripts to reproduce the results are available [here](data).


## <a name="future"></a>Future Work
- We plan to add support for diploid haplotype reconstruction.
- Improve scaling to larger pangenome graphs.


## <a name="pub"></a>Publications
- **Ghanshyam Chandra, Md Helal Hossen, Stephan Scholz, Alexander T Dilthey, Daniel Gibney and Chirag Jain**. "[Integer programming framework for pangenome-based genome inference](https://www.biorxiv.org/)". *bioRxiv* 2024.