<div align="center">
  <img src="test/logo/logo_phi.png" alt="PHI Logo" width="200">
</div>

## <div align="center"><span style="color:red;"><b>PHI</b></span> (<span style="color:red;"><b>P</b></span>angenome-based <span style="color:red;"><b>H</b></span>aplotype <span style="color:red;"><b>I</b></span>nference)</div>


## <a name="started"></a>Getting Started

### Prerequisites

#### Before using PHI, ensure you have the following dependencies installed:

1. **GCC 9 or later** - [GCC](https://gcc.gnu.org/)
2. **Zlib** - [zlib](https://zlib.net/)
3. **Gurobi** - [Gurobi](https://www.gurobi.com)

### Get PHI

```bash
git clone https://github.com/gsc74/PHI
cd PHI
make GUROBI_HOME=/path/to/gurobo_home (i.e. /opt/gurobi1101/linux64)

# test run with ILP
./PHI -t32 -q0 -g test/MHC_4.gfa.gz -r test/CHM13_reads.fq.gz -o CHM13.fa

# test run with IQP (default)
./PHI -t32 -q1 -g test/MHC_4.gfa.gz -r test/CHM13_reads.fq.gz -o CHM13.fa

# test run with vcf and reference
./vcf2gfa.py -v test/MHC_4.vcf.gz -r test/MHC-CHM13.0.fa.gz | gzip > test/MHC_4_vcf.gfa.gz
./PHI -t32 -g test/MHC_4_vcf.gfa.gz -r test/CHM13_reads.fq.gz -o CHM13.fa
```

## Description
PHI is a tool designed to reconstruct haploid haplotypes from low-coverage short reads using a haplotype-aware pangenome graph represented as a Directed Acyclic Graph (DAG). Additionally, a VCF file and a reference genome against which the VCF was built can be used with the script `vcf2gfa.py` to generate a pangenome graph as input. PHI uses short reads to reconstruct haploid haplotypes through two methods:

1. **Integer Linear Programming (ILP)**: Enabled by passing the `-q0` flag. This uses an ILP-based formulation.
2. **Integer Quadratic Programming (IQP)**: Enabled by passing the `-q1` flag, this is the default method and uses an IQP-based formulation.

The details of these formulations are described in our [paper](#publications).

## Results
(Results will be added later)

## Future Work
1. We plan to add support for diploid haplotype reconstruction.

## Publications
(To be added)