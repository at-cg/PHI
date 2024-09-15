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

# test run
./PHI -t32 -g test/MHC_4.gfa.gz -r test/CHM13_reads.fq -o CHM13.fa
```

## Description
PHI is a tool designed to reconstruct haploid haplotypes from low coverage short-reads.

## Benchmark
(To be added).