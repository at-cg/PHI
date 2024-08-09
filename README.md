# <div align="center">AlphaASM</div>

## <div align="center">Haplotype Reconstruction from Low Coverage Reads with Pangenome Graphs</div>

## <a name="started"></a>Getting Started

### Prerequisites

#### Before using AlphaASM, ensure you have the following dependencies installed:

1. **GCC 13 or later** - [GCC](https://gcc.gnu.org/)
2. **Zlib** - [zlib](https://zlib.net/)
3. **Jemalloc** - [Jemalloc](https://github.com/jemalloc/jemalloc)
4. **Gurobi** - [Gurobi](https://www.gurobi.com)

### Get AlphaASM

```bash
git clone https://github.com/gsc74/AlphaASM
cd AlphaASM
make GUROBI_HOME=/path/to/gurobo_home (i.e. /opt/gurobi1101/linux64)

# test run
./AlphaASM -t32 -g test/MHC-CHM13.0.gfa.gz -r test/CHM13_test.fa -o CHM13.fa
```

## Description
AlphaASM is a tool designed to reconstruct haploid haplotypes from low coverage short-reads.

## Benchmark
(To be added).