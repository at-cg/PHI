# <div align="center">AlphaASM</div>

## <div align="center">Haplotype reconstruction from low coverage reads with Pangenome Graphs</div>

## <a name="started"></a>Getting Started

### Get Minichain

```bash
git clone https://github.com/gsc74/AlphaASM
cd AlphaASM && make

# test run
./AlphaASM -t4 -g test/MHC-CHM13.0.gfa.gz -r test/CHM13_3X_HiFi.fq.gz -o out_hap.fa

```