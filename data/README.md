
# Reproduce Results

## 1. Prerequisites

Before starting the process, ensure the following tools and datasets are available:

- **Miniforge**: [Miniforge](https://github.com/conda-forge/miniforge)
- **Required Tools**:
    - edlib-aligner
    - PHI
    - VG
    - PanGenie
    - seqkit
    - bcftools (with plugins)
    - snakemake
    - vcflib
    - Gurobi
    - gfa2gbwt

Ensure all dependencies are installed or available via conda.

---

## 2. Create Conda Environments

To manage dependencies, create separate conda environments for each tool. Here are example commands:

```bash
# Create environment and vcflib
conda create -n vcflib bioconda::vcflib -y

# Create environment for snakemake
conda create -n snakemake -c bioconda -c conda-forge snakemake=8.20.3 -y

# Install other dependencies for postprocessing
conda create -n edlib
conda activate edlib
conda install anconda::pip -y
pip3 install biopython
pip3 install edlib
```

---

## 3. Preprocessing

### Download binary of gfa2gbwt, seqkit and PanGenie as well as export path to .bashrc
```bash
cd ..
./Installdeps
# Add extra/bin and extra/lib to .bashrc
echo 'export PATH="$(pwd)/extra/bin:$PATH"' >> ~/.bashrc
echo 'export LD_LIBRARY_PATH="$(pwd)/extra/lib:$LD_LIBRARY_PATH"' >> ~/.bashrc
# export bcftools plugin
echo 'export BCFTOOLS_PLUGINS=$(pwd)/extra/plugins' >> ~/.bashrc
source ~/.bashrc
cd data

# Install seqkit
wget https://github.com/shenwei356/seqkit/releases/download/v2.8.2/seqkit_linux_amd64.tar.gz
tar -xvf seqkit_linux_amd64.tar.gz
mv seqkit ../extra/bin

# Install PanGenie
git clone https://github.com/eblerjana/pangenie.git  
cd pangenie  
conda env create -f environment.yml  
conda activate pangenie
conda install conda-forge::cereal -y
mkdir build; cd build; cmake .. ; make -j4
cp src/PanGenie ../../../extra/bin
cp src/PanGenie-index ../../../extra/bin
cp src/libPanGenieLib.so ../../../extra/lib
cd ../..
rm -rf pangenie
```
Note: Please remove the paths from `.bashrc` after the reproducing of the results are finished.

### Download MHC haplotypes, build graph and subsample reads

```bash
python3 preprocess.py -t16
cd Ground_truth
gunzip *.gz
cd ..
```

---

## 4. Running ILP and IQP

For ILP and IQP, execute the following scripts:

```bash
# ILP execution (-b is batch size)
python3 run_batch_5_milp.py -b 2
python3 run_batch_6_milp.py -b 2


# ILP execution (no relaxation)
python3 run_batch_5.py -b 2
python3 run_batch_6.py -b 2

# IQP execution (-b is batch size)
python3 run_batch_3_miqp.py -b 2
python3 run_batch_4_miqp.py -b 2

# IQP execution (no relaxation)
python3 run_batch_3.py -b 2
python3 run_batch_4.py -b 2
```

---

## 5. Running Progressive Imputation

For progressive imputation, execute the following scripts:

```bash
python3 run_batch_9.py -b 2 # 13 haps
python3 run_batch_10.py -b 2 # 25 haps
python3 run_batch_11.py -b 2 # 49 haps
python3 run_batch_12.py -b 2 # 7 haps
python3 run_batch_13.py -b 2 # 3 haps
```

---

## 6. Run VG and PanGenie

```bash
# Running VG (-b is batch size)
python3 run_VG.py -b 2

# Running PanGenie
python3 run_PG.py -b 2
```

---

## 7. Postprocessing

After the runs are completed, extract relevant metrics:

```bash
# ILP 
python3 postprocessing_3_MILP.py 

# ILP (no relaxation)
python3 postprocessing_3.py 

# IQP 
python3 postprocessing_2_MIQP.py 

# IQP (no relaxation)
python3 postprocessing_2.py

# Progressive imputation
python3 postprocessing_9.py # 13 haps
python3 postprocessing_10.py # 25 haps
python3 postprocessing_11.py # 49 haps
python3 postprocessing_12.py # 7 haps
python3 postprocessing_13.py # 3 haps

# VG
python3 postprocessing_VG.py

# PanGenie
python3 postprocessing_PG.py
```
---