#!/usr/bin/bash

hal2vg data/MHC-49_MC_out/*.hal --hdf5InMemory --chop 30 --progress --ignoreGenomes _MINIGRAPH_,Anc0 > MHC_49-MC_30.pg && vg ids --sort MHC_49-MC_30.pg > temp && mv temp MHC_49-MC_30.pg && vg convert -f MHC_49-MC_30.pg > data/MHC_49-MC_30.gfa
rm MHC_49-MC_30.pg

# Prepare-vcf-MC
git clone https://github.com/eblerjana/genotyping-pipelines.git
cd genotyping-pipelines
cd prepare-vcf-MC
cp ../../config.yaml .
cp ../../data/MHC-49_MC_out/MHC-49-MC.raw.vcf.gz .
cp ../../data/MHC-49_MC_out/MHC-49-MC.raw.vcf.gz.tbi .
cp ../../data/MHC-49_MC_out/MHC-49-MC.gfa.gz .
gunzip MHC-49-MC.gfa.gz
cat MHC-49-MC.gfa | grep "W" | awk 'BEGIN {OFS="\t"}; {print $2, $3}' > sample-info.tsv # get sample-info.tsv

source ${HOME}/.bashrc
conda activate snakemake

# use snakemake to prepare the VCF
snakemake -j 32

# get filtered VCF
cp results/vcf/MC/MC_filtered_ids.vcf .

# VCF file needs to be transformed such that no variants overlap each other using vcfbub
vcfbub -l 0 -r 100000 -i MC_filtered_ids.vcf > MHC_49-MC.vcf
# normalize the VCF file
bcftools norm -m -any MHC_49-MC.vcf | bgzip > MHC_49-MC_norm.vcf.gz && tabix -p vcf MHC_49-MC_norm.vcf.gz

# copy VCF files to the data folder
cp MHC_49-MC_norm.vcf.gz ../../data/
cp MHC_49-MC_norm.vcf.gz.tbi ../../data/
cp MHC_49-MC.vcf ../../data/