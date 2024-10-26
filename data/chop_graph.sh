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

cd ../../ # data folder

source ~/.bashrc

VCF=data/MHC_49-MC.vcf
REF=data/hprc_haps/MHC-CHM13.0.fa
awk 'BEGIN {OFS="\t"} {if ($1 == "0") $1 = "CHM13#0"; print $0}' $VCF | bgzip  > MHC_49-MC_.vcf.gz && tabix -f -p vcf MHC_49-MC_.vcf.gz
awk '/^>/ {$0=">CHM13#0"} {print}'  $REF > REF.fasta
samtools faidx REF.fasta

vg construct -v MHC_49-MC_.vcf.gz -r REF.fasta -aS -p > X.xg
vg gbwt -x X.xg -v MHC_49-MC_.vcf.gz -o X.gbwt
vg gbwt -x X.xg -E -o X.paths.gbwt
vg gbwt -m X.gbwt X.paths.gbwt -o X.combined.gbwt
vg gbwt -x X.xg X.combined.gbwt --gbz-format -g X.gbz
ignore_genomes=$(python get_ids.py 23) # 23/24 samples are ignored, 1 sample is used
ignore_genomes_2=$(python get_ids_2.py "$ignore_genomes" 21) # 21/23 samples are ignored, 1 + 2 samples are used
ignore_genomes_3=$(python get_ids_2.py "$ignore_genomes_2" 18) # 18/21 samples are ignored, 1 + 2 + 3 samples are used
ignore_genomes_4=$(python get_ids_2.py "$ignore_genomes_3" 12) # 12/18 samples are ignored, 1 + 2 + 3 + 6 samples are used
vg gbwt -x X.xg X.combined.gbwt ${ignore_genomes} --gbz-format -g X3.gbz
ignore_genomes=$(python get_ids.py 21)
vg gbwt -x X.xg X.combined.gbwt ${ignore_genomes_2} --gbz-format -g X7.gbz
ignore_genomes=$(python get_ids.py 18)
vg gbwt -x X.xg X.combined.gbwt ${ignore_genomes_3} --gbz-format -g X13.gbz
ignore_genomes=$(python get_ids.py 12)
vg gbwt -x X.xg X.combined.gbwt ${ignore_genomes_4} --gbz-format -g X25.gbz
gfa2gbwt -d X -p -m 30 && mv X.gfa data/MHC_49-MC_30_2.gfa
gfa2gbwt -d X3 -p -m 30 && mv X3.gfa data/MHC_3-MC_30_2.gfa
gfa2gbwt -d X7 -p -m 30 && mv X7.gfa data/MHC_7-MC_30_2.gfa
gfa2gbwt -d X13 -p -m 30 && mv X13.gfa data/MHC_13-MC_30_2.gfa
gfa2gbwt -d X25 -p -m 30 && mv X25.gfa data/MHC_25-MC_30_2.gfa
cp X.gbz data/MHC_49-MC_30.gbz

rm REF.fasta*
rm MHC_49-MC_.vcf.gz*
rm X.*
rm X3.*
rm X7.*
rm X13.*
rm X25.*