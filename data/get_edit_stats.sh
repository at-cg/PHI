#!/usr/bin/bash

for file in data/hprc_haps/*.fa
do
    edlib-aligner Ground_truth/APD.fasta $file > ${file}_temp && awk 'NR==11 {print "'${file}': ", $2}' ${file}_temp && rm ${file}_temp &

done > APD.txt
wait

for file in data/hprc_haps/*.fa
do
    edlib-aligner Ground_truth/DBB.fasta $file > ${file}_temp && awk 'NR==11 {print "'${file}': ", $2}' ${file}_temp && rm ${file}_temp &

done > DBB.txt
wait

for file in data/hprc_haps/*.fa
do
    edlib-aligner Ground_truth/MANN.fasta $file > ${file}_temp && awk 'NR==11 {print "'${file}': ", $2}' ${file}_temp && rm ${file}_temp &

done > MANN.txt
wait

for file in data/hprc_haps/*.fa
do
    edlib-aligner Ground_truth/QBL.fasta $file > ${file}_temp && awk 'NR==11 {print "'${file}': ", $2}' ${file}_temp && rm ${file}_temp &

done > QBL.txt
wait

for file in data/hprc_haps/*.fa
do
    edlib-aligner Ground_truth/SSTO.fasta $file > ${file}_temp && awk 'NR==11 {print "'${file}': ", $2}' ${file}_temp && rm ${file}_temp &

done > SSTO.txt
wait

awk '{x=$NF; sum+=x; if(NR==1||x<min) min=x; if(x>max) max=x} END {print "Mean: " sum/NR, "Max: " max, "Min: " min}' APD.txt
awk '{x=$NF; sum+=x; if(NR==1||x<min) min=x; if(x>max) max=x} END {print "Mean: " sum/NR, "Max: " max, "Min: " min}' DBB.txt
awk '{x=$NF; sum+=x; if(NR==1||x<min) min=x; if(x>max) max=x} END {print "Mean: " sum/NR, "Max: " max, "Min: " min}' MANN.txt
awk '{x=$NF; sum+=x; if(NR==1||x<min) min=x; if(x>max) max=x} END {print "Mean: " sum/NR, "Max: " max, "Min: " min}' QBL.txt
awk '{x=$NF; sum+=x; if(NR==1||x<min) min=x; if(x>max) max=x} END {print "Mean: " sum/NR, "Max: " max, "Min: " min}' SSTO.txt
 
