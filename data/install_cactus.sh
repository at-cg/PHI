#!/usr/bin/bash

echo "Installing Cactus v2.9.0..."

# This script installs the Cactus framework on the local machine.
wget -O data/cactus-bin-v2.9.0.tar.gz https://github.com/ComparativeGenomicsToolkit/cactus/releases/download/v2.9.0/cactus-bin-v2.9.0.tar.gz
cd data
tar -xzf cactus-bin-v2.9.0.tar.gz
cd cactus-bin-v2.9.0

# create venv
virtualenv -p python3 venv-cactus-v2.9.0
printf "export PATH=$(pwd)/bin:\$PATH\nexport PYTHONPATH=$(pwd)/lib:\$PYTHONPATH\n" >> venv-cactus-v2.9.0/bin/activate
source venv-cactus-v2.9.0/bin/activate
python3 -m pip install -U setuptools pip wheel
python3 -m pip install -U .
python3 -m pip install -U -r ./toil-requirement.txt

# install additional dependencies
cd bin && for i in wigToBigWig faToTwoBit bedToBigBed bigBedToBed axtChain pslPosTarget bedSort hgGcPercent mafToBigMaf hgLoadMafSummary hgLoadChain; do wget -q http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/${i}; chmod +x ${i}; done
conda activate vcflib
ln -s $(which vcfwave) vcfwave
cd ..
# install pangenometools
mkdir -p build-tools
cp ../../downloadPangenomeTools build-tools/downloadPangenomeTools
sh build-tools/downloadPangenomeTools

echo "Cactus v2.9.0 installed successfully."