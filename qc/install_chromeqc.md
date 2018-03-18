# Install ChromeQC

git clone https://github.com/hackseq/2017_project_6 ChromeQC
cd ChromeQC
wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh -O miniconda.sh
bash miniconda.sh -b -p ./miniconda
export PATH=$(pwd)/miniconda/bin:$PATH
conda create -n chromeqc -c bioconda -c conda-forge minimap2 pysam samtools matplotlib matplotlib -y
source activate chromeqc
