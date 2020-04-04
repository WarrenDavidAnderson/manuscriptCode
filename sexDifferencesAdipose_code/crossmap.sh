
cd /nv/vol192/civeleklab/warren/adipose_project/eqtlLD

ml anaconda/5.2.0-py3.6
source activate crossmap

python CrossMap.py bed hg38ToHg19.over.chain.gz bedhg38.bed bedhg19.bed

# To exit the crossmap environment, simply use the following command:
source deactivate

If you wish to install CrossMap for Python 2, you can use the following:
ml anaconda/5.2.0-py2.7
conda create -n crossmap-py2.7 python=2.7
source activate crossmap-py2.7
pip install CrossMap

To activate the new environment:
ml anaconda/5.2.0-py2.7
source activate crossmap-py2.7

