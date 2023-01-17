# UMD database converter

A small tool to reformat gene related mutation data in http://www.umd.be/ .
The original data are in .html table format. This tool can convert it to .csv file, then parse variants from cDNA representation to gene loci representation.

## Files
- `refGene.txt.gz`: Uncompress it to get default gene annotation file. 
- `umd_converter.py`: source code
- `analysis.ipynb`: code in jupyter notebook format


## Requirement
- reference genome: hg38.fa reference sequence(not included)
- environment:
    (1)hgvs: https://github.com/counsyl/hgvs.git
    (2)see `requirement.txt`

## Usage

`python umd_converter.py -i [input] -t [transcript] -r hg38.fa -g ./refGene.txt -s -c`

>example:
`python umd_converter.py -i D:/CODE/BioData/UMD_gene/DMD/DMD.html -t NM_004006 -r D:/CODE/BioData/hg38/hg38.fa -g D:/CODE/hgvs/pyhgvs/data/hg38/refGene.txt -s -c`

check `python umd_converter.py -h` for more info

