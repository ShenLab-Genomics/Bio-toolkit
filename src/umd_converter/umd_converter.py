#!/usr/bin/env python
# coding: utf-8

import pandas as pd
import pyhgvs as hgvs
import pyhgvs.utils as hgvs_utils
import argparse
from pyfaidx import Fasta


class UMD_Converter:
    def __init__(self, transcript, genome, ref_gene) -> None:
        # Initialization
        self.genome = Fasta(genome)
        self.ts = transcript
        # Read RefSeq transcripts into a python dict.
        # Read genome sequence using pyfaidx.
        with open(ref_gene) as infile:
            self.transcripts = hgvs_utils.read_transcripts(infile)

        self.UMD_columns = ['Protein nomenclature', 'cDNA Nomenclature', 'Exon', 'Codon',
                            'Structure', 'HCD', 'Rearrangement', 'Mutation type',
                            'Mutational event', '# records']

    def task_convert(self, table, target_path, filt_snp=True):
        # Provide a callback for fetching a transcript by its name.
        def get_transcript(name):
            return self.transcripts.get(name)
        res_table = pd.DataFrame(
            columns=['chrom', 'offset', 'ref', 'alt'] + self.UMD_columns, index=None)
        cnt = 0
        for idx, row in table.iterrows():
            try:
                ts = self.ts
                cmark = row['cDNA Nomenclature']
                # print(ts + ':' + str(cmark))
                chrom, offset, ref, alt = hgvs.parse_hgvs_name(
                    ts + ':' + str(cmark), self.genome, get_transcript=get_transcript, normalize=True
                )
                # print(chrom,offset,ref,alt)
                # to make output file not too large
                if (len(ref) > 200) or (len(alt) > 200):
                    continue
                if FILT_SNP and ((len(ref) > 1) or (len(alt) > 1)):
                    continue
                if ref == alt:
                    continue
                feature_dict = {}
                for co in self.UMD_columns:
                    feature_dict[co] = row[co]
                res_table.loc[cnt] = pd.Series({'chrom': chrom, 'offset': offset,
                                                'ref': ref, 'alt': alt,
                                                **feature_dict})
                cnt += 1
            except Exception as e:
                print(e)
        return res_table


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='UMD data convert')
    # data source
    parser.add_argument('-i', '--input', type=str, required=True,
                        help='the .html file downloaded from the UMD database')
    # target transcript
    parser.add_argument('-t', '--ts', type=str, required=True,
                        help='the transcript number for the gene')
    # genome:
    parser.add_argument('-r', '--ref', type=str, required=True,
                        help='file path of reference genome (hg38)')
    # ref_gene:
    parser.add_argument('-g', '--gene', type=str, required=True,
                        help='gene transcript reference file. (default:refGene.txt)')
    # convert:
    parser.add_argument('-c', '--convert', action="store_true",
                        help='whether convert cDNA to geneloc')
    # snp:
    parser.add_argument('-s', '--snp',  action="store_true",
                        help='save snp only (abandon any mutation with ref/alt length more than 1)')
    # output:
    parser.add_argument('-o', '--output', type=str,
                        default='./result.csv', help='path for output file')

    args = parser.parse_args()
    html_file = args.input
    TRANSCRIPT = args.ts
    GENOME = args.ref
    REF_GENE = args.gene
    CONVERT = args.convert
    FILT_SNP = args.snp
    OUTPUT_PATH = args.output

    # Load data
    df = pd.read_html(html_file)
    df = pd.DataFrame(df[0])
    df = df.dropna(subset=['Protein nomenclature'], axis=0)
    columns = ['Protein nomenclature', 'cDNA Nomenclature', 'Exon', 'Codon',
               'Structure', 'HCD', 'Rearrangement', 'Mutation type',
               'Mutational event', '# records']
    df = df[columns]
    print('data size:{}'.format(len(df)))
    # process
    convert = UMD_Converter(TRANSCRIPT, GENOME, REF_GENE)
    if CONVERT:
        res = convert.task_convert(df, OUTPUT_PATH, filt_snp=FILT_SNP)
    else:
        res = df
    res.to_csv(OUTPUT_PATH, sep='\t')
    print('done')
