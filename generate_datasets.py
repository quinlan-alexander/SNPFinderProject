#!/usr/bin/env python3
# Name: Quinlan Alexander (qalexand)
# Group Members: None

import os
import arv
import pandas as pd
from find_snp import SNPFinder
from spinner import Spinner

########################################################################
# CommandLine
########################################################################
class CommandLine():
    '''
    Handle the command line, usage and help requests.
    '''

    def __init__(self, inOpts=None):
        '''
        Implement a parser to interpret the command line argv string using argparse.
        '''

        import argparse
        self.parser = argparse.ArgumentParser(
            description='Generate population datasets from 23andme raw text files',
            epilog='',
            add_help=True,  # default is True
            prefix_chars='-',
            usage='%(prog)s [options] -option1[default] <input >output'
            )
        self.parser.add_argument('inFileDir', action='store', help='23ndMe raw data file directory')
        self.parser.add_argument('outFile', action='store', help='output file name of CSV dataset containing all results')
        self.parser.add_argument('-g', '--genes', action='store', default=['MTHFR'], help='genes to query for SNPs')
        self.parser.add_argument('-v', '--version', action='version', version='%(prog)s 0.1')
        if inOpts is None:
            self.args = self.parser.parse_args()
        else:
            self.args = self.parser.parse_args(inOpts)

########################################################################
# Main
# Here is the main program
#
#
########################################################################


def main(inCL=None):
    '''
    Generate some data.
    '''
    if inCL is None:
        cmd_line = CommandLine()
    else:
        cmd_line = CommandLine(inCL)
    print(cmd_line.args)

    directory = os.fsencode(cmd_line.args.inFileDir)

    all_snps_results_df = pd.DataFrame()

    spinner = Spinner()
    spinner.start()

    for file in os.listdir(directory):
        filename = os.fsdecode(file)
        if filename.endswith(".txt"):
            filepath = os.path.join(str(directory.decode()), filename)
            genome = arv.load(filepath)
            genes = list(cmd_line.args.genes.split(','))
            snpFinder = SNPFinder(filename,genome, genes)
            snpFinder.find_snp_by_genes()
            all_snps_results_df = all_snps_results_df.append(snpFinder.relevant_SNPs_df)

    all_snps_results_df.to_csv(cmd_line.args.outFile, index=False)

    spinner.stop()

    samples_df = pd.read_csv('/Users/ba25714/PycharmProjects/final_project/sample_results.csv')

    score_counts_by_gene_df = pd.DataFrame(columns=['gene', 'warning', 'caution', 'ok'])

    for gene in samples_df['gene'].unique():
        warn = samples_df.loc[(samples_df['gene'] == 'gene') & samples_df['score'] <= .05].agg(['count'])['score']
        caution = samples_df.loc[(samples_df['gene'] == 'gene') & (samples_df['score'] > .05) & (samples_df['score'] <= 0.5)].agg(['count'])['score']
        ok = samples_df.loc[(samples_df['gene'] == 'gene') & (samples_df['score'] > .50) & (samples_df['score'] <= 1.0)].agg(['count'])['score']
        score_counts_by_gene_df.append({'gene': gene, 'warning': warn, 'caution': caution, 'ok': ok}, ignore_index=True)

if __name__ == "__main__":
    main()