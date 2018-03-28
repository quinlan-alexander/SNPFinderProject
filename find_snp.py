#!/usr/bin/env python3
# Name: Quinlan Alexander (qalexand)
# Group Members: None

import arv
import requests
from bs4 import BeautifulSoup
import re
import json
import pandas as pd
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
            description='Find SNPs for given 23andMe raw data',
            epilog='',
            add_help=True,  # default is True
            prefix_chars='-',
            usage='%(prog)s [options] -option1[default] <input >output'
            )
        self.parser.add_argument('inFile', action='store', help='23ndMe raw data')
        self.parser.add_argument('outFile', action='store', help='save results to output file')
        self.parser.add_argument('-g', '--genes', action='store', default=['MTHFR'], help='genes to query for SNPs')
        self.parser.add_argument('-v', '--version', action='version', version='%(prog)s 0.1')
        self.parser.add_argument('--generateCSV', action='store_true', help='Output results to CSV file')
        self.parser.add_argument('--color', action='store_true', help='Output results to standard output using color')
        if inOpts is None:
            self.args = self.parser.parse_args()
        else:
            self.args = self.parser.parse_args(inOpts)

class Twenty3andMeData():
    '''Twenty3andMeData is a class to access 23andMe API endpoints'''
    def __init__(self,gene):
        r = requests.get('https://api.23andme.com/3/marker/?gene_name='+gene)
        self.gene_dict = json.loads(r.text)
        self.snp_list = []
        for data in self.gene_dict['data']:
            if not 'i' in data['id']:
                self.snp_list.append(data['id'])

    def get_gene(self):
        return self.gene_dict

    def get_gene_name(self):
        return self.gene_dict['data'][0]['gene_names'][1]

    def get_snp_collection(self):
        return self.snp_list

class MyVariantInfo():

    def __init__(self):
        self.base_complement = {'A': 'T', 'C': 'G', 'T': 'A', 'G': 'C'}
        self.status_code = 404

    def complement(self, sequence):
        '''return complement of alleles'''
        complement = ""
        for c in sequence:
            if c in self.base_complement:
                complement += self.base_complement[c]
            else:
                complement += c
        return complement

    def compute_score(self):
        '''computes a score based on several cached online sources'''
        self.sift_score = 1  # default to no consequence if no SIFT score
        # check dbnsfp first for SIFT score
        if 'dbnsfp' in self.variant_dict['hits'][0]:
            if 'sift' in self.variant_dict['hits'][0]['dbnsfp'].keys():
                if 'converted_rankscore' in self.variant_dict['hits'][0]['dbnsfp']['sift'].keys():
                    self.sift_score = self.variant_dict['hits'][0]['dbnsfp']['sift']['converted_rankscore']
                    return
        # check cadd for SIFT score
        if 'cadd' in self.variant_dict['hits'][0]:
            if 'sift' in self.variant_dict['hits'][0]['cadd'].keys():
                self.sift_score = self.variant_dict['hits'][0]['cadd']['sift']['val']

    def get_data(self, rsid):
        r = requests.get('http://myvariant.info/v1/query?q=' + rsid)
        self.status_code = r.status_code
        if r.status_code == 404:
            print('ERROR: SNP {} does not exist on myvariant.info.'.format(rsid))
            return
        self.variant_dict = json.loads(r.text)
        if len(self.variant_dict['hits']):
            self.common_allele = self.variant_dict['hits'][0]['vcf']['ref']
            self.alt_allele = self.variant_dict['hits'][0]['vcf']['alt']
            self.sift_score = 1  # unknown consequence if no SIFT score
            # if rsid == 'rs4846051':
            #     print('foo')
            self.compute_score()
            self.orientation = 'plus'

    def hasInfo(self):
        return len(self.variant_dict['hits']) is not 0

    def get_common_allele(self, request_orientation=None):
        if request_orientation == None or request_orientation == self.orientation:
            return self.common_allele
        return self.complement(self.common_allele)

    def get_alt_allele(self, request_orientation=None):
        if request_orientation == None or request_orientation == self.orientation:
            return self.alt_allele
        return self.complement(self.alt_allele)

    def get_sift_score(self):
        return self.sift_score

class SNPediaData():

    def __init__(self):
        self.base_complement = {'A': 'T', 'C': 'G', 'T': 'A', 'G': 'C'}

    def complement(self, sequence):
        '''return complement of alleles'''
        complement = ""
        for c in sequence:
            if c in self.base_complement:
                complement += self.base_complement[c]
            else:
                complement += c
        return complement

    def get_data(self, rsid):
        r = requests.get('https://www.snpedia.com/index.php/'+rsid)
        self.status_code = r.status_code
        if r.status_code == 404:
            print('ERROR: SNP {} does not exist on snpedia.'.format(rsid))
            return
        soup = BeautifulSoup(r.text, "html5lib")
        trows = soup('table')

        # table element 1 contains orientation
        p = re.compile('(<td>)(.+)(<\/td>)')
        orientation = str(trows[1])
        orientationMatch = p.findall(orientation)
        if orientationMatch:
            self.orientation = orientationMatch[0][1]
        else:
            self.orientation = 'Error: orientation for SNP '+rsid+' not found'

        # table element 3 contains SNP allele details
        allele_data = soup.find_all("table", class_="sortable smwtable")
        if len(allele_data) == 0:
            return
        p = re.compile('(\([ACGT];[ACGT]\))')
        alleleMatch = p.findall(allele_data[0].text)
        if alleleMatch:
            self.common_allele = alleleMatch[0]
        else:
            self.common_allele = 'Error: common allele for SNP '+rsid+' not found'

    def get_common_allele(self, request_orientation=None):
        if request_orientation == None or request_orientation == self.orientation:
            return self.common_allele
        return self.complement(self.common_allele)

class bcolors:
    OKGREEN = '\033[92m'
    CAUTION = '\033[93m'
    WARNING = '\033[91m'
    ENDC = '\033[0m'

class SNPFinder():

    def __init__(self, subject, genome, genes):
        self.subject = subject
        self.genome = genome
        self.genes = genes
        self.affected = []

    def subject_ethnicity(self):
        if self.genome['rs1426654'].genotype == 'AA':
            return 'light-skinned, European'
        elif self.genome['rs1426654'].genotype == 'AG':
            return 'mixed-skinned, European + (African or Asian)'
        elif self.genome['rs1426654'].genotype == 'GG':
            return 'darker-skinned, Asian or African'

    def find_snp_by_genes(self):
        relevant_SNPs_df = pd.DataFrame(columns=['subject', 'gene', 'rsId', 'score', 'gender', 'ethnicity'])
        for gene in self.genes:
            #print(gene)
            # look up gene on 23andme to get list of SNP rsid
            t3andMe = Twenty3andMeData(gene)
            geneDict = t3andMe.get_gene()
            # for SNP in rsid collection
            for snp in t3andMe.get_snp_collection():
            #   get SNP data from MyVariantInfo (common allele, orientation)
                # print(snp)
                snpData = MyVariantInfo()
                snpData.get_data(snp)
                if snpData.status_code != 404 and snpData.hasInfo():
                    # compare common alleles from MyVariantInfo and 23andme
                    if snp in self.genome:
                        # if snp == 'rs4846051':
                            # print
                            # print('MyVariation common allele genotype {}{}'.format(snpData.common_allele,snpData.common_allele))
                            # print('23andMe genotype {}'.format(self.genome[snp].genotype))
                            # print

                        if snpData.common_allele+snpData.common_allele != self.genome[snp].genotype:
                            # add SNP and score to genome's affected collection
                            #print('{} SNP SIFT score: {}'.format(snp, snpData.get_sift_score()))
                            #self.affected.append(snp
                            gender = "male" if self.genome.y_chromosome else "female"
                            relevant_SNPs_df = relevant_SNPs_df.append({'subject': self.subject,
                                                                        'gene': gene,
                                                                        'rsId': snp,
                                                                        'score': snpData.get_sift_score(),
                                                                        'gender': "male" if self.genome.y_chromosome else "female",
                                                                        'ethnicity' : self.subject_ethnicity()},
                                                                       ignore_index=True)

        self.relevant_SNPs_df = relevant_SNPs_df.sort_values(['score'])

    def print_human_report(self, out_writer):
        out_writer.write('Relevant SNPs for {}\n'.format(self.subject))
        for row in self.relevant_SNPs_df.itertuples():
            out_writer.write('{} {} {}\n'.format(row.rsId,row.gene, row.score))

    def print_color_report(self):
        print('Relevant SNPs for {}'.format(self.subject))
        for row in self.relevant_SNPs_df.itertuples():
            if row.score <= 0.05:
                color = bcolors.WARNING
            elif row.score <= 0.50:
                color = bcolors.CAUTION
            else:
                color = bcolors.OKGREEN
            print(color+'{} {} {}'.format(row.rsId, row.gene, row.score)+bcolors.ENDC)



########################################################################
# Main
# Here is the main program
#
#
########################################################################


def main(inCL=None):
    '''
    Find some SNPs.
    '''
    if inCL is None:
        cmd_line = CommandLine()
    else:
        cmd_line = CommandLine(inCL)
    print(cmd_line.args)

    # read input 23andme file
    if cmd_line.args.inFile != None:
        genome = arv.load(cmd_line.args.inFile)
    else: # read stdin
        cmd_line.args.inFile('')

    genes = list(cmd_line.args.genes.split(','))

    snpFinder = SNPFinder(cmd_line.args.inFile,genome,genes)

    spinner = Spinner()
    spinner.start()

    snpFinder.find_snp_by_genes()

    spinner.stop()

    if cmd_line.args.generateCSV:
        snpFinder.relevant_SNPs_df.to_csv(cmd_line.args.outFile,index=False)
    else:
        output_writer = open(cmd_line.args.outFile, "w")
        snpFinder.print_human_report(output_writer)
        output_writer.close()

    if cmd_line.args.color:
        snpFinder.print_color_report()


if __name__ == "__main__":
    main()
