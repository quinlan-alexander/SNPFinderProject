import arv
import requests
from bs4 import BeautifulSoup
import re
import json

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
        if inOpts is None:
            self.args = self.parser.parse_args()
        else:
            self.args = self.parser.parse_args(inOpts)

class Twenty3andMeData():

    def __init__(self):
        pass

    def get_gene(self, gene):
        r = requests.get('https://api.23andme.com/3/marker/?gene_name='+gene)
        return json.loads(r.text)


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
        soup = BeautifulSoup(r.text, "html5lib")
        trows = soup('table')

        # table element 1 contains orientation
        p = re.compile('(<td>)(.+)(<\/td>)')
        orientation = str(trows[1])
        m = p.findall(orientation)
        if m:
            self.orientation = m[0][1]
        else:
            self.orientation = 'Error: orientation for SNP '+rsid+' not found'

        # table element 3 contains SNP allele details
        allele_data = str(trows[3])
        p = re.compile('(")(Rs\d+)(\([ACGT];[ACGT]\))(")')
        m = p.findall(allele_data)
        if m:
            self.common_allele = m[0][2]
        else:
            self.common_allele = 'Error: common allele for SNP '+rsid+' not found'

    def get_common_allele(self, request_orientation=None):
        if request_orientation == None or request_orientation == self.orientation:
            return self.common_allele

        return self.complement(self.common_allele)

class SNPFinder():

    def __init__(self, genome, genes):
        self.genome = genome
        self.genes = genes


    def find_snp_by_genes(self):
        for gene in self.genes:
            # look up gene on 23andme to get list of SNP rsid
            pass
            # for SNP in rsid collection
            #   get SNP data from SNPedia (common allele, orientation)
            #   get 23andme genome allele from loaded file
            #   if affected allele (23andme.genome(rsid).allele.format(...) != SNPedia.rsid.common_allele('positive')
            #       add rsid to affected collection

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

    snpFinder = SNPFinder(genome,genes)

    # look up gene in 23andme file and get any SNPs for the subject and return array of rsids

    print(genome["rs1801133"])

    r = requests.get('https://www.snpedia.com/index.php/Rs1234')

    # print(genes)

    output_writer = open(cmd_line.args.outFile, "w")

    # close output file
    output_writer.close()


if __name__ == "__main__":
    main()
