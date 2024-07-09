import pandas as pd
from Bio import SeqIO
import pyNGS as pyNGS


def main():
        import sys
        import argparse
	
        parser = argparse.ArgumentParser(description="Process NGS files and output csv.")

        parser.add_argument('input', help="The input NGS file in fastQ format")
        parser.add_argument('refseq', help="The reference nucleotide sequence in fasta format")
        parser.add_argument('primerseq', help="supply the forward primer for determining sequencing offset in fasta format")
        parser.add_argument('output',help="name of output csv file")
        parser.add_argument('-v', '--verbose', action='store_true', help="Increase output verbosity")
        args = parser.parse_args()
        file = getattr(args, 'input')
        refSeq = getattr(args,'refseq')
        primerseq = getattr(args,'primerseq')
        templateSeq = str(next(SeqIO.parse(refSeq,"fasta")).seq)
        fwdprimer = str(next(SeqIO.parse(primerseq,"fasta")).seq)
        df=pyNGS.readNGSfiles(file,templateSeq,fwdprimer)
        df.to_csv(args.output+'.csv')

	
if __name__ == "__main__":
	main()

