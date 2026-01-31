import pyfaidx
import pyfastx
import argparse
from .. import __version__

parser = argparse.ArgumentParser()

parser.add_argument('-v', '--version', action='version', version=f'%(prog)s {__version__}' )

identity_group = parser.add_argument_group(description='Generate candidate eccDNA from rolling circle reads:')
identity_group.add_argument("--fastq", type=str, required=True, help='input reads in fastq format')
identity_group.add_argument("--paf", type=str, required=True, help='input alignments in PAF format')
identity_group.add_argument("--reference", type=str, required=True,
                            help='reference genome sequences in fasta format')
identity_group.add_argument("--output", '-o', type=str, default='.',
                            help='output folder, the default is the current folder')
identity_group.add_argument('--maxOffset', type=int, default=20,
                            help='maximum offset of start/end positions between two sub-reads to be considered as mapping to the same location, default is 20')
identity_group.add_argument('--minMapQual', type=int, default=30,
                            help='minimum mapping quality of sub-reads, default is 30')

polish_group = parser.add_argument_group(description='Get accurate eccDNA location from candidate eccDNA:')
polish_group.add_argument('--fluctuate', type=int, default=20,
                          help='maximum offset between the candidate eccDNA and the candidate eccDNA in a group, default is 20bp')
polish_group.add_argument('--nf', type=int, default=2,
                          help='minimum number of fullpass to support candidate eccDNA, default is 2')
polish_group.add_argument('--cov', type=int, default=1,
                          help='minimum fullpass used for merging a set of candidate eccDNA, default is 1')
polish_group.add_argument('--nc', type=int, default=2,
                          help='minimum number of candidate eccDNA required to consolidate_candidate_groups overlapped candidate eccDNAs, default is 2')

consensus_group = parser.add_argument_group(description='Generate consensus sequences and variants:')
consensus_group.add_argument('--minDP', type=int, default=4, help='minimum depth to call variants, default is 4')
consensus_group.add_argument('--minAF', type=float, default=0.75,
                             help='minimum alternative allele frequency to call variants, default is 0.75')



args = parser.parse_args()

genome = pyfaidx.Fasta(args.reference)
fastq = pyfastx.Fastq(args.fastq)
