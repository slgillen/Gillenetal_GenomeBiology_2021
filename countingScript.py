#!/usr/bin/python
#these scripts were adapted from https://pythonhosted.org/riboplot/ribocount.html

import sys
import pysam
import os
import argparse
import logging
import countingFunctions
import collections

log = logging.getLogger('riboplot')

ch = logging.StreamHandler()
ch.setLevel(logging.INFO)

formatter = logging.Formatter('%(asctime)s - %(module)s %(levelname)s %(message)s', datefmt='%Y-%m-%d %H:%M:%S')
ch.setFormatter(formatter)
log.addHandler(ch)


def create_parser():
    """Argument parser. """
    parser = argparse.ArgumentParser(
        prog='countingScript.py', description='Output read counts for all transcripts')

    # required arguments
    required = parser.add_argument_group('required arguments')
    required.add_argument('-b', '--BAMinput', help='Ribo-Seq alignment file in BAM format', required=True)
    required.add_argument('-f', '--transcriptome_fasta', help='FASTA format file of the transcriptome', required=True)

    # optional arguments
    parser.add_argument('-l', '--read_lengths', help='Read lengths to consider (default: %(default)s). '
                        'Multiple read lengths should be separated by commas. If multiple read lengths '
                        'are specified, corresponding read offsets should also be specified. If you do '
                        'not wish to apply an offset, please input 0 for the corresponding read length',
                        default='0', type=countingFunctions.lengths_offsets)
    parser.add_argument('-s', '--read_offsets', help='Read offsets (default: %(default)s). '
                        'Multiple read offsets should be separated by commas',
                        default='0', type=countingFunctions.lengths_offsets)

    parser.add_argument('-o', '--outpath', help='Files are saved in this directory', default='output')
    parser.add_argument('-of', '--outfile', help='name to save file', default='outputfile')

    return parser


def main(args):
    (BAMinput,fasta_file,read_lengths,read_offsets,outpath,outfile) = \
         (args.BAMinput,args.transcriptome_fasta,args.read_lengths,args.read_offsets,args.outpath,args.outfile)
    
    log.debug('Supplied arguments\n{}'.format(
        '\n'.join(['{:<20}: {}'.format(k, v) for k, v in vars(args).items()])))

    log.info('check 1')

    with pysam.AlignmentFile(BAMinput, 'rb') as b, pysam.FastaFile(fasta_file) as f:
    
        count=0
    
    #create output directories
        if not os.path.exists(outpath):
            os.mkdir(outpath)
        os.chdir(outpath)
        log.info('directory created')

        log.info('Getting RPF counts for all transcripts in FASTA')
        with open(outfile,'aw') as outputfile:
            outputfile.write('Transcript_ID'+'\t'+'Transcript_Sequence'+'\t'+'pos_RPF_counts'+'\t'+'Total_RPF_counts'+'\n')

            for transcript in f.references:
                transcript_sequence=f.fetch(transcript)
                ribo_counts, ribo_reads = countingFunctions.get_RPF_counts(pysamobj=b,transcript_name=transcript,transcript_length=len(transcript_sequence),read_lengths=read_lengths,read_offsets=read_offsets)
                if ribo_reads>1:
                    outputfile.write(str(transcript)+'\t'+str(transcript_sequence)+'\t'+' '.join(map(str,ribo_counts.values()))+'\t'+str(ribo_reads)+'\n')
               
        outputfile.close()
        log.info('RPF counting complete')


def run():
    """Run program"""
    parsed = create_parser()
    args = parsed.parse_args()
    log.info('args')
    print(args)
    main(args)



if __name__ == '__main__':
    run()
