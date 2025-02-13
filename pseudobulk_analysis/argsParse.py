#! /usr/env/python                          #
# CDC group                                 #
# Human Technopole                          #
#############################################

# Argsparse in standard mode for general pipeline
# Import modules
import argparse
import sys

def parse_args():

    """Parse args with argparse"""
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--data', '-d', dest='data', 
                required=True, help='init GEX matrix')
    parser.add_argument('--samples', '-s', dest='samples', 
                required=True, help='.obs column for samples (replicates)')
    parser.add_argument('--annotation', '-a', dest='annot', 
                required=True, help='.obs column for DE')
    parser.add_argument('--out', '-o', dest='outf', 
                required=True, help='Output directory')
    parser.add_argument('--mode', '-m', dest='modality', 
                required=True, help='Choose between wf1: pairwise or wf2: onevsrest')
    parser.add_argument('--replicates', '-r', dest='repl', 
                required=False, action='store_true',
                help='whether to create replicates for samples within pseudobulk')
    parser.add_argument('--filterobscol', '-b', dest='fcolobs', 
                required=False,
                help='column to filter on adata.obs (if any filters required)')
    parser.add_argument('--filtervarcol', '-v', dest='fcolvar', 
                required=False,
                help='column to filter on adata.var (if any filters required)')
    parser.add_argument('--filterval', '-f', dest='flt', 
                required=False,
                help='value to filter on adata.obs column or adata.var (given in --filtervarcol or --filterobscol)')
 
    args = parser.parse_args()

    return args