#! /usr/env/python                          #
# CDC group                                 #
# Human Technopole                          #
#############################################

# Import modules
## python modules
import anndata
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import scanpy as sc
import seaborn as sns
from itertools import *
import warnings
warnings.filterwarnings('ignore')

import os
import sys
import random
import time

# Import inner modules here
#from coreFunctions import * 
import argsParse

# Functions 
def plot_n_counts(df, what, out):
    
    """"""
    plt.figure(figsize=(15,8))
    sns.violinplot(df, x='condition', y=what)
    plt.savefig(os.path.join(out, f'{what}-violin.png'))


def check_directory(dire):

    """Check if a given path is a directory and if not it creates it"""
    if not os.path.isdir(dire):
        try:
            os.makedirs(dire)
        except OSError as e:
            sys.exit("Error creating {} directory. Check path and permissions".format(dire))
        except:
            sys.exit("Unexpected error:", sys.exc_info()[0])
            raise


def pseudobulk_replicates(comp2, func='sum'):

    """ lets shuffle then and generate 2 replicates each
    in this function they also filter donors that are less than n=30 cells from the 
    populations studied (not doing it here)
    here the comparison is not within cell groups but among cell types so we will 
    just continue with the filter we did before and
    we can increase the number of cells required for the two cell 
    groups interrogated after seeing the results """
    pb = []
    for i, donor in enumerate(comp2.obs.sample_name.unique().tolist()):
        print(f"\tProcessing donor {i+1}")
        ## subset donor
        comp_donor = comp2[comp2.obs.sample_name == donor]
        ## create replicates
        indices = list(comp_donor.obs_names)
        random.shuffle(indices)
        ## split this cell names into two groups
        indices = np.array_split(np.array(indices), 2)
        for r, rep in enumerate(indices):
            comp_rep = comp_donor[rep]
   
            ## aggregate
            if func == 'sum':
                df_donor = sc.AnnData(X = comp_rep.X.sum(axis = 0),
                                    var = comp_rep.var)
            elif func == 'mean':
                df_donor = sc.AnnData(X = comp_rep.X.mean(axis = 0),
                                    var = comp_rep.var)
            df_donor.obs['sample_name'] = f'{donor}_{r}'
            df_donor.obs['replicate'] = r
            df_donor.obs.index = df_donor.obs['sample_name']
            df_donor.obs.index.names = ['']
            ## we can also separate the name between the condition and the sample when creating the sample_name field above
            df_donor.obs['condition'] = comp_donor.obs['condition'].unique().tolist()[0]
            pb.append(df_donor)

     
    return pb


def pseudobulk_noreplicates(comp2, fields, func='sum'):
    """Generate pseudobulk counts without replicates"""
    pb = []
    for i, donor in enumerate(comp2.obs.sample_name.unique().tolist()):
        print(f"\tProcessing donor {i+1}")
        ## subset donor
        comp_donor = comp2[comp2.obs.sample_name == donor]
   
        ## aggregate
        if func == 'sum':
            df_donor = sc.AnnData(X = comp_donor.X.sum(axis = 0),
                                var = comp_donor.var) 
        elif func == 'mean':
            df_donor = sc.AnnData(X = comp_donor.X.mean(axis = 0),
                                var = comp_donor.var) 
        df_donor.obs['sample_name'] = donor
        df_donor.obs['sample'] = comp_donor.obs[fields].unique().tolist()[0]
        df_donor.obs.index = df_donor.obs['sample_name']
        df_donor.obs.index.names = ['']
        ## we can also separate the name between the condition and the sample when 
        ## creating the sample_name field above
        df_donor.obs['condition'] = comp_donor.obs['condition'].unique().tolist()[0]
        pb.append(df_donor)
        
    return pb
        

def pseudobulk_prep(df, sample_field, anot_field, condition1, condition2, 
                    outDir, replicates=False):
    """select the conditions (cell types) to perform the differential 
    expression prep tables"""
    folder = os.path.join(outDir, ('{}_{}').format(condition1, condition2).replace(' ',''))
    out = folder
    if not os.path.isdir(folder):
        os.makedirs(folder)

    ## avoid spaces or other problematic characters on annotation fields
    df.obs['condition'] = np.where(df.obs[anot_field] == condition1, condition1.replace(' ', '-'), 
                                   '')
    df.obs['condition'] = np.where(df.obs[anot_field] == condition2, condition2.replace(' ', '-'), 
                                   df.obs['condition'])

    ## we select only cells for the cell groups to be compared 
    dfsel = df[df.obs['condition'] != '']

    ## create a column with sample + condition
    dfsel.obs['sample_name'] = (dfsel.obs[sample_field].astype(str) + '_' + 
                                dfsel.obs['condition'].astype(str))
    if dfsel.obs['sample_name'].nunique() == dfsel.obs[sample_field].nunique():
        sfield = 'sample_name'
    else:
        sfield = sample_field
        
    ## select samples without observations
    samp_ava = dfsel.obs.groupby([sfield, 'condition']).size().reset_index()
    samp_ava.columns = ['sample', 'condition', 'counts']
    ## we will remove any sample with 0 cells either in one condition or the other
    samplestoremove = samp_ava[samp_ava['counts']==0]['sample'].unique().tolist()

    ## store new df with valid samples
    dfsel2 = dfsel[~dfsel.obs[sfield].isin(samplestoremove)]
    used_samples = len(dfsel2.obs[sfield].unique().tolist())

    ## filters regarding min genes and min cells
    sc.pp.filter_cells(dfsel2, min_genes=200)
    sc.pp.filter_genes(dfsel2, min_cells=3)
    ## plot differences between conditions  
    #plot_n_counts(dfsel2.obs, 'n_counts', out)
    #plot_n_counts(dfsel2.obs, 'n_genes', out)

    ## convert condition columns to category
    dfsel2.obs[sample_field] = dfsel2.obs[sample_field].astype("category")
    dfsel2.obs.sample_name = dfsel2.obs.sample_name.astype("category")

    total_cells = len(dfsel2)
    cond_cells = len(dfsel2[dfsel2.obs[anot_field] == condition1])
    outlog = os.path.join(out, 'info.log')
    with open(outlog, 'w') as flog:
        flog.write(f'# n.samples {len(samplestoremove)} are going to be removed')
        flog.write(f'# Using {used_samples} samples')
        flog.write(f'group,nsamples_used,total_cells,group_cells\n{condition1}-{condition2},'
                    f'{used_samples},{total_cells},{cond_cells}\n')
    
    ## pseudobulk table
    if replicates:
        print('Using replicates')
        pb = pseudobulk_replicates(dfsel2)
    else:
        print('Pseudobulk without replicates')
        pb = pseudobulk_noreplicates(dfsel2, sample_field)

    pbs = sc.concat(pb)
    pbs.layers['counts'] = pbs.X.copy()
    pbs.raw = pbs
    pbs.var_names_make_unique()
    ## plot PCA? we would norm before so that is why here we copy again the raw counts
    #pbs.X = pbs.layers['counts'].copy()
    counts = pd.DataFrame(pbs.X, columns = pbs.var_names, index = pbs.obs.index) 
    #need to do this to pass var names
    ## save also only counts of genes > 1000
    countsfilt = counts[counts.columns[counts.sum() > 1000]]
    ## save
    #pbs.write(os.path.join(folder, 'conditions.h5'))
    counts.to_csv(os.path.join(folder, 'counts.csv'), sep="\t")
    countsfilt.to_csv(os.path.join(folder, 'countsabove1000.csv'), sep="\t")
    pbs.obs.to_csv(os.path.join(folder, 'conditions.csv'), 
                   sep="\t", index=False)


def pseudobulk_prep_comp1(dfsel, sample_field, anot_field, condition1, 
                          outDir, replicates=False):
    """select the conditions (cell types) to perform the differential 
    expression prep tables one cell type vs the rest"""
    out = os.path.join(outDir, condition1.replace(' ',''))
    if not os.path.isdir(out):
        os.makedirs(out)
    ## avoid spaces or other problematic characters on annotation fields
    dfsel.obs['condition'] = np.where(dfsel.obs[anot_field] == condition1, 
                                       condition1.replace(' ', '-'), 
                                       'rest')

    ## create a column with sample + condition
    dfsel.obs['sample_name'] = (dfsel.obs[sample_field].astype(str) + '_' + 
                                dfsel.obs['condition'].astype(str))
    ## select samples without observations
    samp_ava = dfsel.obs.groupby([sample_field,'condition']).size().reset_index()
    samp_ava.columns = ['sample', 'condition', 'counts']
    ## we will remove any sample with 0 cells either in one condition or the other
    samplestoremove = samp_ava[samp_ava['counts']==0]['sample'].unique().tolist()

    ## store new df with valid samples
    dfsel2 = dfsel[~dfsel.obs[sample_field].isin(samplestoremove)]
    used_samples = len(dfsel2.obs[sample_field].unique().tolist())
    
    ## filters regarding min genes and min cells
    sc.pp.filter_cells(dfsel2, min_genes=200)
    sc.pp.filter_genes(dfsel2, min_cells=3)
    ## plot differences between conditions  
    #plot_n_counts(dfsel2.obs, 'n_counts', out)
    #plot_n_counts(dfsel2.obs, 'n_genes', out)

    ## convert condition columns to category
    dfsel2.obs[sample_field] = dfsel2.obs[sample_field].astype("category")
    dfsel2.obs.sample_name = dfsel2.obs.sample_name.astype("category")

    total_cells = len(dfsel2)
    cond_cells = len(dfsel2[dfsel2.obs[anot_field] == condition1])
    outlog = os.path.join(out, 'info.log')
    with open(outlog, 'w') as flog:
        flog.write(f'# n.samples {len(samplestoremove)} are going to be removed')
        flog.write(f'# Using {used_samples} samples')
        flog.write(f'group,nsamples,total_cells,group_cells\n{condition1},{used_samples},{total_cells},{cond_cells}\n')
            
    ## pseudobulk table 
    if replicates:
        print('Using replicates')
        pb = pseudobulk_replicates(dfsel2)
    else:
        print('Pseudobulk without replicates')
        pb = pseudobulk_noreplicates(dfsel2, sample_field)

    pbs = sc.concat(pb)
    #pbs.layers['counts'] = pbs.X.copy()
    #pbs.raw = pbs
    pbs.var_names_make_unique()
    ## plot PCA? only norm before PCA
    #pbs.X = pbs.layers['counts'].copy()
    counts = pd.DataFrame(pbs.X, columns = pbs.var_names, index = pbs.obs.index) #need to do this to pass var names
    ## save also only counts of genes > 1000
    countsfilt = counts[counts.columns[counts.sum() > 1000]]
    ## save
    #pbs.write(os.path.join(out, 'conditions.h5'))
    counts.to_csv(os.path.join(out, 'counts.csv'), sep="\t")
    countsfilt.to_csv(os.path.join(out, 'countsabove1000.csv'), sep="\t")
    pbs.obs.to_csv(os.path.join(out, 'conditions.csv'), 
                   sep="\t", index=False)


def pseudobulk_prep_allclusters(df, sample_field, anot_field, outDir, replicates=False):
    
    """prep tables to perform the differential expression across cell types.
    Be aware that the samples column "Sample" has been taken as literal from
    the dataframe used currently, you will need
    to change this in the function if you work with other dataset
    """
    #out = os.path.join(outDir, '')
    ## avoid spaces or other problematic characters on annotation fields
    df.obs['condition'] = df.obs[anot_field]
    ## we select only cells for the cell groups to be compared 
    dfsel = df[~df.obs['condition'].isnull()]

    ## create a column with sample + condition
    dfsel.obs['sample_name'] = (dfsel.obs['Sample'].astype(str) + '_' + 
                                dfsel.obs['condition'].astype(str))
    
    ## filters regarding min genes and min cells
    sc.pp.filter_cells(dfsel, min_genes=200)
    sc.pp.filter_genes(dfsel, min_cells=3)
    ## convert condition columns to category
    dfsel.obs.Sample = dfsel.obs.Sample.astype("category")
    dfsel.obs.sample_name = dfsel.obs.sample_name.astype("category")

    ## pseudobulk table with no replicates
    if replicates:
        print('Using replicates')
        pb = pseudobulk_replicates(dfsel)
    else:
        print('Pseudobulk without replicates')
        pb = pseudobulk_noreplicates(dfsel)
    pbs = sc.concat(pb)
    pbs.layers['counts'] = pbs.X.copy()
    pbs.raw = pbs
    pbs.var_names_make_unique()
    ## plot PCA?
    pbs.X = pbs.layers['counts'].copy()
    counts = pd.DataFrame(pbs.X, columns = pbs.var_names,
                            index=pbs.obs.index) #need to do this to pass var names
    ## save also only counts of genes > 1000
    countsfilt = counts[counts.columns[counts.sum() > 1000]]
    ## save
    folder = os.path.join(outDir, 'edgeRdata')
    if not os.path.isdir(folder):
        os.makedirs(folder)
    #pbs.write_h5ad(os.path.join(folder, 'conditions.h5'))
    counts.to_csv(os.path.join(folder, 'counts.csv'), sep="\t")
    countsfilt.to_csv(os.path.join(folder, 'countsabove1000.csv'), sep="\t")
    pbs.obs.to_csv(os.path.join(folder, 'conditions.csv'), 
                   sep="\t", index=False)


###################################### MAIN ########################################

def main():

    """Generate pseudobulk data for running deseq2 after a single-cell expression matrix.
    """

    # Get arguments and unpack
    args = argsParse.parse_args()
    data, samples, column, outfolder, modality, replicates = (args.data, args.samples, 
                                                                args.annot, 
                                                                args.outf, args.modality,
                                                                args.repl)
    ## create output folder
    check_directory(outfolder)
    ## read anndata object 
    adata = sc.read_h5ad(data)
    ## if filter argument, filter on a column in obs(keep supported)
    if args.fcolobs and args.flt:
        assert isinstance(args.flt, str)
        assert args.fcolobs in adata.obs
        adata = adata[adata.obs[args.fcolobs] == args.flt].copy()
    
    elif args.fcolobs or args.flt:
        if args.fcolvar and args.flt:
            assert args.fcolvar in adata.var
            adata = adata[:, adata.var[args.fcolvar] == args.flt].copy()
        else:
            sys.exit('Not provided both arguments required for filtering')
    else:
        pass

    ## remove samples with nan values
    adata = adata[~adata.obs[samples].isna()]
    ## replace some strange characters in the condition field
    adata.obs[column] = adata.obs[column].str.replace('/', '-')
    ## Define logging
    #log_setup(args)

    # define global variables
    global cdate
    cdate = time.strftime('%d-%m-%Y_%H-%M-%S')

    # we need raw counts so will ensure this at the beginning
    if 'counts' in adata.layers.keys():
        adata.X = adata.layers['counts'].copy()
    else:
        sys.exit('No raw counts stored in adata.layers as "counts"... Exiting')

    ## modalities dict
    dm = {'pairwise' : pseudobulk_prep, 
            'onevsrest' : pseudobulk_prep_comp1}
    if not modality in dm:
        sys.exit('incorrect modality workflow')
    ## choose workflow
    ### quick workaround
    lfields = adata.obs[column].unique().tolist()
    if modality == 'pairwise':
        l = [p for p in combinations(lfields, r=2)]
        for e in l:
            print(f'Generating table for comparison {e}')
            dm[modality](adata, samples, column, *(e), outfolder, 
                                            replicates)
    elif modality == 'onevsrest':
        for i in lfields:
            print(f'Generating table for {i}')
            dm[modality](adata, samples, column, i, outfolder, replicates)
    else:
        sys.exit('Could not select among available workflows. Use --help to check options.')


if __name__ == "__main__":
    main()





