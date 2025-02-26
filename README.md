# pseudobulk DGE analysis

1. Generate pseudobulk counts from scRNA-seq data.

## Usage
```
python pseudobulk_analysis/generate_pseudobulk.py -d anndata -s aggregation -a condition -o output_folder -m mode

```

```
usage: generate_pseudobulk.py [-h] --data DATA --samples SAMPLES --annotation ANNOT --out OUTF --mode MODALITY
                             [--replicates] [--filterobscol FCOLOBS] [--filtervarcol FCOLVAR] [--filterval FLT]```

options:
  -h, --help            show this help message and exit
  --data DATA, -d DATA  init GEX matrix
  --samples SAMPLES, -s SAMPLES
                        .obs column for samples (replicates)
  --annotation ANNOT, -a ANNOT
                        .obs column for DE
  --out OUTF, -o OUTF   Output directory
  --mode MODALITY, -m MODALITY
                        Choose between wf1: pairwise or wf2: onevsrest
  --replicates, -r      whether to create replicates for samples within pseudobulk
  --filterobscol FCOLOBS, -b FCOLOBS
                        column to filter on adata.obs (if any filters required)
  --filtervarcol FCOLVAR, -v FCOLVAR
                        column to filter on adata.var (if any filters required)
  --filterval FLT, -f FLT
                        value to filter on adata.obs column or adata.var (given in --filtervarcol or --filterobscol)`
```

## Input Parameters

| Argument              | Description      |    
| --------------------- | ---------------- |
| `data`            | anndata object .h5ad format        
| `samples`        | .obs column that contains replicate information, this will be used to aggregate counts into pseudobulks |
| `annotation`          | .obs column containing variable used for DGE | 
| `out`            | Desired output directory | 
| `mode`           | two modalities can be chosen for the contrast: **pairwise**; all paired comparisons per category in *--annotation*. **onevsrest**; each category in *--annotation* against all the rest (aggregated) |
| `replicates`     | if included, create pseudoreplicates |
| `filterobscol`   | column in adata.obs to be filtered (if any filters required) |
| `filtervarcol` | column in adata.var to be filtered (if any filters required) |
|`filterval` | category to be selected in adata.obs/adata.var column given in *--filterobscol* or *--filtervarcol* | 

2. Run DeSeq2 differential gene expression test.

Example script for simple contrasts: pseudobulk_analysis/deseq2.R

```
Rscript /group/dominguez/aft/scripts/IEI_popgen/pseudobulk_analysis/deseq2.R input_folder output_folder control_var treatment_var
```


