# gwas101

This repo has some misc examples and helper scripts dealing with genotype data and GWAS.

* [gwas101.ipynb](gwas101.ipynb) downloads a synthetic dataset (``chr21.[bed/bim/fam]``) with genotypes for N=10.000 subjects and M=150.000 variants on chr21, and uses this data to illustrate GWAS and polygenic risk scoring. To synthesize the output phenotype the tool uses [simu_linux](https://github.com/precimed/simu) tools, which is not available for MAC or Windows (if you have time to fix this please let me know!).
``gwas101.ipynb`` also examplifies ``read_bed``, ``write_bed``, ``nancorr`` functions from [gwas_tools.py](gwas_tools.py) - some basic helper functions equivalent to MATLAB scripts (``PlinkRead_binary2.m``,  ``nancorr.m``)
* [sumstats/clean_scz2_gwas.ipynb](sumstats/clean_scz2_gwas.ipyn) and [sumstats/manipulate_gwas_sumstats.ipynb](sumstats/manipulate_gwas_sumstats.ipynb) give few examples of dealing with external summary statistics. 
* [sumstats/sumstats/daner_PGC_SCZ52_0513a.hq2.yaml](sumstats/sumstats/daner_PGC_SCZ52_0513a.hq2.yaml) gives an example of a .yaml file needed for https://github.com/BioPsyk/cleansumstats pipeline.
* [snp_canonic_id.ipynb](snp_canonic_id.ipynb) is not something you want to use, just ignore this!

# openSNP

Another useful resource is https://github.com/ofrei/opensnp, see its README file for more info.
