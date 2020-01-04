#!/bin/bash
##
## Subsets DNA meth data from BRCA cohort in both CpG and patients axes
##
## 09 Dec 2019
## Izaskun Mallona

WD=$HOME/tmp/tcga
# NUM_PATIENTS=50 # not used, is hardcoded
NUM_CPGS=10000

mkdir -p $WD

cd $WD

wget http://gdac.broadinstitute.org/runs/stddata__2016_01_28/data/BRCA/20160128/gdac.broadinstitute.org_BRCA.Merge_methylation__humanmethylation450__jhu_usc_edu__Level_3__within_bioassay_data_set_function__data.Level_3.2016012800.0.0.tar.gz

tar xzvf gdac.broadinstitute.org_BRCA.Merge_methylation__humanmethylation450__jhu_usc_edu__Level_3__within_bioassay_data_set_function__data.Level_3.2016012800.0.0.tar.gz

## probably unnecesary, but `file` says linefeeds are weird
mac2unix gdac.broadinstitute.org_BRCA.Merge_methylation__humanmethylation450__jhu_usc_edu__Level_3__within_bioassay_data_set_function__data.Level_3.2016012800.0.0/BRCA.methylation__humanmethylation450__jhu_usc_edu__Level_3__within_bioassay_data_set_function__data.data.txt


## check number of records, e.g. whow many columns are there

awk '{print NF}' gdac.broadinstitute.org_BRCA.Merge_methylation__humanmethylation450__jhu_usc_edu__Level_3__within_bioassay_data_set_function__data.Level_3.2016012800.0.0/BRCA.methylation__humanmethylation450__jhu_usc_edu__Level_3__within_bioassay_data_set_function__data.data.txt | head

## get every fourth sample; note that beta values repeat every 4 columns
## also append the probe, chr, start and genename (e.g. extracted from columns
# belonging to the first sample

cut  -f"1,3,4,5,"$(seq -s, 6 16 3541) gdac.broadinstitute.org_BRCA.Merge_methylation__humanmethylation450__jhu_usc_edu__Level_3__within_bioassay_data_set_function__data.Level_3.2016012800.0.0/BRCA.methylation__humanmethylation450__jhu_usc_edu__Level_3__within_bioassay_data_set_function__data.data.txt | head -3

## add header of those samples

cut  -f"1,3,4,5,"$(seq -s, 6 16 3541) gdac.broadinstitute.org_BRCA.Merge_methylation__humanmethylation450__jhu_usc_edu__Level_3__within_bioassay_data_set_function__data.Level_3.2016012800.0.0/BRCA.methylation__humanmethylation450__jhu_usc_edu__Level_3__within_bioassay_data_set_function__data.data.txt | head -2  > randomized_meth_brca.txt

## randomize and get some CpGs from these; append to the hader
cut  -f"1,3,4,5,"$(seq -s, 6 16 3541) gdac.broadinstitute.org_BRCA.Merge_methylation__humanmethylation450__jhu_usc_edu__Level_3__within_bioassay_data_set_function__data.Level_3.2016012800.0.0/BRCA.methylation__humanmethylation450__jhu_usc_edu__Level_3__within_bioassay_data_set_function__data.data.txt | \
    tail -n +2 | \
    shuf -n $NUM_CPGS >> randomized_meth_brca.txt

# mind that the second file will is a header as well! to read into R try

# fh <- readLines('randomized_meth_brca.txt')
# d <- read.table(textConnection(fh[-2]), header = TRUE, sep = "\t")
# colnames(d)[1:4] <- c('probe', 'gene_symbol', 'chr', 'start')
