
WD=$HOME/tmp/tcga
NUM_PATIENTS=50
NUM_CPGS=10000

mkdir -p $WD

cd $WD


wget http://gdac.broadinstitute.org/runs/stddata__2016_01_28/data/PRAD/20160128/gdac.broadinstitute.org_PRAD.Merge_methylation__humanmethylation450__jhu_usc_edu__Level_3__within_bioassay_data_set_function__data.Level_3.2016012800.0.0.tar.gz &

wget http://gdac.broadinstitute.org/runs/stddata__2016_01_28/data/COAD/20160128/gdac.broadinstitute.org_COAD.Merge_methylation__humanmethylation450__jhu_usc_edu__Level_3__within_bioassay_data_set_function__data.Level_3.2016012800.0.0.tar.gz &

wget http://gdac.broadinstitute.org/runs/stddata__2016_01_28/data/LUAD/20160128/gdac.broadinstitute.org_LUAD.Merge_methylation__humanmethylation450__jhu_usc_edu__Level_3__within_bioassay_data_set_function__data.Level_3.2016012800.0.0.tar.gz

for cohort in PRAD COAD LUAD
do
    echo $cohort

    tar xzvf gdac*"$cohort"*tar.gz

    ## probably unnecesary, but `file` says linefeeds are weird
    ## mac2unix gdac.broadinstitute.org_BRCA.Merge_methylation__humanmethylation450__jhu_usc_edu__Level_3__within_bioassay_data_set_function__data.Level_3.2016012800.0.0/BRCA.methylation__humanmethylation450__jhu_usc_edu__Level_3__within_bioassay_data_set_function__data.data.txt


    ## check number of records, e.g. whow many columns are there

    awk '{print NF}' gdac*"$cohort"*/*humanmethylation450*__data.data.txt| head

    ## get every fourth sample; note that beta values repeat every 4 columns
    ## also append the probe, chr, start and genename (e.g. extracted from columns
    # belonging to the first sample

    cut  -f"1,3,4,5,"$(seq -s, 6 16 3541) gdac*"$cohort"*/*humanmethylation450*__data.data.txt| head -3

    ## add header of those samples

    cut  -f"1,3,4,5,"$(seq -s, 6 16 3541) gdac*"$cohort"*/*humanmethylation450*__data.data.txt| head -2  > randomized_meth_"$cohort".txt

    ## randomize and get some CpGs from these; append to the hader
    cut  -f"1,3,4,5,"$(seq -s, 6 16 3541) gdac*"$cohort"*/*humanmethylation450*__data.data.txt| \
	tail -n +2 | \
	shuf -n $NUM_CPGS >> randomized_meth_"$cohort".txt
done

gzip randomized*txt
