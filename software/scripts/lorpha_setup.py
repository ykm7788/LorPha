import argparse
import sys
import pandas as pd

def main():
    parser = argparse.ArgumentParser(
        description="Generate LorPha analysis pipeline shell script",
        formatter_class=argparse.RawTextHelpFormatter
    )

    parser.add_argument('-r', required=True, metavar='R',
                        help='REQUIRED: Path to rawdata.txt (sample information table)')
    parser.add_argument('-p', required=True, metavar='P',
                        help='REQUIRED: Path to parameter.txt (software parameters)')
    parser.add_argument('-t', required=True, metavar='T',
                        help='REQUIRED: Number of threads')
    parser.add_argument('-d', required=True, metavar='D',
                        help='REQUIRED: Path to pre-installed databases')
    parser.add_argument('-o', required=True, metavar='O',
                        help='REQUIRED: Output directory ')

    
    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)

    args = parser.parse_args()

    pd.set_option('display.max_columns', None)
    pd.set_option('display.max_rows', None)
    pd.set_option('max_colwidth', None)

    rawdata = pd.read_table(args.r, header=0, dtype=str)
    parameter = pd.read_table(args.p, header=0, dtype=str)

    rawdata_name1 = rawdata['Contig'].str.rsplit('/', n=1).str[-1]
    rawdata_name2 = rawdata_name1.str.rsplit('.', n=1).str[0]

    job1_phabox = f"mkdir -p {args.o}/analysis_contig/{rawdata['SampleID']}/phabox ; "
    job1_genomad = f"mkdir -p {args.o}/analysis_contig/{rawdata['SampleID']}/genomad ; "
    job1_vhost = f"mkdir -p {args.o}/analysis_contig/{rawdata['SampleID']}/vhost ; "

    job1_mkdir = f"""
#### make necessary directories ####
mkdir -p {args.o}/analysis_votu/checkv
mkdir -p {args.o}/analysis_votu/vclust
mkdir -p {args.o}/analysis_votu/dbcan
mkdir -p {args.o}/analysis_votu/VIBRANT
mkdir -p {args.o}/analysis_votu/PhaStyle
mkdir -p {args.o}/analysis_votu/coverm
mkdir -p {args.o}/LorPha_results/Coverm
mkdir -p {args.o}/LorPha_results/GeneFunction
mkdir -p {args.o}/analysis_contig/mmseq/tax
mkdir -p {args.o}/analysis_contig/mmseq/tax_rep
mkdir -p {args.o}/analysis_contig/mmseq/rep
mkdir -p {args.o}/analysis_contig/mmseq/tax_data
mkdir -p {args.o}/analysis_contig/mmseq/tax_out
mkdir -p {args.o}/analysis_contig/mmseq/tax_coverm
mkdir -p {args.o}/analysis_contig/mmseq/cluster
mkdir -p {args.o}/analysis_contig/mmseq/tmp
mkdir -p {args.o}/analysis_contig/defensefinder
mkdir -p {args.o}/analysis_contig/coverm
"""

    job1 = job1_mkdir + "\n" + job1_phabox.to_string(index=False, header=False) + "\n" + job1_genomad.to_string(index=False, header=False) + "\n" + job1_vhost.to_string(index=False, header=False) + "\n\n"

    job2_phabox = f"pixi run --environment phabox phabox2 --task end_to_end --dbdir {args.d}/phabox --outpth {args.o}/analysis_contig/{rawdata['SampleID']}/phabox --contigs {rawdata['Contig']} --threads {args.t} {parameter.iloc[0,1]} ; "
    job2_genomad = f"pixi run --environment genomad genomad end-to-end --cleanup {rawdata['Contig']} --threads {args.t} {parameter.iloc[1,1]} {args.o}/analysis_contig/{rawdata['SampleID']}/genomad {args.d}/genomad ; "
    job2_lexicmap1 = f"pixi run lexicmap search -d {args.d}/vhost {rawdata['Contig']} -o {args.o}/analysis_contig/{rawdata['SampleID']}/vhost/contig.vhost.tsv {parameter.iloc[2,1]} --threads {args.t} ; "
    job2_lexicmap2 = f"pixi run find_identical {args.o}/analysis_contig/{rawdata['SampleID']}/vhost/contig.vhost.tsv {args.d}/vhost/virus.seqID.txt -format 1 -site1 4 -site2 1 -type 3 > {args.o}/analysis_contig/{rawdata['SampleID']}/vhost/contig.vhost.result.txt ;"
    job2_merge = f"pixi run merge.phabox_genomad_vhost_contamination -p {args.o}/analysis_contig/{rawdata['SampleID']}/phabox/final_prediction/final_prediction_summary.tsv -g {args.o}/analysis_contig/{rawdata['SampleID']}/genomad/{rawdata_name2}_summary/{rawdata_name2}_virus_summary.tsv -v {args.o}/analysis_contig/{rawdata['SampleID']}/vhost/contig.vhost.result.txt -o {args.o}/analysis_contig/{rawdata['SampleID']} ; "
    job2_seq = f"pixi run seqtk subseq {args.o}/analysis_contig/{rawdata['SampleID']}/phabox/filtered_contigs.fa {args.o}/analysis_contig/{rawdata['SampleID']}/virus.contig.txt > {args.o}/analysis_contig/{rawdata['SampleID']}/virus.pre_check.fa ; "
    job2_combine = f"cat {args.o}/analysis_contig/*/virus.pre_check.fa > {args.o}/analysis_votu/all.virus.pre_check.fa ;"
    job2_checkv = f"pixi run checkv end_to_end -d {args.d}/checkv {args.o}/analysis_votu/all.virus.pre_check.fa {args.o}/analysis_votu/checkv -t {args.t} ; "
    job2_checkv_combine = f"cat {args.o}/analysis_votu/checkv/proviruses.fna {args.o}/analysis_votu/checkv/viruses.fna > {args.o}/analysis_votu/checkv/virus.all.fna ; "
    job2_filter1 = f"pixi run filter_script -a T -b 8 -c Not-determined -d 2 -i {args.o}/analysis_votu/checkv/quality_summary.tsv -o {args.o}/analysis_votu/checkv/filter.virus.tsv ; "
    job2_filter2 = f"cut -f 1 {args.o}/analysis_votu/checkv/filter.virus.tsv > {args.o}/analysis_votu/checkv/filter.contig.tsv ; "
    job2_filter3 = f"pixi run seqtk subseq {args.o}/analysis_votu/checkv/virus.all.fna {args.o}/analysis_votu/checkv/filter.contig.tsv > {args.o}/analysis_votu/all.virus.post_check.fa ; "
    job2_vclust1 = f"pixi run vclust prefilter -i {args.o}/analysis_votu/all.virus.post_check.fa -o {args.o}/analysis_votu/vclust/fltr.txt {parameter.iloc[3,1]} -t {args.t} ; "
    job2_vclust2 = f"pixi run vclust align -i {args.o}/analysis_votu/all.virus.post_check.fa -o {args.o}/analysis_votu/vclust/ani.tsv --filter {args.o}/analysis_votu/vclust/fltr.txt -t {args.t} ; "
    job2_vclust3 = f"pixi run vclust cluster -i {args.o}/analysis_votu/vclust/ani.tsv -o {args.o}/analysis_votu/vclust/clusters.rep.tsv --ids {args.o}/analysis_votu/vclust/ani.ids.tsv {parameter.iloc[4,1]} ; "
    job2_votu1 = f"cut -f 2 {args.o}/analysis_votu/vclust/clusters.rep.tsv | sort | uniq > {args.o}/analysis_votu/uniq.votu.txt ; "
    job2_votu2 = f"pixi run seqtk subseq {args.o}/analysis_votu/all.virus.post_check.fa {args.o}/analysis_votu/uniq.votu.txt > {args.o}/analysis_votu/votu.fa ; "
    job2_votu3 = f"cat {args.o}/analysis_contig/*/virus.txt > {args.o}/analysis_votu/all.virus.profile.txt ;"
    job2_votu4 = f"pixi run merge.virus_profile.checkv.votu -p {args.o}/analysis_votu/all.virus.profile.txt -c {args.o}/analysis_votu/checkv/quality_summary.tsv -v {args.o}/analysis_votu/uniq.votu.txt -o{args.o}/analysis_votu ; "

    job2_votu = "\n#### extract vOTUs ####\n"
    job2 = job2_votu + job2_phabox.to_string(index=False, header=False) + "\n" + job2_genomad.to_string(index=False, header=False) + "\n" + job2_lexicmap1.to_string(index=False, header=False) + "\n" + job2_lexicmap2.to_string(index=False, header=False) + "\n" + job2_merge.to_string(index=False, header=False) + "\n" + job2_seq.to_string(index=False, header=False) + "\n" + job2_combine + "\n" + job2_checkv + "\n" + job2_checkv_combine + "\n" + job2_filter1 + "\n" + job2_filter2 + "\n" + job2_filter3 + "\n" + job2_vclust1 + "\n" + job2_vclust2 + "\n" + job2_vclust3 + "\n" + job2_votu1 + "\n" + job2_votu2 + "\n" + job2_votu3 + "\n" + job2_votu4 + "\n\n"

    job3_cluster = f"pixi run mmseqs easy-cluster {rawdata['Contig']} {args.o}/analysis_contig/mmseq/cluster/{rawdata['SampleID']}.cluster {args.o}/analysis_contig/mmseq/tmp {parameter.iloc[6,1]} --threads {args.t} ; "
    job3_createdb = f"pixi run mmseqs createdb {args.o}/analysis_contig/mmseq/cluster/{rawdata['SampleID']}.cluster_rep_seq.fasta {args.o}/analysis_contig/mmseq/cluster/{rawdata['SampleID']}.cluster_rep_seq.db ; "
    job3_taxonomy = f"pixi run mmseqs taxonomy {args.o}/analysis_contig/mmseq/cluster/{rawdata['SampleID']}.cluster_rep_seq.db {args.d}/gtdb/GTDB {args.o}/analysis_contig/mmseq/tax/{rawdata['SampleID']}.result {args.o}/analysis_contig/mmseq/tmp {parameter.iloc[7,1]} --threads {args.t} ; "
    job3_createtsv = f"pixi run mmseqs createtsv {args.o}/analysis_contig/mmseq/cluster/{rawdata['SampleID']}.cluster_rep_seq.db {args.o}/analysis_contig/mmseq/tax/{rawdata['SampleID']}.result {args.o}/analysis_contig/mmseq/tax_out/{rawdata['SampleID']}.tax.txt --threads {args.t} ; "
    job3_rep = f"cut -f 1 {args.o}/analysis_contig/mmseq/cluster/{rawdata['SampleID']}.cluster_cluster.tsv | sort | uniq > {args.o}/analysis_contig/mmseq/rep/{rawdata['SampleID']}.cluster.txt ; "
    job3_taxid = f"pixi run find_identical {args.o}/analysis_contig/mmseq/rep/{rawdata['SampleID']}.cluster.txt {args.o}/analysis_contig/mmseq/tax_out/{rawdata['SampleID']}.tax.txt -format 1 -site1 1 -site2 1 -type 3 > {args.o}/analysis_contig/mmseq/tax_rep/{rawdata['SampleID']}.rep.tax.txt ; "
    job3_repfa = f"pixi run seqtk subseq {rawdata['Contig']} {args.o}/analysis_contig/mmseq/rep/{rawdata['SampleID']}.cluster.txt > {args.o}/analysis_contig/mmseq/rep/{rawdata['SampleID']}.cluster.fa ; "

    job3_mmseq = "\n#### extract prokaryotes ####\n"
    job3 = job3_mmseq + job3_cluster.to_string(index=False, header=False) + "\n" + job3_createdb.to_string(index=False, header=False) + "\n" + job3_taxonomy.to_string(index=False, header=False) + "\n" + job3_createtsv.to_string(index=False, header=False) + "\n" + job3_rep.to_string(index=False, header=False) + "\n" + job3_taxid.to_string(index=False, header=False) + "\n" + job3_repfa.to_string(index=False, header=False) + "\n\n"

    job4_covermContig = f"pixi run coverm contig --coupled {rawdata['Cleandata_R1']} {rawdata['Cleandata_R2']} --reference {rawdata['Contig']} {parameter.iloc[8,1]} -o {args.o}/analysis_contig/coverm/{rawdata['SampleID']}.coverm.txt --no-zeros -t {args.t} ; "
    job4_ContigID = f"pixi run merge.tax_coverm -t {args.o}/analysis_contig/mmseq/tax_rep/{rawdata['SampleID']}.rep.tax.txt -c {args.o}/analysis_contig/coverm/{rawdata['SampleID']}.coverm.txt -o {args.o}/analysis_contig/mmseq/tax_coverm/{rawdata['SampleID']}.tax.coverm.txt ; "
    job4_bacteria_merger = f"pixi run bacteria_merger_multi -d {args.o}/analysis_contig/mmseq/tax_coverm/ -o {args.o}/LorPha_results/Coverm ; "
    job4_covermVotu = f"pixi run coverm contig --coupled {rawdata['Cleandata_R1']} {rawdata['Cleandata_R2']} --reference {args.o}/analysis_votu/votu.fa {parameter.iloc[8,1]} -o {args.o}/analysis_votu/coverm/{rawdata['SampleID']}.votu.coverm.txt --no-zeros -t {args.t} ; "
    job4_votu_merger = f"pixi run votu_merger_multi -d {args.o}/analysis_votu/coverm/ -o {args.o}/LorPha_results/Coverm ; "

    job4_coverm = "\n#### calculate abundances ####\n"
    job4 = job4_coverm + job4_covermContig.to_string(index=False, header=False) + "\n" + job4_ContigID.to_string(index=False, header=False) + "\n" + job4_bacteria_merger + "\n" + job4_covermVotu.to_string(index=False, header=False) + "\n" + job4_votu_merger + "\n\n"

    job5_dbcan1 = f"pixi run --environment dbcan3 run_dbcan easy_substrate --db_dir {args.d}/dbcan --output_dir {args.o}/analysis_votu/dbcan --mode meta --input_raw_data {args.o}/analysis_votu/votu.fa --threads {args.t} {parameter.iloc[5,1]} ; "
    job5_dbcan2 = f"cp {args.o}/analysis_votu/dbcan/overview.tsv {args.o}/LorPha_results/GeneFunction/dbcan.txt ; "
    job5_PhaStyle1 = f"pixi run PhaStyle --fastain {args.o}/analysis_votu/votu.fa --out {args.o}/analysis_votu/PhaStyle/PhaStyle.txt --ftmodel {args.d}/PhaStyle --num-cores {args.t} ; "
    job5_PhaStyle2 = f"pixi run merge.votu_profile.phastyle -v {args.o}/analysis_votu/votu.profile.txt -p {args.o}/analysis_votu/PhaStyle/PhaStyle.txt -o {args.o}/LorPha_results ;"
    job5_PhaStyle3 = f"pixi run filter.votu_meta -i {args.o}/LorPha_results/votu.meta.raw.txt -o {args.o}/LorPha_results ;"
    job5_vibrant1 = f"pixi run --environment vibrant VIBRANT_run.py -i {args.o}/analysis_votu/votu.fa -no_plot -virome -t {args.t} -folder {args.o}/analysis_votu/VIBRANT -d {args.d}/vibrant ; "
    job5_vibrant2 = f"cp {args.o}/analysis_votu/VIBRANT/VIBRANT_votu/VIBRANT_results_votu/VIBRANT_annotations_votu.tsv {args.o}/LorPha_results/GeneFunction/VIBRANT_annotation.txt ; "
    job5_vibrant3 = f"cp {args.o}/analysis_votu/VIBRANT/VIBRANT_votu/VIBRANT_results_votu/VIBRANT_AMG_individuals_votu.tsv {args.o}/LorPha_results/GeneFunction/VIBRANT_potentialAMG.txt ; "
    job5_defensefinder1 = f"pixi run --environment defensefinder defense-finder run {args.o}/analysis_contig/mmseq/rep/{rawdata['SampleID']}.cluster.fa -o {args.o}/analysis_contig/defensefinder --models-dir {args.d}/defensefinder/ -w {args.t} ; "
    job5_defensefinder2 = f"pixi run merge.defensefinder_tax -d {args.o}/analysis_contig/defensefinder -t {args.o}/analysis_contig/mmseq/tax_out/ -o {args.o}/LorPha_results/GeneFunction ; "

    job5_function = "\n#### PhaStyle & gene functions ####\n"
    job5 = job5_function + job5_dbcan1 + "\n" + job5_dbcan2 + "\n" + job5_PhaStyle1 + "\n" + job5_PhaStyle2 + "\n" + job5_PhaStyle3 + "\n" + job5_vibrant1 + "\n" + job5_vibrant2 + "\n" + job5_vibrant3 + "\n" + job5_defensefinder1.to_string(index=False, header=False) + "\n" + job5_defensefinder2 + "\n\n"

    job_all = job1 + job2 + job3 + job4 + job5
    out_job = f"{args.o}/job.sh"

    with open(out_job, 'w') as f:
        f.write(job_all)

if __name__ == '__main__':
    main()
