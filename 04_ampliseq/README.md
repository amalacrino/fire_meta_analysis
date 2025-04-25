# ampliseq

The code below runsr `nf-core/ampliseq` and requires:
* `mafft`
* `FastTree`

`$INDIR` is the path to the directory where data was stored after running `cutadapt`.

`$OUTDIR` is the path to the output directory.

`$REFDIR` is the path to the directory where databases for taxonomical assignment are stored (i.e., SILVA and UNITE).

## Run ampliseq for 16S

```bash
nextflow run nf-core/ampliseq -r 2.10.0 -profile singularity  \
--input_folder $INDIR \
--outdir $OUTDIR \
--extension "/*_{1,2}.fastq.gz" \
--multiple_sequencing_runs \
--dada_ref_tax_custom $REFDIR/silva_nr99_v138.1_wSpecies_train_set.fa.gz \
--dada_ref_tax_custom_sp $REFDIR/silva_species_assignment_v138.1.fa.gz \
--trunclenf 250 \
--trunclenr 250 \
--skip_dada_quality \
--skip_cutadapt \
--skip_qiime \
--skip_barplot \
--skip_abundance_tables \
--skip_alpha_rarefaction \
--skip_diversity_indices \
--skip_ancom \
--ignore_empty_input_files \
--ignore_failed_trimming \
--ignore_failed_filtering

cd $OUTDIR

mafft --thread $NTHREADS ASV_seqs.fasta > asv_aligned.fasta

FastTree -gtr -nt < asv_aligned.fasta > tree.tre
```

## Run ampliseq for ITS

```bash
nextflow run nf-core/ampliseq -r 2.10.0 -profile singularity  \
--input_folder $INDIR \
--outdir $OUTDIR \
--extension "/*_{1,2}.fastq.gz" \
--multiple_sequencing_runs \
--illumina_pe_its \
--dada_ref_tax_custom $REFDIR/sh_general_release_dynamic_s_19.02.2025.fasta \
--skip_dada_addspecies \
--trunclenf 250 \
--trunclenr 250 \
--skip_dada_quality \
--skip_cutadapt \
--skip_qiime \
--skip_barplot \
--skip_abundance_tables \
--skip_alpha_rarefaction \
--skip_diversity_indices \
--skip_ancom \
--ignore_empty_input_files \
--ignore_failed_trimming \
--ignore_failed_filtering

cd $OUTDIR

mafft --thread $NTHREADS ASV_seqs.fasta > asv_aligned.fasta

FastTree -gtr -nt < asv_aligned.fasta > tree.tre
```
