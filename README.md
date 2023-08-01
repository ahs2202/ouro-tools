# Scarab-Count (Single-cell Catch-All Analysis Toolkit)
 a python package for analyzing every genomic regions of single-cell RNA/ATAC-seq data. 


![scarab-logo](doc/img/scarab_logo.png)



### Installation 

Currently, Scarab-Count is <u>***unpublished***</u>, and maintained in the private repo. 

In order to install Scarab-Count, run the following commands in bash shell.

```bash
git clone https://github.com/ahs2202/scarab-count.git
cd scarab-count
pip install .
```



scarab-count can be used in command line, in a Python script, or in an interactive Python interpreter.

##### Bash shell

To print the command line usage from the bash shell, please type the following command.

```bash
scarab_count count -h
```



##### Python script

```python
#run scarab-count v0.2.2
import scarab_count # import scarab-count
if __name__ == '__main__' : # protect the entry point (important!)
    # retrieve path of the input bam files
    df_file = scarab_count.pd.read_csv( '/path/to/tsv/file/containing/bam/file/paths', sep = '\t' )
    # retrieve bam files and compose a list of output folders
    l_path_file_bam_input = df_file.path.values
    l_path_folder_output = [ f'{path_file}.scarab-count.{scarab_count.__version__}.output/' for path_file in l_path_file_bam_input ]

    # run scarab-count
    scarab_count.count( 
        path_folder_ref = '/path/where/reference/will/be/saved/and/loaded/Homo_sapiens.GRCh38.105.v0.2.2/', 
        l_path_file_bam_input = l_path_file_bam_input,
        l_path_folder_output = l_path_folder_output,
        l_str_mode_scarab_count = [ 'gex' ],
    )
```

##### IPython environment

Detailed arguments and the docstring can be printed in IPython environment (Jupyter Notebook).

```python
scarab_count.count?
```



### Building an index for scarab-count

Scarab-Count utilizes <u>genome, transcriptome, and gene annotations</u> to assign reads to **genes, isoforms, and genomic bins (tiles across the genome)**. The index building process is automatic; <u>there is no needs to run a separate command in order to build the index</u>. Once Scarab-Count processes these information before analyzing an input BAM file(s), the program saves an index in order to load the information much faster next time.

We recommends using <u>***Ensembl*** reference genome, transcriptome, and gene annotations of the same version</u> (release number). 



##### Pre-build index files 

pre-built index can be downloaded using the following links (should be extracted to a folder using **tar -xf** command):

<u>*Human (GRCh38, Ensembl version 105)*</u> : https://www.dropbox.com/s/8agizrykiorpnag/Homo_sapiens.GRCh38.105.v0.1.1.tar?dl=0

*<u>Mouse (GRCm38, Ensembl version 102)</u>* : https://www.dropbox.com/s/ehen2xg0bmc573g/Mus_musculus.GRCm38.102.v0.1.1.tar?dl=0



##### Building index from scratch

A Scarab-Count index can be built on-the-fly from the input genome, transcriptome, and gene annotation files. For example, below are the list of files that were used for the pre-built Scarab-Count index "<u>*[Human (GRCh38, Ensembl version 105)](https://www.dropbox.com/s/8agizrykiorpnag/Homo_sapiens.GRCh38.105.v0.1.1.tar?dl=0)*</u>".



*required annotations* (*Ensemble version 105*):

* **path_file_fa_genome** : https://ftp.ensembl.org/pub/release-105/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
  * A genome FASTA file. Either gzipped or plain FASTA file can be accepted.
* **path_file_gtf_genome** :  https://ftp.ensembl.org/pub/release-105/gtf/homo_sapiens/Homo_sapiens.GRCh38.105.gtf.gz
  * A GTF file. Either gzipped or plain GTF file can be accepted. Currently GFF3 format files are not supported.
  * Following arguments can be used to set attribute names for identifying gene and transcript annotations in its attributes column.
    * str_name_gtf_attr_for_id_gene : (default: '**gene_id**')
    * str_name_gtf_attr_for_name_gene : (default: '**gene_name**')
    * str_name_gtf_attr_for_id_transcript : (default: '**transcript_id**')
    * str_name_gtf_attr_for_name_transcript : (default: '**transcript_name**')
  * An example of GTF annotation file for gene annotations:

```
1	ensembl_havana	gene	1211340	1214153	.	-	.	gene_id "ENSG00000186827"; gene_version "11"; gene_name "TNFRSF4"; gene_source "ensembl_havana"; gene_biotype "protein_coding";
1	ensembl_havana	transcript	1211340	1214153	.	-	.	gene_id "ENSG00000186827"; gene_version "11"; transcript_id "ENST00000379236"; transcript_version "4"; gene_name "TNFRSF4"; gene_source "ensembl_havana"; gene_biotype "protein_coding"; transcript_name "TNFRSF4-201"; transcript_source "ensembl_havana"; transcript_biotype "protein_coding"; tag "CCDS"; ccds_id "CCDS11"; tag "basic"; transcript_support_level "1 (assigned to previous version 3)";
1	ensembl_havana	exon	1213983	1214153	.	-	.	gene_id "ENSG00000186827"; gene_version "11"; transcript_id "ENST00000379236"; transcript_version "4"; exon_number "1"; gene_name "TNFRSF4"; gene_source "ensembl_havana"; gene_biotype "protein_coding"; transcript_name "TNFRSF4-201"; transcript_source "ensembl_havana"; transcript_biotype "protein_coding"; tag "CCDS"; ccds_id "CCDS11"; exon_id "ENSE00001832731"; exon_version "2"; tag "basic"; transcript_support_level "1 (assigned to previous version 3)";
```

* **path_file_fa_transcriptome** : https://ftp.ensembl.org/pub/release-105/fasta/homo_sapiens/cdna/Homo_sapiens.GRCh38.cdna.all.fa.gz
  * A transcriptome FASTA file. Either gzipped or plain FASTA file can be accepted.



*optional annotations:*

* **path_file_tsv_repeatmasker_ucsc** : [Table Browser (ucsc.edu)](https://genome.ucsc.edu/cgi-bin/hgTables?hgsid=1576143313_LetmEyQf9yggiQJAXajCua4TGOGl&clade=mammal&org=Human&db=hg38&hgta_group=rep&hgta_track=knownGene&hgta_table=0&hgta_regionType=genome&position=chr2%3A25%2C160%2C915-25%2C168%2C903&hgta_outputType=primaryTable&hgta_outFileName=GRCh38_RepeatMasker.tsv.gz) [click "get output" to download the annotation]

  * repeat masker annotations from the UCSC Table Browser

* **path_file_gff_regulatory_element** : https://ftp.ensembl.org/pub/current_regulation/homo_sapiens/homo_sapiens.GRCh38.Regulatory_Build.regulatory_features.20221007.gff.gz

  * The latest regulatory build from **Ensembl**.
  * Annotations from other sources, or custom annotations can be used. Currently only the GFF3 file format is supported (with **.gff** extension).
  * The following argument can be used to set the attribute name for identifying regulatory region
    * str_name_gff_attr_id_regulatory_element : (default: '**ID**')

  * An example of GFF annotation file for regulatory elements:

```
18	Regulatory_Build	enhancer	35116801	35120999	.	.	.	ID=enhancer:ENSR00000572865;bound_end=35120999;bound_start=35116801;description=Predicted enhancer region;feature_type=Enhancer
8	Regulatory_Build	TF_binding_site	37967115	37967453	.	.	.	ID=TF_binding_site:ENSR00001137252;bound_end=37967531;bound_start=37966339;description=Transcription factor binding site;feature_typ
6	Regulatory_Build	enhancer	90249202	90257999	.	.	.	ID=enhancer:ENSR00000798348;bound_end=90257999;bound_start=90249202;description=Predicted enhancer region;feature_type=Enhancer
3	Regulatory_Build	CTCF_binding_site	57689401	57689600	.	.	.	ID=CTCF_binding_site:ENSR00000687477;bound_end=57689600;bound_start=57689401;description=CTCF binding site;feature_type=CTCF
```





### Running the Scarab-Count in a Python script (with examples)

```python
#run scarab-count v0.2.2
#pipeline v0.0.20230110
import scarab_count # import scarab-count
if __name__ == '__main__' :
    # setting
    scarab_count.logger.info( f"pipeline started" )

    # set paths of the reference files
    path_folder_ref = '/path/where/reference/will/be/saved/and/loaded/Homo_sapiens.GRCh38.105.v0.2.2/'
    path_fa_for_cram = '/home/shared/cellranger_reference/refdata-gex-GRCh38-2020-A/fasta/genome.fa' 
    path_file_vcf_for_filtering_variant = '/home/shared/database/ncbi/clinvar/clinvar_20221231.vcf.gz'

    # retrieve path of the input bam files
    df_file = scarab_count.pd.read_csv( '/path/to/tsv/file/containing/bam/file/paths', sep = '\t' )

    # retrieve bam files and compose a list of output folders
    l_path_file_bam_input = df_file.path.values
    l_path_folder_output = [ f'{path_file}.scarab-count.{scarab_count.__version__}.output/' for path_file in l_path_file_bam_input ]

    # run scarab-count
    scidx = scarab_count.count( 
        # scidx = scidx, # previously loaded index can be re-used 
        flag_usage_from_command_line_interface = False, 
        path_folder_ref = path_folder_ref, 
        l_path_file_bam_input = l_path_file_bam_input,
        l_path_folder_output = l_path_folder_output,
        l_str_mode_scarab_count = [ 'gex' ],
        int_num_sam_records_for_each_chunk = 100000, 
        int_num_samples_analyzed_concurrently = 4, # starting 4 pipelines, each pipeline processing a single BAM file at a time independently from other pipelines.
        n_threads = 64, # total number of CPUs used for runs. for example. when 'int_num_samples_analyzed_concurrently' is 4 and 'n_threads' is 64, the number of CPUs for analyzing each sample will be 16, and 4 samples will be analyzed concurrently at the same time.
        verbose = True, # for debugging purpose, we recommend to turn on the verbose setting.
        path_file_fa_for_cram = path_fa_for_cram, # if input files are CRAM files, reference FASTA file for the CRAM files should be given as an input.
        path_file_vcf_for_filtering_variant = path_file_vcf_for_filtering_variant,
        flag_skip_read_analysis_summary_output_bam_file = False, # output annotated BAM file
        flag_skip_read_analysis_summary_output_tsv_file = True, # does not generate a TSV file containing the analysis result for each individual read (useful for debugging).
        int_bp_for_bins = 100, # the size of the genomic bin (genomic tiling)
        flag_include_strand_specific_counts = True, # count reads in a strand-specific manner. 
        flag_output_variant_information_with_annotations = True, # assign reads to each annotation with variant information. 
        dict_num_manager_processes_for_each_data_object = { 
            'dict_it_promoter' : 0,
            'dict_t_splice_junc_to_info_genome' : 0,
            'dict_it_exon' : 0,
            'dict_it_exon_transcriptome' : 3, # 4.5 GB
            'dict_it_splice_junc_transcriptome' : 3, # 4.5 GB
            'dict_it_rpmk' : 4, # each 12 GB
            'dict_it_reg' : 5, # 2 GB
            'dict_fa_transcriptome' : 3, # 1.5 GB
        },
        l_seqname_to_skip = [ 'MT' ], # skip mitochondrial sequences
    )


```



---------------

Scarab-Count was developed by Hyunsu An at Gwangju Institute of Science and Technology under the supervision of Professor JiHwan Park. 

© 2022 Functional Genomics Lab, Gwangju Institute of Science and Technology
