####### This is the repository for code used for alignment and variant detection for APL sequencing project.

###### Parameters used for alignment and variant detection are hardcoded into the script. May have to edit accordingly.

#### Alignment:
`fqtobam.py sample_id sample_name sample_library fq1.gz fq2.gz`

#### Variant calling using VarScan2:
`bash run_varScan.sh normal.bam tumor.bam basename`