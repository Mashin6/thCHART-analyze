# thCHART-analyze
Simple pipeline for processing of sequencing data from thCHART, ChIP-seq, APEX-ChIP and similar experiments


This is a script to process .fastq sequencing files from thCHART, CHART or ChIP-seq experiments and genereate fold-enrichment tracks in bigWig format. 

Script must be started from folder with all .fastq files. Sample names must match filename format: `NAME_R1.fastq NAME_R2.fastq`

Usage: `chart-analyze.sh -g [genome version] -r [read type] [-u] -i [input[,input2]] -s [sample1[,sample2...]]`  
e.g.   `chart-analyze.sh -g dm6 -r PE -i MM200201_01 -s MM200201_02,MM200201_03,MM200201_04`


```
parameters:      -s|--samples           - comma-separated list of sample name(s)
                 -i|--input             - input sample name
                 -g|--genome            - name of genome to use [dm6, dm3, hg38, hg19, hg18, mm9, mm10]
                 -r|--reads             - type of reads [PE, SE]
                 -u|--uniqalign         - keep only uniquely aligned reads by applying mapping quality filter in samtools -q 2
                 -h|--help              - this message 
```
