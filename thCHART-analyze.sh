#!/bin/bash
#    File name: thCHART-analyze.sh
#    Author: Martin Machyna
#    Email: machyna@gmail.com
#    Date created: 10/9/2020
#    Date last modified: 16/9/2021
#    Version: 1.0.3
#    License: GPLv3

# Pipeline for creating fold enrichemnt tracks from thCHART genomic sequencing data 
# 
# Folder must contain sequencing data in format NAME_R1.fastq NAME_R2.fastq
# Following software must be installed an accesible form PATH: 
#   Cutadapt >1.7, Bowtie 2.2.9, MACS 2.1.0, SAMtools 1.4, BEDTools 2.26.0, GNU parallel
#   bedGraphToBigWig (http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/bedGraphToBigWig)


#SBATCH --partition=general
#SBATCH --job-name=thCHART-analyze
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=119G 
#SBATCH --mail-type=ALL
#SBATCH --mail-user=user@email.com


function PrintUsage {
    echo "Script for analyzing thCHART, CHART or ChIP-seq data.
Script must be started from folder with all .fastq files. Sample names must match filename format: NAME_R1.fastq NAME_R2.fastq

Usage: thCHART-analyze.sh -g [genome version] -r [read type] [-u] -i [input[,input2]] -s [sample1[,sample2...]]
e.g.   thCHART-analyze.sh -g dm6 -r PE -i MM200201_01 -s MM200201_02,MM200201_03,MM200201_04
    "

    echo "parameters:      -s|--samples           - comma-separated list of sample name(s)
                 -i|--input             - input sample name or comma-separated list of input names
                 -g|--genome            - name of genome to use [dm6, dm3, hg38, hg19, hg18, mm9, mm10] (default: dm6)
                 -r|--reads             - type of reads [PE, SE]
                 -u|--uniqalign         - keep only uniquely aligned reads by applying mapping quality filter in samtools -q 2 (default: off)
                 -h|--help              - this message "

    exit
}

# Set common path
    fastuniq_dir='/home/path/to/fastuniq'
    bdgToBw_dir='/home/path/to/bedtools'



# Check if parameters are set
    if [ -z $1 ]; then
        PrintUsage
        exit
    fi

    # Read parameters
    POSITIONAL=()
    while [[ $# -gt 0 ]]; do
        case $1 in

            -s|--samples)
                samples=($(echo "$2" | sed 's/,/\n/g'))
                shift # past argument
                shift # past value
            ;;
            -i|--input)
                inputs=($(echo "$2" | sed 's/,/\n/g'))
                shift # past argument
                shift # past value
            ;;
            -g|--genome)
                genome="$2"
                shift # past argument
                shift # past value
            ;;
            -r|--reads)
                if [[ $2 = "PE" ]] || [[ $2 = "SE" ]]; then reads="$2"
                else
                    echo "Error: Read type must be PE or SE only"
                    PrintUsage
                    exit 1
                fi
                shift # past argument
                shift # past value
            ;;
            -u|--uniqalign)
                unique="TRUE"
                shift
            ;;
            -h|--help)
                PrintUsage
                exit 1
            ;;
            *)    # unknown option
                POSITIONAL+=("$1") # save it in an array for later
                shift # past argument
            ;;
        esac
    done
    set -- "${POSITIONAL[@]}" # restore positional parameters

# Defaults
unique=${unique:-"FALSE"}
genome=${genome:-"dm6"}

# Check
    if [ ${#inputs[@]} -eq 0 ]; then echo "You must specify at least 1 input dataset"; PrintUsage; exit 1; fi
    if [ ${#samples[@]} -eq 0 ]; then echo "You must specify at least 1 sample"; PrintUsage; exit 1; fi

# Set genome specific parameters
# Extract genome type
    if [ $genome = "dm6" ]; then
        bowtie_index="/home/mm2594/software/genomes/bowtie2/dm6"
        # chrom_sizes="/home/mm2594/software/bedtools/dm6.chrom.sizes"
        macs_gen="dm"
    elif [ $genome = "dm3" ]; then
        bowtie_index="/home/mm2594/software/genomes/bowtie2/dm3"
        # chrom_sizes="/home/mm2594/software/bedtools/dm3.chrom.sizes"
        macs_gen="dm"
    elif [ $genome = "mm9" ]; then
        bowtie_index="/home/mm2594/software/genomes/bowtie2/mm9"
        # chrom_sizes="/home/mm2594/software/bedtools/mm9.chrom.sizes"
        macs_gen="mm"
    elif [ $genome = "mm10" ]; then
        bowtie_index="/home/mm2594/software/genomes/bowtie2/mm10"
        # chrom_sizes="/home/mm2594/software/bedtools/mm10.chrom.sizes"
        macs_gen="mm"
    elif [ $genome = "hg38" ]; then
        bowtie_index="/home/mm2594/software/genomes/bowtie2/GRCh38"
        #chrom_sizes="/home/mm2594/software/bedtools/hg38.chrom.sizes"
        chrom_sizes=
        macs_gen="hs"
    elif [ $genome = "hg19" ]; then
        bowtie_index="/home/mm2594/software/genomes/bowtie2/hg19"
        # chrom_sizes="/home/mm2594/software/bedtools/hg19.chrom.sizes"
        macs_gen="hs"
    elif [ $genome = "hg18" ]; then
        bowtie_index="/home/mm2594/software/genomes/bowtie2/hg18"
        # chrom_sizes="/home/mm2594/software/bedtools/hg18.chrom.sizes"s
        macs_gen="hs"
    else
        echo "Unknown genome specified"; PrintUsage; exit 1;
    fi

echo "*** Pipeline started"
echo $(date)
echo
echo "Parameters used:"
echo "                  samples:" $(echo ${samples[@]} | sed 's| |,|g')
echo "                  input:"  $(echo ${inputs[@]} | sed 's| |,|g')
echo "                  unique: ${unique}"
echo "                  genome: ${genome}"
echo "                  read-type: ${reads}"


### Data processing

    # ### Filtering for duplicate reads   
    #     echo "* Running Fastuniq"
    #     for sample in ${samples[@]} ${inputs[@]}; do
    #         ${fastuniq_dir}/fastuniq \
    #             -i <(echo "$sample"_R1.fastq; echo "$sample"_R2.fastq) \
    #             -o "$sample"_R1.u.fastq \
    #             -p "$sample"_R2.u.fastq
    #     done

    module load parallel/20210222-GCCcore-10.2.0
       
    ### Trim Illumina adaptor sequences with cutadapt
        module load cutadapt/3.2-GCCcore-10.2.0-Python-3.8.6

        if [ $reads = "PE" ]; then

            for sample in ${samples[@]} ${inputs[@]}; do

                echo "* Running cutadapt in PE reads mode for sample $sample"

                cutadapt \
                        -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC \
                        -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT \
                        -m 3 \
                        --cores=20 \
                        -o "$sample"_R1.t.fastq \
                        -p "$sample"_R2.t.fastq \
                        "$sample"_R1.fastq "$sample"_R2.fastq

                #rm "$sample"_R{1,2}.u.fastq
            done

        else 

            for sample in ${samples[@]} ${inputs[@]}; do

        
                echo "* Running cutadapt in SE reads mode for sample $sample"
       
                cutadapt \
                        -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC \
                        -m 3 \
                        --cores=20 \
                        -o "$sample"_R1.t.fastq \
                        "$sample"_R1.fastq
            done

        fi




    ### Align wiht Bowtie 2
        module load Bowtie2/2.4.2-GCCcore-10.2.0

        if [ $reads = "PE" ]; then
            for sample in ${samples[@]} ${inputs[@]}; do

                echo "* Running alignment with Bowtie2 for sample $sample"

                bowtie2 -p 20 \
                        -x "$bowtie_index" \
                        -1 "$sample"_R1.t.fastq \
                        -2 "$sample"_R2.t.fastq \
                        -S "$sample".sam

                rm "$sample"_R{1,2}.t.fastq 
            done


            # Set parameter for MACS2
            macs_param="BAMPE"

        else 

            for sample in ${samples[@]} ${inputs[@]}; do

                echo "* Running alignment with Bowtie2 for sample $sample"
                
                bowtie2 -p 20 \
                        -x "$bowtie_index" \
                        -U "$sample"_R1.t.fastq \
                        -S "$sample".sam

                rm "$sample"_R1.t.fastq 
            done


            # Set parameter for MACS2
            macs_param="BAM"

        fi



    ### Filter for non-unique alignments and convert to .bam 
        module load SAMtools/1.12-GCCcore-10.2.0

        for sample in ${samples[@]} ${inputs[@]}; do

            echo "* Converting .sam to .bam for sample $sample"

            samtools view -b \
                $( if [ $unique = 'TRUE' ]; then echo "-q 2"; fi ) \
                -@ 20 \
                -o "$sample".bam "$sample".sam
                
            # sort BAM
            samtools sort -@ 20 "$sample".bam -o "$sample"_sort.bam
        
            # remove .sam and old .bam
            rm "$sample".sam
            rm "$sample".bam
            
            # rename sorted bam
            mv "$sample"_sort.bam "$sample".bam
            

        done    



    ### Call peaks with MACS2
        module load MACS2/2.2.7.1-foss-2020b-Python-3.8.6

        echo "* Calling peaks with MACS2"
            
        parallel -j 20  macs2 callpeak \
                            -t {1}.bam \
                            -c ${inputs[@]/%/.bam} \
                            -n {1} \
                            -f "$macs_param" \
                            -g "$macs_gen" \
                            -B \
                            ::: ${samples[@]}



    ### Create fold enrichment track
        echo "* Creating fold enrichemnt tracks"
  
        parallel -j 20 macs2 bdgcmp \
                        -t {1}_treat_pileup.bdg \
                        -c {1}_control_lambda.bdg \
                        -o {1}_FE.bdg \
                        -m FE \
                        ::: ${samples[@]}


    ### Sort BedGraph files uppercase letter before lowercase
        echo "* Sorting bedgraph files"

        parallel -j 20 "LC_COLLATE=C sort -k1,1 -k2,2n {1}_FE.bdg > {1}_FE_sort.bdg" ::: ${samples[@]}


    
    ### Convert .bedgraph to .bigWig
        echo "* Converting .bedgraoh to .bigwig"

        if [ -z $chrom_sizes ]; then
                samtools view -H ${samples}.bam \
                    | awk -v OFS="\t" ' $1 ~ /^@SQ/ {split($2, chr, ":")
                                                     split($3, size, ":")
                                                     print chr[2], size[2]}' > "$genome".chrom.sizes

                chrom_sizes="$genome".chrom.sizes
        fi

            
        parallel rm {1}_FE.bdg ::: ${samples[@]}
        parallel mv {1}_FE_sort.bdg {1}_FE.bdg ::: ${samples[@]}
        

        parallel -j 20 ${bdgToBw_dir}/bedGraphToBigWig \
                                            {1}_FE.bdg \
                                            ${chrom_sizes} \
                                            {1}_FE.bigWig \
                                            ::: ${samples[@]}

 
        # Remove bedgraphs
            rm *_FE.bdg
            rm *_treat_pileup.bdg
            rm *_treat_pileup_sort.bdg
            rm *_control_lambda.bdg

        
    ### Create raw pileup coverage tracks
        module load BEDTools/2.30.0-GCCcore-10.2.0

        echo "* Creating pileup coverage tracks"
        # Create .bedgraph coverage file from .bam file
                
        parallel -j 20 'bedtools genomecov -bg \
                                           $( if [ "$reads" = "PE" ]; then echo "-pc"; fi ) \
                                           -ibam {1}.bam  > {1}_pileup.bdg' ::: ${samples[@]} ${inputs[@]}


        # Sort
        parallel -j 20 "LC_COLLATE=C sort -k1,1 -k2,2n {1}_pileup.bdg > {1}_pileup_sort.bdg" ::: ${samples[@]} ${inputs[@]}

        # Convert to bigWig

        parallel rm {1}_pileup.bdg ::: ${samples[@]} ${inputs[@]}
        parallel mv {1}_pileup_sort.bdg {1}_pileup.bdg ::: ${samples[@]} ${inputs[@]}


        parallel -j 20 ${bdgToBw_dir}/bedGraphToBigWig \
                                        {1}_pileup.bdg  \
                                        ${chrom_sizes} \
                                        {1}_pileup.bigWig \
                                        ::: ${samples[@]} ${inputs[@]}


        # Remove bedgraphs
            rm *_pileup.bdg

echo '*** Pipeline finished'
printf 'Runtime: %02dh:%02dm:%02ds\n' $(($SECONDS/3600)) $(($SECONDS%3600/60)) $(($SECONDS%60))




