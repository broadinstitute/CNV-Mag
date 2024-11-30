version 1.0

workflow CNV_Mag {
    input{
        String sampleName
        String cnvProfiler_Docker = "us.gcr.io/tag-public/covprofileviz:0.0.4"
        File cramOrBamFile
        File cramOrBamIndexFile
        File genomeBoundaryFile = "gs://fc-325cb421-bf1a-4e99-b50c-3f785d6b994a/genomeBoundries/Homo_sapiens_assembly38.genome"
        File? cnvBedFile
        Array[String]? cnvIntervals
        Boolean heterozygosityCheck = false
        File? hardFilteredVcfFile
    }

    if (defined(cnvIntervals)) {
        call CreateBedFromIntervals {
            input:
                cnvIntervals = cnvIntervals,
                cnvProfiler_Docker = cnvProfiler_Docker
        }
    }

    File cnvBedFile = select_first([cnvBedFile, CreateBedFromIntervals.output_bed])

    call GetPaddedCnvBed {
        input:
            cnvBedFile = cnvBedFile,
            genomeBoundaryFile = genomeBoundaryFile
    }
    call SamtoolsDepth {
        input:
            sampleName = sampleName,
            alignedBam = cramOrBamFile,
            alignedBai = cramOrBamIndexFile,
            target_bed = GetPaddedCnvBed.paddedCnvBed

    }
#    call cnvDepthProfiler {
#        input:
#            sampleName = sampleName,
#            depthProfile = SamtoolsDepth.depth_profile,
#            cnvBedFile = cnvBedFile,
#            cnvProfiler_Docker = cnvProfiler_Docker,
#            PaddedcnvBedFile = GetPaddedCnvBed.paddedCnvBed
#    }
    if (heterozygosityCheck) {
        call HeterozygosityCheck {
            input:
                sampleName = sampleName,
                hardFilteredVcfFile = hardFilteredVcfFile,
                cnvBedFile = cnvBedFile,
                cnvProfiler_Docker = cnvProfiler_Docker
        }
    }
    output {
    File mapq0_depth_profile = SamtoolsDepth.mapq0_depth_profile
    File mapq20_depth_profile = SamtoolsDepth.mapq20_depth_profile
    Array[File]? heterozygosity_plot = HeterozygosityCheck.heterozygosity_plot
    }
    meta {
        description: "This workflow takes a BAM or CRAM file and a CNV bed file as input and generates a coverage profile for the CNV regions in the bed file. Optionally, it can also generate a heterozygosity plot using a hard-filtered VCF file."
        author: "Yueyao Gao"
        email: "tag@broadinstitute.org"
        }
}


task CreateBedFromIntervals {
    input {
        Array[String]? cnvIntervals
        String cnvProfiler_Docker
        Int mem_gb = 1
        Int cpu = 1
        Int disk_size_gb = 10
    }
    command <<<
        # Write the CNV intervals to a file
        cnvIntervals=(~{sep=" " cnvIntervals})
        for interval in "${cnvIntervals[@]}"; do
            echo $interval >> cnv_intervals.txt
        done

        # Create a bed file from the CNV intervals
        source activate env_viz
        python3 <<CODE

        with open('cnv_intervals.txt', 'r') as f:
                cnvIntervals = f.readlines()

        with open('cnv_intervals.bed', 'a') as f:
            for interval in cnvIntervals:
                chr = interval.split(':')[0]
                start = interval.split(':')[1].split('-')[0]
                end = interval.split(':')[1].split('-')[1]
                f.write(f"{chr}\t{start}\t{end}" + '\n')
        CODE
    >>>
    runtime {
        docker: cnvProfiler_Docker
        cpu: cpu
        memory: mem_gb + " GB"
        disks: "local-disk " + disk_size_gb + " HDD"
    }
    output {
        File output_bed = "cnv_intervals.bed"
    }
}

task GetPaddedCnvBed {
    input {
        File cnvBedFile
        File genomeBoundaryFile
        String bedtools_docker = "us.gcr.io/broad-dsde-methods/vcfeval_docker:v1.0"
        Int mem_gb = 4
        Int cpu = 1
        Int disk_size_gb = 10
    }

    command <<<
        # Create a padded CNV bed file
        # Extend the CNV regions by 200% of the interval size on each side
        bedtools slop -i ~{cnvBedFile} -g ~{genomeBoundaryFile} -b 2 -pct > padded_cnv.bed
    >>>
    runtime {
        docker: bedtools_docker
        cpu: cpu
        memory: mem_gb + " GB"
        disks: "local-disk " + disk_size_gb + " HDD"
    }
    output {
        File paddedCnvBed = "padded_cnv.bed"
    }
}

task SamtoolsDepth {
        input {
            String sampleName
            File alignedBam
            File alignedBai
            File target_bed
            File HG1bam = "gs://fc-a76d0374-93e7-4c1a-8302-2a88079b480d/DRAGEN_4.3.6_NIST_default/NA12878_HG001_1000ng_3_NVX/NA12878_HG001_1000ng_3_NVX.cram"
            File HG1bai = "gs://fc-a76d0374-93e7-4c1a-8302-2a88079b480d/DRAGEN_4.3.6_NIST_default/NA12878_HG001_1000ng_3_NVX/NA12878_HG001_1000ng_3_NVX.cram.crai"
            File HG2bam = "gs://fc-a76d0374-93e7-4c1a-8302-2a88079b480d/DRAGEN_4.3.6_NIST_default/NA24385_HG002_1_NVX/NA24385_HG002_1_NVX.cram"
            File HG2bai = "gs://fc-a76d0374-93e7-4c1a-8302-2a88079b480d/DRAGEN_4.3.6_NIST_default/NA24385_HG002_1_NVX/NA24385_HG002_1_NVX.cram.crai"
            Int mem_gb = 64
            Int cpu = 8
            Int disk_size_gb = 500
            Boolean use_ssd = true
            String samtools_docker = "euformatics/samtools:1.20"
    }
    command <<<
        # Create output directory
        mkdir output

        # Run samtools depth to get MAPQ20 depth
        # Counting fragments instead of reads using -s option
        for mq in 0 20; do
            samtools depth \
            -@ ~{cpu} \
            -b ~{target_bed} \
            --min-BQ 20 \
            --min-MQ ${mq} \
            -s \
            ~{alignedBam} \
            ~{HG1bam} \
            ~{HG2bam} \
            -o output/~{sampleName}_MAPQ${mq}_samtools.depth;
        done

    >>>
    output {
        File mapq0_depth_profile = "output/~{sampleName}_MAPQ0_samtools.depth"
        File mapq20_depth_profile = "output/~{sampleName}_MAPQ20_samtools.depth"
    }
    runtime {
        memory: mem_gb * 1000 + " MB"
        cpu: cpu
        docker: samtools_docker
        disks: "local-disk " + disk_size_gb + if use_ssd then " SSD" else " HDD"
        preemptible: 0
        maxRetries: 3
    }
}

#task cnvDepthProfiler{
#    input {
#            String sampleName
#            String cnvProfiler_Docker
#            File depthProfile
#            File cnvBedFile
#            File PaddedcnvBedFile
#            Int intervalPadding = 0
#            Int mem_gb = 64
#            Int cpu = 8
#            Int preemptible = 0
#            Int disk_size_gb = 500
#            Int maxRetries = 1
#            Boolean use_ssd = true
#        }
#        command <<<
#            set -e
#            mkdir output
#
#            # Run the coverage profile visualization script
#            conda run --no-capture-output \
#            -n env_viz \
#            python3 /BaseImage/CovProfileViz/scripts/CNV_Depth_Profiler.py \
#            -c ~{depthProfile} \
#            -b ~{cnvBedFile} \
#            -n ~{sampleName} \
#            -pb ~{PaddedcnvBedFile} \
#            -p ~{intervalPadding} \
#            -o output
#
#        >>>
#        output {
#            Array[File] cnv_depth_profile = glob("output/*png")
#        }
#        runtime {
#            memory: mem_gb + " GB"
#            cpu: cpu
#            docker: cnvProfiler_Docker
#            disks: "local-disk " + disk_size_gb + if use_ssd then " SSD" else " HDD"
#            preemptible: preemptible
#            maxRetries: maxRetries
#        }
#}

task HeterozygosityCheck{
    input {
        String sampleName
        String cnvProfiler_Docker
        File HG1_vcf_path = "gs://fc-a76d0374-93e7-4c1a-8302-2a88079b480d/DRAGEN_4.3.6_NIST_default/NA12878_HG001_1000ng_3_NVX/NA12878_HG001_1000ng_3_NVX.hard-filtered.vcf.gz"
        File HG2_vcf_path = "gs://fc-a76d0374-93e7-4c1a-8302-2a88079b480d/DRAGEN_4.3.6_NIST_default/NA24385_HG002_1_NVX/NA24385_HG002_1_NVX.hard-filtered.vcf.gz"
        File? hardFilteredVcfFile
        File cnvBedFile
        Int mem_gb = 64
        Int cpu = 8
        Int preemptible = 0
        Int disk_size_gb = 500
        Int maxRetries = 1
        Boolean use_ssd = true
    }
    command <<<
        set -e
        mkdir output

        # Run the coverage profile visualization script
        conda run --no-capture-output \
        -n env_viz \
        python3 /BaseImage/CovProfileViz/scripts/CNV_SNP_HET_Profiler.py \
        -v1 ~{hardFilteredVcfFile} \
        -v2 ~{HG1_vcf_path} \
        -v3 ~{HG2_vcf_path} \
        -b ~{cnvBedFile} \
        -n1 ~{sampleName} \
        -n2 HG001 \
        -n3 HG002  \
        -o output

    >>>
    output {
        Array[File] heterozygosity_plot = glob("output/*png")
    }
    runtime {
        memory: mem_gb + " GB"
        cpu: cpu
        docker: cnvProfiler_Docker
        disks: "local-disk " + disk_size_gb + if use_ssd then " SSD" else " HDD"
        preemptible: preemptible
        maxRetries: maxRetries
    }
}







