version 1.0

workflow CNV_Mag {
    input{
        String sampleName
        String dockerImage = "us.gcr.io/tag-public/cnv-mag:v0.1"
        File cramOrBamFile
        File cramOrBamIndexFile
        String refGenome = "hg38"
        File? cnvBedFile
        Array[String]? cnvIntervals
        File hardFilteredVcfFile
    }

    if (defined(cnvIntervals)) {
        call CreateBedFromIntervals {
            input:
                cnvIntervals = cnvIntervals,
                dockerImage = dockerImage
        }
    }

    File cnvBedFile = select_first([cnvBedFile, CreateBedFromIntervals.output_bed])

    call GetPaddedCnvBed {
        input:
            cnvBedFile = cnvBedFile,
            refGenome = refGenome,
            dockerImage = dockerImage
    }

    call SamtoolsDepth {
        input:
            sampleName = sampleName,
            alignedBam = cramOrBamFile,
            alignedBai = cramOrBamIndexFile,
            target_bed = GetPaddedCnvBed.paddedCnvBed
    }
    call MagDepth {
        input:
            sampleName = sampleName,
            dockerImage = dockerImage,
            mapq20_depth_profile = SamtoolsDepth.mapq20_depth_profile,
            mapq0_depth_profile = SamtoolsDepth.mapq0_depth_profile,
            refGenome = refGenome,
            cnvBedFile = cnvBedFile,
            PaddedcnvBedFile = GetPaddedCnvBed.paddedCnvBed
    }

    call MagSNP {
        input:
            sampleName = sampleName,
            hardFilteredVcfFile = hardFilteredVcfFile,
            cnvBedFile = cnvBedFile,
            dockerImage = dockerImage
    }

    output {
    File mapq0_depth_profile = SamtoolsDepth.mapq0_depth_profile
    File mapq20_depth_profile = SamtoolsDepth.mapq20_depth_profile
    Array[File] magDepthPlots = MagDepth.magDepthPlots
    Array[File] magSNPPlots = MagSNP.magSNPPlots
    }
    meta {
        description: "CNV-Mag: A Tool to Visualize CNV Events in WGS Data"
        author: "Yueyao Gao"
        email: "gaoyueya@broadinstitute.org"
        }
}


task CreateBedFromIntervals {
    input {
        Array[String]? cnvIntervals
        String dockerImage
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
        source activate CNV-Mag
        python3 <<CODE

        with open('cnv_intervals.txt', 'r') as f:
                cnvIntervals = [line.strip() for line in f if line.strip()]

        with open('cnv_intervals.bed', 'a') as f:
            for i, interval in enumerate(cnvIntervals):
                chr = interval.split(':')[0]
                start = interval.split(':')[1].split('-')[0]
                end = interval.split(':')[1].split('-')[1]
                # Write without a newline for the last item
                if i == len(cnvIntervals) - 1:
                    f.write(f"{chr}\t{start}\t{end}")
                else:
                    f.write(f"{chr}\t{start}\t{end}\n")
        CODE
    >>>
    runtime {
        docker: dockerImage
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
        String refGenome
        String dockerImage
        Int padpct = 2 # Percentage to pad the CNV regions
        Int mem_gb = 4
        Int cpu = 1
        Int disk_size_gb = 10
    }

    command <<<
        if [[ ~{refGenome} == "hg19" ]]; then
            genomeBoundaryFile="/BaseImage/MagRef/Homo_sapiens_assembly19.genome"
        elif [[ ~{refGenome} == "hg38" ]]; then
            genomeBoundaryFile="/BaseImage/MagRef/Homo_sapiens_assembly38.genome"
        else
            echo "Reference genome $refGenome not supported"
            exit 1
        fi

        # Create a padded CNV bed file
        # Extend the CNV regions by 200% of the interval size on each side
        bedtools slop -i ~{cnvBedFile} -g ${genomeBoundaryFile} -b ~{padpct} -pct > padded_cnv.bed
    >>>
    runtime {
        docker: dockerImage
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
            File HG001Bam = "gs://fc-a76d0374-93e7-4c1a-8302-2a88079b480d/DRAGEN_4.3.6_CNV_Mag_resource/HG001.bam"
            File HG001Bai = "gs://fc-a76d0374-93e7-4c1a-8302-2a88079b480d/DRAGEN_4.3.6_CNV_Mag_resource/HG001.bai"
            File HG002Bam = "gs://fc-a76d0374-93e7-4c1a-8302-2a88079b480d/DRAGEN_4.3.6_CNV_Mag_resource/HG002.bam"
            File HG002Bai = "gs://fc-a76d0374-93e7-4c1a-8302-2a88079b480d/DRAGEN_4.3.6_CNV_Mag_resource/HG002.bai"
            Int minBQ = 20
            Int mem_gb = 64
            Int cpu = 8
            Int disk_size_gb = 500
            Boolean use_ssd = true
            String samtoolsDocker = "euformatics/samtools:1.20"
    }
    command <<<
        # Create output directory
        mkdir output

        # Run samtools depth to get MAPQ20 depth & MAPQ0 depth
        # Counting fragments instead of reads using -s option
        for mq in 0 20; do
            samtools depth \
            -@ ~{cpu} \
            -b ~{target_bed} \
            --min-BQ ~{minBQ} \
            --min-MQ ${mq} \
            -s \
            ~{alignedBam} \
            ~{HG001Bam} \
            ~{HG002Bam} \
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
        docker: samtoolsDocker
        disks: "local-disk " + disk_size_gb + if use_ssd then " SSD" else " HDD"
        preemptible: 0
        maxRetries: 3
    }
}

task MagDepth{
    input {
            String sampleName
            String dockerImage
            File mapq0_depth_profile
            File mapq20_depth_profile
            File cnvBedFile
            File PaddedcnvBedFile
            String refGenome
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

            # Run the MagDepth script
            conda run --no-capture-output \
            -n CNV-Mag \
            python3 /BaseImage/CNV-Mag/MagDepth.py \
            --maq20 ~{mapq20_depth_profile} \
            --maq0 ~{mapq0_depth_profile} \
            --bed ~{cnvBedFile} \
            --padded_bed ~{PaddedcnvBedFile} \
            --ref_genome ~{refGenome} \
            -o output

        >>>
        output {
            Array[File] magDepthPlots = glob("output/*png")
        }
        runtime {
            memory: mem_gb + " GB"
            cpu: cpu
            docker: dockerImage
            disks: "local-disk " + disk_size_gb + if use_ssd then " SSD" else " HDD"
            preemptible: preemptible
            maxRetries: maxRetries
        }
}

task MagSNP{
    input {
        String sampleName
        String dockerImage
        File hardFilteredVcfFile
        File HG001FilteredVcfFile = "gs://fc-a76d0374-93e7-4c1a-8302-2a88079b480d/DRAGEN_4.3.6_NIST_default/NA12878_HG001_1000ng_3_NVX/NA12878_HG001_1000ng_3_NVX.hard-filtered.vcf.gz"
        File HG002FilteredVcfFile = "gs://fc-a76d0374-93e7-4c1a-8302-2a88079b480d/DRAGEN_4.3.6_NIST_default/NA24385_HG002_1_NVX/NA24385_HG002_1_NVX.hard-filtered.vcf.gz"
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
        -n CNV-Mag \
        python3 /BaseImage/CNV-Mag/MagSNP.py \
        -v1 ~{hardFilteredVcfFile} \
        -v2 ~{HG001FilteredVcfFile} \
        -v3 ~{HG002FilteredVcfFile} \
        -b ~{cnvBedFile} \
        -n1 ~{sampleName} \
        -n2 HG001 \
        -n3 HG002  \
        -o output

    >>>
    output {
        Array[File] magSNPPlots = glob("output/*png")
    }
    runtime {
        memory: mem_gb + " GB"
        cpu: cpu
        docker: dockerImage
        disks: "local-disk " + disk_size_gb + if use_ssd then " SSD" else " HDD"
        preemptible: preemptible
        maxRetries: maxRetries
    }
}







