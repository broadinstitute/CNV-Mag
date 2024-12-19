# CNV-Mag
A tool to visualize CNV Events in WGS Data. CNV-Mag is taking use of BAM and short variant VCF files to visualize CNV events in WGS data in simple and direct way. It can help you to verify if a CNV event is real or not.
This tool is part of the Broad Clinical Labs’s CNV reporting workflow.

## Background
Verify if a CNV event is real or not is a challenging task. Loading aligned BAM file and VCF into IGV and check 
the read depth, soft clip, and SNP Allele Fraction (AF) is a common way to verify CNV events. However, it 
is not easy to check the read depth in IGV, especially when the breakpoint is not accurate and even an 
experienced Bioinformatician may make mistakes. CNV-Mag is designed to help you to visualize CNV events in a
simple and direct way. It will show you the read depth, and SNP AF of the CNV region in a single figure.

The soft-clip information is not included in the current version of CNV-Mag. I am working on it and will
release it in the future.

## Workflow
![image](misc/CNV-Mag TubeMap.png)
There are two main workflows in CNV-Mag:
- Mag-SNP
  - Visualize all the PASS SNP AF of the CNV region.
![image](misc/SNP_AF_Interval_chr15_73536771-81697074.png)
- Mag-Depth
  - Visualize the read depth of the CNV region.
![image](misc/SDSM-WZ_MagDepth_chr15_73536771-81697074.png)

## Running CNV-Mag on Terra
You can run CNV-Mag on Terra. The workflow is available at https://dockstore.org/workflows/github.com/broadinstitute/CNV-Mag/CNV-Mag:dev?tab=info.
You can deploy the workflow to your workspace and run it on Terra.

## Runtime specification
CNV-Mag is part of Broad Clinical Lab CNV reporting workflow. The run cost is less than $2 per run. 
