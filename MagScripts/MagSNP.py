import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib.ticker as ticker
from pysam import VariantFile
import logging
import sys
import argparse

argparser = argparse.ArgumentParser(description='MagSNP: Visualize PASS SNP AF distribution within the target intervals')
argparser.add_argument('-v1', '--vcf1', help='Path to VCF file 1', required=True)
argparser.add_argument('-v2', '--vcf2', help='Path to VCF file 2', required=True)
argparser.add_argument('-v3', '--vcf3', help='Path to VCF file 3', required=True)
argparser.add_argument('-b', '--bed', help='bed file with CNV interval list', required=True)
argparser.add_argument('-n1', '--name1', help='name of VCF file 1', required=False)
argparser.add_argument('-n2', '--name2', help='name of VCF file 2', required=False)
argparser.add_argument('-n3', '--name3', help='name of VCF file 3', required=False)
argparser.add_argument('-o', '--output', help='output directory path', default='.')
args = argparser.parse_args()

vcf1_path = args.vcf1
vcf2_path = args.vcf2
vcf3_path = args.vcf3
output_dir = args.output

# print the input arguments
print(f"VCF file 1: {vcf1_path}")
print(f"VCF file 2: {vcf2_path}")
print(f"VCF file 3: {vcf3_path}")

# Set up logging
logging.basicConfig(
    stream=sys.stdout,
    level=logging.DEBUG,
    format="%(levelname)s\t%(asctime)s\t%(name)s\t%(message)s",
    datefmt="%Y-%m-%d %H:%M:%S")
log = logging.getLogger("Create SNP AF Distribution Viz for CNV-Mag")

def annotate_af(vcf_path, interval):
    interval_chr = interval.split(":")[0]
    interval_start = int(interval.split(":")[1].split("-")[0])
    interval_end = int(interval.split(":")[1].split("-")[1])

    vcf = VariantFile(vcf_path)
    sample_name = list(vcf.header.samples)[0]
    pos_af_list = []
    for rec in vcf.fetch(contig=interval_chr, start=interval_start, end=interval_end):
        if rec.filter.keys() == ['PASS']:
            af = np.round(rec.samples[sample_name]["AF"][0], 2)
            pos_af_list.append((rec.pos, af))

    if len(pos_af_list) == 0:
        raise ValueError(f"No SNPs on {interval_chr} wihtin {interval_start} and {interval_end}")

    interval_pass_snp_df = pd.DataFrame(pos_af_list, columns=['POS', 'ALLELE_FRACTION'])
    return interval_pass_snp_df

def examine_interval_with_snp(vcf1_path: str,vcf2_path: str,interval_list: list, vcf3_path: str, dragen_call=None, vcf1_name="", vcf2_name="HG001", vcf3_name="HG002"):
    for interval in interval_list:
        try:
            # Get components of interval
            interval_chr = interval.split(":")[0]
            interval_pos = int(interval.split(":")[1].split("-")[0])
            interval_end = int(interval.split(":")[1].split("-")[1])

            # select SNPs in the interval
            print(f"VCF1: {vcf1_name}")
            vcf1_interval_snp_df = annotate_af(vcf1_path, interval)
            print(f"PASS SNP Count within {interval}: {vcf1_interval_snp_df.shape[0]}\n")
            print(f"VCF2: {vcf2_name}")
            vcf2_interval_snp_df = annotate_af(vcf2_path, interval)
            print(f"PASS SNP Count within {interval}: {vcf2_interval_snp_df.shape[0]}\n")
            print(f"VCF3: {vcf3_name}")
            vcf3_interval_snp_df = annotate_af(vcf3_path, interval)
            print(f"PASS SNP Count within {interval}: {vcf3_interval_snp_df.shape[0]}\n")


            # Alpha for scatter plot
            if vcf1_interval_snp_df.shape[0] <= 1000:
                alpha = 1
            else:
                alpha = 0.1

            f, axs = plt.subplots(nrows=2, ncols=3, figsize=(15, 5), sharex='col', sharey='row')
            # Plotting VCF1
            sns.scatterplot(data=vcf1_interval_snp_df, y='POS', x='ALLELE_FRACTION', ax=axs[0,0], alpha=alpha,color='black',s=1)
            sns.histplot(data=vcf1_interval_snp_df, x='ALLELE_FRACTION', bins=50, ax=axs[1,0], facecolor='black', linewidth=1, alpha=0.75)
            axs[0, 0].yaxis.set_major_formatter(ticker.FuncFormatter(lambda x, _: '{:,.0f}'.format(x)))
            axs[0, 0].set_ylabel("POS", fontsize=14, fontweight='bold')
            axs[1, 0].set_xlim(0, 1.1)
            axs[1, 0].set_xlabel("Allele Fraction", fontsize=14, fontweight='bold')
            axs[1, 0].set_ylabel("Count", fontsize=14, fontweight='bold')
            axs[1, 0].set_xticks(np.arange(0, 1.1, 0.1))
            interval_af_1_df = vcf1_interval_snp_df[vcf1_interval_snp_df['ALLELE_FRACTION']==1].shape[0]
            percent_interval_vcf1_af_1 = (interval_af_1_df/vcf1_interval_snp_df.shape[0])*100
            axs[1, 0].text(0.1, 0.8, f"Total SNPs: {vcf1_interval_snp_df.shape[0]:,}\nAF=1: {percent_interval_vcf1_af_1:.2f}%", fontsize=10, fontweight='bold', transform=axs[1,0].transAxes, color='blue')
            axs[0, 0].set_title(vcf1_name, fontsize=14, fontweight='bold')

            # Plotting VCF2
            sns.scatterplot(data=vcf2_interval_snp_df, y='POS', x='ALLELE_FRACTION', ax=axs[0,1], alpha=alpha,color='black',s=1)
            sns.histplot(data=vcf2_interval_snp_df, x='ALLELE_FRACTION', bins=50, ax=axs[1,1], facecolor='black', linewidth=1, alpha=0.75)
            axs[0, 1].yaxis.set_major_formatter(ticker.FuncFormatter(lambda x, _: '{:,.0f}'.format(x)))
            axs[0, 1].set_title(vcf2_name, fontsize=14, fontweight='bold')
            axs[1, 1].set_xlim(0, 1.1)
            axs[1, 1].set_xlabel("Allele Fraction", fontsize=14, fontweight='bold')
            axs[1, 1].set_xticks(np.arange(0, 1.1, 0.1))
            interval_af_1_df = vcf2_interval_snp_df[vcf2_interval_snp_df['ALLELE_FRACTION']==1].shape[0]
            percent_interval_vcf2_af_1 = (interval_af_1_df/vcf2_interval_snp_df.shape[0])*100
            axs[1, 1].text(0.1, 0.8, f"Total SNPs: {vcf2_interval_snp_df.shape[0]:,}\nAF=1: {percent_interval_vcf2_af_1:.2f}%", fontsize=10, fontweight='bold', transform=axs[1,1].transAxes, color='blue')

            # Plotting VCF3
            sns.scatterplot(data=vcf3_interval_snp_df, y='POS', x='ALLELE_FRACTION', ax=axs[0,2], alpha=alpha,color='black',s=1)
            sns.histplot(data=vcf3_interval_snp_df, x='ALLELE_FRACTION', bins=50, ax=axs[1,2], facecolor='black', linewidth=1, alpha=0.75)
            axs[0, 2].yaxis.set_major_formatter(ticker.FuncFormatter(lambda x, _: '{:,.0f}'.format(x)))
            axs[0, 2].set_title(vcf3_name, fontsize=14, fontweight='bold')
            axs[1, 2].set_xlim(0, 1.1)
            axs[1, 2].set_xlabel("Allele Fraction", fontsize=14, fontweight='bold')
            axs[1, 2].set_xticks(np.arange(0, 1.1, 0.1))
            interval_af_1_df = vcf3_interval_snp_df[vcf3_interval_snp_df['ALLELE_FRACTION']==1].shape[0]
            percent_interval_vcf3_af_1 = (interval_af_1_df/vcf3_interval_snp_df.shape[0])*100
            axs[1, 2].text(0.1, 0.8, f"Total SNPs: {vcf3_interval_snp_df.shape[0]:,}\nAF=1: {percent_interval_vcf3_af_1:.2f}%", fontsize=10, fontweight='bold', transform=axs[1,2].transAxes, color='blue')

            if dragen_call:
                for i in [0, 1]:
                    axs[0, i].axhline(y=dragen_call[0], color='red', linestyle='--', label='DRAGEN CNV Call')
                    axs[0, i].axhline(y=dragen_call[1], color='red', linestyle='--')
                    axs[0, i].legend()
            fmt_interval = f"{interval_chr}:{interval_pos:,}-{interval_end:,}"
            plt.suptitle(t=f"SNP AF at Interval {fmt_interval}", fontsize=16)

            if output_dir.endswith('/'):
                plt.savefig(f"{output_dir}SNP_AF_Interval_{interval.replace(':','_')}.png", dpi=300)
            else:
                plt.savefig(f"{output_dir}/SNP_AF_Interval_{interval.replace(':','_')}.png", dpi=300)
        except Exception as e:
            print(f"Error processing interval {interval}: {e}")
            continue

# Initialize CNV interval list
cnv_interval_list = []
with open(args.bed, 'r') as f:
    f = [line for line in f if line.strip()]
    for line in f:
        chr = line.strip().split('\t')[0]
        start = line.strip().split('\t')[1]
        end = line.strip().split('\t')[2]
        print(f'{chr}:{start}-{end}')
        cnv_interval_list.append(f'{chr}:{start}-{end}')

if args.name1 and args.name2 and args.name3:
    vcf1_name = args.name1
    vcf2_name = args.name2
    vcf3_name = args.name3
else:
    vcf1_name = args.vcf1.split('/')[-1].split('.')[0]
    vcf2_name = args.vcf2.split('/')[-1].split('.')[0]
    vcf3_name = args.vcf3.split('/')[-1].split('.')[0]

examine_interval_with_snp(vcf1_path=vcf1_path,vcf2_path=vcf2_path, vcf3_path=vcf3_path, interval_list=cnv_interval_list, dragen_call=None, vcf1_name=vcf1_name, vcf2_name=vcf2_name, vcf3_name=vcf3_name)
