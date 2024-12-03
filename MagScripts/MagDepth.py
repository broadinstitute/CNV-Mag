import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
from matplotlib.lines import Line2D
import matplotlib.patches as patches
import seaborn as sns
import argparse
import logging
import sys

argparser = argparse.ArgumentParser(prog='MagDepth', description='MagDepth: Visualize read depth distribution near target intervals')
argparser.add_argument('-mq20', '--maq20', help='Path to samtools MapQ20 depth file', required=True)
argparser.add_argument('-mq0', '--maq0', help='Path to samtools MapQ0 depth file', required=True)
argparser.add_argument('-r', '--ref_genome', help='Reference genome. Either "hg38" or "hg19"', required=True)
argparser.add_argument('-b', '--bed', help='BED file with CNV interval list', required=True)
argparser.add_argument('-pb', '--padded_bed', help='Padded BED file of input CNV interval list', required=True)
argparser.add_argument('-o', '--output', help='(Optional) Output directory path', default='.')


args = argparser.parse_args()

samtools_mq20_depth_path = args.maq20
samtools_mq0_depth_path = args.maq0
ref_genome = args.ref_genome
bed_path = args.bed
padded_cnv_bed_path = args.padded_bed
output_dir = args.output

# print the input arguments
print(f"Samtools MAPQ20 depth file: {samtools_mq20_depth_path}")
print(f"Samtools MAPQ0 depth file: {samtools_mq0_depth_path}")
print(f"Reference genome: {ref_genome}")
print(f"CNV interval BED file: {bed_path}")
print(f"Padded CNV interval BED file: {padded_cnv_bed_path}")
print(f"Output directory: {output_dir}")


# Set up logging
logging.basicConfig(
    stream=sys.stdout,
    level=logging.DEBUG,
    format="%(levelname)s\t%(asctime)s\t%(name)s\t%(message)s",
    datefmt="%Y-%m-%d %H:%M:%S")
log = logging.getLogger("Create Depth Viz for CNV-Mag")

# Create a dictionary to store chromosome sizes AND a df to store centromere positions
if ref_genome == 'hg38':
    chrom_size_path = "/BaseImage/MagRef/Homo_sapiens_assembly38.genome"
    centromere_bed_path = "/BaseImage/MagRef/Homo_sapiens_assembly38_centromeres.bed"
elif ref_genome == 'hg19':
    chrom_size_path = "/BaseImage/MagRef/Homo_sapiens_assembly19.genome"
    centromere_bed_path = "/BaseImage/MagRef/Homo_sapiens_assembly19_centromeres.bed"
else:
    raise ValueError('INPUT REF GENOME MUST BE EITHER "hg38" or "hg19"')

chrom_size_dict = {}
with open(chrom_size_path,'r') as f:
    for line in f.readlines():
        chr, size = line.split()
        chrom_size_dict[chr] = int(size)

centromere_df = pd.read_csv(filepath_or_buffer=centromere_bed_path, sep='\t', names=['chr', 'start', 'end','name'])


def custom_formatter(chr, x, pos):
    return f'{chr}:{x:,.0f}'


def get_smoothed_depth(interval_depth:pd.DataFrame, bin_size=5000):
    cov_list = interval_depth['cov'].tolist()
    pos_list = interval_depth['pos'].tolist()
    smoothed_cov_dict = {}
    for i in range(0, len(cov_list), bin_size):
        bin_pos = int(np.round(np.mean(pos_list[i:i+bin_size]), decimals=0))
        bin_cov = np.mean(cov_list[i:i+bin_size])
        smoothed_cov_dict[bin_pos] = bin_cov
    smoothed_cov_df = pd.DataFrame(list(smoothed_cov_dict.items()), columns=['pos', 'cov'])
    return smoothed_cov_df

def get_binned_histogram(interval_depth:pd.DataFrame, bin_size=5000, depth_bins=30):
    # Create an empty list to hold binned data
    binned_data = []
    # Bin the data into bins of size bin_size, then sub-bin the data into depth_bins
    for i in range(0, len(interval_depth), bin_size):
        pos_window = interval_depth['pos'][i:i+bin_size]
        cov_window = interval_depth['cov'][i:i+bin_size]
        # If there are fewer than bin_size positions left, stop
        if len(pos_window) == 0:
            break
        # Divide the window into depth_bins
        bin_edges = np.linspace(pos_window.min(), pos_window.max(), depth_bins+1)
        # Digitize the positions into these sub-bins
        pos_digitized = np.digitize(pos_window, bin_edges)
        # For each sub-bin, calculate the mean position and coverage
        for j in range(1, depth_bins+1):
            mask = pos_digitized == j
            if np.any(mask):  # Check if there are any values in this bin
                mean_pos = np.mean(pos_window[mask])  # Mean position within this sub-bin
                mean_cov = np.mean(cov_window[mask])  # Mean coverage within this sub-bin
                binned_data.append((mean_pos, mean_cov))
    # Convert the binned data into a DataFrame
    binned_df = pd.DataFrame(binned_data, columns=['pos', 'cov'])
    return binned_df


# Define the function to get CNV depth profile
def plot_cnv_depth(mapq20_depth_path, mapq0_depth_path, target_intervals, padded_intervals):
    sample_name = mapq0_depth_path.split("/")[-1].split("_")[0]
    mapq0_proband_depths = pd.read_csv(mapq0_depth_path, sep="\t", usecols=[0, 1, 2], names=["contig", "pos", "proband"])['proband'].tolist()
    cov_y_uplimit = np.percentile(mapq0_proband_depths, q=99)*1.5
    sample_names = [sample_name, "HG001", "HG002"]
    for index, interval in enumerate(target_intervals):
        # Get interval start and end
        interval_chr = interval.split(':')[0]
        interval_pos, interval_end = interval.split(':')[1].split('-')
        interval_pos = int(interval_pos)
        interval_end = int(interval_end)
        print("Processing interval:", interval)
        # Padded interval
        padded_interval = padded_intervals[index]
        print("Padded interval:", padded_interval)
        padded_interval_chr = padded_interval.split(':')[0]
        padded_interval_pos, padded_interval_end = padded_interval.split(':')[1].split('-')
        padded_interval_pos = int(padded_interval_pos)
        padded_interval_end = int(padded_interval_end)
        padded_size = padded_interval_end - padded_interval_pos + 1
        # Get bin size and alpha for scatter plot
        bin_size = int(round(padded_size/500, 0))
        print(f"Bin size: {bin_size}")

        # Sanity check
        if interval_chr != padded_interval_chr:
            raise Exception(f"Interval chromosome mismatch: {interval_chr} != {padded_interval_chr}")
        # For each interval, Create a depth plot for proband, HG001, and HG002
        fig, axs = plt.subplots(nrows=4, ncols=1, figsize=(15, 8), height_ratios=[1, 1, 1, 0.1])
        for i in range(2, 5):
            mq20_depth_df = pd.read_csv(mapq20_depth_path, sep="\t", usecols=[0, 1, i], names=["contig", "pos", "cov"])
            mq0_depth_df = pd.read_csv(mapq0_depth_path, sep="\t", usecols=[0, 1, i], names=["contig", "pos", "cov"])
            # Get interval specific depth
            interval_mq20_depth_cov_df = mq20_depth_df[(mq20_depth_df['contig'].astype(str)==str(interval_chr))&(mq20_depth_df['pos']>=padded_interval_pos)&(mq20_depth_df['pos']<=padded_interval_end)]
            interval_mq0_depth_cov_df = mq0_depth_df[(mq0_depth_df['contig'].astype(str)==str(interval_chr))&(mq0_depth_df['pos']>=padded_interval_pos)&(mq0_depth_df['pos']<=padded_interval_end)]
            # Get the smoothed depth profile
            smoothed_mq20_depth_cov_df = get_smoothed_depth(interval_mq20_depth_cov_df, bin_size)
            smoothed_mq0_depth_cov_df = get_smoothed_depth(interval_mq0_depth_cov_df, bin_size)
            # Get the binned histogram
            binned_mq20_depth_cov_df = get_binned_histogram(interval_mq20_depth_cov_df, bin_size=bin_size, depth_bins=30)
            # Plot MQ20 and MQ0 smoothed depth in line plot, binned MQ20 depth, and highlight the CNV interval
            # Plot cov
            plot_index = i-2
            sns.histplot(x='pos', y='cov', data=binned_mq20_depth_cov_df, color='blue', alpha=1, ax=axs[plot_index], bins=30, thresh=2.5)
            sns.lineplot(x='pos', y='cov', data=smoothed_mq0_depth_cov_df, color='silver', ax=axs[plot_index], linewidth=1.5, alpha=0.7)
            sns.lineplot(x='pos', y='cov', data=smoothed_mq20_depth_cov_df, color='black', ax=axs[plot_index], linewidth=1.5)



            axs[plot_index].xaxis.set_major_formatter(mticker.FuncFormatter(lambda x, pos: custom_formatter(interval_chr, x, pos)))
            axs[plot_index].set_ylim(0, cov_y_uplimit)
            axs[plot_index].set_xlim(padded_interval_pos-padded_size*0.05, padded_interval_end+padded_size*0.05)
            axs[plot_index].set_ylabel('Read Depth', fontsize=14)
            axs[plot_index].set_xlabel('')
            axs[plot_index].tick_params(axis='x', labelsize=12)
            axs[plot_index].tick_params(axis='y', labelsize=12)

            axs[plot_index].set_title(sample_names[plot_index], fontsize=14, fontweight='bold')

            # Create custom legend handles
            legend_handles = [
                Line2D(xdata=[0], ydata=[0], color='black', lw=1.5, label=f'MQ20 Depth'),
                Line2D(xdata=[0], ydata=[0], color='silver', lw=1.5, label=f'MQ0 Depth'),
                axs[plot_index].fill_betweenx(y=[0, cov_y_uplimit], x1=interval_pos, x2=interval_end, color='yellow', alpha=0.2, label='Target Interval')
            ]
            axs[plot_index].legend(handles=legend_handles, bbox_to_anchor=(1.01, 0.5), loc='center left', fontsize=12)

        #Plot centromere and chromosome as a bottom reference track
        axs[3].add_patch(patches.Rectangle(xy=(0, 0), width=chrom_size_dict[interval_chr] , height=.5, edgecolor='black', fill=False, label='Chromosome', linewidth=1))
        axs[3].set_xlim(0-chrom_size_dict[interval_chr]*0.005, chrom_size_dict[interval_chr]*1.005)
        axs[3].set_ylim(-0.5, 1)

        for index, row in centromere_df[centromere_df['chr']==interval_chr].iterrows():
            axs[3].add_patch(patches.Rectangle((row['start'], 0), row['end']-row['start'], 0.5, edgecolor='r', facecolor='r', fill=True, linewidth=0, label='Centromere'))

        # Plot the target interval on the chromosome
        axs[3].add_patch(patches.Rectangle((interval_pos, 0), interval_end-interval_pos, 0.5, edgecolor='y', facecolor='y', fill=True, linewidth=1, label='Interval'))
        axs[3].scatter((interval_pos+interval_end)/2, -0.5, color='darkorange', label='indicator', s=100, marker='^')

        # Clean up the chromosome track
        axs[3].set_yticks([])
        axs[3].set_ylabel('')
        axs[3].axis('off')

        # Add text labels (p-arm and q-arm)
        axs[3].text(0, 0.8, 'p-arm', fontsize=10, ha='left', color='black')
        axs[3].text(chrom_size_dict[interval_chr], 0.8, 'q-arm', fontsize=10, ha='right', color='black')
        axs[3].text(centromere_df[centromere_df['chr']==interval_chr]['start'].mean(), 0.8, 'centromere', fontsize=10, ha='center', color='red')
        axs[3].text((interval_pos+interval_end)/2, 1, f'Target', fontsize=10, ha='center', color='darkorange', fontweight='bold')

        # Add overall display description
        axs[3].text(chrom_size_dict[interval_chr] / 2, -2, f'Overall display of {interval_chr}', fontsize=12, ha='center', color='black', fontweight='bold')
        fmt_interval = f"{interval_chr}:{interval_pos:,}-{interval_end:,}"
        plt.suptitle(t=f'Read Depth Distribution Near {fmt_interval}', fontsize=16)
        plt.tight_layout()
        # Save the plot in the output directory
        if output_dir.endswith('/'):
            plt.savefig(fname=f"{output_dir}{sample_name}_MagDepth_{interval.replace(':','_')}.png", dpi=600)
        else:
            plt.savefig(fname=f"{output_dir}/{sample_name}_MagDepth_{interval.replace(':','_')}.png", dpi=600)


# Initialize CNV interval list
cnv_interval_list = []
with open(bed_path, 'r') as f:
    for line in f:
        chrom = line.strip().split('\t')[0]
        start = line.strip().split('\t')[1]
        end = line.strip().split('\t')[2]
        print(f'{chrom}:{start}-{end}')
        cnv_interval_list.append(f'{chrom}:{start}-{end}')



padded_interval_list = []
with open(padded_cnv_bed_path, 'r') as f:
    for line in f:
        chrom = line.strip().split('\t')[0]
        start = line.strip().split('\t')[1]
        end = line.strip().split('\t')[2]
        print(f'{chrom}:{start}-{end}')
        padded_interval_list.append(f'{chrom}:{start}-{end}')

# Get CNV depth profile
plot_cnv_depth(samtools_mq20_depth_path, samtools_mq0_depth_path, cnv_interval_list, padded_interval_list)

# TODO: Edit Logging Part of MagDepth.py