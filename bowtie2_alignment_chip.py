#bowtie2_cutnrun.py
from datetime import datetime
import os
import glob
import shutil
import sys
from subprocess import check_output

startTime = datetime.now()

# Change directory within container
nav1 = os.path.isdir('/data/trim_out')
nav2 = os.path.isdir('/data/fastq')
if nav1 == True:
	os.chdir('/data/trim_out')
	print 'Current working directory is:' +  os.getcwd()
	print '\n'
elif nav2 == True:
	os.chdir('/data/fastq')
	print 'Current working directory is:' +  os.getcwd()
	print '\n'
else:
	os.chdir('/data')
	print 'Current working directory is:' +  os.getcwd()
	print '\n'

#prepare input files
input_files = sorted(glob.glob('*trimmed*'))

if len(input_files) == 0:
	print 'No quality-controlled input files from trim_galore, checking input folder for fastqc output...'
	print '\n'
	input_files = sorted(glob.glob('*.fastq'))
	if len(input_files) == 0:
		print 'No uncompressed fastq files detected, looking for fastq.gz...'
		print '\n'
		input_files = sorted(glob.glob('*.fastq.gz'))
		if len(input_files) == 0:
			print 'No valid input files detected. Exiting...'	
			sys.exit()
		else:
			print 'Decompressing .gz files'
			#create fastq filename from .fastq.gz
			temp_list = [f.replace('.gz', '') for f in input_files]
			#run system command to decompress file with zcat
			for i in range(len(input_files)):
				temp_str= 'zcat ' + input_files[i] + ' > ' + temp_list[i]
				check_output(temp_str, shell=True)
			input_files = temp_list



print 'Files for single end alignment:'
print '\n'.join(input_files)
print '\n'

#detect spike-in genome preferences
fly_spike = os.path.isfile('/genomes/spike/dm6/Sequence/Bowtie2Index/genome.1.bt2')

yeast_spike = os.path.isfile('/genomes/spike/sacCer3/Sequence/Bowtie2Index/genome.1.bt2')

#make sam output names
sam_names = []
spike_sam_names = []

for i in range(len(input_files)):
	sam_name = input_files[i].split('.f')[0] + '.sam'
	sam_names.append(sam_name)

print 'Output sam filenames for data species:'
print '\n'.join(sam_names)
print '\n'


if fly_spike == True and yeast_spike == True:
	print 'Two different spike in genomes detected. Exiting...'

elif fly_spike == True and yeast_spike != True:
	spike_species = 'dm6'
	print 'Spike-in chromatin from Drosophila'
	for i in range(len(sam_names)):
		temp_name = sam_names[i] + '.dm6'
		spike_sam_names.append(temp_name)
	print 'Output sam filenames for spike species:'
	print '\n'.join(spike_sam_names)
	print '\n'
	print 'Data will be aligned to hg19 and ' + spike_species

elif fly_spike != True and yeast_spike == True:
	spike_species = 'sacCer3'
	print 'Spike-in chromatin from Saccharomyces'
	for i in range(len(sam_names)):
		temp_name = sam_names[i] + '.sacCer3'
		spike_sam_names.append(temp_name)
	print 'Output sam filenames for spike species:'
	print '\n'.join(spike_sam_names)
	print '\n'
	print 'Data will be aligned to hg19 and ' + spike_species

else:
	spike_species = 'null'
	print 'Data will be aligned to hg19'
	print 'No spike-in alignment will be performed.'
	print '\n'
	print '\n'

#define bowtie2 index for hg19
hg19_index = '/genomes/test/Sequence/Bowtie2Index/genome'
print 'Bowtie index for data files found at:'
print hg19_index
print '\n'

#run bowtie for hg19
print 'Running bowtie2 alignment for hg19'
print '\n'
for i in range(len(input_files)):
	print 'count = ' + str(i)
	print '\n'
	#create string for system command
	temp_str = 'bowtie2 --sensitive --no-unal --phred33 -q -I 10 -X 500 --threads 16' \
	+ ' -x ' + hg19_index + ' -U ' + input_files[i] + ' >> ' + sam_names[i]

	print temp_str

	check_output(temp_str, shell=True)
	print '\n'

#define bowtie2 index for spike_species
if spike_species == 'dm6': 
	spike_index = '/genomes/spike/dm6/Sequence/Bowtie2Index/genome'
	print 'Bowtie index for spike files found at:'
	print spike_index
	print '\n'

elif spike_species == 'sacCer3':
	spike_index = '/genomes/spike/sacCer3/Sequence/Bowtie2Index/genome'
	print 'Bowtie index for spike files found at:'
	print spike_index
	print '\n'

elif spike_species not in ['dm6', 'sacCer3']:
	print 'Bowtie index unavailable for spike species or no spike specied index detected.'
	print '\n'

#run bowtie for spike_species
if spike_species in ['dm6', 'sacCer3']:
	print 'Running bowtie2 alignment for ' + spike_species
	print '\n'
	for i in range(len(input_files)):
		print 'count = ' +str(i)
		print '\n'
		#create string for system command
		temp_str = 'bowtie2 --very-sensitive --no-unal --phred33 -q -I 10 -X 700 --threads 16' \
		+ ' -x ' + spike_index + ' -U ' + input_files[i] + ' >> ' + spike_sam_names[i]

		print temp_str

		check_output(temp_str, shell=True)
		print '\n'

#make new directory for output
os.mkdir('/data/sams')
print 'Current working directory is:' +  os.getcwd()
print '\n'

#copy files to output folder
output_dir = '/data/sams/'
print 'Moving sam files to output folder'
print '\n'
for i in range(len(sam_names)):
	shutil.move(sam_names[i], output_dir)

if spike_species in ['dm6', 'sacCer3']:
	for i in range(len(spike_sam_names)):
		shutil.move(spike_sam_names[i], output_dir)


print 'Alignment Runtime (hh:mm:ss): ' + str(datetime.now() - startTime)
print '\n'

###SAM conversion to bam, bedgraph, and BigWig

os.mkdir('/data/bams')
os.chdir('/data/sams')
print 'Current working directory is:' +  os.getcwd()
print '\n'

import pybedtools
from pybedtools import BedTool
import pandas as pd
from pybedtools.helpers import chromsizes
from pybedtools.contrib.bigwig import bedgraph_to_bigwig
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

startTime = datetime.now()

#######################################################
## spike in or average normalization with file generation#########
#######################################################

if spike_species in ['dm6', 'sacCer3']:	
	print 'Converting sams to bams, bedgraph, and BigWig with spike-in normalization.'
	print '\n'
	datafiles = sorted(glob.glob('*.sam'))
	spikefiles = sorted(glob.glob('*.' + spike_species))
	
	print '\n'
	print 'Data files loaded:'
	print '\n'.join(datafiles)
	print '\n'
	print 'Spike files loaded:'
	print '\n'.join(spikefiles)
	print '\n'

	##check equal spike and data files
	if len(datafiles) == 0:
		print 'No data files imported. Exiting...'
		sys.exit()
	
	if len(datafiles) != len(spikefiles):
		print 'Unequal number of data and spike files. Exiting...'
		sys.exit()

	for i in range(len(datafiles)):
		if datafiles[i] not in spikefiles[i]:
			print 'Unmatched pairs and spike files.'
			print datafiles[i] + ' not paired with ' + spikefiles[i]
			print 'Exiting...'
			sys.exit()

else:
	print 'Converting sams to bams, bedgraph, and BigWig without spike-in normalization.'
	print '\n'
	datafiles = sorted(glob.glob('*.sam'))

	print '\n'
	print 'Data files loaded:'
	print '\n'.join(datafiles)
	print '\n'
	print 'No spike-in normalization'
	print '\n'

##convert to bam format
print 'Converting to bam format'
print '\n'
print '\n'

bam_names = [f.replace('sam', 'bam') for f in datafiles]

print '\n'
print '\n'.join(bam_names)
print '\n'
print '\n'
print 'SAM to BAM'
print '\n'

bam_string = []

for i in range(len(bam_names)):
        bam_string.append('samtools view -b -S ' + datafiles[i] + ' > /data/bams/' + bam_names[i])

for item in bam_string:
        check_output(item, shell = True)

datafiles = bam_names

##Sort and index
print '\n'
print 'Sorting bams'
print '\n'

sorted_bam_names = []
for i in range(len(bam_names)):
	sorted_bam_name = 'sorted.' + bam_names[i]
	sorted_bam_names.append(sorted_bam_name)

bam_names_sorted = [f.replace('.bam', '') for f in sorted_bam_names]

sort_string = []

for i in range(len(bam_names)):
		sort_string.append('samtools sort /data/bams/' + bam_names[i] + ' /data/bams/' + bam_names_sorted[i])

for item in sort_string:
		check_output(item, shell = True)

print '\n'
print 'Bam files to index:'
print '\n'.join(bam_names_sorted)
print '\n'


print 'Indexing bams'
print '\n'             
os.chdir('/data/bams')
print 'Current working directory is:' +  os.getcwd()

index_string = []

for i in range(len(bam_names_sorted)):
        index_string.append('samtools index ' + sorted_bam_names[i])

for item in index_string:
		check_output(item, shell = True)

sorted_files = sorted(glob.glob('sorted*'))
os.mkdir('/data/bams/sorted_bams')
output_dir = '/data/bams/sorted_bams'
for i in range(len(sorted_files)):
	shutil.move(sorted_files[i], output_dir)


##generate bed files from bam files
print 'Generating bed files representing whole insert from paired end reads in the data files'
print '\n'

print 'Current working directory is:' +  os.getcwd()
print '\n'

#generate bed file names
bed_names = [f.replace('bam', 'bed') for f in datafiles]
stripped_names = [f.replace('.bam', '_stripped.bed') for f in datafiles]

#generate bed files with bam_to_bed tool (makes bed12 format)
for i in range(len(datafiles)):
	temp_bed = BedTool(datafiles[i]).bam_to_bed(bedpe=False).to_dataframe()

	temp_bed.to_csv(bed_names[i], sep="\t", header = False, index = False)

	temp_bed_stripped = temp_bed.iloc[:,[0,1,2]].sort_values(by = ['chrom', 'start', 'end'])

	temp_bed_stripped['length'] = temp_bed_stripped['end'] - temp_bed_stripped['start']

	temp_bed_stripped.to_csv(stripped_names[i], sep="\t", header = False, index = False)

print 'Finished generating bed files:'
print '\n'
print 'whole insert bed files:' + '\n' + '\n'.join(bed_names)
print '\n'

#generate normalized bedgraphs
print 'Current working directory is:' +  os.getcwd()
print '\n'

#generate bedgraph names
bg_names = [f.replace('bed', 'bg') for f in bed_names]

if spike_species in ['dm6', 'sacCer3']:
	print 'Generating spike-normalized bedgraphs for all the bed files'
	print '\n'
	#calculate number of reads in spike files
	#spikecount is a list of reads populated by the 'samtools view -c' shell command for each spikefile	
	spikecount = []
	spike_string = []

	for item in spikefiles:
		spike_string.append('samtools view -c -S /data/sams/' + item)

	for item in spike_string:
		spikecount.append(int(check_output(item, shell = True))/2)

	print 'spike counts'
	print spikecount
	print '\n'

	#calculate scaling factors
	scaling_factor = []
	for item in spikecount:
		scaling_factor.append(float(10000)/item)

	for i in range(len(spikefiles)):
		count = str(spikecount[i])
		print '\n'
		print 'for ' + spikefiles[i] + ' count:' 
		print count
		print '\n'
		print 'scaling factor for ' + bed_names[i] + ' is:'
		print scaling_factor[i]

else:
	print 'Generating average normalized bedgraphs for all the bed files'
	print '\n'
	#count total number of reads in each bed file (before size selection)
	read_count = []
	for item in stripped_names:
		read_count.append(BedTool(item).count())

	print read_count

	#calculate genome size
	genome_file = chromsizes('hg19')
	DF = pd.DataFrame.from_dict(genome_file, orient='index')
	genome_size = DF[1].sum()
	
	scaling_factor = []
	for i in range(len(read_count)):
		scaling_factor.append(float(genome_size) / read_count[i])

	for i in range(len(bed_names)):
		count = str(read_count[i])	
		print '\n'
		print 'for ' + bed_names[i] + ' read count:'
		print count
		print '\n'
		print 'scaling factor for ' + bed_names[i] + ' is:'
		print scaling_factor[i]

#run bedtools genomecov to generate bedgraph files
for i in range(len(bg_names)):
	BedTool(stripped_names[i]).genome_coverage(bg = True, genome = 'hg19', scale = scaling_factor[i]).moveto(bg_names[i])

print 'Finished generating bedgraph files:'
print '\n'
print 'whole insert bedgraph files:' + '\n' + '\n'.join(bg_names)
print '\n'


##make bigwig files
print 'Current working directory is:' +  os.getcwd()
print '\n'

print 'Generating big_wig files for all the bedgraphs'
print '\n'

#generate bigwig names
bw_names = [f.replace('bg', 'bw') for f in bg_names]

#run bedgraph_to_bigwig tool
for i in range(len(bg_names)):
	bedgraph_to_bigwig(BedTool(bg_names[i]), 'hg19', bw_names[i])

print 'Finished generating bigwig files:'
print '\n'
print 'whole insert bigwig files:' + '\n' + '\n'.join(bw_names)
print '\n'

os.mkdir('/data/bams/beds')
os.mkdir('/data/bams/norm_bedgraphs')
os.mkdir('/data/bams/norm_bigwigs')

output_dir0 = '/data/bams/beds'
for i in range(len(bed_names)):
	shutil.move(bed_names[i], output_dir0)
	shutil.move(stripped_names[i], output_dir0)

print 'Moving bedgraphs to output folder'
print '\n'
output_dir1 = '/data/bams/norm_bedgraphs'
for i in range(len(bg_names)):
	shutil.move(bg_names[i], output_dir1)

print 'Moving bigwigs to output folder'
print '\n'
output_dir2 = '/data/bams/norm_bigwigs'
for i in range(len(bw_names)):
	shutil.move(bw_names[i], output_dir2)

print 'Finished'
print '\n'
print 'Runtime file conversion (hh:mm:ss): ' + str(datetime.now() - startTime)

#run plotFingerprint

startTime = datetime.now()

os.chdir('/data/bams/sorted_bams')
print 'Current working directory is:' +  os.getcwd()

print 'Running plotFingerprint for all bam files in folder'
print '\n'
temp_str = 'plotFingerprint -b *.bam --smartLabels --minFragmentLength 20 -T "Mapping Fingerprints" --plotFile mapqc_fingerprint.png --outQualityMetrics mapqc.tab'

print temp_str

check_output(temp_str, shell=True)
print '\n'

print 'Finished'
print '\n'
print 'Runtime mapping QC (hh:mm:ss): ' + str(datetime.now() - startTime)


	
