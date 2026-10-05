import numpy as np
import pandas as pd
import os, subprocess
import cooler
import cooltools
from multiprocessing import Pool
import bioframe

############################
######## Parameters ########
############################

# os.chdir() calls to move to the appropriate dirs; for HEK see os.chdir below
#os.chdir('/oak/stanford/groups/tttt/users/achuthan/research/plate-c_granule_plate/merge_240419a/coolfiles/balanced_2')  # granule
#os.chdir('/oak/stanford/groups/tttt/users/achuthan/research/plate-c_ipsc/merged_ARVA2_deep_240928a/coolfiles/balanced_2') # ipsc

# HEK
# Balance the following two directories separately
os.chdir('/oak/stanford/groups/tttt/users/achuthan/research/plate-c_HEK/primary_hek/coolfiles/balanced_2')  
cool_files_dir = '/oak/stanford/groups/tttt/users/achuthan/research/plate-c_HEK/primary_hek/coolfiles/unbalanced/' # Should end with "/"

#os.chdir('/oak/stanford/groups/tttt/users/achuthan/research/plate-c_HEK/merge_hek/coolfiles/balanced_2/')  
#cool_files_dir = '/oak/stanford/groups/tttt/users/achuthan/research/plate-c_HEK/merge_hek/coolfiles/unbalanced/' # Should end with "/"

# granule
#cool_files_dir = '/oak/stanford/groups/tttt/users/achuthan/research/plate-c_granule_plate/merge_240419a/coolfiles/unbalanced/' # Should end with "/"
# ipsc
#cool_files_dir = '/oak/stanford/groups/tttt/users/achuthan/research/plate-c_ipsc/merged_ARVA2_deep_240928a/coolfiles/unbalanced/' # Should end with "/"

ref_genome_ID = 'hg19'
resolution = 1000000
resolution_str = '1Mb'
ref_genome = '/oak/stanford/groups/tttt/users/achuthan/tools/cooltools/ref_genomes/'+ref_genome_ID+'.fa'

# For balancing - default parameters
balance_params = {'mad_max':5,      
                  'min_nnz':10,      
                  'min_count':0}

# For gc cov
bins = pd.read_csv("/oak/stanford/groups/tttt/users/achuthan/tools/cooltools/bins_sample_"+resolution_str+"_"+ref_genome_ID+".tsv", sep="\t")

# Download ref genome if it doesn't exist
if not os.path.isfile(ref_genome):
    ## note downloading a ~1Gb file can take a minute
    subprocess.call('wget --progress=bar:force:noscroll https://hgdownload.cse.ucsc.edu/goldenpath/mm10/bigZips/mm10.fa.gz', shell=True)
    subprocess.call('gunzip mm10.fa.gz', shell=True)

cool_files = [cool_files_dir+i for i in os.listdir(cool_files_dir) if i.endswith('.mcool')] 
cool_files[0]

### Check total number of files
print(len(cool_files))

############################
############ Run ###########
############################

######### IMPORTANT NOTE #########
# Modify function definition first in the important note section there
##################################

# The parallel run box
args = [(i, cool_files[i], resolution, resolution_str, balance_params, ref_genome_ID) for i in range(len(cool_files))]

# Use multiprocessing to parallelize the loop
with Pool() as pool:
    pool.map(process_cool_file, args)

def process_cool_file(args):
    file_idx, cool_file, resolution, resolution_str, balance_params, ref_genome_ID = args
    print('Starting '+str(file_idx)+'\n')
    clr = cooler.Cooler(cool_file + '::resolutions/' + str(resolution))
    
    balanced_weights = cooler.balance.balance_cooler(
        clr, 
        mad_max=balance_params['mad_max'],  
        min_nnz=balance_params['min_nnz'],    
        min_count=balance_params['min_count']
    )
    
    bins = clr.bins()[:]
    bins['weight'] = balanced_weights[0]

    ######### IMPORTANT NOTE #########
    # The following line exists only for hg19 (HEK, iPSC) and not for mm10 (granule)
    if ref_genome_ID[0] == "hg19":
        bins['chrom'] = 'chr' + bins['chrom'].astype(str)
    ##################################
    
    bins = bins[bins['chrom'] != 'chrMT']
    bins = bins[bins['chrom'] != 'chrM']
    
    balanced_file = os.path.basename(cool_file).rstrip('.mcool') + '_balanced_' + resolution_str + '.cool'
    
    cooler.create_cooler(
        balanced_file,
        pixels=clr.pixels()[:], 
        bins=bins
    )
    print('Ending '+str(file_idx)+'\n')

############################
####### Make gc_cov ########
############################

########### IMPORTANT NOTE ###########
########### This cell allows to see the bins variable and gc_cov variable. The chroms and bin divisions should match!
bins
ref_genome_bioframe = bioframe.load_fasta(ref_genome);
# Ensure, before going to the next command that the right bins variable exists (matching resolution, resolution_str)
gc_cov = bioframe.frac_gc(bins[['chrom', 'start', 'end']], ref_genome_bioframe)
gc_cov.to_csv('/oak/stanford/groups/tttt/users/achuthan/tools/cooltools/ref_genomes/'+ref_genome_ID+'_gc_cov_'+resolution_str+'.tsv',index=False,sep='\t')
display(gc_cov)
gc_cov = pd.read_csv('/oak/stanford/groups/tttt/users/achuthan/tools/cooltools/ref_genomes/'+ref_genome_ID+'_gc_cov_'+resolution_str+'.tsv',sep='\t')
