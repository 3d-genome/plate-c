import sys
import numpy as np
import pandas as pd
import matplotlib
from matplotlib import pyplot as plt
from matplotlib import gridspec
import matplotlib.colors as mcolors
from matplotlib.colors import Normalize
import warnings
warnings.filterwarnings("ignore", category=FutureWarning)
pd.options.mode.chained_assignment = None  # default='warn'

# Make sure to run this from the hicstraw environment so there's no import conflicts.
# Note that this script does not take into account width of bars, so the start and end positions should be multiples of 1Mb to match the bedgraph file.

def plot_ab_tracks_heatmap(bedgraph_filepath, chr, start_position, end_position, output_png_path,
                   ylim_a=10.0, ylim_b=-10.0, input_hic_shape=200, width_scale_factor=30, height_scale_factor=10, normalization=False):
    """
    Plots the A/B compartments over a 1D heatmap-like track.
    
    Parameters:
    bedgraph_filepath (str): File path to the bedgraph of interest. 
    chr (str): Chromosome to visualize, matching bedgraph headers, i.e. chr11. 
    start_position (int): Genomic start position to visualize. Make sure this is a multiple of 1Mb.
    end_position (int): Genomic end position to visualize. Make sure this is a multiple of 1Mb.
    output_png_path (str): File path to image that you want to save.
    ylim_a (float): Sets y-axis limits for A compartment (positive values) such that the color mapping is normalized to this.
    ylim_b (float): Sets y-axis limits for B compartment (negative values) such that the color mapping is normalized to this.
    input_hic_shape (int): Shape of the input Hi-C array, same as in plot_matrix in hicstraw. Used for image width.
    width_scale_factor (int): Scaling factor for image width. Defaults to 30 but may want higher scaling factors for longer ranges.
    height_scale_factor (int): Scaling factor for image height. Defaults to 10 but may want higher scaling factors for longer ranges.
    normalization (boolean): Whether to apply normalization by the maximum value so the plot scales -1 to 1.
    """
    # Load the BEDGraph file
    bedgraph_data = pd.read_csv(bedgraph_filepath, sep='\t', header=None, names=['chrom', 'start', 'end', 'value'])
    
    # Filter the data for the specified chromosome and section
    filtered_data = bedgraph_data[(bedgraph_data['chrom'] == chr) &
                                (bedgraph_data['start'] < end_position) & 
                                (bedgraph_data['end'] > start_position)]
    
    # Add in data for regions with no data (bedgraph data starts at 2.5Mb)
    resolution = 500000 # 1Mb resolution
    
    new_bedgraph_data = pd.DataFrame(columns=['chrom', 'start', 'end', 'value'])
    for current_position in range(start_position, end_position, resolution):
        row_exists = ((filtered_data['start'] == current_position) & (filtered_data['end'] == current_position+resolution)).any() # check if entry exists
        if (not row_exists):
            # check if either the start or end matches a known start or end
            adjacent_left_row = ((filtered_data['start'] == current_position)).any()
            adjacent_right_row = ((filtered_data['end'] == current_position+resolution)).any()
            if (adjacent_left_row):
                adjacent_left_value = filtered_data[filtered_data['start'] == current_position]['value']
                new_row = pd.DataFrame({'chrom': chr, 
                                        'start': current_position, 
                                        'end': current_position + resolution, 
                                        'value': float(adjacent_left_value.values[0])}, index = [0])
            elif (adjacent_right_row): 
                adjacent_right_value = filtered_data[filtered_data['end'] == current_position+resolution]['value']
                new_row = pd.DataFrame({'chrom': chr, 
                                        'start': current_position, 
                                        'end': current_position + resolution, 
                                        'value': float(adjacent_right_value.values[0])}, index = [0])
            else:
                new_row = pd.DataFrame({'chrom': chr, 
                                        'start': current_position, 
                                        'end': current_position + resolution, 
                                        'value': 0.0}, index = [0])
            new_bedgraph_data = pd.concat([new_bedgraph_data, new_row], ignore_index=True)
            
    # Reshape for plotting
    data = new_bedgraph_data['value']
    data_array = data.to_numpy()
    reshaped_data = data_array.reshape(1, data_array.shape[0])

    # Create a plot
    fig_width = input_hic_shape/30 # Size of input array divided by scaling factor
    fig_height = fig_width/height_scale_factor
    dpi = 300 # Default publication dpi
    plt.figure(figsize=(fig_width, fig_height), dpi=dpi)
    
    norm = Normalize(vmin=ylim_b, vmax=ylim_a)
    
    # Create a 1D heatmap
    plt.imshow(reshaped_data, 
               aspect='auto', 
               cmap= mcolors.LinearSegmentedColormap.from_list("custom_cmap", [(0.8,0,0.8), (1, 1, 1), (0, 0.8, 0)], N=256),
              norm=norm)
    
    plt.axis('off')
    plt.savefig(output_png_path, bbox_inches='tight', pad_inches=0, format='png')

if __name__ == "__main__":
    # Check if correct number of arguments provided
    if len(sys.argv) < 5:
        print("Usage: python plot_ab_tracks_heatmap.py <bedgraph_filepath> <chr> <start_position> <end_position> <output_png_path> <optional: ylim_a> <optional: ylim_b> <optional: input_hic_shape> <optional: width_scale_factor> <optional:height_scale_factor> <optional: normalization> \n")
        sys.exit(1)
    
    # Get arguments from command line
    bedgraph_filepath = sys.argv[1]
    chr = sys.argv[2]
    start_position = int(sys.argv[3])
    end_position = int(sys.argv[4])
    output_png_path = sys.argv[5]
    ylim_a = 1
    ylim_b = -1
    input_hic_shape = 200
    width_scale_factor = 30
    height_scale_factor = 10
    normalization = False
    if (len(sys.argv) > 6):
    	ylim_a = float(sys.argv[6])
    if (len(sys.argv) > 7):
    	ylim_b = float(sys.argv[7])
    if (len(sys.argv) > 8):
        input_hic_shape = int(sys.argv[8])
    if (len(sys.argv) > 9):
        width_scale_factor = int(sys.argv[9])
    if (len(sys.argv) > 10):
        height_scale_factor = int(sys.argv[10])
    if (len(sys.argv) > 11):
    	normalization = bool(sys.argv[11])
    
    # Call the function
    plot_ab_tracks_heatmap(bedgraph_filepath, chr, start_position, end_position, output_png_path,
                   ylim_a, ylim_b, input_hic_shape, width_scale_factor, height_scale_factor, normalization)
