%% load scA/B matrix
% for Dip-C mouse in vivo experiments

clear
clc

folder = '/Users/tanlongzhi/Research/dip-c/brain2/summary/colors';

filename = 'tsa_in_vivo_50k_240927a.cpg_b1m.color2s.txt'; % primary in vivo experiment
% filename = 'tsa_injections_240829b.cpg_b1m.color2s.txt'; % different dosing schemes
% filename = 'tsa_time_course_50k_241016a.cpg_b1m.color2s.txt'; % single-dose time course

% load data
temp_data = importdata([folder, '/', filename]);
ref_names = temp_data.data(:, 1);
ref_loci = temp_data.data(:, 2);
cell_names = temp_data.textdata(1, 3:end);
color_data = temp_data.data(:, 3:end);
clear temp_data;

% reformat cell names
for i = 1:size(cell_names, 2)
    cell_name_data = strsplit(cell_names{i}, {'/', '.'});
    cell_names{i} = cell_name_data{end-2};
end

% filter cells
cell_ids = 1:size(cell_names, 2);
cell_names = cell_names(cell_ids);
color_data = color_data(:, cell_ids);

% filter (and sort) genomic loci
locus_ids = 1:size(ref_loci, 1);
locus_ids = find(min(color_data, [], 2) >= 0 & ref_names < 20); % mouse autosome, no missing values for for all cells
[C,I]=sort(ref_names(locus_ids));
locus_ids = locus_ids(I);
ref_names = ref_names(locus_ids);
ref_loci = ref_loci(locus_ids);
color_data = color_data(locus_ids, :);

% find chromosome starts
chr_start_ids = [];
chr_start_names = [];
for i = 1:size(ref_loci, 1)
    if i == 1 || ref_names(i - 1) ~= ref_names(i)
        chr_start_ids = [chr_start_ids; i];
        chr_start_names = [chr_start_names; ref_names(i)];
    end
end

% % plot raw scA/B matrix
% figure();
% imagesc(color_data);
% %caxis([0.005, 0.02]);
% %caxis([0, 1]);
% %caxis([0,100]);
% %caxis([-1, 10]);
% xticks(1:size(cell_names, 2));
% xticklabels(cell_names);
% xtickangle(90);
% yticks(chr_start_ids);
% yticklabels(chr_start_names);
% set(gca, 'ticklabelinterpreter', 'none');

% rank-normalize each cell to 0 to 1
for i = 1:size(color_data, 2)
    disp(['normalizing cell ', num2str(i), ': ', cell_names{i}]);
    rows = find(color_data(:, i) >= 0);
    cell_data = color_data(rows, i);
    % rank normalize, ok to have missing data
    color_data(rows, i) = (tiedrank(cell_data) - 1)./(size(cell_data, 1) - 1);
end

% % plot normalized scA/B matrix
% figure();
% imagesc(color_data);
% %imagesc(color_data(I2,I1));
% caxis([0, 1]);
% xticks(1:size(cell_names, 2));
% xticklabels(cell_names);
% xtickangle(90);
% yticks(chr_start_ids);
% yticklabels(chr_start_names);
% set(gca, 'ticklabelinterpreter', 'none');


% PCA

% specify which genomic loci
pca_rows = 1:size(color_data,1);
% specific which cells
pca_cols = 1:size(color_data,2);
% run PCA
[coeff, score, latent, tsquared, explained, mu] = pca(color_data(pca_rows, pca_cols)');

% project all data onto the PCA of the specified cells and genomic loci
score=(color_data' - ones(size(color_data'))*diag(mu))/coeff';
raw_score=score;
disp('pca done');

% pick PC space
pc_space = raw_score(:,1:20); % first 20 PCs for Dip-C experiments


% UMAP with PCA initialization
umap_score = run_umap(pc_space,'init',raw_score(:,1:2),'randomize','false','method','c','verbose','text');
score = umap_score;
disp('UMAP done');

% initialize colors for plotting
scatter_colors=0;

%% plot UMAP

plot_filename = 'output.pdf';

% plot and dot sizes
plot_width = 2;plot_height = 2;dot_size = 2; % TSA primary plot 241115a
% plot_width = 3;plot_height = 3;dot_size = 2; % TSA injections plot 240908a
% plot_width = 2.5;plot_height = 2.5;dot_size = 2; % TSA time course plot 241117a

% set dot size
scatter_sizes=zeros(size(color_data,2),1)+dot_size;

% plot
h = figure('color', 'w');
set(h, 'visible', 'off');

hold on;

% color by mean scA/B of a list of genomic regions

title('scA/B');
color_row_info = [];

% specify list of genomic regions
% color_row_info = load('tsa_primary_granule_up.b1m.txt', '-ascii'); c_min=mean(scatter_colors)-0.04; c_max=mean(scatter_colors)+0.04;
% color_row_info = load('tsa_primary_granule_down.b1m.txt', '-ascii'); c_min=mean(scatter_colors)-0.04; c_max=mean(scatter_colors)+0.04;

% calculate mean scA/B of list
color_rows = [];
for i = 1:size(color_row_info, 1)
    %disp(color_row_info(i, :));
    color_row = find(ref_names == color_row_info(i, 1) & ref_loci == color_row_info(i, 2));
    if size(color_row, 1) == 1
        color_rows(end+1, :) = color_row;
    end
end
scatter_colors=mean(color_data(color_rows, :),1);


% alternative: color cell by loading a vector

title('Treatment in Primary In Vivo Experiment'); scatter_colors = 1.*contains(cell_names,'dmso_replicate')+2.*contains(cell_names,'tsa_replicate'); c_min=-0.1; c_max=2;
% title('Different Dosing Schemes'); scatter_colors = 1.*contains(cell_names,'vvvv')+2.*contains(cell_names,'dvvv')+3.*contains(cell_names,'ddvv')+4.*contains(cell_names,'dddv')+5.*contains(cell_names,'vvdd')+6.*contains(cell_names,'vddd')+7.*contains(cell_names,'dddd'); c_min=-0.2; c_max=7;
% title('Single-dose Time Course'); scatter_colors = 1.*contains(cell_names,'dmso_01h')+2.*contains(cell_names,'dmso_02h')+3.*contains(cell_names,'dmso_04h')+4.*contains(cell_names,'dmso_24h')+5.*contains(cell_names,'tsa_01h')+6.*contains(cell_names,'tsa_02h')+7.*contains(cell_names,'tsa_04h')+8.*contains(cell_names,'tsa_24h'); c_min=-0.2; c_max=8;
% title('Replicates for Primary Experiment'); scatter_colors = 1.*contains(cell_names,'dmso_replicate_01')+2.*contains(cell_names,'dmso_replicate_02')+3.*contains(cell_names,'dmso_replicate_03')+4.*contains(cell_names,'tsa_replicate_01')+5.*contains(cell_names,'tsa_replicate_02')+6.*contains(cell_names,'tsa_replicate_03'); c_min=1; c_max=6;

% title('Cell Type in Primary In Vivo Experiment'); scatter_colors = readcell('/Users/tanlongzhi/Research/dip-c/brain2/aux_data/ward6_tsa_primary_241115a_merged.txt','Delimiter','\t');scatter_colors=cell2mat(scatter_colors(:,2));c_min=1;c_max=6;
% title('Cell Type in Different Dosing Schemes'); scatter_colors = readcell('/Users/tanlongzhi/Research/dip-c/brain2/aux_data/ward17_tsa_time_course_241117a_merged.txt','Delimiter','\t');scatter_colors=cell2mat(scatter_colors(:,2));c_min=1;c_max=6;
% title('Cell Type in Single-dose Time Course'); scatter_colors = readcell('/Users/tanlongzhi/Research/dip-c/brain2/aux_data/ward25_tsa_injections_241117a_merged.txt','Delimiter','\t');scatter_colors=cell2mat(scatter_colors(:,2));c_min=1;c_max=6;


% specify which cells to plot
rows = 1:size(color_data,2);
rows = randperm(size(score,1)); % randomly order to avoid overlaps

% main plotting command
scatter(score(rows, 1), score(rows, 2), scatter_sizes(rows), scatter_colors(rows), 'filled');

% specific color ranges
caxis([c_min,c_max]);
% colorbar;

% custom color map (green to magenta)
low_color = [1,0,1];high_color = [0,1,0];
num_colors = 256;
low_to_high = [];
for i=1:3
    low_to_high = [low_to_high, linspace(low_color(i), high_color(i), num_colors)'];
end

% custom color map (others)
cluster_colors=flipud(parula(256));cluster_colors=cluster_colors(floor(256*0.05):end,:); % parula yellow: low, blue: high
% cluster_colors=[252,204,255; 217,166,255; 200,150,100; 71,191,141; 138,231,255; 255,82,0; 186,186,186]./255; % cell type clusters 240512a
% cluster_colors=[209,209,209; 199,199,199; 179,179,179; 159,159,159; 255,187,99; 255,157,31; 239,134,0; 171,117,48]./255; % tsa in vivo time course 241117a
% cluster_colors=[0.75,0.75,0.75;0.65,0.65,0.65;0.55,0.55,0.55;0.45,0.45,0.45;0.6,0,0;0.8,0,0;1,0,0]; % TSA injections 240908a

% choose a color map
colormap(cluster_colors);
% colormap(low_to_high);


% format
box on;
axis equal;

xlim([min(score(rows,1)), max(score(rows,1))]+(max(score(rows,1)) - min(score(rows,1)))*0.05*[-1, 1]);
ylim([min(score(rows,2)), max(score(rows,2))]+(max(score(rows,2)) - min(score(rows,2)))*0.05*[-1, 1]);
xticks([]);
yticks([]);

ax = gca;
ax.XColor = [0,0,0];
ax.YColor = [0,0,0];

% print
set(gcf,...
'paperunits', 'inches',...
'papersize', [plot_width, plot_height],...
'paperposition', [0, 0, plot_width, plot_height],...
'units', 'inches',...
'position', [0, 0, plot_width, plot_height]);
print('-dpdf', '-painters', [folder, '/', plot_filename]);
close();
disp('printed');


%% hierachical clustering with Wald's method

distance_cutoff=6; % for TSA primary 241115a
% distance_cutoff=17; % for TSA time course 241115a
% distance_cutoff=25; % for TSA injections 241117a


linkage_data = linkage(pc_space, 'ward', 'euclidean');
cluster_data = cluster(linkage_data,'CutOff',distance_cutoff,'Criterion','distance');

temp_cluster_data=cluster_data;
temp_linkage_data=linkage_data;
temp_cluster_ids=1:size(cell_names,2);

tmp = temp_cluster_data;

% manually and annotate merge clusters

% % TSA primary 241115a
% temp_cluster_data=temp_cluster_data.*0+6; % 0 = unknown
% temp_cluster_data(tmp==9|tmp==27|tmp==2|tmp==14|tmp==29|tmp==13|tmp==21)=1; % 1 = cerebellar granule cells
% temp_cluster_data(tmp==16|tmp==25|tmp==22|tmp==8|tmp==6|tmp==12|tmp==7|tmp==3)=2; % 2 = other neurons
% temp_cluster_data(tmp==4)=3; % 3 = astrocytes
% temp_cluster_data(tmp==28|tmp==11)=4; % 4 = oligodendrocytes
% temp_cluster_data(tmp==1|tmp==20|tmp==10)=5; % 5 = microglia

% % TSA time course 241117a
% temp_cluster_data=temp_cluster_data.*0+6;
% temp_cluster_data(tmp==11|tmp==3|tmp==6)=1;
% temp_cluster_data(tmp==9|tmp==4|tmp==2|tmp==10)=2;
% temp_cluster_data(tmp==13)=3;
% temp_cluster_data(tmp==1|tmp==7)=4;
% temp_cluster_data(tmp==8)=5; % microglia

% % TSA injections 241117a
% temp_cluster_data=temp_cluster_data.*0+6;
% temp_cluster_data(tmp==3|tmp==4)=1;
% temp_cluster_data(tmp==8|tmp==6)=2;
% temp_cluster_data(tmp==7)=3;
% temp_cluster_data(tmp==1)=4;
% temp_cluster_data(tmp==5)=5; % microglia

cluster_data=temp_cluster_data;

% write to file
outfile=fopen('ward.txt','w');

for i = 1:size(cell_names,2)
    %fprintf(outfile, [num2str(cluster_data(i)),'\n']);
    fprintf(outfile, [cell_names{i},'\t',num2str(cluster_data(i)),'\n']);
end


fclose(outfile);
disp('written');

% visualize cluster
figure('color','w');
subplot('position',[0 0.5 0.5 0.5]);
[H,~,outperm]=dendrogram(temp_linkage_data,0,'ColorThreshold',distance_cutoff);
set(H,'LineWidth',2);
axis off;

subplot('position',[0 0 0.5 0.5]);
plot(temp_cluster_data(outperm),'k.');
for cluster_id=unique(temp_cluster_data)'
    sub_rows=find(temp_cluster_data(outperm)==cluster_id);
    text(mean(sub_rows), cluster_id+0.5, num2str(cluster_id), 'fontsize', 10);
end
xlim([1,size(temp_cluster_data,1)]);
ylim([0,max(temp_cluster_data)+1]);
yticks(1:100);
grid on;

subplot('position',[0.5 0 0.5 1]);
rows=temp_cluster_ids;
scatter(score(rows, 1), score(rows, 2), 100, temp_cluster_data,'filled');
% text(score(rows, 1), score(rows, 2), num2str(temp_cluster_data), 'fontsize', 10);
for cluster_id=unique(temp_cluster_data)'
    sub_rows=find(temp_cluster_data==cluster_id);
    text(mean(score(rows(sub_rows), 1),1), mean(score(rows(sub_rows), 2),1), num2str(cluster_id), 'fontsize', 10);
end
axis equal;
axis off;


%% find differential scA/B regions (can have missing values)

clc

% significance thresholds
min_diff = 0.02; max_fdr = 0.01; % for scA/B, TSA primary 241116a

% initial axis limits
x_lims = [-1,1]*0.12; y_lims = [0,20]; % time course
x_lims = [-1,1]*0.18; y_lims = [0,30]; % primary

plot_width = 0.7; plot_height = 0.7; dot_size=3;
min_num_1 = 1;
min_num_2 = 1;

% load cell types
cell_type_data = readcell('/Users/tanlongzhi/Research/dip-c/brain2/aux_data/ward6_tsa_primary_241115a_merged.txt','Delimiter','\t'); cell_type_data=cell2mat(cell_type_data(:,2));
% cell_type_data = readcell('/Users/tanlongzhi/Research/dip-c/brain2/aux_data/ward25_tsa_injections_241117a_merged.txt','Delimiter','\t'); cell_type_data=cell2mat(cell_type_data(:,2));
% cell_type_data = readcell('/Users/tanlongzhi/Research/dip-c/brain2/aux_data/ward17_tsa_time_course_241117a_merged.txt','Delimiter','\t'); cell_type_data=cell2mat(cell_type_data(:,2));

% specify cells to compare
output_prefix='tsa_primary_granule';cell_ids_1=find((cell_type_data==1)&contains(cell_names,'tsa_')');cell_ids_2=find((cell_type_data==1)&contains(cell_names,'dmso_')');

% output_prefix='tsa_injections_granule';cell_ids_1=find((cell_type_data==1)&contains(cell_names,'d_replicate')');cell_ids_2=find((cell_type_data==1)&contains(cell_names,'v_replicate')');
% output_prefix='tsa_injections_other_neuron';cell_ids_1=find((cell_type_data==2)&contains(cell_names,'d_replicate')');cell_ids_2=find((cell_type_data==2)&contains(cell_names,'v_replicate')');

% output_prefix='tsa_time_course_granule_1_2_4h';cell_ids_1=find((cell_type_data==1)&contains(cell_names,'tsa_0')');cell_ids_2=find((cell_type_data==1)&contains(cell_names,'dmso_0')');
% output_prefix='tsa_time_course_other_neuron_1_2_4h';cell_ids_1=find((cell_type_data==2)&contains(cell_names,'tsa_0')');cell_ids_2=find((cell_type_data==2)&contains(cell_names,'dmso_0')');
% output_prefix='tsa_time_course_granule_1h';cell_ids_1=find((cell_type_data==1)&contains(cell_names,'tsa_01h')');cell_ids_2=find((cell_type_data==1)&contains(cell_names,'dmso_01h')');
% output_prefix='tsa_time_course_granule_2h';cell_ids_1=find((cell_type_data==1)&contains(cell_names,'tsa_02h')');cell_ids_2=find((cell_type_data==1)&contains(cell_names,'dmso_02h')');
% output_prefix='tsa_time_course_granule_4h';cell_ids_1=find((cell_type_data==1)&contains(cell_names,'tsa_04h')');cell_ids_2=find((cell_type_data==1)&contains(cell_names,'dmso_04h')');


disp(['cell type 1 (treatment): ',num2str(length(cell_ids_1)),' cells']);
disp(['cell type 2 (control): ',num2str(length(cell_ids_2)),' cells']);

% get mean values
% for missing value (-1)
nan_color_data = color_data;
nan_color_data(color_data<0) = nan;
mean_colors_1 = nanmean(nan_color_data(:, cell_ids_1), 2);
mean_colors_2 = nanmean(nan_color_data(:, cell_ids_2), 2);
mean_colors_diff = mean_colors_1 - mean_colors_2;

% perform tests
p_values = [];
for i = 1:size(color_data,1)
    if sum(color_data(i, cell_ids_1)>=0) >= min_num_1 & sum(color_data(i, cell_ids_2)>=0) >= min_num_2
        [h, p] = ttest2(nan_color_data(i, cell_ids_1),nan_color_data(i, cell_ids_2)); % t test, equal variance
    else
        p = nan; % flag for too many missing data
    end
    p_values = [p_values; p];
end
FDR = mafdr(p_values,'BHFDR',true);

% rescale axis if needed
if max(-log10(p_values)) > y_lims(2)
    y_lims(2) = max(-log10(p_values));
    disp(['rescaled y axis to :', num2str(y_lims(2))]);
end
if max(abs(mean_colors_diff)) > x_lims(2)
    x_lims = [-1,1] * max(abs(mean_colors_diff));
    disp(['rescaled x axis to:', num2str(x_lims(2))]);
end

% sort and calculate significance
significant_rows_up = find(FDR<max_fdr & mean_colors_diff > min_diff);
significant_rows_down = find(FDR<max_fdr & mean_colors_diff < -min_diff);
disp(['significant loci with FDR < ', num2str(max_fdr), ', min difference > ', num2str(min_diff), ' (up): ', num2str(length(significant_rows_up)), '  loci']);
disp(['significant loci with FDR < ', num2str(max_fdr), ', min difference < ', num2str(-min_diff), ' (down): ', num2str(length(significant_rows_down)), '  loci']);

% plot
h = figure('color', 'w');
set(h, 'visible', 'off');

scatter(mean_colors_diff (p_values <=1), -log10(p_values(p_values <=1)),dot_size,[1,1,1]*0.7,'filled');
hold on;

red_rows = [significant_rows_up; significant_rows_down];
scatter(mean_colors_diff(red_rows), -log10(p_values(red_rows)),dot_size,[1,0,0],'filled');

xlim(x_lims);
ylim(y_lims);

axis off

set(gcf,...
'paperunits','inches',...
    'papersize',[plot_width,plot_height],...
'paperposition',[0,0,plot_width,plot_height],...
'units','inches',...
'position',[0,0,plot_width,plot_height],...
'inverthardcopy','off');
set(gca, 'units', 'inches', 'position', [0, 0, plot_width, plot_height]); % axis full
print('-dpng', '-image','-r600', [output_prefix,'.png']);
display('Done!');
close();

% write output
bin_size =1e6; % 1Mb bins

outfile=fopen([output_prefix,'.bed'],'w');
outfile2_up=fopen([output_prefix,'_up.b1m.txt'],'w');
outfile2_down=fopen([output_prefix,'_down.b1m.txt'],'w');
outfile3=fopen([output_prefix,'.diff.txt'],'w');
outfile4=fopen([output_prefix,'.all_b1m_diff.txt'],'w');
outfile5=fopen([output_prefix,'.supp_table.txt'],'w');
outfile6=fopen([output_prefix,'.signed_significance.bedgraph'],'w');

% optional: write protein-coding genes in each genomic region
gene_midpoint_data=importdata('/Users/tanlongzhi/R/plate-c/gene_position.mm10.vM25.midpoint_matlab_protein_coding.txt');
gene_midpoint_chr_data=gene_midpoint_data.data(:,1);
gene_midpoint_loci_data=gene_midpoint_data.data(:,2);
gene_midpoint_name_data=gene_midpoint_data.textdata;


fprintf(outfile3,'#chr\tstart\tend\tmean_1\tmean_2\tmean_diff\t-log10(p_value)\n');
fprintf(outfile5,'#chr\tstart\tend\tscab_1\tscab_2\tscab_diff\tminus_log10_p\tminus_log10_fdr\tgenes\n');
for i = 1:length(red_rows)
    row = red_rows(i);
    
    % optional: find genes
    gene_rows = find(gene_midpoint_chr_data==ref_names(row) & abs(gene_midpoint_loci_data-ref_loci(row))<=0.5*bin_size);
    gene_string = strjoin(gene_midpoint_name_data(gene_rows),',');
    
    fprintf(outfile, ['chr',num2str(ref_names(row)),'\t',num2str(ref_loci(row)-0.5*bin_size),'\t',num2str(ref_loci(row)+0.5*bin_size),'\n']);
    if mean_colors_diff(row) > 0
        fprintf(outfile2_up, [num2str(ref_names(row)),'\t',num2str(ref_loci(row)),'\n']);
    elseif mean_colors_diff(row) < 0
        fprintf(outfile2_down, [num2str(ref_names(row)),'\t',num2str(ref_loci(row)),'\n']);
    end
    fprintf(outfile3, [num2str(ref_names(row)),'\t',num2str(ref_loci(row)),'\t',num2str(mean_colors_1(row)),'\t',num2str(mean_colors_2(row)),'\t',num2str(mean_colors_diff(row)),'\t',num2str(-log10(p_values(row))),'\n']);
    fprintf(outfile5, ['chr',num2str(ref_names(row)),'\t',num2str(ref_loci(row)-0.5*bin_size),'\t',num2str(ref_loci(row)+0.5*bin_size),'\t',num2str(mean_colors_1(row)),'\t',num2str(mean_colors_2(row)),'\t',num2str(mean_colors_diff(row)),'\t',num2str(-log10(p_values(row))),'\t',num2str(-log10(FDR(row))),'\t',gene_string,'\n']);
end

for row = 1:length(mean_colors_diff)
    fprintf(outfile4, [num2str(ref_names(row)),'\t',num2str(ref_loci(row)),'\t',num2str(mean_colors_diff(row)),'\n']);
    fprintf(outfile6, ['chr',num2str(ref_names(row)),'\t',num2str(ref_loci(row)-0.5*bin_size),'\t',num2str(ref_loci(row)+0.5*bin_size),'\t',num2str(sign(mean_colors_1(row)-mean_colors_2(row)).*-log10(FDR(row))),'\n']);
end

fclose(outfile);
fclose(outfile2_up);
fclose(outfile2_down);
fclose(outfile3);
fclose(outfile4);
fclose(outfile5);
fclose(outfile6);
disp('written');


%% plot scA/B mean +/- SEM by category, TSA primary

clc

spread_width=0.8;

y_lims=[0,1];

y_lims=[0.395,0.435];
% y_lims=[0.455,0.495];
% y_lims=[0.35,0.39];
% y_lims=[0.54,0.58];
% 
% y_lims=[0.395,0.465];
% 
y_lims=[0.47,0.53];
% y_lims=[0.465,0.525];
% y_lims=[0.485,0.545];
% y_lims=[0.385,0.445];

y_lims=[0.47,0.51]; % primary granule up, in injections
y_lims=[0.49,0.53]; % primary granule down, in injections


y_lims=[0.64,0.67]; %
y_lims=[0.60,0.63]; %
y_lims=[0.65,0.68]; %
y_lims=[0.60,0.63]; %

y_lims=[0.69,0.72]; %
y_lims=[0.64,0.67]; %
y_lims=[0.66,0.69]; %
y_lims=[0.653,0.673]; %

y_lims=[0.54,0.6]; %



plotwidth=1;plotheight=0.8;% for plots

barthickness=0.5;
boundary_amount=1e-5;
plot_filename='types.pdf';

h=figure('color','none');
set(h,'Visible','off');

% cluster colors
cluster_colors = zeros([10,3]);

% load cell type assignments


% for TSA primary
cell_type_data = readcell('/Users/tanlongzhi/Research/dip-c/brain2/aux_data/ward6_tsa_primary_241115a_merged.txt','Delimiter','\t'); cell_type_data=cell2mat(cell_type_data(:,2));
replicate_data = contains(cell_names,'dmso_replicate_01')'*1  + contains(cell_names,'dmso_replicate_02')'*2 + contains(cell_names,'dmso_replicate_03')'*3 + contains(cell_names,'tsa_replicate_01')'*4  + contains(cell_names,'tsa_replicate_02')'*5 + contains(cell_names,'tsa_replicate_03')'*6;
cluster_data = replicate_data.* (cell_type_data==1); % granule
baseline_data = cluster_data<=3 & cluster_data>=1;



% print cluster info
tabulate(cluster_data);


classifications=cluster_data;
clusters_to_plot=1:max(classifications);
centers_to_plot=1:max(classifications);
widths_to_plot=spread_width*ones([max(classifications),1]);

% plot baseline level
baseline_color = mean(scatter_colors(find(baseline_data)));
line([0,max(classifications)+1],[1,1]*baseline_color,'color',[1,1,1]*0.7,'linewidth',barthickness);

for i=1:size(clusters_to_plot,2)
    %classification_name = classification_names{i};
    %disp(classification_name);
    cluster_to_plot=clusters_to_plot(i);
    plot_color=[0,0,0];
    disp(cluster_to_plot);
    rows=find(classifications==cluster_to_plot);
    data_to_plot=scatter_colors(rows)';
    disp(['  n = ',num2str(length(rows))]);
    disp(['  mean = ',num2str(mean(data_to_plot))]);
    disp(['  std = ',num2str(std(data_to_plot))]);

    
    if length(rows)<=2
        continue
    end
    
    % plot mean
    middle=mean(data_to_plot);
    %middle=median(data_to_plot);
    line(centers_to_plot(i)+[-1,1]*widths_to_plot(i)/2,[1,1]*middle,'color',[1,1,1]*0,'linewidth',barthickness);

    % plot errorbar
    line(centers_to_plot(i)*[1,1],middle+std(data_to_plot)./sqrt(length(data_to_plot))*[-1, 1],'color',[1,1,1]*0,'linewidth',barthickness);

    
    hold on;
end

box on;

xlim([0.5,max(cluster_data)+0.5]);
ylim(y_lims);
ax=gca;
ax.XTick=[];
ax.XColor = [0,0,0];
ax.YColor = [0,0,0];
ax.YTick=[];
 

set(gcf,...
'paperunits','inches',...
'papersize',[plotwidth,plotheight],...
'paperposition',[0,0,plotwidth,plotheight],...
'units','inches',...
'position',[0,0,plotwidth,plotheight],...
'inverthardcopy','off');
% set(gca, 'units', 'inches', 'position', [0, 0, plotwidth, plotheight]); % axis full
print('-dpdf', '-painters', [folder, '/', plot_filename]);
display('Done!');
close();


% get p-value info from single cells
[h,p]=ttest2(scatter_colors(find(cluster_data<=3 & cluster_data>=1)), scatter_colors(find(cluster_data<=6 & cluster_data>=4)));
mean(scatter_colors(find(cluster_data<=6 & cluster_data>=4)))-mean(scatter_colors(find(cluster_data<=3 & cluster_data>=1)))
disp(p);

% get p-value info from pseudo-bulks
[h,p]=ttest2([mean(scatter_colors(find(cluster_data==1))),mean(scatter_colors(find(cluster_data==2))),mean(scatter_colors(find(cluster_data==3)))], [mean(scatter_colors(find(cluster_data==4))),mean(scatter_colors(find(cluster_data==5))),mean(scatter_colors(find(cluster_data==6)))]);
disp(p);


%% plot scA/B mean +/- SEM by category, TSA injections

clc

spread_width=0.8;

y_lims=[0,1];

y_lims=[0.395,0.435];
% y_lims=[0.455,0.495];
% y_lims=[0.35,0.39];
% y_lims=[0.54,0.58];
% 
% y_lims=[0.395,0.465];
% 
y_lims=[0.47,0.53];
% y_lims=[0.465,0.525];
% y_lims=[0.485,0.545];
% y_lims=[0.385,0.445];

y_lims=[0.47,0.51]; % primary granule up, in injections
% y_lims=[0.49,0.53]; % primary granule down, in injections


% y_lims=[0.64,0.67]; %
% y_lims=[0.60,0.63]; %

% y_lims=[0.69,0.72]; %

% y_lims=[0.653,0.663]; %


plotwidth=1;plotheight=0.8;% for plots

barthickness=0.5;
boundary_amount=1e-5;
plot_filename='types.pdf';

h=figure('color','none');
set(h,'Visible','off');

% cluster colors
cluster_colors = zeros([10,3]);

% load cell type assignments
% for TSA injections
cell_type_data = readcell('/Users/tanlongzhi/Research/dip-c/brain2/aux_data/ward25_tsa_injections_241117a_merged.txt','Delimiter','\t'); cell_type_data=cell2mat(cell_type_data(:,2));
cluster_data = contains(cell_names,'dddd') + contains(cell_names,'vddd').*2 + contains(cell_names,'vvdd').*3 + contains(cell_names,'dddv').*4 + contains(cell_names,'ddvv').*5 + contains(cell_names,'dvvv').*6 + contains(cell_names,'vvvv').*7;
cluster_data = cluster_data.*(cell_type_data==1)';
baseline_data = contains(cell_names,'vvvv');
baseline_data = baseline_data.*(cell_type_data==1)';

% % for TSA primary
% cell_type_data = readcell('/Users/tanlongzhi/Research/dip-c/brain2/aux_data/ward6_tsa_primary_241115a_merged.txt','Delimiter','\t'); cell_type_data=cell2mat(cell_type_data(:,2));
% replicate_data = contains(cell_names,'dmso_replicate_01')'*1  + contains(cell_names,'dmso_replicate_02')'*2 + contains(cell_names,'dmso_replicate_03')'*3 + contains(cell_names,'tsa_replicate_01')'*4  + contains(cell_names,'tsa_replicate_02')'*5 + contains(cell_names,'tsa_replicate_03')'*6;
% cluster_data = replicate_data + (cell_type_data-1)*6;
% baseline_data = cluster_data<=3;

% % for TSA time course
% cell_type_data = readcell('/Users/tanlongzhi/Research/dip-c/brain2/aux_data/ward17_tsa_time_course_241117a_merged.txt','Delimiter','\t'); cell_type_data=cell2mat(cell_type_data(:,2));
% condition_data = 5.*contains(cell_names,'dmso_01h')+6.*contains(cell_names,'dmso_02h')+7.*contains(cell_names,'dmso_04h')+8.*contains(cell_names,'dmso_24h')+1.*contains(cell_names,'tsa_01h')+2.*contains(cell_names,'tsa_02h')+3.*contains(cell_names,'tsa_04h')+4.*contains(cell_names,'tsa_24h');
% cluster_data = condition_data.*(cell_type_data==2)';
% baseline_data = cluster_data==5;


% print cluster info
tabulate(cluster_data);


classifications=cluster_data;
clusters_to_plot=1:max(classifications);
centers_to_plot=1:max(classifications);
widths_to_plot=spread_width*ones([max(classifications),1]);

% plot baseline level
baseline_color = mean(scatter_colors(find(baseline_data)));
line([0,max(classifications)+1],[1,1]*baseline_color,'color',[1,1,1]*0.7,'linewidth',barthickness);

for i=1:size(clusters_to_plot,2)
    %classification_name = classification_names{i};
    %disp(classification_name);
    cluster_to_plot=clusters_to_plot(i);
    plot_color=[0,0,0];
    disp(cluster_to_plot);
    rows=find(classifications==cluster_to_plot);
    data_to_plot=scatter_colors(rows)';
    disp(['  n = ',num2str(length(rows))]);
    disp(['  mean = ',num2str(mean(data_to_plot))]);
    disp(['  std = ',num2str(std(data_to_plot))]);
    disp(['  median = ',num2str(median(data_to_plot))]);
    disp(['  >=3 = ',num2str(mean(data_to_plot>=3))]);
    disp(['  >=5 = ',num2str(mean(data_to_plot>=5))]);
    
    if length(rows)<=2
        continue
    end
    
    % plot mean
    middle=mean(data_to_plot);
    %middle=median(data_to_plot);
    line(centers_to_plot(i)+[-1,1]*widths_to_plot(i)/2,[1,1]*middle,'color',[1,1,1]*0,'linewidth',barthickness);

    % plot errorbar
    line(centers_to_plot(i)*[1,1],middle+std(data_to_plot)./sqrt(length(data_to_plot))*[-1, 1],'color',[1,1,1]*0,'linewidth',barthickness);

    
    hold on;
end

% t-test
for i=1:6
    disp(i)
    disp(mean(scatter_colors(find(classifications==i)))-baseline_color);
    [h, p] = ttest2(scatter_colors(find(classifications==i)),scatter_colors(find(baseline_data))); % t test
    disp(p);
end

box on;

xlim([0.5,max(cluster_data)+0.5]);
ylim(y_lims);
ax=gca;
ax.XTick=[];
ax.XColor = [0,0,0];
ax.YColor = [0,0,0];
ax.YTick=[];
 

set(gcf,...
'paperunits','inches',...
'papersize',[plotwidth,plotheight],...
'paperposition',[0,0,plotwidth,plotheight],...
'units','inches',...
'position',[0,0,plotwidth,plotheight],...
'inverthardcopy','off');
% set(gca, 'units', 'inches', 'position', [0, 0, plotwidth, plotheight]); % axis full
print('-dpdf', '-painters', [folder, '/', plot_filename]);
display('Done!');
close();

% get p-value info from single cells
[h,p]=ttest2(scatter_colors(find(cluster_data<=3 & cluster_data>=1)), scatter_colors(find(cluster_data<=7 & cluster_data>=4)));
mean(scatter_colors(find(cluster_data<=3 & cluster_data>=1)))- mean(scatter_colors(find(cluster_data<=7 & cluster_data>=4)))
disp(p);

% get p-value info from pseudo-bulks
[h,p]=ttest2([mean(scatter_colors(find(cluster_data==1))),mean(scatter_colors(find(cluster_data==2))),mean(scatter_colors(find(cluster_data==3)))], [mean(scatter_colors(find(cluster_data==4))),mean(scatter_colors(find(cluster_data==5))),mean(scatter_colors(find(cluster_data==6))),mean(scatter_colors(find(cluster_data==7)))]);
disp(p);


%% plot scA/B mean +/- SEM by category, TSA time course

clc

spread_width=0.8;

y_lims=[0,1];

y_lims=[0.395,0.435];
% y_lims=[0.455,0.495];
% y_lims=[0.35,0.39];
% y_lims=[0.54,0.58];
% 
% y_lims=[0.395,0.465];
% 
y_lims=[0.47,0.53];
% y_lims=[0.485,0.545];
% y_lims=[0.385,0.445];

y_lims=[0.47,0.51]; % primary granule up, in injections
y_lims=[0.49,0.53]; % primary granule down, in injections

y_lims=[0.47,0.53]; % primary granule up & down, in time course

% y_lims=[0.66,0.69]; % 

% y_lims=[0.62,0.65]; % 
y_lims=[0.65,0.68]; % 

y_lims=[0.60,0.63]; % 

y_lims=[0.69,0.72]; % 

y_lims=[0.653,0.673]; % 

y_lims=[0.55,0.61]; % 
y_lims=[0.42,0.48]; % 
y_lims=[0.36,0.42]; % 
y_lims=[0.55,0.62]; % 
y_lims=[0.47,0.54]; % 
y_lims=[0.41,0.48]; % 
y_lims=[0.38,0.45]; % 
y_lims=[0.50,0.57]; % 
y_lims=[0.41,0.48]; % 
y_lims=[0.35,0.42]; % 



plotwidth=0.6;plotheight=0.8;% for plots

barthickness=0.5;
boundary_amount=1e-5;
plot_filename='types.pdf';

h=figure('color','none');
set(h,'Visible','off');

% cluster colors
cluster_colors = zeros([10,3]);

% load cell type assignments
% for TSA time course
cell_type_data = readcell('/Users/tanlongzhi/Research/dip-c/brain2/aux_data/ward17_tsa_time_course_241117a_merged.txt','Delimiter','\t'); cell_type_data=cell2mat(cell_type_data(:,2));
condition_data = 5.*contains(cell_names,'dmso_01h')+6.*contains(cell_names,'dmso_02h')+7.*contains(cell_names,'dmso_04h')+8.*contains(cell_names,'dmso_24h')+1.*contains(cell_names,'tsa_01h')+2.*contains(cell_names,'tsa_02h')+3.*contains(cell_names,'tsa_04h')+4.*contains(cell_names,'tsa_24h');
cluster_data = condition_data.*(cell_type_data==2)'; % cell type to plot
baseline_data = cluster_data==5;


% print cluster info
tabulate(cluster_data);


classifications=cluster_data;
clusters_to_plot=1:max(classifications);
centers_to_plot=[1:4,1:4];
widths_to_plot=spread_width*ones([max(classifications),1]);

% plot baseline level
baseline_color = mean(scatter_colors(find(baseline_data)));
line([0,max(classifications)+1],[1,1]*baseline_color,'color',[1,1,1]*0.7,'linewidth',barthickness);

for i=1:size(clusters_to_plot,2)
    %classification_name = classification_names{i};
    %disp(classification_name);
    cluster_to_plot=clusters_to_plot(i);
    plot_color=[0,0,0];
    disp(cluster_to_plot);
    rows=find(classifications==cluster_to_plot);
    data_to_plot=scatter_colors(rows)';
    disp(['  n = ',num2str(length(rows))]);
    disp(['  mean = ',num2str(mean(data_to_plot))]);
    disp(['  std = ',num2str(std(data_to_plot))]);
    
    if length(rows)<=2
        continue
    end
    
    % plot mean
    middle=mean(data_to_plot);
    %middle=median(data_to_plot);
    line(centers_to_plot(i)+[-1,1]*widths_to_plot(i)/2,[1,1]*middle,'color',[1,1,1]*0,'linewidth',barthickness);

    % plot errorbar
    line(centers_to_plot(i)*[1,1],middle+std(data_to_plot)./sqrt(length(data_to_plot))*[-1, 1],'color',[1,1,1]*0,'linewidth',barthickness);

    
    hold on;
end

% t-tests
disp(mean(scatter_colors(find(classifications==1)))-mean(scatter_colors(find(classifications==5))));
[h, p] = ttest2(scatter_colors(find(classifications==1)),scatter_colors(find(classifications==5))); % t test
disp(p);

disp(mean(scatter_colors(find(classifications==2)))-mean(scatter_colors(find(classifications==6))));
[h, p] = ttest2(scatter_colors(find(classifications==2)),scatter_colors(find(classifications==6))); % t test
disp(p);

disp(mean(scatter_colors(find(classifications==3)))-mean(scatter_colors(find(classifications==7))));
[h, p] = ttest2(scatter_colors(find(classifications==3)),scatter_colors(find(classifications==7))); % t test
disp(p);

disp(mean(scatter_colors(find(classifications==4)))-mean(scatter_colors(find(classifications==8))));
[h, p] = ttest2(scatter_colors(find(classifications==4)),scatter_colors(find(classifications==8))); % t test
disp(p);

box on;

xlim([0.5,max(cluster_data)/2+0.5]);
ylim(y_lims);
ax=gca;
ax.XTick=[];
ax.XColor = [0,0,0];
ax.YColor = [0,0,0];
ax.YTick=[];
 

set(gcf,...
'paperunits','inches',...
'papersize',[plotwidth,plotheight],...
'paperposition',[0,0,plotwidth,plotheight],...
'units','inches',...
'position',[0,0,plotwidth,plotheight],...
'inverthardcopy','off');
% set(gca, 'units', 'inches', 'position', [0, 0, plotwidth, plotheight]); % axis full
print('-dpdf', '-painters', [folder, '/', plot_filename]);
display('Done!');
close();


% get p-value info from pseudo-bulks (1,2,4h)
[h,p]=ttest2([mean(scatter_colors(find(cluster_data==1))),mean(scatter_colors(find(cluster_data==2))),mean(scatter_colors(find(cluster_data==3)))], [mean(scatter_colors(find(cluster_data==4))),mean(scatter_colors(find(cluster_data==5))),mean(scatter_colors(find(cluster_data==6)))]);
disp(p);

