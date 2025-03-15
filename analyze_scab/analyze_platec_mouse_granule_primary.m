%% load scA/B matrix
% for Plate-C mouse cultured cerebellar granule cells primary experiment

clear
clc

folder = '/Users/tanlongzhi/Research/dip-c/brain2/summary/colors';

filename = 'granule_200k_241003a.cpg_b1m.color2s.txt';

% load treatment data
% treatment_data=readcell('/Users/tanlongzhi/Research/dip-c/brain2/aux_data/treatments_granule_200k_241003a.txt','delimiter','');
treatment_data=readcell('/Users/tanlongzhi/Research/dip-c/brain2/aux_data/treatments_granule_200k_241003a_dmso_numbered.txt','delimiter','');

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
pc_space = raw_score(:,1:20); % first 10 PCs for Plate-C experiment

% initialize colors for plotting
scatter_colors=0;

%% plot PCA

plot_filename = 'output.pdf';

% plot and dot sizes
plot_width = 2.8;plot_height = 2.8;dot_size = 6; % plot 240928a
% plot_width = 1.8;plot_height = 1.8;dot_size = 3; % plot 241007a ward clustering plot

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

title('Drug Type'); drug_type_data=readcell('/Users/tanlongzhi/Research/dip-c/brain2/aux_data/drug_types_granule_200k_241003a.txt','Delimiter',''); scatter_colors = 1*contains(drug_type_data,'control') + 2*contains(drug_type_data,'HDAC inhibitor') + 3*contains(drug_type_data,'HDM inhibitor') + 4*contains(drug_type_data,'HAT inhibitor') + 5*contains(drug_type_data,'HMT inhibitor') + 6*contains(drug_type_data,'DNMT inhibitor') + 7*contains(drug_type_data,'BET inhibitor'); c_min=0; c_max=7;
% title('Ward Clusters'); cell_type_data=readcell('/Users/tanlongzhi/Research/dip-c/brain2/aux_data/ward3.3_granule_drug_241007b_manual_merge.txt','Delimiter','\t'); cell_type_data=cell2mat(cell_type_data(:,2)); scatter_colors = cell_type_data; c_min=0; c_max=5;

% title('Compartment Strength'); scatter_colors = importdata('/Users/tanlongzhi/Research/plate-c/manuscript/trainee_data/compartment_strength_granule/AABBABBAvalues_1Mb_granule_extent10_50groups_calculation.well_only.tsv'); scatter_colors = log2(scatter_colors.data);  c_min=0; c_max=2;
% title('% Inter-chromosomal'); scatter_colors = importdata('/Users/tanlongzhi/Research/dip-c/brain2/aux_data/contact_stats_granule_200k_241003a.txt'); scatter_colors = 100-scatter_colors.data(:,4);  c_min=0; c_max=50;


% specific which cells to plot
rows = 1:size(color_data,2);
rows = randperm(size(score,1)); % randomly order to avoid overlaps

% plot which 2 PCs
cols = [1,2]; % PC 1 vs. 2
% cols = [1,3]; % PC 1 vs. 3
% cols = [1,4]; % PC 1 vs. 4

% main plotting command
% scatter(score(rows, cols(1)), score(rows, cols(2)), scatter_sizes(rows), scatter_colors(rows), 'filled'); % PC 1 vs 2

% optional: label drug names
% text(score(rows, cols(1)), score(rows, cols(2)), strcat(treatment_data(rows)), 'fontsize', 1, 'interpreter', 'none');

% alternative plotting command: label significant compounds as solid, others as circles
drug_names_as_dots = readcell('/Users/tanlongzhi/Research/dip-c/brain2/aux_data/chen-qin_granule_241003a_significant_drugs.txt','Delimiter','');
as_dots = [];
for i = 1:size(color_data, 2)
    as_dots(i) = max(ismember(drug_names_as_dots, treatment_data{i}));
end
rows_as_dots = find(as_dots);
rows = randperm(length(rows_as_dots)); % randomly order to avoid overlaps
scatter(score(rows_as_dots(rows), cols(1)), score(rows_as_dots(rows), cols(2)), scatter_sizes(rows_as_dots(rows)), scatter_colors(rows_as_dots(rows)), 'filled');
rows_as_circles = find(~as_dots);
rows = randperm(length(rows_as_circles)); % randomly order to avoid overlaps
scatter(score(rows_as_circles(rows), cols(1)), score(rows_as_circles(rows), cols(2)), scatter_sizes(rows_as_circles(rows))*0.8, scatter_colors(rows_as_circles(rows)),'linewidth',0.2);

% specify color ranges
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
cluster_colors=[0,0,0; 179,179,179; 255,157,31; 47,212,0; 0,178,244; 19,0,233; 225,116,255; 162,0,180]./255; % drug type 240429a: 0=other, 1=control, 2=HDACi, 3=HDMi, 4=HATi, 5=HMTi, 6=DNMTi, 7=BETi
% cluster_colors=[0.85,0.85,0.85; 1,0,0; 0.8,0,0; 0.6,0,0; 0.9,0.75,0; 1,0.1,1]; % ward clustering 241007a


% choose a color map
colormap(cluster_colors);
% colormap(low_to_high);


% format
box on;
axis equal;

xlim([min(score(rows,cols(1))), max(score(rows,cols(1)))]+0.5*[-1, 1]); ylim([min(score(rows,cols(2))), max(score(rows,cols(2)))]+0.5*[-1, 1]); % fixed padding size
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

distance_cutoff=3.3; % granule clusters 241007b


linkage_data = linkage(pc_space, 'ward', 'euclidean');
cluster_data = cluster(linkage_data,'CutOff',distance_cutoff,'Criterion','distance');

temp_cluster_data=cluster_data;
temp_linkage_data=linkage_data;
temp_cluster_ids=1:size(cell_names,2);

tmp = temp_cluster_data;

% granule drug 241007a
temp_cluster_data=temp_cluster_data.*0;
temp_cluster_data(tmp==9|tmp==23)=1; % HDAC strong
temp_cluster_data(tmp==14)=2; % HDAC medium
temp_cluster_data(tmp==7|tmp==21)=3; % HDAC weak
temp_cluster_data(tmp==13|tmp==16|tmp==2)=4;
temp_cluster_data(tmp==1)=5;

cluster_data=temp_cluster_data;

% write to file
outfile=fopen('output.txt','w');

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

%% plot clustering tree


% optimized leaf ordering, can take ~0.5 hours
outperm=optimalleaforder(temp_linkage_data,pdist(pc_space),'Criteria','group');
disp('leaf order optimized');


plot_width=3; plot_height=3; plot_filename='tree.pdf';

figure('color','w');

% plot the tree
H=dendrogram(temp_linkage_data,0,'Reorder',outperm);


%set(H,'LineWidth',2);
set(H,'Color',[0,0,0]);
axis off;

% print
set(gcf,...
'paperunits', 'inches',...
'papersize', [plot_width, plot_height],...
'paperposition', [0, 0, plot_width, plot_height],...
'units', 'inches',...
'position', [0, 0, plot_width, plot_height]);
set(gca,'position', [0, 0, 1, 1]);
print('-dpdf', '-painters', [folder, '/', plot_filename]);
close();
disp('printed');


%% find differential scA/B regions (can have missing values)

clc

% significance thresholds
min_diff = 0.02; max_fdr = 0.01; % for scA/B

% initial axis limits
x_lims = [-1,1]*0.12; y_lims = [0,20]; % time course
x_lims = [-1,1]*0.18; y_lims = [0,30]; % primary

plot_width = 0.7; plot_height = 0.7; dot_size=3;
min_num_1 = 1;
min_num_2 = 1;

% load Ward clusters
cell_type_data=readcell('/Users/tanlongzhi/Research/dip-c/brain2/aux_data/ward3.3_granule_drug_241007b_manual_merge.txt','Delimiter','\t'); cell_type_data=cell2mat(cell_type_data(:,2));

% specify cells to compare
% cell_ids_1 is sample, cell_ids_2 is control

% for ward
output_prefix='granule_drug_cluster_1_vs_control_241009a';cell_ids_1=find(cell_type_data==1&~strcmp(treatment_data,'DMSO'));cell_ids_2 = find(strcmp(treatment_data,'DMSO'));
output_prefix='granule_drug_cluster_2_vs_control_241009a';cell_ids_1=find(cell_type_data==2&~strcmp(treatment_data,'DMSO'));cell_ids_2 = find(strcmp(treatment_data,'DMSO'));
output_prefix='granule_drug_cluster_3_vs_control_241009a';cell_ids_1=find(cell_type_data==3&~strcmp(treatment_data,'DMSO'));cell_ids_2 = find(strcmp(treatment_data,'DMSO'));
output_prefix='granule_drug_cluster_4_vs_control_241009a';cell_ids_1=find(cell_type_data==4&~strcmp(treatment_data,'DMSO'));cell_ids_2 = find(strcmp(treatment_data,'DMSO'));
output_prefix='granule_drug_cluster_5_vs_control_241009a';cell_ids_1=find(cell_type_data==5&~strcmp(treatment_data,'DMSO'));cell_ids_2 = find(strcmp(treatment_data,'DMSO'));


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

% division line

% plot([0,0],y_lims,'color',[0,0,0],'linewidth',0.5);
disp(['mean: ',num2str(mean(mean_colors_diff(red_rows)))]);

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
bin_size =1e6;

outfile=fopen([output_prefix,'.bed'],'w');
outfile2_up=fopen([output_prefix,'_up.b1m.txt'],'w');
outfile2_down=fopen([output_prefix,'_down.b1m.txt'],'w');
outfile3=fopen([output_prefix,'.diff.txt'],'w');
outfile4=fopen([output_prefix,'.all_b1m_diff.txt'],'w');
outfile5=fopen([output_prefix,'.supp_table.txt'],'w');
outfile6=fopen([output_prefix,'.signed_significance.bedgraph'],'w');

% optional: find genes
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


%% find treatments that significantly changed genome-wide scA/B profiles by Chen-Qin tests


% significance threshold
chen_qin_test_stat_threshold = 9.7979;


data_matrix = color_data';

control_drug_name = 'DMSO';
control_drug_name_prefix = 'DMSO';

unique_drug_names = unique(treatment_data); % for treating numbered DMSO as "null drugs"
% unique_drug_names = unique_drug_names(find(~strcmp(unique_drug_names, control_drug_name)));

% optional: annotate compound type in the output file
drug_type_data=readcell('/Users/tanlongzhi/Research/dip-c/brain2/aux_data/drug_types_granule_200k_241003a.txt','Delimiter',''); 

chen_qin_test_stats=[];
mean_distances=[];

% write to file
outfile=fopen('chen-qin.txt','w');
fprintf(outfile, '#drug\tdrug type\tsignificance\tnumber of replicates\tChen-Qin test statistic (T/sigma)\tChen-Qin T\tChen-Qin sigma\t||delta scA/B||\n');

for i = 1:length(unique_drug_names) 
    unique_drug_name = unique_drug_names{i};
    disp(unique_drug_name);
    drug_cell_ids = find(strcmp(treatment_data,unique_drug_name));
    control_cell_ids = find(startsWith(treatment_data,control_drug_name_prefix));

    control_cell_ids = setdiff(control_cell_ids, drug_cell_ids); % remove self for DMSO as "null drugs"

    % find drug type
    drug_type = drug_type_data{drug_cell_ids(1)};


    if (length(control_cell_ids))==0
        continue
    end
    
    drug_cell_ids = drug_cell_ids(randperm(length(drug_cell_ids)));
    control_cell_ids = control_cell_ids(randperm(length(control_cell_ids)));

    % calculate statistics
    x_1 = data_matrix(drug_cell_ids,:);
    x_2 = data_matrix(control_cell_ids,:);
    n_1 = length(drug_cell_ids);
    n_2 = length(control_cell_ids);
    
    mean_distance = norm(mean(x_1,1) - mean(x_2,1));
    
    t_value = mean(squareform(tril(x_1*x_1',-1)),'all') + mean(squareform(tril(x_2*x_2',-1)),'all') - 2 * mean(x_1*x_2','all');
    
    sigma_square_1 = 0;
    for j = 1:n_1
        for k = (j+1):n_1
            mean_excluding_j_k = mean(x_1(setdiff(1:n_1,[j,k]),:),1);
            sigma_square_1 = sigma_square_1 + x_1(j,:)*(x_1(k,:)-mean_excluding_j_k)'*x_1(k,:)*(x_1(j,:)-mean_excluding_j_k)';
        end
    end
    sigma_square_1 = sigma_square_1 * 2/n_1/(n_1-1);

    sigma_square_2 = 0;
    for j = 1:n_2
        for k = (j+1):n_2
            mean_excluding_j_k = mean(x_2(setdiff(1:n_2,[j,k]),:),1);
            sigma_square_2 = sigma_square_2 + x_2(j,:)*(x_2(k,:)-mean_excluding_j_k)'*x_2(k,:)*(x_2(j,:)-mean_excluding_j_k)';
        end
    end
    sigma_square_2 = sigma_square_2 * 2/n_2/(n_2-1);

    sigma_square_1_2 = 0;
    for j = 1:n_1
        for k = 1:n_2
            mean_excluding_j = mean(x_1(setdiff(1:n_1,j),:),1);
            mean_excluding_k = mean(x_2(setdiff(1:n_2,k),:),1);
            sigma_square_1_2 = sigma_square_1_2 + x_1(j,:)*(x_2(k,:)-mean_excluding_k)'*x_2(k,:)*(x_1(j,:)-mean_excluding_j)';
        end
    end
    sigma_square_1_2 = sigma_square_1_2 /n_1/n_2;

    sigma_square = 2/n_1/(n_1-1)*sigma_square_1 + 2/n_2/(n_2-1)*sigma_square_2 + 4/n_1/n_2*sigma_square_1_2;
    
    chen_qin_test_stat = t_value./sqrt(sigma_square);
    

    chen_qin_test_stats(i) = chen_qin_test_stat;
    mean_distances(i) = mean_distance;

    is_significant = 'no';
    if chen_qin_test_stat>=chen_qin_test_stat_threshold
        is_significant = 'yes';
    end

    % write output
    fprintf(outfile, [unique_drug_name, '\t', drug_type,'\t',is_significant,'\t',num2str(n_1),'\t',num2str(chen_qin_test_stat),'\t',num2str(t_value),'\t',num2str(sqrt(sigma_square)),'\t',num2str(mean_distance),'\n']);

end

fclose(outfile);
disp('written');


%% plot a histogram of the Chen-Qin test statistic

plot_filename='chen-qin.pdf';

plot_width = 4;plot_height = 0.6;dot_size = 6; % plot 241005a


% plot
h = figure('color', 'w');
set(h, 'visible', 'off');

histogram(chen_qin_test_stats(startsWith(unique_drug_names,control_drug_name_prefix)),-20:2:120,'facealpha',1,'facecolor', [179,179,179]/255,'edgecolor','none');
histogram(chen_qin_test_stats(~startsWith(unique_drug_names,control_drug_name_prefix)),-20:2:120,'facealpha',1,'facecolor', [128,128,128]/255 ,'edgecolor','none');

% format


box off;


xlim([-20,120]);
ylim([0,13]); % control
ylim([0,28]); % drugs
xticks(-200:20:200);
yticks([]);
xticklabels([]);
yticklabels([]);



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
