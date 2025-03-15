%% plot pca for a specific batch (1 or 2) from deseq2 output
clear
clc

data=importdata('/Users/tanlongzhi/R/plate-c/rna_ci994_241007a_deseq2_matrix.csv');

expression_data=data.data;
gene_names=data.textdata(2:end,1);
sample_names=split(data.textdata(1,1),',');

% create metadata
num_groups = 10;
condition_data = contains(sample_names,'DIV1-Start')*1 + contains(sample_names,'DIV2-Vehicle')*2 +contains(sample_names,'DIV4-Vehicle')*3 +  contains(sample_names,'DIV2-CI-994')*4 +contains(sample_names,'DIV4-CI-994')*5;
batch_data = contains(sample_names,'Batch1')*1 + contains(sample_names,'Batch2')*2;
condition_batch_data = condition_data + (num_groups/2)*(batch_data - 1);

% plot parameters
plot_width = 2.4;plot_height = 1.3;dot_size = 10; padding_size=10; % plot 241012a pca

% plot a specific batch
plot_filename = 'pca.pdf';
rows = find(batch_data==1);
%rows = find(batch_data==2);

rows = rows(randperm(length(rows))); % randomly permutate to avoid overlap

% pca
[coeff, score, latent, tsquared, explained, mu] = pca(expression_data(:,rows)');
cols=[1,2];

% calculate average coordinates for each condition
mean_score = [];
for i = 1:(num_groups/2)
    mean_score(i,:) = mean(score(condition_data(rows)==i,:));
end

h = figure('color', 'w');
set(h, 'visible', 'off');

hold on;

% plot dots
scatter(score(:,cols(1)), score(:,cols(2)), dot_size, 'k', 'filled'); % PC 1 vs 2

% format
box on;
axis equal;

xlim([min(score(:,cols(1))), max(score(:,cols(1)))]+padding_size*[-1, 1]); ylim([min(score(:,cols(2))), max(score(:,cols(2)))]+padding_size*[-1, 1]);
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
print('-dpdf', '-painters', plot_filename);
close();
disp('printed');


%% plot a gene's trajectory 

clear
clc

axis_limits =5;
axis_limits =2;
gene_to_plot='Gabra2';
gene_to_plot='Syt5';
gene_to_plot='Bdnf';
gene_to_plot='Kif18a';
gene_to_plot='Sox11';
gene_to_plot='Ccdc184';
gene_to_plot='Ppp1r14a';
gene_to_plot='Serpine1';
gene_to_plot='Stmn4';
gene_to_plot='Fgfr3';
gene_to_plot='Fgf11';
gene_to_plot='Tubb3';
gene_to_plot='Stc2';
gene_to_plot='Selenop';
gene_to_plot='Rbm20';
gene_to_plot='Nes';
gene_to_plot='Nhlh1';
gene_to_plot='Nhlh2';
gene_to_plot='Fos';
gene_to_plot='Ccdc184';
gene_to_plot='Syt7';
gene_to_plot='Ash2l';
gene_to_plot='Rflnb';
gene_to_plot='Dusp15';
gene_to_plot='Stac2';
gene_to_plot='Doc2b';
gene_to_plot='Necab1';
gene_to_plot='Rflnb';
gene_to_plot='Sema3f';
gene_to_plot='Nes';
gene_to_plot='Fgfr3';
gene_to_plot='Syt13';
gene_to_plot='Sema3c';
gene_to_plot='Ash2l';
gene_to_plot='Gabra2';
gene_to_plot='Bdnf';
gene_to_plot='Gabra6';
% gene_to_plot='Atoh1';
% gene_to_plot='Snap25';
% gene_to_plot='Nhlh1';

data=importdata('/Users/tanlongzhi/R/plate-c/rna_ci994_241007a_deseq2_matrix.csv');

expression_data=data.data;
gene_names=data.textdata(2:end,1);
sample_names=split(data.textdata(1,1),',');

% create metadata
num_groups = 10;
condition_data = contains(sample_names,'DIV1-Start')*1 + contains(sample_names,'DIV2-Vehicle')*2 +contains(sample_names,'DIV4-Vehicle')*3 +  contains(sample_names,'DIV2-CI-994')*4 +contains(sample_names,'DIV4-CI-994')*5;
batch_data = contains(sample_names,'Batch1')*1 + contains(sample_names,'Batch2')*2;
condition_batch_data = condition_data + (num_groups/2)*(batch_data - 1);

% calculate each gene's trajectory, merging batches
mean_expression_data_by_group = zeros(length(gene_names),num_groups/2);
sem_expression_data_by_group = zeros(length(gene_names),num_groups/2);
for i=1:(num_groups/2)
    mean_expression_data_by_group(:,i) = mean(expression_data(:,condition_data==i),2);
    sem_expression_data_by_group(:,i) = std(expression_data(:,condition_data==i),0,2)./sqrt(size(expression_data(:,condition_data==i),2));
end

% plot parameters
plot_width = 0.5;plot_height = 1;dot_size = 7; errorbar_size=3;

% plot
plot_filename = [gene_to_plot,'.pdf'];

gene_to_plot_id = find(strcmp(gene_names,gene_to_plot));
baseline_data = mean_expression_data_by_group(gene_to_plot_id,1);

h = figure('color', 'w');
set(h, 'visible', 'off');

plot([1,4],[0,0],'k-');

hold on;

% plot lines
errorbar([1,2,4],mean_expression_data_by_group(gene_to_plot_id,[1,4,5])-baseline_data,sem_expression_data_by_group(gene_to_plot_id,[1,2,3]),'Color',[255,157,31]/255,'Marker','.','MarkerSize',dot_size,'CapSize',errorbar_size);
errorbar([1,2,4],mean_expression_data_by_group(gene_to_plot_id,[1,2,3])-baseline_data,sem_expression_data_by_group(gene_to_plot_id,[1,2,3]),'Color',[179,179,179]/255,'Marker','.','MarkerSize',dot_size,'CapSize',errorbar_size);

% format
box on;

xlim([1,4]);
ylim([-1,1]*axis_limits);
xticks([1,2,4]);
yticks([]);
xticklabels([]);
yticklabels([]);
title(gene_to_plot);

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
print('-dpdf', '-painters', plot_filename);
close();
disp('printed');


%% clustering genes with k-means

% center the matrix by subtracting the baseline
data_to_cluster = expression_data - mean(expression_data(:,condition_data==1),2);

% k-means clustering
num_clusters = 12;
disp('k means clustering starting');
[gene_cluster_data, cluster_centroids, sum_point_to_centroid_distances, point_to_centroid_distances] = kmeans(data_to_cluster, num_clusters,'MaxIter',1e3, 'Replicates',1e3);
disp('k means clustering done');
tabulate(gene_cluster_data);

% reorder clusters by average changes in CI-994
mean_cluster_centroids_by_group = zeros(num_clusters,num_groups/2);
for i=1:(num_groups/2)
    mean_cluster_centroids_by_group(:,i) = mean(cluster_centroids(:,condition_data==i),2);
end
mean_treatment_changes =  mean(mean_cluster_centroids_by_group(:,[4,5]),2) - mean(mean_cluster_centroids_by_group(:,[2,3]),2);
[C,outperm] = sort(mean_treatment_changes);

% write clustering results
outfile=fopen('kmeans.txt','w');
fprintf(outfile,'#gene\tcluster');
for j=1:num_clusters
    fprintf(outfile,['\tdistance to cluster ',num2str(j),' centroid']);
end
fprintf(outfile,'\n');
for j=1:num_clusters
    gene_ids = find(gene_cluster_data == outperm(j));
    for i=1:length(gene_ids)
        fprintf(outfile,[gene_names{gene_ids(i)},'\t',num2str(j)]);
        for k=1:num_clusters
            fprintf(outfile,['\t',num2str(point_to_centroid_distances(gene_ids(i),outperm(k)))]);
        end
        fprintf(outfile,'\n');
    end
end

fclose(outfile);
disp('written');

% create heatmap by reordering
rows = [];
cols = [];
for j=1:num_clusters
    gene_ids = find(gene_cluster_data == outperm(j));
    gene_ids = gene_ids(randperm(length(gene_ids)));
    rows = [rows; gene_ids];
end
for i=1:num_groups
    sample_ids = find(condition_batch_data == i);
    sample_ids = sample_ids(randperm(length(sample_ids)));
    cols = [cols; sample_ids];
end

% output heatmap
heatmap_data = data_to_cluster(rows, cols);
% c_min = -1; c_max = 1;
c_min = -2; c_max = 2;
x_rescale = 1;
y_rescale = 50;

% custom color map
low_color = [0,0,1];mid_color=[1,1,1];high_color = [1,0,0];
num_colors = 256;
low_to_high = [];
for i=1:3
    low_to_high = [low_to_high, [linspace(low_color(i), mid_color(i), num_colors/2), linspace(mid_color(i), high_color(i), num_colors/2)]'];
end

% plot heatmap
color_map = low_to_high;
final_value = heatmap_data-c_min;
value_range = c_max - c_min;
final_value = round(max(min(final_value / value_range, 1), 0) * (num_colors - 1) + 1);
final_value = imresize(final_value,[size(final_value, 1)*x_rescale,size(final_value, 2)*y_rescale], 'nearest');
imwrite(final_value, color_map, 'output.png');
disp('printed');


%% plot mean trajectory of genes in each k-means cluster

% calculate each cluster's trajectory, merge batches
mean_expression_data_by_group = zeros(num_clusters,num_groups/2);
for i=1:(num_groups/2)
    for j=1:num_clusters
        mean_expression_data_by_group(j,i) = mean(expression_data(gene_cluster_data == outperm(j),condition_data==i),'all');
    end
end
% plot parameters
plot_width = 0.5;plot_height = 0.5;line_width = 1; axis_limits=4; % plot 241012a pca

% plot
for j=1:num_clusters
plot_filename = ['kmeans_',num2str(j),'.pdf'];

baseline_data = mean_expression_data_by_group(j,1);

h = figure('color', 'w');
set(h, 'visible', 'off');

hold on;

% plot lines
plot([1,2,4],mean_expression_data_by_group(j,[1,4,5])-baseline_data,'Color',[255,157,31]/255,'LineWidth',line_width);
plot([1,2,4],mean_expression_data_by_group(j,[1,2,3])-baseline_data,'Color',[179,179,179]/255,'LineWidth',line_width);

plot([1,4],[0,0],'k-');

% format
box on;

xlim([1,4]);
ylim(mean(mean_expression_data_by_group(j,:)-baseline_data)+[-0.5,0.5]*axis_limits);
xticks([]);
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
print('-dpdf', '-painters', plot_filename);
close();
disp('printed');
end


%% plot mean trajectory of a list of genes

% read gene list
gene_list_to_plot=readcell('/Users/tanlongzhi/R/plate-c/rna_ci994_241007a_deseq2_log2fc_vehicle_div2_vs_div1.genes_down.txt');
% gene_list_to_plot=readcell('/Users/tanlongzhi/R/plate-c/rna_ci994_241007a_deseq2_log2fc_vehicle_div2_vs_div1.genes_up.txt');
% gene_list_to_plot=readcell('/Users/tanlongzhi/R/plate-c/rna_ci994_241007a_deseq2_log2fc_ci994_vs_vehicle_div2.genes_down.txt');
% gene_list_to_plot=readcell('/Users/tanlongzhi/R/plate-c/rna_ci994_241007a_deseq2_log2fc_ci994_vs_vehicle_div2.genes_up.txt');

gene_list_to_plot=readcell('/Users/tanlongzhi/R/plate-c/rna_ci994_241007a_deseq2_log2fc_vehicle_div4_vs_div1.genes_down.txt');
gene_list_to_plot=readcell('/Users/tanlongzhi/R/plate-c/rna_ci994_241007a_deseq2_log2fc_vehicle_div4_vs_div1.genes_up.txt');
gene_list_to_plot=readcell('/Users/tanlongzhi/R/plate-c/rna_ci994_241007a_deseq2_log2fc_ci994_vs_vehicle_div4.genes_down.txt');
gene_list_to_plot=readcell('/Users/tanlongzhi/R/plate-c/rna_ci994_241007a_deseq2_log2fc_ci994_vs_vehicle_div4.genes_up.txt');

% calculate the gene list's mean trajectory, merge batches
mean_expression_data_by_group = zeros(num_groups/2,1);
sem_expression_data_by_group = zeros(num_groups/2,1);
for i=1:(num_groups/2)
    mean_expression_data_by_group(i) = mean(expression_data(ismember(gene_names, gene_list_to_plot),condition_data==i),'all');
    sem_expression_data_by_group = zeros(length(gene_names),num_groups/2);
end
% plot parameters
plot_width = 0.5;plot_height = 0.5;line_width = 1; axis_limits=2; % plot 241012a pca
plot_width = 0.5;plot_height = 0.5;dot_size = 7; errorbar_size=3; axis_limits=2; % plot 241012a pca

% plot
plot_filename = 'trajectory.pdf';

baseline_data = mean_expression_data_by_group(1);

h = figure('color', 'w');
set(h, 'visible', 'off');

hold on;

% plot lines with sem
errorbar([1,2,4],mean_expression_data_by_group([1,4,5])-baseline_data,sem_expression_data_by_group([1,2,3]),'Color',[255,157,31]/255,'Marker','.','MarkerSize',dot_size,'CapSize',errorbar_size);
errorbar([1,2,4],mean_expression_data_by_group([1,2,3])-baseline_data,sem_expression_data_by_group([1,2,3]),'Color',[179,179,179]/255,'Marker','.','MarkerSize',dot_size,'CapSize',errorbar_size);

plot([1,4],[0,0],'k-');


% format
box on;

xlim([1,4]);
ylim(mean(mean_expression_data_by_group-baseline_data)+[-0.5,0.5]*axis_limits);
xticks([]);
yticks([]);
xticklabels([]);
yticklabels([]);

ax = gca;
ax.XColor = [0,0,0];
ax.YColor = [0,0,0];

% title(num2str(sum(gene_cluster_data == outperm(j))));

% print
set(gcf,...
'paperunits', 'inches',...
'papersize', [plot_width, plot_height],...
'paperposition', [0, 0, plot_width, plot_height],...
'units', 'inches',...
'position', [0, 0, plot_width, plot_height]);
print('-dpdf', '-painters', plot_filename);
close();
disp('printed');
