%% plot delta scA/B comparison between 2 sets of differential scA/B analysis

clear
clc

% set x and y axis limits
% axis_limits = 0.2;
axis_limits = 0.15; % 250219a, Fig 4
% axis_limits = 0.13; % 250220a
% axis_limits = 0.21; % 250220b

% load 2 sets of delta scA/B profiles
raw_gene_data_1= load('tsa_time_course_granule_1_2_4h.all_b1m_diff.txt','-ascii'); % in vivo HDACi
raw_gene_data_2= load('fast_tsa_1_2_4h.all_b1m_diff.txt','-ascii'); % in vitro HDACi
% raw_gene_data_2 = load('granule_5_vs_4.all_b1m_diff.txt','-ascii');
% raw_gene_data_2 = load('granule_4_vs_3.all_b1m_diff.txt','-ascii');
% raw_gene_data_2 = load('granule_3_vs_2.all_b1m_diff.txt','-ascii');
% raw_gene_data_2 = load('granule_2_vs_1.all_b1m_diff.txt','-ascii');
% raw_gene_data_1 = load('granule_5_vs_3.all_b1m_diff.txt','-ascii');
% raw_gene_data_1 = load('granule_4_vs_2.all_b1m_diff.txt','-ascii');
% raw_gene_data_1 = load('granule_3_vs_1.all_b1m_diff.txt','-ascii');
% raw_gene_data_1 = load('granule_5_vs_2.all_b1m_diff.txt','-ascii');
% raw_gene_data_1 = load('granule_4_vs_1.all_b1m_diff.txt','-ascii');
% raw_gene_data_1 = load('granule_5_vs_1.all_b1m_diff.txt','-ascii');
% raw_gene_data_2 = load('tsa_injections_granule.all_b1m_diff.txt','-ascii');
% raw_gene_data_1 = load('tsa_primary_granule.all_b1m_diff.txt','-ascii');
% raw_gene_data_2 = load('tsa_time_course_other_neuron_1_2_4h.all_b1m_diff.txt','-ascii');
% raw_gene_data_1 = load('tsa_time_course_granule_1h.all_b1m_diff.txt','-ascii');
% raw_gene_data_2 = load('tsa_time_course_granule_2h.all_b1m_diff.txt','-ascii');
 

% find common genomic loci (rows) between the 2 files
[C, locus_indices_1, locus_indices_2] = intersect([raw_gene_data_1(:,1), raw_gene_data_1(:,2)-500e3, raw_gene_data_1(:,2)+500e3], [raw_gene_data_2(:,1), raw_gene_data_2(:,2)-500e3, raw_gene_data_2(:,2)+500e3], 'rows');
gene_data_1 = raw_gene_data_1(locus_indices_1,3);
gene_data_2 = raw_gene_data_2(locus_indices_2,3);

% calculate correlation
disp(['Pearson correlation coefficient: ',num2str(corr(gene_data_1, gene_data_2))]);

% calculate slope
[coeff, score, latent, tsquared, explained, mu] = pca([gene_data_1, gene_data_2]);
disp(['Total least squares linear regression slope: ',num2str(coeff(2,1)./coeff(1,1))]);

% plot parameters
plot_width = 0.6;plot_height = 0.6;dot_size = 0.5; % plot 241014a FC vs FC
plot_width = 0.8;plot_height = 0.8;dot_size = 0.5; % plot 250219a FC vs FC

% plot
plot_filename = 'versus.png';

h = figure('color', 'w');
%set(h, 'visible', 'off');

hold on;

rows=randperm(length(gene_data_1));

% plot dots
scatter(gene_data_1(rows), gene_data_2(rows), dot_size,[1,1,1]*0.5, 'filled');

% format
box on;
axis equal;
axis off;


xlim([-1,1]*axis_limits);
ylim([-1,1]*axis_limits);
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
set(gca,...
'position', [0, 0, 1, 1]);
print('-dpng', '-image','-r600', plot_filename);
% exportgraphics(h,plot_filename);
close();
disp('printed');


%% plot correlation matrix heatmap between many conditions

clear
clc

data=load('/Users/tanlongzhi/Research/plate-c/manuscript/tables/correlation_250221a.tsv','-ascii');

plot_width = 0.8;plot_height = 0.8;color_limit=0.8; % plot 250219a FC vs FC

% plot
plot_filename = 'correlation.pdf';

% plot
h = figure('color', 'w');
%set(h, 'visible', 'off');

hold on;

% plot
imagesc(data);

% red white blue
num_colors=256;
num_colors_one_side=(num_colors+1)/2;
t=linspace(0,1,num_colors_one_side)';
redwhiteblue=[[1,0,0].*(1-t)+[1,1,1].*t;[1,1,1].*(1-t(2:end))+[0,0,1].*t(2:end),];

% for fine cell types
% cluster_colors = parula(max(scatter_colors));
cluster_colors = flipud(redwhiteblue);

% set color map
colormap(cluster_colors);

% format
box on;
axis off;
axis equal;
caxis([-1,1]*color_limit);

xlim([0.5,size(data,1)+0.5]); ylim([0.5,size(data,2)+0.5]);
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
set(gca,'position', [0, 0, 1, 1]);
print('-dpdf','-painters', plot_filename);
close();
disp('printed');
