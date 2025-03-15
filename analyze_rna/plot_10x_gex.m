%% plot UMAP from seurat data

clear
clc

% load umap coordinates
score = load('/Users/tanlongzhi/R/plate-c/umap_unintegrated.txt', '-ascii');
score = load('/Users/tanlongzhi/R/plate-c/umap.txt', '-ascii');
% score = load('/Users/tanlongzhi/R/plate-c/umap_harmony.txt', '-ascii');

% load cell colors
% cell types
scatter_colors = load('/Users/tanlongzhi/R/plate-c/cell_types.txt', '-ascii'); c_min = 1; c_max = 20;
% scatter_colors = load('/Users/tanlongzhi/R/plate-c/cell_types_split.txt', '-ascii'); c_min = 14; c_max = 17; % 15 is progentior, 16 is mature

% treatments
% scatter_colors = load('/Users/tanlongzhi/R/plate-c/orig.ident.txt', '-ascii'); c_min = 1; c_max = 4;

% expression of a specific gene
% scatter_colors = load('/Users/tanlongzhi/R/plate-c/expression_Gabra6.txt', '-ascii'); c_min = 0; c_max = 3;
% scatter_colors = load('/Users/tanlongzhi/R/plate-c/expression_Foxp2.txt', '-ascii'); c_min = 0; c_max = 3;
% scatter_colors = load('/Users/tanlongzhi/R/plate-c/expression_Rasgef1b.txt', '-ascii'); c_min = 0; c_max = 3;
% scatter_colors = load('/Users/tanlongzhi/R/plate-c/expression_Cntnap2.txt', '-ascii'); c_min = 0; c_max = 3;
% scatter_colors = load('/Users/tanlongzhi/R/plate-c/expression_Scn2a.txt', '-ascii'); c_min = 0; c_max = 3;
% scatter_colors = load('/Users/tanlongzhi/R/plate-c/expression_Trp73.txt', '-ascii'); c_min = 0; c_max = 3;
% scatter_colors = load('/Users/tanlongzhi/R/plate-c/expression_Cux2.txt', '-ascii'); c_min = 0; c_max = 3;
% scatter_colors = load('/Users/tanlongzhi/R/plate-c/expression_Rorb.txt', '-ascii'); c_min = 0; c_max = 3;
% scatter_colors = load('/Users/tanlongzhi/R/plate-c/expression_Tle4.txt', '-ascii'); c_min = 0; c_max = 3;
% scatter_colors = load('/Users/tanlongzhi/R/plate-c/expression_Nr4a2.txt', '-ascii'); c_min = 0; c_max = 3;
% scatter_colors = load('/Users/tanlongzhi/R/plate-c/expression_Satb2.txt', '-ascii'); c_min = 0; c_max = 3;

% plot parameters
plot_filename = 'output.png';
plot_width = 2; dot_size = 1; padding_size=0.5; % 250214a umap big
plot_height = 2; dot_size = 1; padding_size=0.5; % 250214b umap big
plot_height = 2; dot_size = 0.15; padding_size=0.5; % 250305a umap for cerebellar granule cells
plot_height = 1; dot_size = 0.25; padding_size=0.5; % 250305a small umap for other markers

% plot
h = figure('color', 'w');
set(h, 'visible', 'off');

hold on;

% cells to plot
rows = randperm(size(score,1)); % randomly order to avoid overlaps
% plot_height = plot_width./(max(score(rows,1))-min(score(rows,1))+2.*padding_size).*(max(score(rows,2))-min(score(rows,2))+2.*padding_size);
plot_width = plot_height.*(max(score(rows,1))-min(score(rows,1))+2.*padding_size)./(max(score(rows,2))-min(score(rows,2))+2.*padding_size);

% plot
scatter(score(rows, 1), score(rows, 2), dot_size, scatter_colors(rows), 'filled');

% custom gene expression colors
low_color = [1,1,1]*0.8;high_color = [1,0,0];
num_colors = 256;
cluster_colors = [];
for i=1:3
    cluster_colors = [cluster_colors, linspace(low_color(i), high_color(i), num_colors)'];
end

% color for cell types
cluster_colors = [...
    175,119,255; ... % cerebellar granule cell/progenitor
    140,196,255; ... % olfactory ensheathing cell
    254,82,0; ... % microglia
    204,0,5; ... % endothelial cell
    227,143,131; ... % choroid plexus epithelial cell
    255,191,197; ... % vascular smooth muscle cell/pericyte
    255,148,128; ... % vascular and leptomeningeal cell
    209,149,151; ... % arachnoid barrier cell
    22,112,64; ... % ependymal cell
    138,231,255; ... % oligodendrocyte progenitor
    19,125,255; ... % oligodendrocyte
    71,191,141; ... % astrocyte
    214,98,162; ... % purkinje cell
    255,132,0; ... % medium spiny neuron
    255,216,125; ... % cerebellar interneuron
    109,107,255; ... % hippocampal granule cell
    175,125,75; ... % other excitatory neuron
    214,200,81; ... % other inhibitory neuron
    179,172,113; ... % unknown neuron
    186,186,186; ... % unknown
    ]/255;

% color for treatments
% cluster_colors = [199,199,199;179,179,179;255,157,31;239,134,0]/255;

% color for granule cell progenitor vs mature
% cluster_colors = [186,186,186;216,166,255;140,77,254;186,186,186]/255;

% set color map
colormap(cluster_colors);

% format
box on;
axis off;
axis equal;
caxis([c_min,c_max]);

xlim([min(score(rows,1)), max(score(rows,1))]+padding_size*[-1, 1]); ylim([min(score(rows,2)), max(score(rows,2))]+padding_size*[-1, 1]);
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
print('-dpng', '-image','-r600', plot_filename);
close();
disp('printed');

%% plot psuedobulk comparison

clear
clc

% load data to plot
data_all_genes = load('/Users/tanlongzhi/R/plate-c/pseudobulk_all.txt', '-ascii');
data_degs = load('/Users/tanlongzhi/R/plate-c/pseudobulk_all_degs.txt', '-ascii');

% load gene names
metadata_degs = [readtable('/Users/tanlongzhi/R/plate-c/degs_up.all.csv'); readtable('/Users/tanlongzhi/R/plate-c/degs_down.all.csv')];

disp(max(data_all_genes,[],'all'));

% choose the 2 samples to plot against each other
sample_1 = 1; sample_2 = 3;
% sample_1 = 2; sample_2 = 4;

plot_filename = 'output.png';
plot_width = 0.8; dot_size = 20; axis_limit = 5; % 250215a TSA vs DMSO
plot_height = plot_width;

% plot
h = figure('color', 'w');
set(h, 'visible', 'off');

hold on;

% plot
scatter(data_all_genes(:,sample_1), data_all_genes(:,sample_2), dot_size*0.25, [1,1,1]*0.5, 'filled');
scatter(data_degs(:,sample_1), data_degs(:,sample_2), dot_size, 'r', 'filled');

% label gene names
text(data_degs(:,sample_1), data_degs(:,sample_2), metadata_degs.Var1); 


% format
box on;
axis off;
axis equal;

xlim([0,axis_limit]); ylim([0,axis_limit]);
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
print('-dpng', '-image','-r600', plot_filename);
close();
disp('printed');


%% plot log2fc comparison between 2 cell types
% note: x axis is cell type 2, y axis is cell type 1

clear
clc

% load data to plot

cell_type_1='cerebellar_granule_cell_and_progenitor';
cell_type_2='medium_spiny_neuron';
% cell_type_2='purkinje_cell';
% cell_type_2='cerebellar_interneuron';
% cell_type_2='differentiation';


data_all_genes = readtable(['/Users/tanlongzhi/R/plate-c/degs.between_cell_types.',cell_type_1,'.vs.',cell_type_2,'.csv']);
degs_1 = [readtable(['/Users/tanlongzhi/R/plate-c/degs_up.',cell_type_1,'.csv']); readtable(['/Users/tanlongzhi/R/plate-c/degs_down.',cell_type_1,'.csv'])];
degs_2 = [readtable(['/Users/tanlongzhi/R/plate-c/degs_up.',cell_type_2,'.csv']); readtable(['/Users/tanlongzhi/R/plate-c/degs_down.',cell_type_2,'.csv'])];

deg_names_both = intersect(degs_1.gene, degs_2.gene);
deg_names_one = setdiff(union(degs_1.gene, degs_2.gene), deg_names_both);

deg_ids_both = find(ismember(data_all_genes.gene, deg_names_both));
deg_ids_one = find(ismember(data_all_genes.gene, deg_names_one));

% find
disp(max(data_all_genes{:,[4,9]},[],'all'));
disp(min(data_all_genes{:,[4,9]},[],'all'));

% load gene names
plot_filename = 'output.png';
plot_width = 1.5; plot_height = 1.5; dot_size = 2; axis_limit_x = 6; axis_limit_y = 3;

% set dot sizes and colors
scatter_colors = zeros(height(data_all_genes),1);
scatter_colors(deg_ids_both)=1;
scatter_colors(deg_ids_one)=2;


% plot
h = figure('color', 'w');
set(h, 'visible', 'off');

hold on;

% plot
rows = randperm(height(data_all_genes)); % randomly order to avoid overlaps

scatter(data_all_genes{rows,9}, data_all_genes{rows,4}, dot_size, scatter_colors(rows), 'filled','o'); % x=cell_type_2, y=cell_type_1

% colormap([0.5,0.5,0.5;1,0,0;21/255,152/255,255/255]);
colormap([0.5,0.5,0.5;1,0,0;0,0,1]);

% label gene names
text(data_all_genes{deg_ids_both,9}, data_all_genes{deg_ids_both,4}, data_all_genes.gene(deg_ids_both),'color','r'); 
text(data_all_genes{deg_ids_one,9}, data_all_genes{deg_ids_one,4}, data_all_genes.gene(deg_ids_one),'color','b'); 


% format
box on;
axis off;
% axis equal;

xlim([-1,1]*axis_limit_x); ylim([-1,1]*axis_limit_y);caxis([0,3]);
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
print('-dpng', '-image','-r600', plot_filename);
% close();
disp('printed');


%% plot mean expression level (and SEM) of a specific gene across reactions and cell types
% left to right: granule, MSN, purkinje, interneuron

clear
clc

% format 
plotwidth=1.3;plotheight=0.8; spread_width = 0.8; barthickness=0.5;

% cell types
cell_type_data = load('/Users/tanlongzhi/R/plate-c/cell_types.txt', '-ascii'); % 1 to 20

% treatment
treatment_data = load('/Users/tanlongzhi/R/plate-c/orig.ident.txt', '-ascii'); % 1 to 4

% gene expression
scatter_colors = load('/Users/tanlongzhi/R/plate-c/expression_Gabra6.txt', '-ascii'); y_lims=[0,5];
% scatter_colors = load('/Users/tanlongzhi/R/plate-c/expression_Pde7b.txt', '-ascii'); y_lims=[0,5];
% scatter_colors = load('/Users/tanlongzhi/R/plate-c/expression_Slc1a6.txt', '-ascii'); y_lims=[0,5];
% scatter_colors = load('/Users/tanlongzhi/R/plate-c/expression_Gad2.txt', '-ascii'); y_lims=[0,5];
% scatter_colors = load('/Users/tanlongzhi/R/plate-c/expression_Itpr1.txt', '-ascii'); y_lims=[0,5];
% scatter_colors = load('/Users/tanlongzhi/R/plate-c/expression_Dscaml1.txt', '-ascii'); y_lims=[0,5];
% scatter_colors = load('/Users/tanlongzhi/R/plate-c/expression_Prkcg.txt', '-ascii'); y_lims=[0,5];
% scatter_colors = load('/Users/tanlongzhi/R/plate-c/expression_Inpp5a.txt', '-ascii'); y_lims=[0,5];


classifications = cell_type_data.*4 + treatment_data;

clusters_to_plot = [1*4+(1:4), 14*4+(1:4), 13*4+(1:4), 15*4+(1:4)];
widths_to_plot=spread_width*ones([max(clusters_to_plot),1]);
colors_to_plot = [199,199,199;179,179,179;255,157,31;239,134,0]/255;
colors_to_plot = repmat(colors_to_plot,[4,1]);


h=figure('color','none');
set(h,'Visible','off');

for i=1:size(clusters_to_plot,2)
    cluster_to_plot=clusters_to_plot(i);
    plot_color=[0,0,0];
    disp(cluster_to_plot);
    rows=find(classifications==cluster_to_plot);
    data_to_plot=scatter_colors(rows)';
    
    if length(rows)<=2
        continue
    end
    
    % plot mean
    line(i+[-1,1]*widths_to_plot(i)/2,[1,1]*mean(data_to_plot,'all'),'color',colors_to_plot(i,:),'linewidth',barthickness*2);

    % plot errorbar
    line(i*[1,1],mean(data_to_plot,'all')+std(data_to_plot)./sqrt(length(data_to_plot))*[-1, 1],'color',colors_to_plot(i,:),'linewidth',barthickness);

    
    hold on;
end

% add division lines
line([1,1]*4.5,y_lims,'color',[0,0,0],'linewidth',barthickness);
line([1,1]*8.5,y_lims,'color',[0,0,0],'linewidth',barthickness);
line([1,1]*12.5,y_lims,'color',[0,0,0],'linewidth',barthickness);


box on;

xlim([0.5,length(clusters_to_plot)+0.5]);
ylim(y_lims);
ax=gca;
ax.XTick=[];
ax.XColor = [0,0,0];
ax.YColor = [0,0,0];
ax.YTick=[0:5];
ax.YTickLabels=[];
 

set(gcf,...
'paperunits','inches',...
'papersize',[plotwidth,plotheight],...
'paperposition',[0,0,plotwidth,plotheight],...
'units','inches',...
'position',[0,0,plotwidth,plotheight],...
'inverthardcopy','off');
% set(gca, 'units', 'inches', 'position', [0, 0, plotwidth, plotheight]); % axis full
print('-dpdf', '-painters', 'types.pdf');
display('Done!');
close();

%% plot mean expression level (and SEM) of a specific gene across reactions and in more mature vs. less mature cerebeller granule cells

clear
clc

% format 
plotwidth=0.7;plotheight=0.8; spread_width = 0.8; barthickness=0.5;

% cell types
cell_type_data = load('/Users/tanlongzhi/R/plate-c/cell_types_split.txt', '-ascii'); % 1 to 20

% treatment
treatment_data = load('/Users/tanlongzhi/R/plate-c/orig.ident.txt', '-ascii'); % 1 to 4

% gene expression
scatter_colors = load('/Users/tanlongzhi/R/plate-c/expression_Gabra6.txt', '-ascii'); y_lims=[0,5];
scatter_colors = load('/Users/tanlongzhi/R/plate-c/expression_Scn1a.txt', '-ascii'); y_lims=[0,5];
scatter_colors = load('/Users/tanlongzhi/R/plate-c/expression_Slit2.txt', '-ascii'); y_lims=[0,5];
scatter_colors = load('/Users/tanlongzhi/R/plate-c/expression_Gpc3.txt', '-ascii'); y_lims=[0,5];
scatter_colors = load('/Users/tanlongzhi/R/plate-c/expression_Setbp1.txt', '-ascii'); y_lims=[0,5];
scatter_colors = load('/Users/tanlongzhi/R/plate-c/expression_Cog7.txt', '-ascii'); y_lims=[0,5];
scatter_colors = load('/Users/tanlongzhi/R/plate-c/expression_Dpf3.txt', '-ascii'); y_lims=[0,5];
scatter_colors = load('/Users/tanlongzhi/R/plate-c/expression_Tmem108.txt', '-ascii'); y_lims=[0,5];
scatter_colors = load('/Users/tanlongzhi/R/plate-c/expression_Fgfr1.txt', '-ascii'); y_lims=[0,5];
scatter_colors = load('/Users/tanlongzhi/R/plate-c/expression_Htr2c.txt', '-ascii'); y_lims=[0,5];
scatter_colors = load('/Users/tanlongzhi/R/plate-c/expression_Prdm10.txt', '-ascii'); y_lims=[0,5];
scatter_colors = load('/Users/tanlongzhi/R/plate-c/expression_Crppa.txt', '-ascii'); y_lims=[0,5];

classifications = cell_type_data.*4 + treatment_data;

clusters_to_plot = [15*4+(1:4), 16*4+(1:4)];
widths_to_plot=spread_width*ones([max(clusters_to_plot),1]);
colors_to_plot = [199,199,199;179,179,179;255,157,31;239,134,0]/255;
colors_to_plot = repmat(colors_to_plot,[2,1]);


h=figure('color','none');
set(h,'Visible','off');

for i=1:size(clusters_to_plot,2)
    cluster_to_plot=clusters_to_plot(i);
    plot_color=[0,0,0];
    disp(cluster_to_plot);
    rows=find(classifications==cluster_to_plot);
    data_to_plot=scatter_colors(rows)';
    
    if length(rows)<=2
        continue
    end
    
    % plot mean
    line(i+[-1,1]*widths_to_plot(i)/2,[1,1]*mean(data_to_plot,'all'),'color',colors_to_plot(i,:),'linewidth',barthickness*2);

    % plot errorbar
    line(i*[1,1],mean(data_to_plot,'all')+std(data_to_plot)./sqrt(length(data_to_plot))*[-1, 1],'color',colors_to_plot(i,:),'linewidth',barthickness);

    
    hold on;
end

% add division lines
line([1,1]*4.5,y_lims,'color',[0,0,0],'linewidth',barthickness);
line([1,1]*8.5,y_lims,'color',[0,0,0],'linewidth',barthickness);
line([1,1]*12.5,y_lims,'color',[0,0,0],'linewidth',barthickness);


box on;

xlim([0.5,length(clusters_to_plot)+0.5]);
ylim(y_lims);
ax=gca;
ax.XTick=[];
ax.XColor = [0,0,0];
ax.YColor = [0,0,0];
ax.YTick=[0:5];
ax.YTickLabels=[];
 

set(gcf,...
'paperunits','inches',...
'papersize',[plotwidth,plotheight],...
'paperposition',[0,0,plotwidth,plotheight],...
'units','inches',...
'position',[0,0,plotwidth,plotheight],...
'inverthardcopy','off');
% set(gca, 'units', 'inches', 'position', [0, 0, plotwidth, plotheight]); % axis full
print('-dpdf', '-painters', 'types.pdf');
display('Done!');
close();

