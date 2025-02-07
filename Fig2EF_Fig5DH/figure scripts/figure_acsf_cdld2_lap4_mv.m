%% CdCl2 + L-AP4 membrane potential experiments (Figure 2E)

% set paths
addpath(genpath('C:\Users\nickf\Dropbox\Norway\mGluR Paper NF\code'))
dataDir = 'C:\Users\nickf\Dropbox\Norway\mGluR Paper NF\experiments\CdCl2+L-AP4 ephys';
fileName = fullfile(dataDir,'data.mat');

% load data
load(fileName)

% organize data
group1 = vertcat(data.combined.no_drug.Vrest_median);
group2 = vertcat(data.combined.CdCl2.Vrest_median);
group3 = vertcat(data.combined.CdCl2_LAP4.Vrest_median);

% compute stats
stats12 = compute_signrank_stats(group1, group2, 'group 1 & 2');
group12PV = stats12.p;
stats23 = compute_signrank_stats(group2, group3, 'group 2 & 3');
group23PV = stats23.p;
stats13 = compute_signrank_stats(group1, group3, 'group 1 & 3');
group13PV = stats13.p;

% compute mean and sem 
asem = nanstd(group1)/sqrt(size(group1,1)); ya = nanmean(group1); 
asemb = nanstd(group2)/sqrt(size(group2,1)); yb = nanmean(group2); 
asemc = nanstd(group3)/sqrt(size(group3,1)); yc =nanmean(group3);

% set figure options 
ylab = 'Membrane Potential (mV)';
xlabels = {'ACSF','CdCl_2','CdCl_2+ L-AP4'};
xlabels = cellfun(@(x) strrep(x,' ','\newline'), xlabels,'UniformOutput',false);
errorBarMarker = '+'; 
errorBarColor = 'k';
errorBarLineW = 1.5;

% initialize figure
hf = init_figure(9, 9);
ax1 = gca;

% plot groups
hold on
x1 = ones(length(group1),1);
x2 = ones(length(group2),1) + 1;
x3 = ones(length(group3),1) + 2;
for i = 1:length(group1)
   plot([x1(i),x2(i),x3(i)],[group1(i),group2(i),group3(i)], '--o', 'color',[0.5,0.5,0.5],'MarkerEdgeColor', [0.5,0.5,0.5], 'linewidth', 0.5,'MarkerFaceColor',[0.5,0.5,0.5])
end 

% error bars
e1 = errorbar(1,ya,asem); e1.Marker = errorBarMarker; e1.LineWidth = errorBarLineW; e1.Color = errorBarColor;
e2 = errorbar(2,yb,asemb); e2.Marker = errorBarMarker; e2.LineWidth = errorBarLineW; e2.Color = errorBarColor;
e3 = errorbar(3,yc,asemc); e3.Marker = errorBarMarker; e3.LineWidth = errorBarLineW; e3.Color = errorBarColor;

% add axis labels
ax1.XTickLabel = xlabels;
ylabel(ylab)

% format axis 
set(ax1,'XTick',([1 2 3]));
set(ax1,'XLim',[0.5 3.5])
set(ax1, 'YLim', [-90, -30]);
axis square

% add p-values 
yBuffer1 = 1;
yBuffer2 = 5;
maxY = ceil(max([group1;group2;group3]));
add_significance_line(hf,ax1,1,2, group12PV,'height',maxY+yBuffer1)
add_significance_line(hf,ax1,2,3, group23PV,'height',maxY+yBuffer1)
add_significance_line(hf,ax1,1,3, group13PV,'height',maxY+yBuffer2)