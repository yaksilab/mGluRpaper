%% CPPG membrane potential experiments (Figure 2F)

% set paths
addpath(genpath('C:\Users\nickf\Dropbox\Norway\mGluR Paper NF\code'))
dataDir = 'C:\Users\nickf\Dropbox\Norway\mGluR Paper NF\experiments\ACSF+CPPG ephys';
fileName = fullfile(dataDir,'data.mat');

% plotting options 
subtractInputCurrent = false; % option to subtract the input current 
figSize = [9,9]; % cm (width x height)
yLim = []; % [-90, -30]
group1Color = [0.5,0.5,0.5];
group2Color = [0.94, 0.29, 0.60];
markerSize = 15; 
errMarkerType = '+';                 
errLineWidth = 1.5;                  
errColor = 'k';               
groupLabels = {'ACSF', 'CPPG'};     % X-axis labels
yAxisLabel = 'Membrane Potential (mV)'; % Y-axis label

% load data
load(fileName)

% organize cppg data
drugVrest = {data.cppg.vRest};
drugVnorm = {data.cppg.vNorm};
drugVmed = cellfun(@median, drugVrest);
drugVmedNorm = cellfun(@median, drugVnorm);

% organize acsf data
controlVrest = {data.acsf.vRest};
controlVnorm = {data.acsf.vNorm};
Vmed = cellfun(@median, controlVrest);
VmedNorm = cellfun(@median, controlVnorm);

% select which condition (input R)
if subtractInputCurrent
    group1 = VmedNorm; group2 = drugVmedNorm;
else 
    group1 = Vmed; group2 = drugVmed; 
end 

% Wilcoxon rank-sum test
stats = compute_ranksum_stats(group1, group2, 'cppg & acsf');
p = stats.p;

% Compute mean and SEM
y1 = mean(group1); sem1 = std(group1) / sqrt(numel(group1));
y2 = mean(group2); sem2 = std(group2) / sqrt(numel(group2));

% Create figure
hf = init_figure(figSize(1), figSize(2));
hold on;

% Plot
scatter(ones(length(group1),1),group1,markerSize,group1Color,"filled")
scatter(ones(length(group2),1)+1,group2,markerSize,group2Color,"filled")

% Error bars
errorbar(1, y1, sem1, errColor, 'Marker', errMarkerType, 'LineWidth', errLineWidth);
errorbar(2, y2, sem2, errColor, 'Marker', errMarkerType, 'LineWidth', errLineWidth);

% add p-value
yBuffer1 = 2;
maxY = ceil(max([group1,group2]));
add_significance_line(hf,gca,1,2, p,'height',maxY+yBuffer1)

% Formatting
set(gca, 'XTick', [1, 2], 'XTickLabels', groupLabels);
ylabel(yAxisLabel);
set(gca,'XLim',[0 3])
if ~isempty(yLim), set(gca, 'YLim', yLim); end 
axis square