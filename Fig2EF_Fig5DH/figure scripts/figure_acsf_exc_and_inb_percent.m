%% ACSF percent of excitation and inhibition (Figure 5E)


%% Load data 

% set paths
addpath(genpath('C:\Users\nickf\Dropbox\Norway\mGluR Paper NF\code'))
acsfDir =  'C:\Users\nickf\Dropbox\Norway\mGluR Paper NF\experiments\ACSF stim imaging';

% load data
acsf = load(fullfile(acsfDir,'data.mat')); 

% organize data 
perE = [acsf.responses.perE] .* 100;
perI = [acsf.responses.perI] .* 100;

%% Plotting 

% plotting options 
figSize = [9,9]; % cm (width x height)
exColor = [0.8500 0.3250 0.0980];
inColor = [0 0.4470 0.7410];
barWidth = 0.6;
barLineW =  0.5;
barAlpha = 0.25;
errorColor = 'k';
errorLineW = 1;
groupLabels = {'Excitation', 'Inhibition'}; % X-axis labels
xAxisLabel = 'ACSF'; % X-axis label
yAxisLabel = '% Neurons'; % Y-axis label

% Create figure
hf = init_figure(figSize(1), figSize(2));

% aCSF excitation
plot_density_scatter(perE, 1, exColor,'rand');
bar(1, mean(perE), barWidth, 'facecolor', exColor, 'edgecolor', exColor, 'linew', barLineW, 'facealpha', barAlpha);
sem1 = std(perE) / sqrt(length(perE));
e1 = errorbar(1, mean(perE), sem1, sem1);
e1.Color = errorColor; e1.LineWidth = errorLineW;

% aCSF inhibition
plot_density_scatter(perI, 2, inColor,'rand');
bar(2, mean(perI), barWidth, 'facecolor', inColor, 'edgecolor', inColor, 'linew', barLineW, 'facealpha', barAlpha);
sem2 = std(perI) / sqrt(length(perI));
e2 = errorbar(2, mean(perI), sem2, sem2);
e2.Color = errorColor; e2.LineWidth = errorLineW;

% stats 
stats = compute_signrank_stats(perE, perI, 'acsf exc & inb percent');
p12 = stats.p;

% add p-value
yBuffer1 = 1;
maxY = ceil(max([perE,perI]));
add_significance_line(hf,gca,1,2, p12,'height',maxY+yBuffer1)

% Formatting 
set(gca, 'XTick', 1:2, 'XTickLabels', groupLabels);
xlabel(xAxisLabel)
ylabel(yAxisLabel);
set(gca,'XLim',[0.5 2.5])
axis square