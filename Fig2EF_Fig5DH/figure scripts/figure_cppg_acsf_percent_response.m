%% ACSF + CPPG percent of responsive cells (Figure 5G/H)


%% Load data 

% set paths
addpath(genpath('C:\Users\nickf\Dropbox\Norway\mGluR Paper NF\code'))
acsfDir =  'C:\Users\nickf\Dropbox\Norway\mGluR Paper NF\experiments\ACSF stim imaging';
cppgDir = 'C:\Users\nickf\Dropbox\Norway\mGluR Paper NF\experiments\CPPG stim imaging';

% load data
acsf = load(fullfile(acsfDir,'data.mat')); 
cppg = load(fullfile(cppgDir,'data.mat'));

% organize data 
perE = [acsf.responses.perE] .* 100;
perI = [acsf.responses.perI] .* 100;
cppgPerE = [cppg.responses.perE] .* 100;
cppgPerI = [cppg.responses.perI] .* 100;

%% Plotting 

% plotting options 
figSize = [9,9]; % cm (width x height)
acsfColor = [0.5,0.5,0.5];
cppgColor = [0.94, 0.29, 0.60];
barWidth = 0.6;
barLineW =  0.5;
barAlpha = 0.25;
errorColor = 'k';
errorLineW = 1;
groupLabels = {'ACSF', 'CPPG','ACSF', 'CPPG'}; % X-axis labels
xAxisLabel = 'Excitation        Inhibition'; % X-axis label
yAxisLabel = '% Neurons'; % Y-axis label

% Create figure
hf = init_figure(figSize(1), figSize(2));

% aCSF excitation
plot_density_scatter(perE, 1, acsfColor,'rand');
bar(1, mean(perE), barWidth, 'facecolor', acsfColor, 'edgecolor', acsfColor, 'linew', barLineW, 'facealpha', barAlpha);
sem1 = std(perE) / sqrt(length(perE));
e1 = errorbar(1, mean(perE), sem1, sem1);
e1.Color = errorColor; e1.LineWidth = errorLineW;

% CPPG excitation
plot_density_scatter(cppgPerE, 2, cppgColor,'rand');
bar(2, mean(cppgPerE), barWidth, 'facecolor', cppgColor, 'edgecolor', cppgColor, 'linew', barLineW, 'facealpha', barAlpha);
sem2 = std(cppgPerE) / sqrt(length(cppgPerE));
e2 = errorbar(2, mean(cppgPerE), sem2, sem2);
e2.Color = errorColor; e2.LineWidth = errorLineW;

% aCSF inhibition
plot_density_scatter(perI, 3, acsfColor,'rand');
bar(3, mean(perI), barWidth, 'facecolor', acsfColor, 'edgecolor', acsfColor, 'linew', barLineW, 'facealpha', barAlpha);
sem3 = std(perI) / sqrt(length(perI));
e3 = errorbar(3, mean(perI), sem3, sem3);
e3.Color = errorColor; e3.LineWidth = errorLineW;

% CPPG inhibition
plot_density_scatter(cppgPerI, 4, cppgColor,'rand');
bar(4, mean(cppgPerI), barWidth, 'facecolor', cppgColor, 'edgecolor', cppgColor, 'linew', barLineW, 'facealpha', barAlpha);
sem4 = std(cppgPerI) / sqrt(length(cppgPerI));
e4 = errorbar(4, mean(cppgPerI), sem4, sem4);
e4.Color = errorColor; e4.LineWidth = errorLineW;

% stats 
statsExc = compute_ranksum_stats(perE, cppgPerE, 'cppg & acsf exc');
p12 = statsExc.p;
statsInb = compute_ranksum_stats(perI, cppgPerI, 'cppg & acsf inb');
p34 = statsInb.p;

% add p-value
yBuffer1 = 1;
maxY = ceil(max([perE,cppgPerE,perI,cppgPerI]));
add_significance_line(hf,gca,1,2, p12,'height',maxY+yBuffer1)
add_significance_line(hf,gca,3,4, p34,'height',maxY+yBuffer1)

% Formatting 
set(gca, 'XTick', 1:4, 'XTickLabels', groupLabels);
ylabel(yAxisLabel);
set(gca,'XLim',[0.5 4.5])
axis square
xlabel(xAxisLabel)