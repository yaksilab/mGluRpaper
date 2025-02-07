%% ACSF Stim responses (Figure 5D)

%% Load data

% set paths
addpath(genpath('C:\Users\nickf\Dropbox\Norway\mGluR Paper NF\code'))
dataDir =  'C:\Users\nickf\Dropbox\Norway\mGluR Paper NF\experiments\ACSF stim imaging';

% load data
data = load(fullfile(dataDir,'data.mat'));

preBuff = 20;
postBuff = 80;
timeS = (-preBuff:postBuff-1)/10;

%% Plotting

% plotting options
figSize = [18,18]; % cm (width x height)
patchColor = [0.6350 0.0780 0.1840];
exColor = [0.8500 0.3250 0.0980];
inColor = [0 0.4470 0.7410];
nrColor = [0.5,0.5,0.5];
respLineW = 1;
patchLineW = 2;

% Create figure
hf = init_figure(figSize(1), figSize(2));
tile = tiledlayout('flow','TileSpacing','none','Padding','tight');
for iT = 1:size(data.responses,2)
    nexttile, hold on
   
    % get indices 
    nrIdx = find([data.responses(iT).nrCellIdx]);
    inIdx = find([data.responses(iT).inCellIdx]);
    exIdx = find([data.responses(iT).exCellIdx]);

    % no response
    for ii = 1:length(nrIdx)
        tmpTrace = mean(squeeze(data.responses(iT).cellTraces(:,:,nrIdx(ii))),2);
        plot(timeS,tmpTrace,'Color',nrColor,'LineWidth',respLineW)
        xticks([]), yticks([])
    end

    % inhibited
    for ii = 1:length(inIdx)
        tmpTrace = mean(squeeze(data.responses(iT).cellTraces(:,:,inIdx(ii))),2);
        plot(timeS,tmpTrace,'Color',inColor,'LineWidth',respLineW)
        xticks([]), yticks([])
    end

    % excited
    for ii = 1:length(exIdx)
        tmpTrace = mean(squeeze(data.responses(iT).cellTraces(:,:,exIdx(ii))),2);
        plot(timeS,tmpTrace,'Color',exColor,'LineWidth',respLineW)
        xticks([]), yticks([])
    end

    % patched cell
    tmpTrace = mean(data.responses(iT).patchResp,2);
    yyaxis right
    plot(timeS,tmpTrace,'Color',patchColor,'LineWidth',patchLineW)
    xticks([]), yticks([])
    ax = gca;
    align_yyaxis_zero(ax)
end


