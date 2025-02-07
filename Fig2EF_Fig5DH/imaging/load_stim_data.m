function data = load_stim_data(dataDir,expName,fileName)

% load data
% dataDir = 'E:\data\Kavli\exps\';
% expName = 'Exp. 2 connectivity spont & stim';
% expName = 'Exp. 6 CPPG connectivity';
% fileName = 'recordings_stim.mat';

% load exp info struct
load(fullfile(dataDir,expName,fileName));

% combine data
data = struct;  tic
for iCell = 1:size(recordings_stim,2)
    
    % load data
    cellName = recordings_stim(iCell).recordName;
    cellName = erase(cellName,'Stim');
    load(fullfile(dataDir,expName,'cells',cellName,'expDataStim.mat'))
    
     % select calcium data
    dataTmp = expData.calciumSigs;
    dataTmp = permute(reshape(dataTmp,size(dataTmp,1),[],length(recordings_stim(iCell).sweeps)),[2,1,3]);
    
    data(iCell).patchTrace = squeeze(dataTmp(:,1,:));
    data(iCell).allTrace = dataTmp(:,1:end,:);
    data(iCell).calTimeS = expData.imagingData.times;
    
    % change the stim times to account for the data lost because of the light pulse
    pulseBuffer = (length(expData.imagingData(1).times)-size(dataTmp,1))/2;
    data(iCell).calStims = [expData.imagingData(1).stims] - pulseBuffer;
    
    % select ephys data
    data(iCell).volt = expData.ephysData(end).voltage;
    data(iCell).voltTimeS = expData.ephysData.times;
    data(iCell).voltStims = expData.ephysData.stims;
    
    % select distances
    data(iCell).dist = expData.cellDist;
    clear expData
    
end
toc 