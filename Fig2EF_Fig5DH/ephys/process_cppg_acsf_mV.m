%% Process CPPG & aCSF membrane potential experiments 

% set paths
addpath(genpath('C:\Users\nickf\Dropbox\Norway\mGluR Paper NF\code'))
outputDir = 'C:\Users\nickf\Dropbox\Norway\mGluR Paper NF\experiments\ACSF+CPPG ephys';
acsfPath1 = 'E:\data\Kavli\experiments\Exp. 2 connectivity spont & stim\cells\'; % local drive
acsfPath2 = 'E:\data\Kavli\experiments\Exp. 2 connectivity spont & stim\recordings.mat';
cppgPath1 = 'E:\data\Kavli\experiments\Exp. 6 CPPG connectivity\cells\';
cppgPath2 = 'E:\data\Kavli\experiments\Exp. 6 CPPG connectivity\recordings.mat';

% process data 
buffer = 2000; % buffer before and after stim
data.acsf = processExperiment(acsfPath1, acsfPath2, buffer);
data.cppg = processExperiment(cppgPath1, cppgPath2, buffer);

% save 
save(fullfile(outputDir,'data.mat'),"data","-v7.3")

%% processExperiment
function experimentData = processExperiment(dataInfoLoc, recordingFile, buffer)
% Load experiment data
load(recordingFile, 'recordings');

% Identify stimulation experiments
expIdx = logical(vertcat(recordings.expType));
stimExps = recordings(expIdx);

% Initialize data structure
experimentData = struct( ...
    'expName', [], ...
    'voltages', [], ...
    'rInputIdx', [], ...
    'stimIdx', [], ...
    'vRest', [], ...
    'vNorm', [], ...
    'iHold', [] ...
    );

% Process each experiment
for stimI = 1:numel(stimExps)
    expName = erase(stimExps(stimI).recordName, 'Stim');
    expDataFile = fullfile(dataInfoLoc, expName, 'expDataStim.mat');
    load(expDataFile, 'expData');

    % Extract IR test indices and stim indices
    rInputIdx = 1:expData.ephysData(end).RinputTest(end) + buffer;
    stimIdx = getStimIndices(expData.ephysData(end), buffer);

    % Extract voltage and apply NaN mask
    voltages = vertcat(expData.ephysData(end).voltage);
    voltages([rInputIdx, stimIdx]) = NaN;

    % Compute resting potential
    vRest = nanmean(voltages);
    iHold = vertcat(expData.ephysData.Ihold);
    vNorm = vRest' - iHold;

    % Store 
    experimentData(stimI).expName = expName;
    experimentData(stimI).voltages = voltages;
    experimentData(stimI).rInputIdx = rInputIdx;
    experimentData(stimI).stimIdx = stimIdx;
    experimentData(stimI).vRest = vRest;
    experimentData(stimI).vNorm = vNorm;
    experimentData(stimI).iHold = iHold;
end
end

%% getStimIndices
function stimIdx = getStimIndices(ephysData, buffer)
stimIdx = [];
for stim = ephysData.stims
    stimIdx = [stimIdx, (stim - buffer):(stim + ephysData.stimDuration + buffer)];
end
end
