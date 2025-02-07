
%% Set up 

% set paths 
addpath(genpath('C:\Users\nickf\Dropbox\Norway\mGluR Paper NF\code'))
dataDir = 'E:\data\Kavli\experiments\'; % local drive

% set params 
pix2micron = 100/154;
preBuff = 20;
postBuff = 80;

%% ACSF data 
expName = 'Exp. 2 connectivity spont & stim'; 
outputDir = 'C:\Users\nickf\Dropbox\Norway\mGluR Paper NF\experiments\ACSF stim imaging';
fileName = 'recordings_stim.mat';
trialData = load_stim_data(dataDir,expName,fileName); % load data 
responses = get_stim_responses(trialData,preBuff,postBuff,pix2micron); % get cell responses 
if ~isfolder(outputDir), mkdir(outputDir), end % check for folder 
save(fullfile(outputDir,'data.mat'),"trialData","responses","-v7.3") % save data 
clearvars expName outputDir fileName trialData responses % clean workspace

%% CPPG
expName = 'Exp. 6 CPPG connectivity';
outputDir = 'C:\Users\nickf\Dropbox\Norway\mGluR Paper NF\experiments\CPPG stim imaging';
fileName = 'recordings_stim.mat';
trialData = load_stim_data(dataDir,expName,fileName); % load data 
responses = get_stim_responses(trialData,preBuff,postBuff,pix2micron); % get cell responses 
if ~isfolder(outputDir), mkdir(outputDir), end % check for folder 
save(fullfile(outputDir,'data.mat'),"trialData","responses","-v7.3") % save data 
clearvars expName outputDir fileName trialData responses % clean workspace