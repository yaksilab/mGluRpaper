function [data] = AO_loadCSVfileZantiks(data_path, cfg, T)
%AO_loadCSVfileZantiks - One line description of what the function or script performs (H1 line)
%   Author: Anna Maria Ostenrath
%   Optional file header info (to give more details about the function than in the H1 line)
%   Optional file header info (to give more details about the function than in the H1 line)
%
%   Syntax:
%       output = function(input1, input2)
%       output = function(input1, input2, input3)
%
%   Description:
%       function() - description
%    
%   Inputs:
%       input1 - Description
%       input2 - Description
%       input3 - Description
%
%   Outputs:
%       output1 - Description
%       output2 - Description
%
%   Examples: 
%       Line 1 of example
%       Line 2 of example
%       Line 3 of example
%
%   Other m-files required: none
%   Subfunctions: none
%   MAT-files required: none
%
%   See also: OTHER_FUNCTION_NAME1,  OTHER_FUNCTION_NAME2
%   Author: Anna Maria Ostenrath 
%   Date : March 2024	

%% load the files
stk_files = dir(fullfile(data_path, '*.csv*' ));
% in the software released 2022, there are 2 files, one with x,y position
% and one with distance

%% IMPORT XY DATA========================================================
% load the x,y data. The file is labelled "_XY.csv"
stk_files_XY = stk_files(contains({stk_files.name}, '-XY'));

%% Import the data

% set the number of columns to import in the variable 'NumColumn'
% number of column to import is 1 + 2* the number of Arena, because column 1 is time in sec, then there are x and y for each arena
NumColumn=1+2*cfg.NumberArena ;
opts = delimitedTextImportOptions("NumVariables", NumColumn);

% Specify range and delimiter
opts.DataLines = [2, Inf]; %the line 1 is the text on the CSV file. We only import from line 2
opts.Delimiter = ",";

% Specify column names and types
opts.VariableTypes = repmat("double",1,NumColumn);

% Specify file level properties
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";

% Import the data
temp = readtable([data_path stk_files_XY.name], opts);

%% Convert to output type
temp = table2array(temp);

%% output the data in a cell array where each fish is a cell and the data for each fish is a structure array
data=cell(cfg.NumberArena,1);

for i=1:cfg.NumberArena
data{i,1}.time=temp(:,1);
data{i,1}.x=temp(:,i*2);
data{i,1}.y=temp(:,i*2+1);
end


%% Clear temporary variables
clear opts temp NumColumn

%% IMPORT DISTANCE AND METADATA DATA========================================================
% load the distance data,  the file is labelled ".csv" and do not contain
% XY
stk_files_dist = stk_files(~contains({stk_files.name}, '-XY'));

%% Import the data

% set the number of columns to import in the variable 'NumColumn'
% number of column to import is 4 plus the number of arena
NumColumn=4+cfg.NumberArena ;
opts = delimitedTextImportOptions("NumVariables", NumColumn);

% Specify range and delimiter
opts.DataLines = [5, Inf]; %the line 1 is the text on the CSV file. We only import from line 2
opts.Delimiter = ",";

% Specify column names and types
opts.VariableNames = ["RUNTIME", "TEMPERATURE", "LIGHT", repmat("well",1,cfg.NumberArena)];
opts.VariableTypes = ["double", "double", "string", repmat("double",1,cfg.NumberArena)];

% Specify file level properties
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";

% Import the distance data
distance = readtable([data_path stk_files_dist.name], opts);

% Clear temporary variables
clear opts
%% Import the metadata of the experiments, aka date, machine ID

fileID = fopen([data_path stk_files_dist.name]);
delimiter = ',';
formatSpec = '%*q%*s%s%s%*[^\n]'; 
temp = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'TextType', 'string',  'ReturnOnError', false);
metadata.date=temp{1,1}(1);
metadata.apparatus=temp{1, 2}(2);
metadata.unitID=temp{1, 2}(3);

%% output the data in a cell array where each fish is a cell and the data for each fish is a structure array

for i=1:cfg.NumberArena
data{i,1}.binTime=table2array(distance(:,1));
data{i,1}.binTemperature=table2array(distance(:,2));
data{i,1}.binDistance=table2array(distance(:,3+i));
data{i,1}.binStimulus=distance(:,3);
data{i,1}.metadata=metadata;
data{i,1}.ArenaNumber=i;
end

%% Adding some extra data from the table
% so now I need to find the right fish in the table T
% I need to identify the right index in the table 
pathparts = strsplit(data_path,filesep); %%%% changed this for the new exp!! 
% exp_name = pathparts{5};
exp_name = pathparts{4};
cur_exp_no = str2num(exp_name([5]));
% if length(exp_name) == 12
%     cur_exp_no = str2num(exp_name(6));
% else
%     cur_exp_no = str2num(exp_name([6 7]));
% 
% end
% if length(exp_name) == 12
%     cur_exp_no = str2num(exp_name(6));
% else
%     cur_exp_no = str2num(exp_name([6 7]));
% 
% end
table_idx = find(T.NumberOfExperiment == cur_exp_no);
% table_idx = find(contains(T.NumberOfExperiment, exp_name(1)));


for no_fish=1:cfg.NumberArena
    data{no_fish,1}.group = T.Group(table_idx(no_fish)); 
    data{no_fish,1}.genotype = T.genotype(table_idx(no_fish)); 
    data{no_fish,1}.stable = T.TrackingStable_unstable(table_idx(no_fish)); 
    data{no_fish,1}.fishNo = T.Fish(table_idx(no_fish)); 
    data{no_fish,1}.stableVib = T.Stable_V(table_idx(no_fish)); 
end

%% Stimulus train should also be recorded

% Light stimulus 
LightInfoMatrix = []; %define empty matrix
LightInfoMatrix = data{1,1}.binStimulus; %extract light information from cell for one fish (sufficient because light stimuli identical for each fish)
LightInfoMatrix(end,:) = []; %delete last row, contains meta data which are not needed
LightInfoMatrix = table2array(LightInfoMatrix); %convert table to array to be able to use strcmp in the following 

ChangeLightsOn  = diff(strcmp(LightInfoMatrix, 'LIGHT_ON')); %logical indexing when lights are on =1, lights off =0 and compare in same step the difference between cells, so +1 is onset of light, -1 offset of light
LightsOn = find(ChangeLightsOn == 1); %indices at which second light is turned on
LightsOn = LightsOn+1; %bc diff(strcmp..) made vector one cell shorter
LightsOff = find(ChangeLightsOn == -1); %indices at which second light is turned off
LightsOff = LightsOff+1;


% Vibration stim ... they are all called vibr although some are just light
ChangeVibr = diff(strcmp(LightInfoMatrix, 'VIBRATION'));
all_vibr = find(ChangeVibr == 1);
vibr_conNames = {'Vib', 'Light', 'VibLight', 'Light'}; 
vibr_conTrials = [[1:10]; [11:20]; [21:30]; [31:40]];

% save onset and offset in each fish
for no_fish=1:cfg.NumberArena

    data{no_fish,1}.LDSstimuliOnset=LightsOn;
    data{no_fish,1}.LDSstimuliOffset=LightsOff;

    data{no_fish,1}.VibstimuliOnset=all_vibr;
    data{no_fish,1}.VibrConNames=vibr_conNames;
    data{no_fish,1}.VibrConTrials=vibr_conTrials;
   
end

%% Velo
% Rebin the distance 
velBin = 0.5; 
% len_velo_trace = length(fish_speed(new_time>floor(stim_times(1))-duration & new_time<floor(stim_times(1))+duration));  % double check why you floor it

for no_fish=1:length(data) %parfor
    % Load variables
    t           = data{no_fish,1}.time;
    dt          = diff(data{no_fish, 1}.time); %diff(all_fish{fish, 1}.time);
%     s           = all_fish{fish,1}.binDistance;
    n_bins      = floor(max(data{no_fish,1}.time)/velBin);
    
    % calculate the distance 
    s = sqrt(diff(data{no_fish, 1}.x).^2 + diff(data{no_fish, 1}.y).^2);
    % Calculate speed over time and delta-time
    V1 = nan(1,length(dt));
    for j=1:length(dt)-1
        V1(j)= (s(j+1))/dt(j);
    end
    
    %Calculate velocity per second
    start_tim = t(1);
    bV_temp=nan(1,n_bins);
%     disp(ROI)
    new_time = [];
    for i=1:n_bins
%         disp(i)
        new_time = [new_time; start_tim + (i-1)*velBin];
        try
            bV_temp(i)= sum(s(t>(i-1)*velBin & t<=i*velBin))/...
                        sum(dt(t>(i-1)*velBin & t<=i*velBin));
        catch
            disp('Here is weird thing with dt and length of t')
            disp(no_fish)
            disp(i)
            bV_temp(i) = nan; 
        end
    end
    % Save data inn cell array
    data{no_fish,1}.speed_over_time                   =   V1;
    if velBin < 1
        velo_string = num2str(velBin); 
        velo_string(2) = '_';
        data{no_fish,1}.(['binnedVel_' velo_string])  =   bV_temp;
        data{no_fish,1}.(['speed_over_time_' velo_string])  =   V1;
        data{no_fish,1}.(['new_time_' velo_string])  =   new_time;
    else
        data{no_fish,1}.(['binnedVel_' num2str(velBin)])  =   bV_temp;
        data{no_fish,1}.(['speed_over_time_' num2str(velBin)]) =   V1;
        data{no_fish,1}.(['new_time_'  num2str(velBin)])  =   new_time;
    end

    % and then i could also already do the new stimulus onsets...
    for trial = 1:length(data{no_fish,1}.LDSstimuliOnset)
        new_on = find(new_time < data{no_fish,1}.binTime(data{no_fish,1}.LDSstimuliOnset(trial)));
        new_on = new_on(end); 
        data{no_fish,1}.LDS_stimOnset_bin(trial,1) = new_on; 


    end

    for trial = 1:length(data{no_fish,1}.LDSstimuliOffset)
        new_on = find(new_time < data{no_fish,1}.binTime(data{no_fish,1}.LDSstimuliOffset(trial)));
        new_on = new_on(end); 
        data{no_fish,1}.LDS_stimOffset_bin(trial,1) = new_on; 


    end

     for trial = 1:length(data{no_fish,1}.VibstimuliOnset)
        new_on = find(new_time < data{no_fish,1}.binTime(data{no_fish,1}.VibstimuliOnset(trial)));
        new_on = new_on(end); 
        data{no_fish,1}.Vib_stimOn_bin(trial,1) = new_on; 


    end


end


end