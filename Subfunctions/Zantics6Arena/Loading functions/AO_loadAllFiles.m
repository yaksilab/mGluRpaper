function   all_fish = AO_loadAllFiles(source, timebin, folder_path_save, fish_name,conversion_factor,date_exp)
%AO_loadAllFiles - Loads the results from all experiments at once
%   Author: Anna Maria Ostenrath
%
%   Syntax:
%       all_fish = AO_loadAllFiles(source, timebin, folder_path_save, fish_name,conversion_factor,date_exp)
%
%   Description:
%       AO_loadAllFiles() - description
%    
%   Inputs:
%       source - Path to all the experimental folders
%       timebin - Timebin for the velocity
%       folder_path_save - Path to where the data should be saved in
%       fish_name - Name of the exp
%       conversion_factor - Factor to convert from pixel to mm
%       date_exp - Date when the exp was performed
%
%   Outputs:
%       all_fish - Cell array with struct for each fish
%
%   Examples: 
%       Line 1 of example
%       Line 2 of example
%       Line 3 of example
%
%   Other m-files required: AO_aj_loadFiles, AO_aj_loaddata, AO_aj_velocityCalc
%       AO_aj_pdfPosition, AO_aj_distCalculation
%   Subfunctions: none
%   MAT-files required: none
%
%   See also: OTHER_FUNCTION_NAME1,  OTHER_FUNCTION_NAME2
%   Author: Anna Maria Ostenrath 
%   Date : Nov 2022	

%% Loading the data over multiple experiments
[all_fish] = AO_aj_loadFiles(source,timebin,folder_path_save,...
             fish_name,conversion_factor,date_exp);

nROIs = size(all_fish, 1);

%% Calculation - Mean velocity
% make this that it loops over the different time bins 
tic
[all_fish]=AO_aj_velocityCalc(all_fish,nROIs,timebin,...
           folder_path_save,fish_name,conversion_factor);

clc;
disp('Jo, I am done.')
toc 

tic
[all_fish]=AO_aj_velocityCalc(all_fish,nROIs,1,...
           folder_path_save,fish_name,conversion_factor);

clc;
disp('Jo, I am done.')
toc 
%% Burst duration
% framerate = 15; 
% [all_fish] = aj_burstDuration_20211206(all_fish,nROIs,timebin,framerate);
% clc; 
% disp('Jo, I am done.')

%% Calculation - Sum of distance
% timebin=2;                           % Binlength(s) for calculation of velocity etc.
% timebin=1;
% timebin=0.5;
% timebin=0.1;
[all_fish]=AO_aj_distCalculation(all_fish,nROIs,timebin,...
           folder_path_save,fish_name,date_exp,conversion_factor);
toc
[all_fish]=AO_aj_distCalculation(all_fish,nROIs,1,...
           folder_path_save,fish_name,date_exp,conversion_factor);
% save([folder_path_save,filesep,fish_name,...
%     date_exp '_data.mat'],'all_fish','-v7.3');

clc;
disp('Jo, I am done.')

%% Calculation - probability distribution of position of the fish
% [all_fish] = aj_pdfPosition(all_fish,nROIs,idxStimulus,folder_path_save,...
%     fish_name,date_exp,figures_subfolder,group1,group2,...
%     name_group1,name_group2);
% do calculation for novel tank here too
duration = 30; 
[all_fish] = AO_aj_pdfPosition(all_fish,nROIs,duration); 
clc;
disp('Jo, I am done.')
% toc
%Then we save everything
save([folder_path_save,filesep,fish_name,...
    date_exp '_data.mat'],'all_fish','-v7.3');


end