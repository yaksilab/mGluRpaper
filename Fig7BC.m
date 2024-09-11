%Fig 7BC

%code us based on the zantics-6-arena codes 

%%%%% YOUR DATA HERE
source='W:\Anna\Behaviour_fromBram\Zantics\mGluR6a_NT_Vibration\';               % folder where data for experiment is found
folder_path_save=['W:\Anna\Behaviour_fromBram\Zantics\mGluR6a_NT_Vibration\'];   % folder for saving .mat file with analyzed data
fish_name='mGluR_NT_Vib';                                                % name of dataset (e.g. name of fish line + length of recording)
date_exp='_2022';                                                           % date of recording(s)                                                                 % number of experiment

table_path = dir(fullfile(source, '*.xlsx'));                            % information tabel about the experiments
T = readtable(fullfile(source,table_path(1).name));
%% Constant variables
mkdir(folder_path_save,'\Figures');                                         % makes a 'Figures' directory under the 'folder_path' directory (gets error message if already exists)
figures_subfolder='Figures';
conversion_factor=0.26458;                                                  % to convert from pixels to mm
timebin=0.5;                                                                % Binlength(s) for calculation of velocity etc.

%% Loading or creating the dataset
% locate if data file already exists
data_file = dir([source filesep '*data.mat']);
if isempty(data_file)
    all_fish = AO_loadAllFiles(source, timebin, folder_path_save,...
             fish_name,conversion_factor,date_exp);
         
else 
    % loading the file that is already there
    load(fullfile(source, data_file(1).name))
end

%% Additional variables 
% now we want to make some additional variables to make the plotting and
% everything easier later on

groups = {group1; group2;group3}; 
names = {name_group1;  name_group2; name_group3}; 
n_groups = 3;
cmap_hex = ['18191A'; '23395D'; '8B0000'; 'FFA500']; 
% 
RGB = hex2rgb(cmap_hex);
cmap = RGB/255;

%% Mechanical stimulations
% in my case I am switching to a slightly different group of fish using
% vibration as an indicator
group1 = []; %find(T.Group == 1); 
group2 = []; %find(T.Group == 2);
group3 = []; %find(T.Group == 3);
for fish = 1: size(all_fish,1)
%     if ismember(all_fish{fish, 1}.exp, proper_exp)
        if all_fish{fish, 1}.vibration == 1
%             if all_fish{fish, 1}.selected == 1
                if all_fish{fish, 1}.group == 1
%                     group1 = [group1; all_fish{fish, 1}.realNum];
                    group1 = [group1; fish];
                elseif all_fish{fish, 1}.group == 2
%                     group2 = [group2; all_fish{fish, 1}.realNum];
                    group2 = [group2; fish];

                elseif all_fish{fish, 1}.group == 3
%                     group3 = [group3; all_fish{fish, 1}.realNum];
                    group3 = [group3; fish];

                end
%             end
        end
%     end
end

% here you can also use the flexible code for y and velo with different
% time interval and obviously a "true" for stimulus 
[combined_y_pos_vibr, p_values_vibr] = AO_average_y_trace_flexibel(all_fish, n_groups, groups, names, cmap, folder_path_save, figures_subfolder, 1200,2590, 0, 1, 60)
[collected_velo_vibr] = AO_average_velo_trace_felxibel(all_fish, n_groups, groups, names, cmap, folder_path_save, figures_subfolder, timebin, 1200,2590, 1);

% Now we want to calculate the y_position around the stimulus as well as
% the change in y pos
duration = 30; % duration for before and after the stimulus to cut out from the trace
baseline_dur = 2; % length of the baseline before the stimulus for the change 
y_position_groups = {};
y_pos_change_groups = {};
stim_indicator = {}; % should always be the same but just in case 
for group = 1:n_groups
    [collected_yposition_1,change_yposition_1, marker_stim] = AO_plot_ypos_average_and_change(all_fish, duration, groups{group}, names{group}, folder_path_save, figures_subfolder, 0, baseline_dur);
    y_position_groups{group} = collected_yposition_1;
    y_pos_change_groups{group} = change_yposition_1;
    stim_indicator{group} = marker_stim;
end

% now let's plot the change and the y pos
AO_combined_traces_plot(all_fish, n_groups, groups, names, y_pos_change_groups, cmap, folder_path_save, figures_subfolder, 'Y Pos change', stim_indicator{1}, 'Change_y_pos_groups')
AO_combined_traces_plot(all_fish, n_groups, groups, names, y_position_groups, cmap, folder_path_save, figures_subfolder, 'Y Pos [mm]', stim_indicator{1}, 'Y_pos_groups')


% the same for velocity

timebin = 0.5; % timebin used for velo
duration = 30; % duration for before and after the stimulus to cut out from the trace
baseline_dur = 2; % length of the baseline before the stimulus for the change 
velo_groups = {};
velo_change_groups = {};
stim_indicator = {}; % should always be the same but just in case 
for group = 1:n_groups
    [collected_velo_1,change_velo_1, collected_stim] = AO_plot_velo_average_and_change(all_fish, groups{group}, names{group}, duration, timebin, folder_path_save, figures_subfolder,0, baseline_dur);% avg traces for each group for y and velo
    velo_groups{group} = collected_velo_1;
    velo_change_groups{group} = change_velo_1;
    stim_indicator{group} = collected_stim;
end

% something is weird with the x axis here!!! 
% now let's plot the change and the y pos
AO_combined_traces_plot(all_fish, n_groups, groups, names, velo_change_groups, cmap, folder_path_save, figures_subfolder, 'Velo 0.5s change', stim_indicator{1}(1), 'Change_Velo_groups')
AO_combined_traces_plot(all_fish, n_groups, groups, names, velo_groups, cmap, folder_path_save, figures_subfolder, 'Velo 0.5s binned', stim_indicator{1}(1), 'Velo_groups')
