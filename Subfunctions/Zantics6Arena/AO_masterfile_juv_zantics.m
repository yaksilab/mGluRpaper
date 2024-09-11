%% Masterfile 6 Arena Zantics
%
%   Author: Anna Maria Ostenrath
%   Date : November 2022	

%% Variables to be adjusted for each recording
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

%% Upcoming 
% seperate the fish into your groups
% plotting NTT and vibrations 
% heatmaps%& balbla
% Changes for PDF calculations, how to load metadata tables 
% bla bla
%% Seperate animals into different groups 
% This is an example of 3 different groups 
% you might need to add the select column in both the Table and the loading
% function
name_group1 = 'Wild Type'; 
name_group2 = 'Heterozygous';
name_group3 = 'Homozygous'; 

% In this example NTT and vibration was split into two columns in the
% metatable and loaded as such
% here I seperate the groups according to NT 
% proper_exp = 6:26; 
group1 = []; %find(T.Group == 1); 
group2 = []; %find(T.Group == 2);
group3 = []; %find(T.Group == 3);
for fish = 1: size(all_fish,1)
%     if ismember(all_fish{fish, 1}.exp, proper_exp)
        if all_fish{fish, 1}.NTT == 1
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

%for some codes you still need the total number of fish, nROIS 
nROIs = size(all_fish,1);

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

%% NTT
% Now let's plot the y pos and velo for NTT 
[combined_y_pos, p_values] = AO_average_y_trace_flexibel(all_fish, n_groups, groups, names, cmap, folder_path_save, figures_subfolder, 1,1200, 0, 0, 60)
[collected_velo] = AO_average_velo_trace_felxibel(all_fish, n_groups, groups, names, cmap, folder_path_save, figures_subfolder, timebin, 1,20*60, 0);

% Heatmaps 
single_plotting = 0; % ploting individual fish of each group

for group = 1:n_groups
    AO_ploting_noveltank_heatmaps_over_fish(all_fish, groups{group}, names{group}, folder_path_save, figures_subfolder, single_plotting)
end

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

% now let's plot the change and the y poseee
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

%% How to use the new PDF plot 000

% This is the new function: [group_pdf] = AO_pdf_felxible(all_fish, group, names, folder_path_save, figures_subfolder, start_points,end_points, single_plot, save_name, plotting)
% It will create the heatmaps for you and save them in the group_pdf
% variable

% example NTT 
% making two time points 1: min 1-2 and 2: min 19-20
start_points = [1*60 19*60]; 
end_points = [2*60 20*60]; 
single_plot = 0; 
save_name = {'Early NTT', 'Late NTT'};
plotting = 1; % yes i want to plot the avg over the group
seconds = 1; % yes my values are in seconds
[wt_group_pdf] = AO_pdf_felxible(all_fish, group1, name_group1, folder_path_save, figures_subfolder, start_points,end_points, single_plot, save_name, plotting, seconds);

% example of looping over all the groups
groups = {group1; group2;group3}; 
names = {name_group1;  name_group2; name_group3}; 
n_groups = 3;

pdfs_per_group = {}; 
start_points = [1*60 19*60]; 
end_points = [2*60 20*60]; 
single_plot = 0; 
save_name = {'Early NTT', 'Late NTT'};
plotting = 1; % yes i want to plot the avg over the group
seconds = 1; % yes my values are in seconds
for group = 1:n_groups
    [pdfs_per_group{group}] = AO_pdf_felxible(all_fish, groups{group}, names{group}, folder_path_save, figures_subfolder, start_points,end_points, single_plot, save_name, plotting, seconds);
end

% example for a stimulus train will follow I might adapt some things then
%% New plot heatmaps for speed 
mkdir(folder_path_save,'\TestJune'); 
save_path = fullfile(folder_path_save, '\TestJune')
% assuming we already got the cutout traces we could just simply plot them
% as a heatmap trial by trial... 
% but maybe we also want to just see the whole trace ? 

% lets just start with the whole trace 
% this variable already has each one: collected_velo_vibr
start_point =  1200; 
end_point =  2590; 
n_bins = floor(max(all_fish{groups{1}(1),1}.t)/timebin);
new_time = linspace(0, max(all_fish{groups{1}(1),1}.t), n_bins);
start_time = find(new_time >=start_point);
start_time = start_time(1);
end_time = find(new_time >=end_point);
end_time = end_time(1); 
x_values = new_time(start_time:end_time); %1:size(collected_yposition_upper,1);

stim_times = all_fish{groups{1}(1), 1}.timeStampStim(find(~contains(all_fish{groups{1}(1), 1}.stimInfo, "VIB_0")));
stim_train = []; 
for stim = 1:length(stim_times)
    idx = find(stim_times(stim)< x_values); 
    stim_train = [stim_train; idx(1)];
end

figure('units','centimeters','Position',[2 2 28 15])
for group = 1:n_groups
    subplot(3,1,group)
    imagesc(collected_velo_vibr{group,1}(:,groups{group})') % plotting only the fish for each group
    colormap(flipud(hot))
    caxis([0 5])
    title([names{group}, ' Speed'])
        hold on
    for stim = 1:length(stim_times)
        xline(stim_train(stim), 'Color', 'b')
    end
    hold off
    colorbar
end

figure
for group = 1:n_groups
    subplot(3,1,group)
%     imagesc(squeeze(nanmean(velo_change_groups{1,group}(:,:,:),3))') % plotting only the fish for each group
    imagesc(squeeze(nanmean(velo_groups{1,group}(:,:,:),3))') % plotting only the fish for each group

    colormap jet %(flipud(hot))
    caxis([0 3])
    title([names{group}, ' Speed'])
        hold on
    for stim = 1:length(stim_times)
        xline(stim_train(stim), 'Color', 'b')
    end
    hold off
    colorbar
end
%% Heatmaps for each trial for each fish? 

% y_position_groups y_pos_change_groups
figure
for group = 1:n_groups
    subplot(3,1,group)
    imagesc(squeeze(nanmean(y_pos_change_groups{1,group}(:, :,:),2))') 
%     imagesc(squeeze(nanmean(y_position_groups{1,group}(:, :,:),2))')
    hold on
    xline(stim_indicator{1,group}, 'LineWidth', 2, 'Color', 'k', 'LineStyle', '--')
    title([names{group}, ' Y pos change'])
    caxis([-10 10])
    colorbar
    colormap(colormap_blueblackred)
%     for fish = 1:size(groups{group}, 1)
%         figure
%         imagesc(squeeze(y_pos_change_groups{1,group}(:, fish,:))') 
%         hold on
%         colormap(colormap_blueblackred)
%         xline(stim_indicator{1,group}, 'LineWidth', 2, 'Color', 'y', 'LineStyle', '--')
%         caxis([-5 5])
%         colorbar
%     end
end

% maybe if I bin the y pos I could also have a similar plot? 


figure('units','centimeters','Position',[2 2 28 15])
for group = 1:n_groups
    subplot(3,1,group)
    imagesc(combined_y_pos_vibr{group,1}(:,groups{group})') % plotting only the fish for each group
    colormap(flipud(hot))
%     caxis([0 5])
    title([names{group}, ' Y Pos'])
        hold on
%     for stim = 1:length(stim_times)
%         xline(stim_times(stim), 'Color', 'b')
%     end
    hold off
    colorbar
end

%% How about probability of response? 
n_trials = 15; 
baseline_ped = [1:55]; 
stim_perdiod = [59:65];
% what if I use velocity as a measure? 
perc_resp_trial_group = nan(n_groups, n_trials); %change the hard coding here
perc_resp_trial_group_fish = {}; 
for group = 1:n_groups
    resp_fish_per_trial = zeros(15, size(velo_groups{1,group},2)); 
    for fish = 1:size(velo_groups{1,group},2)
        cur_fish = squeeze(velo_groups{1,group}(:,fish,:)); 
        for trial = 1: n_trials 
           pre_mean_velo =  nanmean(cur_fish(1:baseline_end,trial),1);
           std_velo = nanstd(cur_fish(1:baseline_end,trial),0,1); 
           post = nanmean(cur_fish(stim_perdiod,trial),1); % dont make it hard_coded
           if post >  pre_mean_velo%+2*std_velo
               resp_fish_per_trial(trial,fish) = resp_fish_per_trial(trial,fish) + 1; 
           end
        end
        
    end
    perc_resp_trial_group_fish{group} = resp_fish_per_trial; 
    perc_resp_trial_group(group,:) = sum(resp_fish_per_trial,2)/size(velo_groups{1,group},2); 
end

figure
plots = []; 
for group = 1:n_groups
    subplot(2,3,[1 2 3])
    hold on
    p1 = plot(perc_resp_trial_group(group,:), 'Color', cmap(group,:), 'LineWidth', 2)
    title('Perc resp per trial')
    
    plots = [plots , p1];
    subplot(2,3,group+3)
    imagesc(perc_resp_trial_group_fish{group}')
    colormap hot
    title(names{group})
    
end
subplot(2,3,[1 2 3])
legend(plots, names)

saveas(gcf, fullfile(save_path, 'Perc_resp_pertrial_group.png'))

%% So now do the freezing 
%def... freezing is the immobility for more than 5 s... so it is 10 frames
%for my vibration period 

% i can use the binned distance... and sum it with a moving window maybe to
% see how many times it has been freezing? 

D = all_fish{22, 1}.dt; %all_fish{20, 1}.calcBinnedDistance_1; %all_fish{20, 1}.distance % 
X = all_fish{22, 1}.x; 
Y = all_fish{22, 1}.y; 

[freezing_percentage,freezing_time]=me_freezing_fromFP(X,Y,D,[2 2]) 

%% I want to split my vibrations into early, middle and late 

% first define which stimuli are which "condition"
% we make a new variable called conIdx that will store this
no_con = 3; % in how many conditions do you want to split? 
trials_per_con = 5; %how many trials will there be per condition
total_no_trials = 15; % how many trials do we have in total? 

trial_list = [1:total_no_trials]; 
conIdx = [1:5; 6:10; 11:15]; 

% here we define the name of the conditions
conNames =  {'Early', 'Middle', 'Late'}; 

figure
for con = 1:no_con
    subplot(1,3,con)
    hold on
    plots = [];
    for group = 1:n_groups
        p = plot(mean(mean(y_position_groups{1,group}(:,:,conIdx(con,:)),3),2), 'Color', cmap(group,:))
        plots = [plots, p];
    end
    title(conNames{con})
    legend(plots, names)
end
