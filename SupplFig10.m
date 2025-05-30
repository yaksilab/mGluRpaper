% Supplementarz Figure 10

%% General Info 

% Load in the data by dropping in the file or using load("YOUR URL")
metadata.data_save='X:\anna\Manuscript\Data and matlab\FigData\Fig7'; %replace the save path for your device 
load("X:\anna\code\Repositories\mGluRpaper\Subfunctions\beachVibes.mat") % replace with your link

group_names = {'NO', 'Con','CPPG'}; %change the group names 
% group_names = {'Wt', 'Het','Hom'}; %change the group names depending on the experiment you look at

load("X:\anna\code\Repositories\mGluRpaper\Subfunctions\beachVibes.mat") % replace with your link


%% Now I want to make my group variables 
groups_LDS = cell(size(metadata.GroupName,1),1); % this is for the LDS
groups_Vib = cell(size(metadata.GroupName,1),1); % this is for the startle resp
for fish = 1:size(all_fish,1)
    if all_fish{fish, 1}.group ~= 0
        if all_fish{fish,1}.stable == 1
            groups_LDS{all_fish{fish, 1}.group,1}(end +1) = fish; 
        end
        if all_fish{fish,1}.stableVib == 1
           groups_Vib{all_fish{fish, 1}.group,1}(end +1) = fish; 
    
        end
    end

end
no_group = size(metadata.GroupName,1);

group_NTT = groups_LDS; 
group_NTT{2,1}(4) = [];
group_NTT{3,1}(18) = [];

%% Plotting background 
cmap_wt = ['E4E6EB'; 'B0B3B8'; '18191A']; 
cmap_het = ['00FFFF'; '40E0D0'; '008080'];
cmap_hom = ['FF007F'; 'FF007F'; 'A94064'];
% cmap_hom = ['F89B29'; 'FF0F7B'; 'FF0F7B']; cmap_hom = ['FC8EAC';
% 'DE5D83'; 'FF0F7B']; ['F89B29'; 'FC5552'; 'FF0F7B'];
% cmap_hom = ['F89B29'; 'FF0F7B'; 'FF0F7B'];
RGB = hex2rgb(cmap_wt);
map_bef = RGB/255;

RGB = hex2rgb(cmap_het);
map_con = RGB/255;

RGB = hex2rgb(cmap_hom);
map_dru = RGB/255;
all_cmap = {map_bef, map_con, map_dru};


% if you have more groups, you should add another color. 
short_cmap = [map_bef(2,:); map_con(2,:); map_dru(2,:)] ;


%%

%% Looking at the baseline period 

% plotting each individual fish to see if there are tracking errors
addition = 20; 
collect_dist_NT = cell(no_group,1);
for group = 1:no_group
    curr_dist = [];
    figure
    title(metadata.GroupName(group,:))
    hold on
    for fish =1:size(group_NTT{group},2)
        cur_fish = group_NTT{group}(fish);
        noveltank_period = 1:600 %all_fish{cur_fish, 1}.LDSstimuliOffset(1)-20;

        curr_dist = [curr_dist, all_fish{cur_fish, 1}.binDistance(noveltank_period)];

        plot(all_fish{cur_fish, 1}.binDistance(noveltank_period) + addition*(fish-1))

    end
    collect_dist_NT{group,1} = curr_dist;
    set(gcf,'units','centimeters','Position',[2 2 10 20])
    saveas(gcf, fullfile(metadata.data_save, [metadata.GroupName(group,:) '_individual_traces.png']))
    saveas(gcf, fullfile(metadata.data_save, [metadata.GroupName(group,:) '_individual_traces.png']))
    saveas(gcf, fullfile(metadata.data_save, [metadata.GroupName(group,:) '_individual_traces.png']))
end

%% Binning 
total_length = size(collect_dist_NT{group,1},1);
% lets do binn of three seconds 
binsize = 5;
bin_num = total_length/binsize; 


collect_dist_NT_binned = cell(no_group,1);
for group = 1:no_group
    new_dist_all = [];
    for fish = 1:size(group_NTT{group},2)
        cur_dist = collect_dist_NT{group,1}(:,fish); 
        little_dist = [];
        bin_start = 1; 
        bin_end = bin_start + binsize-1;
        for bin = 1:bin_num
            little_dist = [little_dist, mean(cur_dist(bin_start:bin_end))];
            bin_start = bin_end+1;
            bin_end = bin_end +binsize; 
        
        end
        new_dist_all = [new_dist_all; little_dist];
    end
    collect_dist_NT_binned{group,1} = new_dist_all;

end

plplpl = [];
figure
time = [1:size(collect_dist_NT_binned{1,1}, 2)]; 
for group = 2:no_group
    hold on

    H1 = shadedErrorBar(time, squeeze(mean(collect_dist_NT_binned{group,1}, 1)), std(collect_dist_NT_binned{group,1},[],1)/sqrt(size(collect_dist_NT_binned{group,1},1)))
    H1.mainLine.LineWidth = 2;
    H1.mainLine.Color = short_cmap(group,:);
    H1.patch.FaceColor = short_cmap(group,:); 
    H1.patch.EdgeColor= short_cmap(group,:);

    plplpl = [plplpl, H1.mainLine]; 
end
title('Avg Distance Baseline')
% legend(plplpl, metadata.GroupName)
ylabel('distance (mm)')
xlabel('time (s)')
title('Avg Distance Baseline')
set(gcf,'units','centimeters','Position',[2 2 20 8])
ylim([0 6])
tit = 'Avg Distance Baseline'
% saveas(gcf, fullfile(metadata.data_save, [tit, '_tracesdru2.png']))
% saveas(gcf, fullfile(metadata.data_save, [tit, '_tracesdru2.svg']))
% saveas(gcf, fullfile(metadata.data_save, [tit, '_tracesdru2.fig']))

saveas(gcf, fullfile(metadata.data_save, [tit, '_tracesbinned5s.png']))
saveas(gcf, fullfile(metadata.data_save, [tit, '_tracesbinned5s.svg']))
saveas(gcf, fullfile(metadata.data_save, [tit, '_tracesbinned5s.fig']))

%% Let's look at the baseline in two steps and see if one group recovers quicker

% period1 = 1:120; % this is the two minutes
% period2 = 481:600; 

period1 = 1:180; % this is the three minutes 
period2 = 421:600; 

collect_dist_period1 = cell(no_group,1);
collect_dist_period2 = cell(no_group,1);
collect_dist_difference = cell(no_group,1);

for group = 1:no_group
    curr_dist1 = [];
    curr_dist2 = [];
    diff = [];
    for fish =1:size(group_NTT{group},2)
        cur_fish = group_NTT{group}(fish);
        % noveltank_period = 1:all_fish{cur_fish, 1}.LDSstimuliOffset(1)-20;

        curr_dist1 = [curr_dist1, all_fish{cur_fish, 1}.binDistance(period1)];
        curr_dist2 = [curr_dist2, all_fish{cur_fish, 1}.binDistance(period2)];

        diff = [diff, mean(all_fish{cur_fish, 1}.binDistance(period2)) - mean(all_fish{cur_fish, 1}.binDistance(period1))]; 


    end
    collect_dist_period1{group,1} = curr_dist1;
    collect_dist_period2{group,1} = curr_dist2;
    collect_dist_difference{group,1} = diff;

end

group_oder = [];
combined_data = [];
group_avg = []; 
sems = [];

% I need to organise the data so it can be plotted as a beeswarm with error
% bars
% no_group = size(collect_dist_NT,1); 
x_spots = [1:no_group];
avg_per_fish = cell(no_group,1);
for group = 1:no_group 

    avg_per_fish{group,1} = collect_dist_difference{group,1};;

    group_oder = [group_oder; ones(size(avg_per_fish{group,1},2),1)*x_spots(1,group)]; %ones(size(off_data,2),1)*group
    combined_data = [combined_data; avg_per_fish{group,1}'];
    group_avg = [group_avg; mean(mean(avg_per_fish{group,1},2),1)];
    sems = [sems; squeeze(nanstd(avg_per_fish{group,1},0,2)/sqrt(size(avg_per_fish{group,1},2)))];

end

tit = 'Change in Avg Distance'
figure
x = beeswarm(group_oder,combined_data, 'colormap', short_cmap)
title(tit)
hold on
er = errorbar([x_spots(1,:)],group_avg, sems)
er.Color = [0 0 0];                            
er.LineStyle = 'none';
set(gcf,'units','centimeters','Position',[2 2 10 10])
legend(metadata.GroupName)
ylabel('Change in Avg Distance')
xticks(unique(group_oder))
xticklabels(metadata.GroupName)
xlim([0 no_group+1])

saveas(gcf, fullfile(metadata.data_save, [tit, '_beswarmplot.png']))
saveas(gcf, fullfile(metadata.data_save, [tit, '_beswarmplot.svg']))
saveas(gcf, fullfile(metadata.data_save, [tit, '_beswarmplot.fig']))

% simple stats
p_values_diff = [];
if no_group == 3
    [p_1, h_1] = quick_statistic(mean(avg_per_fish{1,1},1), mean(avg_per_fish{2,1},1));
    [p_2, h_2] = quick_statistic(mean(avg_per_fish{1,1},1), mean(avg_per_fish{3,1},1));
    [p_3, h_3] = quick_statistic(mean(avg_per_fish{2,1},1), mean(avg_per_fish{3,1},1))

    p_values_diff = [p_1; p_2; p_3]; 
end

save(fullfile(metadata.data_save, 'p_values_diff_NT.mat'), 'p_values_diff')



%% Vibration experiment 
% create heatmap for all groups
collect_dist_vib = cell(no_group,1);
collect_stim_onsets_vib = cell(no_group,1);
for group = 1:no_group
    curr_dist = [];
    cur_stimuli = [];
    for fish =1:size(groups_Vib{group},2)
        cur_fish = groups_Vib{group}(fish);
        period_vib = all_fish{cur_fish, 1}.VibstimuliOnset(1)-60:all_fish{cur_fish, 1}.VibstimuliOnset(40)+60;
        % period_vib = (all_fish{cur_fish, 1}.VibstimuliOnset(1)-60)*2:(all_fish{cur_fish, 1}.VibstimuliOnset(40)+60)*2;

        curr_dist = [curr_dist, all_fish{cur_fish, 1}.binDistance(period_vib)];
        % curr_dist = [curr_dist; all_fish{cur_fish, 1}.binnedVel_0_5(period_vib)];
        cur_stimuli = [cur_stimuli, all_fish{cur_fish, 1}.VibstimuliOnset-(all_fish{cur_fish, 1}.VibstimuliOnset(1)-60)]; % I need to move the onsets to the new timeline as I cut it out... 
        
    end
    collect_dist_vib{group,1} = curr_dist;
    collect_stim_onsets_vib{group,1} = cur_stimuli; 
end
vib_stimuli_trigger = collect_stim_onsets_vib{1,1}(:,1);

tit = 'Vibr'
AO_plotavgeragesVib(collect_dist_vib, vib_stimuli_trigger, metadata, short_cmap, tit)


%%
no_con = 4;
vibr_conTrials = [[1:10]; [11:20]; [21:30]; [31:40]];
vibr_conNames = {'Vib', 'Light Bef', 'VibLight', 'Light Aft'};
vibr_data_per_stim = cell(3,1); 
vibr_change_per_stim = cell(3,1); 

%all_fish{6, 1}.VibrConTrials  all_fish{6, 1}.VibrConNames all_fish{6, 1}.VibstimuliOnset
baseline = 5; %how much before the stimulus
% baseline = 10; % for the other bin
for group = 1:no_group
    for fish = 1:size(groups_Vib{group},2)
        
        curr_fish = groups_Vib{group}(fish);
        tiral_length = length(all_fish{curr_fish, 1}.VibstimuliOnset(1)-baseline:all_fish{curr_fish, 1}.VibstimuliOnset(1+1)-baseline);

        % tiral_length = 120; %length(all_fish{curr_fish, 1}.Vib_stimOn_bin(1)-baseline:all_fish{curr_fish, 1}.Vib_stimOn_bin(1+1)-baseline);


        % vib_per_tiral = nan(tiral_length,size(all_fish{curr_fish, 1}.Vib_stimOn_bin,1));
        % vibr_change = nan(tiral_length,size(all_fish{curr_fish, 1}.Vib_stimOn_bin,1));
        % for trials = 1:size(all_fish{curr_fish, 1}.Vib_stimOn_bin,1)
        %     start_tim = all_fish{curr_fish, 1}.Vib_stimOn_bin(trials)-baseline; 
        %     end_tim = all_fish{curr_fish, 1}.Vib_stimOn_bin(trials)+120-baseline;
        %     trial_data = all_fish{curr_fish, 1}.binnedVel_0_5(start_tim:end_tim); 
        %     % time_intervals = sort(cat(1,all_fish{curr_fish, 1}.LDSstimuliOnset,all_fish{curr_fish, 1}.LDSstimuliOffset));
        %     % time_intervals = [all_fish{curr_fish, 1}.LDSstimuliOffset(1)-300,time_intervals',all_fish{curr_fish, 1}.LDSstimuliOffset(11)+300];
        %     vib_per_tiral(:,trials) = trial_data(1:tiral_length);
        %     vibr_change(:,trials) = abs(trial_data(1:tiral_length)) - abs(mean(trial_data(1:baseline)));
        % end

        vib_per_tiral = nan(tiral_length,size(all_fish{curr_fish, 1}.VibstimuliOnset,1));
        vibr_change = nan(tiral_length,size(all_fish{curr_fish, 1}.VibstimuliOnset,1));
        for trials = 1:size(all_fish{curr_fish, 1}.VibstimuliOnset,1)
            start_tim = all_fish{curr_fish, 1}.VibstimuliOnset(trials)-baseline; 
            end_tim = all_fish{curr_fish, 1}.VibstimuliOnset(trials)+60-baseline;
            trial_data = all_fish{curr_fish, 1}.binDistance(start_tim:end_tim); 
            % time_intervals = sort(cat(1,all_fish{curr_fish, 1}.LDSstimuliOnset,all_fish{curr_fish, 1}.LDSstimuliOffset));
            % time_intervals = [all_fish{curr_fish, 1}.LDSstimuliOffset(1)-300,time_intervals',all_fish{curr_fish, 1}.LDSstimuliOffset(11)+300];
            vib_per_tiral(:,trials) = trial_data;
            vibr_change(:,trials) = abs(trial_data) - abs(mean(trial_data(1:baseline)));
        end

       
        if fish == 1
            gr_togeth = vib_per_tiral;
            gr_togeth_change = vibr_change;

        else
            gr_togeth = cat(3, gr_togeth, vib_per_tiral);
            gr_togeth_change = cat(3, gr_togeth_change, vibr_change);
        end
        
    end
    vibr_data_per_stim{group,1} = gr_togeth;
    vibr_change_per_stim{group,1} = gr_togeth_change;
end

[p_vals_per_con_chnage, p_vals_per_con, p_values_vibr] = AO_plotperstimandbee_Vib(vibr_data_per_stim, vibr_change_per_stim, short_cmap, metadata, baseline)

