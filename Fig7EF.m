%% Fig7EF

% 
%% General Info

% Data paths 
cfg.data_path ='W:\Anna\mGluR LDS\mGluRBehav\'; %folder where data for experiment is found and saved
cfg.stk_files   = dir(fullfile(cfg.data_path, '*Exp*')); 

% ENTER VARIABLEs FROM EXPERIMENTS
cfg.NumberArena=6; %enter the numbers of arenas used in your aquisition
cfg.numStimuli = 10; %number of stimuli used in protocol  

cfg.GroupName=char('Wt','Het','Hom'); %add group names %BT:indicates days of treatment
group_names = {'Wt','Het','Hom'};
% Table path
table_path = dir(fullfile(cfg.data_path, '*.xlsx')); %define path of excel table with infos about experiments
T= readtable([table_path(1).folder filesep table_path(1).name]); %read table 
disp('Table read.')

clear table_path
cfg.data_save=[cfg.data_path 'Analyzed' filesep];
load("X:\anna\code\Repositories\Anna-Code-Collection\everyday functions\beachVibes.mat")

%% Loading the data 
% Identify all folders where each single experiment is found (1 plate) and load simultaneously both XY and distance
% data into one cell per fish

tic
[all_fish]=[];
for i=24:size(cfg.stk_files,1)
data_path_temp=[cfg.stk_files(i).folder filesep cfg.stk_files(i).name filesep] %load data from subfolders
[temp] = AO_loadCSVfileZantiks(data_path_temp, cfg, T); %restructure data to have all information from one fish in one cell
[all_fish] = cat(1, all_fish, temp); %make variable with information of all fish together
data_path_temp=[];
temp=[];
end
toc

clear data_path_temp temp

%% Save data, specify experiment
%make a folder to save data
mkdir([cfg.data_path 'Analyzed' filesep]); %make a new folder to save analyzed data
cfg.data_save=[cfg.data_path 'Analyzed' filesep]; %put info into structure
save([cfg.data_save 'all_fish_data.mat'] , 'all_fish', 'cfg', '-v7.3'); %save variable 'all_fish' and 'cfg' = metadata



%% Now I want to make my group variables 
groups_LDS = cell(size(cfg.GroupName,1),1); % this is for the LDS
groups_Vib = cell(size(cfg.GroupName,1),1); % this is for the startle resp
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
no_group = size(cfg.GroupName,1);

%% create heatmap for all groups
collect_dist_LDS = cell(no_group,1);
for group = 1:no_group
    curr_dist = [];
    for fish =1:size(groups_LDS{group},2)
        cur_fish = groups_LDS{group}(fish);
        period_LDS = all_fish{cur_fish, 1}.LDSstimuliOffset(1)-300:all_fish{cur_fish, 1}.LDSstimuliOffset(11)+300;
%         period_LDS = all_fish{cur_fish, 1}.LDS_stimOffset_bin(1)-300:all_fish{cur_fish, 1}.LDS_stimOffset_bin(11)+300;

        curr_dist = [curr_dist, all_fish{cur_fish, 1}.binDistance(period_LDS)];
%         curr_dist = [curr_dist; all_fish{cur_fish, 1}.binnedVel_0_5(period_LDS)];
        
    end
    collect_dist_LDS{group,1} = curr_dist;
    
end
figure('units','pixel','Position',[100 100 1200 1000])
%for i=1 %for WT group only
for i=1:no_group  %for all groups
        subplot(no_group,1,i), imagesc(collect_dist_LDS{i,1}') %select fish from certain group - encoded in T.Groupcounter
        colormap (flipud (hot))
        colorbar
        title (cfg.GroupName(i,:)) %change according to group name
        ylabel('Fish number')
        xlabel('time (s)')
        box ('off')
        set(gca,'TickDir','out')
        caxis([0 10])
        xline(all_fish{(groups_LDS{group}(1)), 1}.LDSstimuliOnset, '--r', 'LineWidth', 2) %vertical line for stimulus onset
        xline(all_fish{(groups_LDS{group}(1)), 1}.LDSstimuliOffset, '--k', 'LineWidth', 2) %vertical line for stimulus onset
end
sgtitle(num2str(cfg.data_path))
saveas(gcf, fullfile(cfg.data_save, ['Heatmap_exp.png']))
saveas(gcf, fullfile(cfg.data_save, ['Heatmap_exp.svg']))
% Plot curves for averaged activity with SEM

% clear BinBinDistance temp
% time=all_fish{1, 1}.binTime(1:end-1); %define x-axis
% %define nice green
% %nicegreen=[153,204,153]/255;%divide by 255 to normalize the values
% col=char('b','c','m');
% 
% binSize=10; %how many seconds to bin. The users can change the bin size.
% binTime=(1:binSize:floor(length(BinDistance)/binSize)*binSize)/60; %binTime in minute
% for ii=1:floor(length(BinDistance)/binSize)
% BinBinDistance(:,ii)=sum(BinDistance(:,ii*binSize-binSize+1:ii*binSize),2);
% end

figure('units','pixel','Position',[100 100 1500 600])
%for i=1    %for WT group only
for i=1:no_group     %for all groups
    
    temp=collect_dist_LDS{i,1}';
    shadedErrorBar(period_LDS, mean(temp,1), std(temp,[],1)/sqrt(size(temp,1)),'lineProps',col(i)), hold on
    ylabel('distance (mm)')
    xlabel('time (min)')
    set(gca,'TickDir','out')
    xline(all_fish{(groups_LDS{group}(1)), 1}.LDSstimuliOnset, '--r')
    xline(all_fish{(groups_LDS{group}(1)), 1}.LDSstimuliOffset, '--k')
   temp=[];
end
hold on;
sgtitle(num2str(cfg.data_path))
legend ("","","",cfg.GroupName(1,:),"","","","","","","","","","","","","","",cfg.GroupName(2,:),"","","","","","","","","","","","","","",cfg.GroupName(3,:),'Light ON',"","","","",'Light OFF',"","","","","")
%legend ("","","",cfg.GroupName(1,:),'Light ON',"","","","",'Light OFF',"","","","","")    %for WT group only
saveas(gcf, fullfile(cfg.data_save, ['Traces_exp.png']))
saveas(gcf, fullfile(cfg.data_save, ['Traces_exp.svg']))
%% Plotting averages per group (Binned distance) (from_ME)

num_intervals = length(time_intervals) - 1;
% 
for i=1:no_group
    avg_distance{i,1} = mean(collect_dist_LDS{i,1},2);
end
% avg_distance = cellfun(@mean,collect_dist_LDS,'UniformOutput',false);

% Initialize a cell array to store divided data
divided_avg_distance = cell(1, no_group);

% Iterate over each cell in avg_distance
for i = 1:no_group
    data = avg_distance{i}; % Get the data from the current cell
    divided_data_cell = cell(1, num_intervals);
    
    % Divide data into intervals
    for j = 1:num_intervals
        start_idx = 1+ 300*(j-1);%time_intervals(j);
        end_idx = 300 + 300*(j-1); %time_intervals(j+1);
        
        interval_data = data(start_idx:end_idx);
        divided_data_cell{j} = interval_data;
    end
    divided_avg_distance{i} = divided_data_cell;
end

average_function = @(inner_cell) mean(inner_cell);
avg_distance2 = cellfun(@(inner_cell) cellfun(average_function, inner_cell), divided_avg_distance, 'UniformOutput', false);

col = char('b','c','m');

figure,
hold on
all_plot = [];
for i=1:no_group
    p1 = scatter(x_values, avg_distance2{1,i}, 'filled', 'MarkerFaceColor', col(i), 'MarkerEdgeColor', col(i));
    plot(x_values, avg_distance2{1,i}, 'LineWidth', 1.5, 'Color', col(i));
    xlabel('Time interval')
    ylabel('Average Distance')
    axis([0 23 0 10])
    all_plot = [all_plot, p1];
    title('Average Binned Distance');
end
p2 = xline(1.5:2:21.5, '--k');
p3 = xline(2.5:2:20.5, '--r');
all_plot = [all_plot, p2(1), p3(1)];
legend(all_plot, {cfg.GroupName(1,:),cfg.GroupName(2,:), cfg.GroupName(3,:),'Light OFF','Light ON'})
saveas(gcf, fullfile(cfg.data_save, ['Avg_binned_distance.png']))
saveas(gcf, fullfile(cfg.data_save, ['Avg_binned_distance.svg']))


%% Extract all stimuli per fish from all distance data and plot

% a cell array with each cell being one fish and each row in the matrix
% being one stimuli
% be aware each datapoint is one second in BinDistance
% this will only work for an experiment with 5 stimui per interstimulus
% interval (isi)
All_fish_stimuli=[];
LDS_data_per_stim = cell(3,1); 
LDS_data_per_stim_change = cell(3,1);
baseline = 60; %how much before the stimulus
use_change = 1;
for group = 1:no_group
    for fish = 1:size(groups_LDS{group},2)
        
        curr_fish = groups_LDS{group}(fish);
        tiral_length = length(all_fish{curr_fish, 1}.LDSstimuliOffset(1)-baseline:all_fish{curr_fish, 1}.LDSstimuliOffset(1+1)+baseline);
        LDS_per_tiral = nan(tiral_length,size(all_fish{curr_fish, 1}.LDSstimuliOffset,1)-1);
        LDS_per_tiral_change = nan(tiral_length,size(all_fish{curr_fish, 1}.LDSstimuliOffset,1)-1);
        for trials = 1:size(all_fish{curr_fish, 1}.LDSstimuliOffset,1)-1
            start_tim = all_fish{curr_fish, 1}.LDSstimuliOffset(trials)-baseline; 
            end_tim = all_fish{curr_fish, 1}.LDSstimuliOffset(trials+1)+baseline;
%             trial_data = all_fish{curr_fish, 1}.binnedVel_0_5(start_tim:end_tim); 
            trial_data = all_fish{curr_fish, 1}.binDistance(start_tim:end_tim); 
            % time_intervals = sort(cat(1,all_fish{curr_fish, 1}.LDSstimuliOnset,all_fish{curr_fish, 1}.LDSstimuliOffset));
            % time_intervals = [all_fish{curr_fish, 1}.LDSstimuliOffset(1)-300,time_intervals',all_fish{curr_fish, 1}.LDSstimuliOffset(11)+300];
            LDS_per_tiral(:,trials) = trial_data;
            
            LDS_per_tiral_change(:,trials) = abs(trial_data) - abs(mean(trial_data(1:baseline)));
        end
        if fish == 1
            gr_togeth = LDS_per_tiral;
            gr_togeth_change = LDS_per_tiral_change;

        else
            gr_togeth = cat(3, gr_togeth, LDS_per_tiral);
            gr_togeth_change = cat(3, gr_togeth_change, LDS_per_tiral_change);
        end
        
    end
    LDS_data_per_stim{group,1} = gr_togeth; 
    LDS_data_per_stim_change{group,1} = gr_togeth_change;
end
if use_change 
    temp = LDS_data_per_stim; 
    LDS_data_per_stim = LDS_data_per_stim_change; 
    name_add = 'Change';
else
    name_add = '';
end
    
figure('units','pixel','Position',[100 100 1200 1000])

for group=1:no_group
   % for fish = 1:size(groups_LDS{group},2)
   % 
   %      curr_fish = groups_LDS{group}(fish);
   %      temp(fish,:) = mean(LDS_data_per_stim{group,1}, 2);
   % end
   % 
   subplot(no_group,1,group), imagesc(squeeze(mean(LDS_data_per_stim{group,1}, 2))')
   %sgtitle(cfg.ISIname(plotISI,:)) %add overall title to subplots
    
    colorbar
    title ([ name_add, ' ', cfg.GroupName(group,:)]) %title for each individual subplot
    ylabel('Fish nr.')
    xlabel('time (s)')
    box ('off')
    set(gca,'TickDir','out')
    if use_change
        caxis([-4 8])
        colormap (beachVibes)
    else
        caxis([0 10])
        colormap (flipud(hot))
    end
    
    xline(baseline, '--k')
    xline(300+baseline, '--k')
    xline(600+baseline, '--k')
    temp=[];
    sgtitle(['average all stimuli for ' num2str(cfg.data_path)])
end
saveas(gcf, fullfile(cfg.data_save, [name_add, '_','Avg_heatmap.png']))
saveas(gcf, fullfile(cfg.data_save, [name_add, '_','Avg_heatmap.svg']))
% plot curves for averaged stimuli with SEM

time=[-baseline:600+baseline]; %define x-axis
col=char('b','c','m');
figure('units','pixel','Position',[100 100 1500 600])
hold on
plpl = [];
for group=1:no_group

    H1 = shadedErrorBar(time, squeeze(mean(mean(LDS_data_per_stim{group,1}, 2),3)), std(mean(LDS_data_per_stim{group,1}, 2),[],3)/sqrt(size(squeeze(mean(LDS_data_per_stim{group,1}, 2)),2)),'lineProps',col(group)), hold on
    H1.mainLine.LineWidth = 2;
    ylabel('distance (mm)')
    xlabel('time (s)')
    set(gca,'TickDir','out')
    xline(0, '--k')
    xline(300, '--r')
   plpl = [plpl, H1.mainLine ]
   %title(cfg.ISIname(plotISI,:))
   sgtitle(['average all stimuli for ' num2str(cfg.data_path) name_add])
%legend("","","",cfg.GroupName(1,:),"","","","","","","","","","","","","","",cfg.GroupName(2,:),"","","","","","","","","","","","","","",cfg.GroupName(3,:),"","","","","","","","","","","","","","",cfg.GroupName(4,:),'Light ON',"","","","",'Light OFF',"","","","","")
end
legend(plpl, cfg.GroupName);

% legend ("","","",cfg.GroupName(1,:),"","","","","",cfg.GroupName(2,:),"","","","","",cfg.GroupName(3,:),'Light OFF','Light ON')

saveas(gcf, fullfile(cfg.data_save, [name_add, '_','Avg_Traces_exp.png']))
saveas(gcf, fullfile(cfg.data_save, [name_add, '_','Avg_Traces_exp.svg']))
x_spots = [[1.5 2 2.5]; [4.5 5 5.5]];
avg_off = cell(3,1);
avg_on = cell(3,1);

time_per = 30;
figure('units','centimeters','Position',[2 2 10 10])
hold on
plplpl = [];
group_oder = [];
combined_data = [];
grouop_avg = []; 
sems = [];
for group = 1:no_group 
    % so i I want to scatter the first 10 s after stim 
    off_data = squeeze(mean(LDS_data_per_stim{group,1}(baseline:baseline+time_per,:,:), 1));
    on_data = squeeze(mean(LDS_data_per_stim{group,1}(300:300+time_per,:,:), 1)); 
    avg_on{group,1} = on_data;
    avg_off{group,1} = off_data;
    
    a = scatter(ones(size(off_data,2),1)*x_spots(1,group), mean(off_data,1), 'filled', 'MarkerFaceColor', col(group), 'MarkerEdgeColor', col(group))
    scatter(ones(size(on_data,2),1)*x_spots(2,group), mean(on_data,1), 'filled', 'MarkerFaceColor', col(group), 'MarkerEdgeColor', col(group))
    plplpl = [plplpl, a];
    er = errorbar([x_spots(:,group)],[mean(mean(off_data,2),1), mean(mean(on_data,2),1)],[squeeze(nanstd(mean(off_data),0,2)/sqrt(size(off_data,2))), squeeze(nanstd(mean(on_data),0,2)/sqrt(size(on_data,2)))]);    
    er.Color = [0 0 0];                            
    er.LineStyle = 'none';
    
    group_oder = [group_oder; ones(size(off_data,2),1)*x_spots(1,group)]; %ones(size(off_data,2),1)*group
    combined_data = [combined_data; mean(off_data,1)'];
    grouop_avg = [grouop_avg; mean(mean(off_data,2),1)];
    sems = [sems; squeeze(nanstd(mean(off_data),0,2)/sqrt(size(off_data,2)))];
end
title([name_add,' with timebin ', num2str(time_per), ' s after stim '])
xticks([2 5]);
xticklabels({'Off trans', 'On trans'})
xlim([0 7])
legend(plplpl, group_names)
ylabel(['Avg binned distance [mm] with timebin ', num2str(time_per) , ' s after stim ' name_add, ' ',])
saveas(gcf, fullfile(cfg.data_save, [name_add, '_',num2str(time_per) '_Avg_transition_scatter.png']))
saveas(gcf, fullfile(cfg.data_save, [name_add, '_',num2str(time_per) '_Avg_transition_scatter.svg']))

figure
x = beeswarm(group_oder,combined_data)
title([name_add,' with timebin ', num2str(time_per), ' s after stim '])
hold on
er = errorbar([x_spots(1,:)],grouop_avg, sems)
er.Color = [0 0 0];                            
er.LineStyle = 'none';
saveas(gcf, fullfile(cfg.data_save, [name_add, '_',num2str(time_per) '_Avg_transition_scatterBEE.png']))
saveas(gcf, fullfile(cfg.data_save, [name_add, '_',num2str(time_per) '_Avg_transition_scatterBEE.svg']))
saveas(gcf, fullfile(cfg.data_save, [name_add, '_',num2str(time_per) '_Avg_transition_scatterBEE.fig']))

[p_1, h_1] = quick_statistic(mean(avg_on{1,1},1), mean(avg_on{2,1},1))
[p_2, h_2] = quick_statistic(mean(avg_on{1,1},1), mean(avg_on{3,1},1))
[p_3, h_3] = quick_statistic(mean(avg_on{2,1},1), mean(avg_on{3,1},1))
pval_LDS.(['timebin_on', num2str(time_per) ]) = [p_1; p_2; p_3]; 
[p_1, h_1] = quick_statistic(mean(avg_off{1,1},1), mean(avg_off{2,1},1))
[p_2, h_2] = quick_statistic(mean(avg_off{1,1},1), mean(avg_off{3,1},1))
[p_3, h_3] = quick_statistic(mean(avg_off{2,1},1), mean(avg_off{3,1},1))
pval_LDS.(['timebin_off', num2str(time_per) ]) = [p_1; p_2; p_3];  
