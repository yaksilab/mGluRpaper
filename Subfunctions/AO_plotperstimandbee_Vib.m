function [p_vals_per_con_change, p_vals_per_con, p_values_vibr] = AO_plotperstimandbee_Vib(vibr_data_per_stim, vibr_change_per_stim, short_cmap, metadata, baseline)
%AO_plotperstimandbee_Vib - plots the average over the diff stimulus
%conditions and compares it between groups
%   Author: Anna Maria Ostenrath
%   Optional file header info (to give more details about the function than in the H1 line)
%   Optional file header info (to give more details about the function than in the H1 line)
%
%   Syntax:
%       output = function(input1, input2)
%       output = function(input1, input2, input3)
%
%   Description:
%       AO_plotperstimandbee_Vib() - plots the average over the diff stimulus
%conditions and compares it between groups
%    
%   Inputs:
%       vibr_data_per_stim - cell array with the average data per fish per
%       stim per group
%       vibr_change_per_stim - cell array with the change per fish per
%       stim per group
%       short_cmap - colors for the traces (rgb)
%       metadata - struct with extra info like savepath and groupnames
%       baseline - integer 
%       input3 - Description
%
%   Outputs:
%       p_vals_per_con - pvalues comparing the different groups only works
%       with 3 groups
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
%   Date : September 2024	


no_con = 1; %4;
vibr_conTrials = [[1:10]; [11:20]; [21:30]; [31:40]];
vibr_conNames = {'Vib', 'Light Bef', 'VibLight', 'Light Aft'};
no_group = size(vibr_data_per_stim,1);
% now we can plot the average for each con 
time = -baseline:60-baseline; 
% time = -baseline:120-baseline; 
figure('units','centimeters','Position',[2 2 30 14])
for con = 1:no_con
    % subplot(2,2,con)
    hold on
    plpl = [];
    for group = 1:no_group
        % plot(mean(mean(vibr_data_per_stim{group,1}(:,vibr_conTrials(con,:),:),2),3), )

        H1 = shadedErrorBar(time(), squeeze(mean(mean(vibr_data_per_stim{group,1}(:,vibr_conTrials(con,:),:),2),3)), std(mean(vibr_data_per_stim{group,1}(:,vibr_conTrials(con,:),:),2),[],3)/sqrt(size(squeeze(mean(vibr_data_per_stim{group,1}(:,vibr_conTrials(con,:),:),2)),2))), hold on
        H1.mainLine.LineWidth = 2;
        H1.mainLine.Color = short_cmap(group,:);
        H1.patch.FaceColor = short_cmap(group,:); 
        H1.patch.EdgeColor= short_cmap(group,:);
        ylabel('distance (mm)')
        xlabel('time (s)')
        title(vibr_conNames{con})
        plpl = [plpl, H1.mainLine];
    end
    xline(0, '--k')
    legend(plpl, metadata.GroupName)
    xlim([-5 10]) % used to be 30
    ylim([0 12])
    
end
sgtitle('Average')
saveas(gcf, fullfile(metadata.data_save, ['Avg_Traces_exp_vib.png']))
saveas(gcf, fullfile(metadata.data_save, ['Avg_Traces_exp_vib.svg']))
saveas(gcf, fullfile(metadata.data_save, ['Avg_Traces_exp_vib.fig']))


figure('units','centimeters','Position',[2 2 30 14])
for con = 1:no_con
    % subplot(2,2,con)
    hold on
    plpl = [];
    for group = 1:no_group
        % plot(mean(mean(vibr_data_per_stim{group,1}(:,vibr_conTrials(con,:),:),2),3), )

        H1 = shadedErrorBar(time(), squeeze(nanmean(nanmean(vibr_change_per_stim{group,1}(:,vibr_conTrials(con,:),:),2),3)), nanstd(nanmean(vibr_change_per_stim{group,1}(:,vibr_conTrials(con,:),:),2),[],3)/sqrt(size(squeeze(mean(vibr_change_per_stim{group,1}(:,vibr_conTrials(con,:),:),2)),2))), hold on
        H1.mainLine.LineWidth = 2;
        H1.mainLine.Color = short_cmap(group,:);
        H1.patch.FaceColor = short_cmap(group,:); 
        H1.patch.EdgeColor= short_cmap(group,:);
        ylabel('chnage distance (mm)')
        xlabel('time (s)')
        title(vibr_conNames{con})
        plpl = [plpl, H1.mainLine];
    end
    xline(0, '--k')
    legend(plpl, metadata.GroupName)
    xlim([-2 5]) %used to be 20
    ylim([-3 12])
    
end
sgtitle('Change')
saveas(gcf, fullfile(metadata.data_save, ['Avg_Traces_exp_vib_change.png']))
saveas(gcf, fullfile(metadata.data_save, ['Avg_Traces_exp_vib_change.svg']))
saveas(gcf, fullfile(metadata.data_save, ['Avg_Traces_exp_vib_change.fig']))

% Now the scatter for the second after maybe? 
x_spots = [[1.5 2 2.5]; [4.5 5 5.5]];
avg_off_vib = cell(3,4);
avg_on_vib = cell(3,4);

figure('units','centimeters','Position',[2 2 30 14])
hold on

for con = 1:no_con
    plplpl = [];
    for group = 1:no_group 
    % so i I want to scatter the first 10 s after stim

        subplot(2,2,con)
        hold on
        
        off_data = squeeze(mean(vibr_change_per_stim{group,1}([8,9,10],vibr_conTrials(con,:),:), 1));
        on_data = squeeze(mean(vibr_change_per_stim{group,1}([7],vibr_conTrials(con,:),:), 1)); 
        % off_data = squeeze(mean(vibr_data_per_stim{group,1}(8,vibr_conTrials(con,:),:), 1));
        % on_data = squeeze(mean(vibr_data_per_stim{group,1}(7,vibr_conTrials(con,:),:), 1)); 
        avg_on_vib{group,con} = on_data;
        avg_off_vib{group,con} = off_data;
        
        a = scatter(ones(size(off_data,2),1)*x_spots(2,group), mean(off_data,1), 'filled', 'MarkerFaceColor', short_cmap(group,:), 'MarkerEdgeColor', short_cmap(group,:))
        scatter(ones(size(on_data,2),1)*x_spots(1,group), mean(on_data,1), 'filled', 'MarkerFaceColor', short_cmap(group,:), 'MarkerEdgeColor', short_cmap(group,:))
        plplpl = [plplpl, a];
        er = errorbar([x_spots([2 1],group)],[mean(mean(off_data,2),1), mean(mean(on_data,2),1)],[squeeze(nanstd(mean(off_data),0,2)/sqrt(size(off_data,2))), squeeze(nanstd(mean(on_data),0,2)/sqrt(size(on_data,2)))]);    
        er.Color = [0 0 0];                            
        er.LineStyle = 'none';
    end
    xticks([2 5]);
    xticklabels({'first s', 'sec s'})
    xlim([0 7])
    legend(plplpl, metadata.GroupName, 'Location', 'Eastoutside')
    ylabel(['Avg binned distance change [mm] ' ])
    title(vibr_conNames{con})
end

saveas(gcf, fullfile(metadata.data_save, ['Change_transition_scatter_vib.png']))
saveas(gcf, fullfile(metadata.data_save, ['Change_transition_scatter_vib.svg']))
saveas(gcf, fullfile(metadata.data_save, ['Change_transition_scatter_vib.fig']))


p_vals_per_con_change = cell(4,2);
if no_group == 3
    for con = 1:no_con
        [p_1, h_1] = quick_statistic(mean(avg_on_vib{1,con},1), mean(avg_on_vib{2,con},1))
        [p_2, h_2] = quick_statistic(mean(avg_on_vib{1,con},1), mean(avg_on_vib{3,con},1))
        [p_3, h_3] = quick_statistic(mean(avg_on_vib{2,con},1), mean(avg_on_vib{3,con},1))
        p_vals_per_con_change{con,1} = [p_1; p_2; p_3]; 
    
        [p_1, h_1] = quick_statistic(mean(avg_off_vib{1,con},1), mean(avg_off_vib{2,con},1))
        [p_2, h_2] = quick_statistic(mean(avg_off_vib{1,con},1), mean(avg_off_vib{3,con},1))
        [p_3, h_3] = quick_statistic(mean(avg_off_vib{2,con},1), mean(avg_off_vib{3,con},1))
        p_vals_per_con_change{con,2} = [p_1; p_2; p_3];
    
    end
end
%%
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
    off_data = squeeze(mean(vibr_change_per_stim{group,1}([8,9,10],vibr_conTrials(1,:),:), 1));
    avg_per_fish{group,1} = mean(off_data);

    group_oder = [group_oder; ones(size(avg_per_fish{group,1},2),1)*x_spots(1,group)]; %ones(size(off_data,2),1)*group
    combined_data = [combined_data; avg_per_fish{group,1}'];
    group_avg = [group_avg; mean(mean(avg_per_fish{group,1},2),1)];
    sems = [sems; squeeze(nanstd(avg_per_fish{group,1},0,2)/sqrt(size(avg_per_fish{group,1},2)))];

end

tit = 'Change after initial'
figure
x = beeswarm(group_oder,combined_data, 'colormap', short_cmap)
title(tit)
hold on
er = errorbar([x_spots(1,:)],group_avg, sems)
er.Color = [0 0 0];                            
er.LineStyle = 'none';
set(gcf,'units','centimeters','Position',[2 2 10 10])
legend(metadata.GroupName, 'Location', 'south')
ylabel('Change in Avg Distance')
xticks(unique(group_oder))
xticklabels(metadata.GroupName)
xlim([0 no_group+1])

saveas(gcf, fullfile(metadata.data_save, [tit, '_beswarmplot.png']))
saveas(gcf, fullfile(metadata.data_save, [tit, '_beswarmplot.svg']))
saveas(gcf, fullfile(metadata.data_save, [tit, '_beswarmplot.fig']))

% simple stats
p_values_vibr_off = [];
if no_group == 3
    [p_1, h_1] = quick_statistic(mean(avg_per_fish{1,1},1), mean(avg_per_fish{2,1},1));
    [p_2, h_2] = quick_statistic(mean(avg_per_fish{1,1},1), mean(avg_per_fish{3,1},1));
    [p_3, h_3] = quick_statistic(mean(avg_per_fish{2,1},1), mean(avg_per_fish{3,1},1))

    p_values_vibr_off = [p_1; p_2; p_3]; 
end
p_values_vibr.p_values_vibr_off = p_values_vibr_off;
% save(fullfile(metadata.data_save, 'p_values_inititalvibr.mat'), 'p_values_vibr')

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
    off_data = squeeze(mean(vibr_change_per_stim{group,1}([7],vibr_conTrials(1,:),:), 1));
    avg_per_fish{group,1} = mean(off_data);

    group_oder = [group_oder; ones(size(avg_per_fish{group,1},2),1)*x_spots(1,group)]; %ones(size(off_data,2),1)*group
    combined_data = [combined_data; avg_per_fish{group,1}'];
    group_avg = [group_avg; mean(mean(avg_per_fish{group,1},2),1)];
    sems = [sems; squeeze(nanstd(avg_per_fish{group,1},0,2)/sqrt(size(avg_per_fish{group,1},2)))];

end

tit = 'Change at initial'
figure
x = beeswarm(group_oder,combined_data, 'colormap', short_cmap)
title(tit)
hold on
er = errorbar([x_spots(1,:)],group_avg, sems)
er.Color = [0 0 0];                            
er.LineStyle = 'none';
set(gcf,'units','centimeters','Position',[2 2 10 10])
legend(metadata.GroupName, 'Location', 'south')
ylabel('Change in Avg Distance')
xticks(unique(group_oder))
xticklabels(metadata.GroupName)
xlim([0 no_group+1])

saveas(gcf, fullfile(metadata.data_save, [tit, '_beswarmplot.png']))
saveas(gcf, fullfile(metadata.data_save, [tit, '_beswarmplot.svg']))
saveas(gcf, fullfile(metadata.data_save, [tit, '_beswarmplot.fig']))

% simple stats
p_values_vibr_on = [];
if no_group == 3
    [p_1, h_1] = quick_statistic(mean(avg_per_fish{1,1},1), mean(avg_per_fish{2,1},1));
    [p_2, h_2] = quick_statistic(mean(avg_per_fish{1,1},1), mean(avg_per_fish{3,1},1));
    [p_3, h_3] = quick_statistic(mean(avg_per_fish{2,1},1), mean(avg_per_fish{3,1},1))

    p_values_vibr_on = [p_1; p_2; p_3]; 
end
p_values_vibr.p_values_vibr_on = p_values_vibr_on;
save(fullfile(metadata.data_save, 'p_values_inititalvibr.mat'), 'p_values_vibr')
%%



% figure('units','centimeters','Position',[2 2 30 14])
% hold on
% avg_off_vib = cell(3,4);
% avg_on_vib = cell(3,4);
% for con = 1:4%no_con
%     plplpl = [];
%     for group = 1:no_group 
%     % so i I want to scatter the first 10 s after stim
% 
%         subplot(2,2,con)
%         hold on
% 
%         % off_data = squeeze(mean(vibr_change_per_stim{group,1}([8,9,10],vibr_conTrials(con,:),:), 1));
%         % on_data = squeeze(mean(vibr_change_per_stim{group,1}([7],vibr_conTrials(con,:),:), 1)); 
%         off_data = squeeze(mean(vibr_data_per_stim{group,1}(8,vibr_conTrials(con,:),:), 1));
%         on_data = squeeze(mean(vibr_data_per_stim{group,1}(7,vibr_conTrials(con,:),:), 1)); 
%         avg_on_vib{group,con} = on_data;
%         avg_off_vib{group,con} = off_data;
% 
%         a = scatter(ones(size(off_data,2),1)*x_spots(2,group), mean(off_data,1), 'filled', 'MarkerFaceColor', short_cmap(group,:), 'MarkerEdgeColor', short_cmap(group,:))
%         scatter(ones(size(on_data,2),1)*x_spots(1,group), mean(on_data,1), 'filled', 'MarkerFaceColor', short_cmap(group,:), 'MarkerEdgeColor', short_cmap(group,:))
%         plplpl = [plplpl, a];
%         er = errorbar([x_spots([2 1],group)],[mean(mean(off_data,2),1), mean(mean(on_data,2),1)],[squeeze(nanstd(mean(off_data),0,2)/sqrt(size(off_data,2))), squeeze(nanstd(mean(on_data),0,2)/sqrt(size(on_data,2)))]);    
%         er.Color = [0 0 0];                            
%         er.LineStyle = 'none';
%     end
%     xticks([2 5]);
%     xticklabels({'first s', 'sec s'})
%     xlim([0 7])
%     legend(plplpl, metadata.GroupName, 'Location', 'Eastoutside')
%     ylabel(['Avg binned distance  [mm] ' ])
%     title(vibr_conNames{con})
% end
% 
% saveas(gcf, fullfile(metadata.data_save, ['Avg_transition_scatter_vib.png']))
% saveas(gcf, fullfile(metadata.data_save, ['Avg_transition_scatter_vib.svg']))
% saveas(gcf, fullfile(metadata.data_save, ['Avg_transition_scatter_vib.fig']))
% 
% 
% p_vals_per_con = cell(4,2);
% if no_group == 3
%     for con = 1:no_con
%         [p_1, h_1] = quick_statistic(mean(avg_on_vib{1,con},1), mean(avg_on_vib{2,con},1))
%         [p_2, h_2] = quick_statistic(mean(avg_on_vib{1,con},1), mean(avg_on_vib{3,con},1))
%         [p_3, h_3] = quick_statistic(mean(avg_on_vib{2,con},1), mean(avg_on_vib{3,con},1))
%         p_vals_per_con{con,1} = [p_1; p_2; p_3]; 
% 
%         [p_1, h_1] = quick_statistic(mean(avg_off_vib{1,con},1), mean(avg_off_vib{2,con},1))
%         [p_2, h_2] = quick_statistic(mean(avg_off_vib{1,con},1), mean(avg_off_vib{3,con},1))
%         [p_3, h_3] = quick_statistic(mean(avg_off_vib{2,con},1), mean(avg_off_vib{3,con},1))
%         p_vals_per_con{con,2} = [p_1; p_2; p_3];
% 
%     end
% end
p_vals_per_con = 0;
end