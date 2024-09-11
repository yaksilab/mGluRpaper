function AO_combined_traces_plot(all_fish, n_groups, groups, names, y_pos_change_groups, cmap, folder_path_save, figures_subfolder, y_label_name, marker_stim, plotname)
%AO_combined_traces_plot - Plots traces in one plot together
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
%   Date : Nov 2022

% there is an issue with one of my groups that is why i need to cut it down
if size(y_pos_change_groups{1},1) ~= size(y_pos_change_groups{2},1)
    disp('Issue with the size of the sec group')
    y_pos_change_groups{2} = y_pos_change_groups{2}(1:end-1,:,:);
end
stim_types = unique(all_fish{groups{1,1}(1), 1}.stimInfo); 
stim_types = stim_types(2:end);
stim_types_indice = all_fish{groups{1,1}(1), 1}.stimInfo(find(~contains(all_fish{groups{1,1}(1), 1}.stimInfo, "VIB_0"))); % this is now the stimuli times for the specific stimulus 

% getting the x axis in time 
if contains(plotname, 'Velo')
    x_values = 1:size(y_pos_change_groups{1},1);
else
    x_values = 1:size(y_pos_change_groups{1},1); 
    x_values = x_values/15;
end
% getting the name of the different vibrations
labels = all_fish{groups{1,1}(1), 1}.stimInfo(find(~contains(all_fish{groups{1,1}(1), 1}.stimInfo, "VIB_0")));%unique(all_fish{1, 1}.stimInfo);
type_labels = unique(labels);
figure('units','centimeters','Position',[2 2 20 8])
for stim = 1:length(stim_types)
    current_stim_ind = find(stim_types_indice == stim_types(stim)); 
%     subplot(1,3,stim)
    hold on
    all_plots = [];
    for group = 1:n_groups
        H2=shadedErrorBar(x_values, squeeze(nanmean(nanmean(y_pos_change_groups{group}(:,:,current_stim_ind),3),2)),squeeze(nanstd(nanmean(y_pos_change_groups{group}(:,:,current_stim_ind),3),0,2)/sqrt(size(groups{group},1))), 'lineProps','r');
        H2.mainLine.Color = cmap(group,:);
        H2.patch.FaceColor = cmap(group,:); 
        H2.mainLine.LineWidth = 1; 
        all_plots = [all_plots; H2.mainLine];
    end
    if contains('Velo', plotname)
        xline(marker_stim, 'Color', 'b' )
    else
        xline(marker_stim/15, 'Color', 'b' )
    end
%     for stimu = 1:length(current_stim_ind)
%         plot(nanmean(collected_velo(:,:,current_stim_ind(stimu)),2))
%     end
    hold off
%     ylim([-10 6])
%     xlim([marker_stim/15-5 60])
    title([type_labels(stim)])
    ylabel([y_label_name])
    xlabel('Time in s')
    hold off
%     legend([H1.mainLine H2.mainLine H3.mainLine H4.mainLine], {name_group1, name_group3, name_group2, name_group4})
    legend(all_plots, names)

end

saveas(gcf, fullfile(folder_path_save, figures_subfolder, [plotname, '.png'])); 
saveas(gcf, fullfile(folder_path_save, figures_subfolder, [plotname, '.svg'])); 

%% Add the scatter plot
starting = find(x_values > 45); 
starting = starting(1);
% starting = 0
group_oder = [];
combined_data = [];
grouop_avg = []; 
sems = [];
no_group = n_groups; 
for group = 1:no_group
    group_oder = [group_oder; ones(size(y_pos_change_groups{group},2),1)*group];
    cur_avg = mean(mean(y_pos_change_groups{group}(starting:end,:,:),3),1);
    grouop_avg = [grouop_avg; mean(cur_avg)];
    combined_data = [combined_data; cur_avg'];
    sems = [sems; squeeze(nanstd(cur_avg,0,2)/sqrt(size(cur_avg,2)))];

end

figure
x = beeswarm(group_oder,combined_data)
title(['45-60 s (15 s after stim)'])
hold on
er = errorbar([1 2 3],grouop_avg, sems)
er.Color = [0 0 0];                            
er.LineStyle = 'none';

[p1, h] = quick_statistic(combined_data(find(group_oder == 1)), combined_data(find(group_oder == 2)))
[p2, h] = quick_statistic(combined_data(find(group_oder == 1)), combined_data(find(group_oder == 3)))
[p3, h] = quick_statistic(combined_data(find(group_oder == 2)), combined_data(find(group_oder == 3)))
saveas(gcf, fullfile(folder_path_save, figures_subfolder, ['beeeswarm', '.png'])); 
saveas(gcf, fullfile(folder_path_save, figures_subfolder, ['beeeswarm', '.svg'])); 

end