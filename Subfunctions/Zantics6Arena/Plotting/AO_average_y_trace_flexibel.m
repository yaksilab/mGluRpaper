function [combined_y_pos, p_values] = AO_average_y_trace_flexibel(all_fish, n_groups, groups, names, cmap, folder_path_save, figures_subfolder, start_point,end_point, seperate_y, stimulus, time_duration)
%AO_average_y_trace_flexibel - Plots the average y traces for each of the
%seperate groups in one plot and then seperately into subplots
%   Author: Anna Maria Ostenrath
%
%   Syntax:
%       combined_y_pos = AO_average_y_trace_flexibel(all_fish, n_groups, groups, names, cmap, folder_path_save, figures_subfolder, start_point,end_point, seperate_y, stim)
%       output = function(input1, input2, input3)
%  
%   Inputs:
%       all_fish - cell array with all the fish data
%       n_groups - number of groups (e.g. 2 for 1 control and 1 treated)
%       groups - cell array with list of fish indices for each group (dim n_groups x 1)
%       names - cell array with names of each group (n_groups x 1 )
%       cmap - colormap with a color (RGB) for each group (n_groups x 3 )
%       folder_path_save, figures_subfolder - The save path and folder to
%       save as png
%       start_point & end_point - start and endpoint of where you want to
%       look in seconds (e.g. novel tank: 1 to 1200 s)
%       seperate_y - logical 1 or 0 if you want to plot the y position seperated into the
%       rows of the arenas
%       stimulus - logical 1 or 0 if you want the stimuli timepoints to be
%       plotted
%       time_duration - Time in s before and after the stimulus or duration
%       of early and late period for NTT
%
%   Outputs:
%       combined_y_pos - cell array (n_groups x 1) with the average y
%       position in the time interval for each group (NaN for fish that don*t belong to that group)
%       p_values - p_values between the groups
%       
%
%   Notes: 
%       The shift of the y pos is hard coded
%       To Do: add plotting for mutlipe control and test groups
%
%   Other m-files required: none
%   Subfunctions: none
%   MAT-files required: none
%
%   See also: AO_noveltank_average_y_trace_groups (made for 4 groups),
%   AO_plot_ypos_average_and_change (plots the average and the change over each stimulus of one conditions)
%   Author: Anna Maria Ostenrath 
%   Date: Nov 2022	

% Here we check which fish was in a lower row because it influences the y
% pos
numExp = size(all_fish,1)/6;
for i=1:numExp
lowerRow (i*3-(3-1):i*3) = [4 5 6]+6*(i-1); % make an array for the lower row of ROIs 
end
t = all_fish{groups{1}(1),1}.t;
x = all_fish{groups{1}(1),1}.x;
y = all_fish{groups{1}(1),1}.y;

start_time = find(t >=start_point);
start_time = start_time(1);
end_time = find(t >=end_point);
end_time = end_time(1);

len_y_trace = size(y(start_time:end_time),1);  % double check why you floor it
collected_y_pos_upper = cell(n_groups,1); 
collected_y_pos_lower = cell(n_groups,1); 
% creating the arrays to store each individual fish in
% it will be nan for the rows of all fish except the ones belonging to the
% group
for group = 1:n_groups
    collected_y_pos_upper{group,1} = nan([len_y_trace, size(all_fish,1)]);
    collected_y_pos_lower{group,1} = nan([len_y_trace, size(all_fish,1)]); 
end
% now looping over all fish %%%%% SHOULD THESE VALUES BE ADAPTED?
for fish = 1:size(all_fish,1)
    for group = 1:n_groups 
       if ismember(fish,groups{group})
           if ismember(fish,lowerRow)
               y = all_fish{fish,1}.y;

               temp_y = -y; % flipping the y position 
               collected_y_pos_lower{group}(:, fish) = temp_y(start_time:end_time)+168; % adding a value to have the "true zero"
           else
               y = all_fish{fish,1}.y;
               temp_y = -y;
               collected_y_pos_upper{group}(:, fish) = temp_y(start_time:end_time)+26; % adding a value to have the "true zero"
           end
       end 
    end
end
% if you want to check if there is a difference between upper and lower row
if seperate_y
    sizes_low = {};
    sizes_up = {};
    for group = 1:n_groups 
        sizes_low{group} = sum(ismember(groups{group},lowerRow));
        sizes_up{group} = length(groups{group})-sizes_low{group}; 

    end


    figure('units','centimeters','Position',[2 2 38 14])
    subplot(2,1,1)
    x_values = t(start_time:end_time); %1:size(collected_yposition_upper,1);
    all_plots = []; 
    hold on 
    %
    for group = 1:n_groups 
        p1 = shadedErrorBar(x_values, squeeze(nanmean(collected_y_pos_upper{group,1},2)),squeeze(nanstd(collected_y_pos_upper{group,1}(:,:),0,2)/sqrt(sizes_up{group})), 'lineProps','k');
        p1.mainLine.Color = cmap(group,:);
        p1.patch.FaceColor = cmap(group,:); 
%         p1.edge.Color = cmap(group,:);
        all_plots = [all_plots, p1]; 

    end
    %%
    if stimulus
        hold on
        for stim = 1:length(stim_times)
            xline(stim_times(stim), 'Color', 'b')
        end
        hold off
    end
    hold off
    title(['Avg Y Upper '])
    % xlim([0 3000])
    ylim([-140 0])
    ylabel('Y position [mm]')
    xlabel('Time [s]')
    legend(names)
    
    subplot(2,1,2)
    all_plots = []; 
    hold on 
    for group = 1:n_groups 
        p1 = shadedErrorBar(x_values, squeeze(nanmean(collected_y_pos_lower{group},2)),squeeze(nanstd(collected_y_pos_lower{group}(:,:),0,2)/sqrt(sizes_low{group})), 'lineProps','k');
        p1.mainLine.Color = cmap(group,:);
        p1.patch.FaceColor = cmap(group,:); 
%         p1.edge.Color = cmap(group,:);
        all_plots = [all_plots, p1]; 

    end
    if stimulus
        hold on
        for stim = 1:length(stim_times)
            xline(stim_times(stim), 'Color', 'b')
        end
        hold off
    end
    hold off
    title(['Avg Y Lower '])
    % xlim([0 3000])
    ylim([-140 0])
    ylabel('Y position [mm]')
    xlabel('Time [s]')
    legend(names)
    saveas(gcf, fullfile(folder_path_save, figures_subfolder, ['Y_pos_comparison_groups_seperated_rows.png']));
    % close; 
end
% now we combine the upper and lower in one
combined_y_pos = {}; 
for group = 1:n_groups 
    combined_y_pos{group,1} = collected_y_pos_upper{group}; 
    combined_y_pos{group,1}(:,lowerRow) = collected_y_pos_lower{group}(:,lowerRow);
end
stim_times = all_fish{groups{1}(1), 1}.timeStampStim(find(~contains(all_fish{groups{1}(1), 1}.stimInfo, "VIB_0")));
%

% plotting the combined ones
figure('units','centimeters','Position',[2 2 30 14])
% subplot(n_groups+1,1,1)
x_values = t(start_time:end_time); %1:size(collected_yposition_upper,1);

hold on
all_plots = [];
for group = 1:n_groups 
    H1=shadedErrorBar(x_values,squeeze(nanmean(combined_y_pos{group},2)),squeeze(nanstd(combined_y_pos{group}(:,:),0,2)/sqrt(size(groups{group},1))), 'lineProps','k');
    H1.mainLine.Color = cmap(group,:);
    H1.patch.FaceColor = cmap(group,:); 
    H1.mainLine.LineWidth = 2; 
    all_plots = [all_plots, H1.mainLine]; 

end

if stimulus
    hold on
    for stim = 1:length(stim_times)
        xline(stim_times(stim), 'Color', 'b')
    end
    hold off
end
hold off
title(['Avg Y'])
% xlim([1300 4200])
ylim([-105 0])
ylabel('Y position [mm]')
xlabel('Time [s]')
legend(all_plots, names) % what about the stim lines all_plots

if stimulus
    saveas(gcf, fullfile(folder_path_save, figures_subfolder, ['Y_pos_comparison_groups_stim_period.png']));
    % close; 
else
%     time_duration = 1*60; % 2 min with 15 frame rate

    early = find(t>20& t<20+floor(time_duration)); % right now the same as the pdfs
    late = find(t>end_point-time_duration & t<=end_point);
    xline(t(early(1)))
    xline(t(early(end)))
    xline(t(late(1)))
    xline(t(late(end)))
    saveas(gcf, fullfile(folder_path_save, figures_subfolder, ['Y_pos_comparison_groups_NTT.png']));
%         close;
end
% now the scatter plots
if stimulus
    disp('Not implemented yet.')
    p_values = [];

else
    % this part is for ploting just the avg y pos early and late
%     time_duration = 1*60; % 30; %  2 min with 15 frame rate

%     early = find(t>60& t<60+floor(time_duration)); 
%     late = find(t>end_point-time_duration & t<=end_point); 
    early = find(t>20& t<20+floor(time_duration)); 
    late = find(t>end_point-time_duration & t<=end_point); 
    % now lets plot the avg position for each fish at that point
    avg_early_1 = nanmean(combined_y_pos{1}(early,:),1);
    avg_early_2 = nanmean(combined_y_pos{2}(early,:),1);
    avg_late_1 = nanmean(combined_y_pos{1}(late,:),1);
    avg_late_2 = nanmean(combined_y_pos{2}(late,:),1);


    sem_Early1 = squeeze(nanstd(avg_early_1,0,2)/sqrt(size(groups{1},1)));
    sem_Early2 = squeeze(nanstd(avg_early_2,0,2)/sqrt(size(groups{2},1)));
    sem_Late1 = squeeze(nanstd(avg_late_1,0,2)/sqrt(size(groups{1},1)));
    sem_Late2 = squeeze(nanstd(avg_late_2,0,2)/sqrt(size(groups{2},1)));
    
    if n_groups == 3
        avg_early_3 = nanmean(combined_y_pos{3}(early,:),1);
        avg_late_3 = nanmean(combined_y_pos{3}(late,:),1);
        sem_Early3 = squeeze(nanstd(avg_early_1,0,2)/sqrt(size(groups{1},1)));
        sem_Late3 = squeeze(nanstd(avg_late_2,0,2)/sqrt(size(groups{2},1)));
    end
    if n_groups == 2
        figure
        a = scatter(ones(length(avg_early_1),1)*0.5, avg_early_1, 'filled', 'MarkerEdgeColor', cmap(1,:), 'MarkerFaceColor', cmap(1,:)); 
        hold on
        b = scatter(ones(length(avg_early_2),1)*1, avg_early_2, 'filled',  'MarkerEdgeColor', cmap(2,:), 'MarkerFaceColor', cmap(2,:));
        scatter(ones(length(avg_late_1),1)*2.5, avg_late_1, 'filled',  'MarkerEdgeColor', cmap(1,:), 'MarkerFaceColor', cmap(1,:));
        scatter(ones(length(avg_late_2),1)*3, avg_late_2, 'filled',  'MarkerEdgeColor', cmap(2,:), 'MarkerFaceColor', cmap(2,:));
            er = errorbar([0.5 1 2.5 3],[nanmean(avg_early_1), nanmean(avg_early_2), nanmean(avg_late_1), nanmean(avg_late_2)],[sem_Early1, sem_Early2,sem_Late1,sem_Late2],[-sem_Early1, -sem_Early2, -sem_Late1,-sem_Late2]);    
        er.Color = [0 0 0];                            
        er.LineStyle = 'none';  
        xlim([0 4])
        ylim([-105 0])
        xticks([1 3])

        xticklabels({'Early','Late'})
        ylabel('Y pos mm');


        [p_early1v2, h_early1v2] = quick_statistic(avg_early_1(:,groups{1,1}), avg_early_2(:,groups{2,1}));
        [p_late1v2, h_late1v2] = quick_statistic(avg_late_1(:,groups{1,1}), avg_late_2(:,groups{2,1}));
            
        p_values = {{p_early1v2}, {p_late1v2}}; % saving the ranksum p-values
        plot([0.5, 1], [nanmean(avg_early_1), nanmean(avg_early_2)], 'Color', 'k')
        plot([2.5, 3], [nanmean(avg_late_1), nanmean(avg_late_2)], 'Color', 'k')
        legend([a b], names)
        
        hold off
        title('Avg Swimming Depth')
    
    elseif n_groups == 3
        figure
        a = scatter(ones(length(avg_early_1),1)*0.5, avg_early_1, 'filled', 'MarkerEdgeColor', cmap(1,:), 'MarkerFaceColor', cmap(1,:)); 
        hold on
        b = scatter(ones(length(avg_early_2),1)*1, avg_early_2, 'filled',  'MarkerEdgeColor', cmap(2,:), 'MarkerFaceColor', cmap(2,:));
        c = scatter(ones(length(avg_early_3),1)*1.5, avg_early_3, 'filled',  'MarkerEdgeColor', cmap(3,:), 'MarkerFaceColor', cmap(3,:));

        scatter(ones(length(avg_late_1),1)*2.5, avg_late_1, 'filled',  'MarkerEdgeColor', cmap(1,:), 'MarkerFaceColor', cmap(1,:));
        scatter(ones(length(avg_late_2),1)*3, avg_late_2, 'filled',  'MarkerEdgeColor', cmap(2,:), 'MarkerFaceColor', cmap(2,:));
        scatter(ones(length(avg_late_3),1)*3.5, avg_late_3, 'filled',  'MarkerEdgeColor', cmap(3,:), 'MarkerFaceColor', cmap(3,:));

        er = errorbar([0.5 1 1.5 2.5 3 3.5],[nanmean(avg_early_1), nanmean(avg_early_2), nanmean(avg_early_3), nanmean(avg_late_1), nanmean(avg_late_2), nanmean(avg_late_3)],[sem_Early1, sem_Early2, sem_Early3, sem_Late1,sem_Late2, sem_Late3],[-sem_Early1, -sem_Early2,-sem_Early3, -sem_Late1,-sem_Late2, -sem_Late3]);    
        er.Color = [0 0 0];                            
        er.LineStyle = 'none';  
        xlim([0 4])
        ylim([-105 0])
        xticks([1 3])

        xticklabels({'Early','Late'})
        ylabel('Y pos mm');


        [p_early1v2, h_early1v2] = quick_statistic(avg_early_1(:,groups{1,1}), avg_early_2(:,groups{2,1}));
        [p_late1v2, h_late1v2] = quick_statistic(avg_late_1(:,groups{1,1}), avg_late_2(:,groups{2,1}));
            
        [p_early1v3, h_early1v3] = quick_statistic(avg_early_1(:,groups{1,1}), avg_early_3(:,groups{3,1}));
        [p_late1v3, h_late1v3] = quick_statistic(avg_late_1(:,groups{1,1}), avg_late_3(:,groups{3,1}));
        [p_early2v3, h_early2v3] = quick_statistic(avg_early_2(:,groups{2,1}), avg_early_3(:,groups{3,1}));
        [p_late2v3, h_late2v3] = quick_statistic(avg_late_2(:,groups{2,1}), avg_late_3(:,groups{3,1}));

        p_values = {{p_early1v2; p_early1v3; p_early2v3}, {p_late1v2; p_late1v3; p_late2v3}}; % saving the ranksum p-values
        
        plot([0.5, 1, 1.5], [nanmean(avg_early_1), nanmean(avg_early_2),nanmean(avg_early_3)], 'Color', 'k')
        plot([2.5, 3, 3.5], [nanmean(avg_late_1), nanmean(avg_late_2),nanmean(avg_late_3)], 'Color', 'k')
        legend([a b c], names, 'Location', 'Northwest')

        hold off
        title('Avg Swimming Depth')

    end

    saveas(gcf, fullfile(folder_path_save, figures_subfolder, ['Y_pos_scatter_plot_NTT.png']));
end
end
