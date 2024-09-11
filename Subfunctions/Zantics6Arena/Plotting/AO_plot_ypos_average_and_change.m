function [collected_yposition,change_yposition, marker_stim] = AO_plot_ypos_average_and_change(all_fish, duration, group, name, folder_path_save, figures_subfolder, single_ploting, baseline_sec)
%AO_plot_ypos_average_and_change - Plotting and collecting the trialwise
%velocity and the change in velo
%   Author: Anna Maria Ostenrath
%
%   Syntax:
%       [collected_yposition,change_yposition, marker_stim] = AO_plot_ypos_average_and_change(all_fish, duration, group, name, folder_path_save, figures_subfolder, single_ploting)
%    
%   Inputs:
%       all_fish - cell array with all the fish data
%       group - list of indices for the fish in this group
%       name - string of the name of the group
%       duration - pre and post duration of stimulus you want to look at
%       folder_path_save, figures_subfolder - The save path and folder to
%       save as png
%       single_ploting - logical 1 or 0 if you want to plot the fish
%       individually 
%       baseline_sec - the baseline duration in s for calculating the
%       change 
%
%   Outputs:
%       collected_yposition - array for each fish with average y pos per trial in the group with dim time x
%       fish x stim
%       change_yposition - array for each fish with change in y pos per trial in the group with dim time x
%       fish x stim
%       marker_stim - point where the stimulus comes in 
%
%   Notes: 
%       The baseline (for change) is hard coded to 2 s before stim onset (baseline_sec)
%       To DO: maybe add colors and make baseline flexibel
%       Line 3 of example
%
%   Other m-files required: none
%   Subfunctions: none
%   MAT-files required: none
%
%   See also: AO_average_y_trace_flexibel (looking at full trace),  OTHER_FUNCTION_NAME2
%   Author: Anna Maria Ostenrath 
%   Date : May 2022	

try
    mkdir(fullfile(folder_path_save, figures_subfolder, name))
catch
end

stim_times = all_fish{group(1), 1}.timeStampStim(find(~contains(all_fish{group(1), 1}.stimInfo, "VIB_0")));

numExp = size(all_fish,1)/6;
for i=1:numExp
lowerRow (i*3-(3-1):i*3) = [4 5 6]+6*(i-1); % make an array for the lower row of ROIs 
end
t = all_fish{group(1), 1}.t;
%            x = all_fish{fish,1}.x;
y = all_fish{group(1), 1}.y;
len_y_trace = size(-y(t>floor(stim_times(end))-duration & t<floor(stim_times(end))+duration),1)-10;  % double check why you floor it
collected_yposition = nan([len_y_trace, size(group,1), length(stim_times)]);
% collected_yposition_lower = nan([len_y_trace, size(group,1), length(stim_times)]); 
change_yposition = nan([len_y_trace, size(group,1), length(stim_times)]);
% baseline_sec = 2; % in s
for stim = 1: length(stim_times)
   
   for fish = 1:size(group,1)

       if ismember(group(fish),lowerRow)
           curr_stim_times = all_fish{group(fish), 1}.timeStampStim(find(~contains(all_fish{group(fish), 1}.stimInfo, "VIB_0")));
           curr_stim_ons = curr_stim_times(stim); 
           t = all_fish{group(fish),1}.t;
    %            x = all_fish{fish,1}.x;
           y = all_fish{group(fish),1}.y;    
           stim_new_time = find(t>floor(curr_stim_ons));
           onset_period = find(t>floor(curr_stim_ons)-duration); 
           marker_stim = stim_new_time(1) - onset_period(1); 
           
           temp_y = -y(t>floor(curr_stim_ons)-duration & t<floor(curr_stim_ons)+duration)+168; 
           collected_yposition(:, fish, stim) = temp_y(1:len_y_trace);
           change_yposition(:, fish, stim) = abs(temp_y(1:len_y_trace)) -abs(nanmean(temp_y(marker_stim-(baseline_sec*15):marker_stim),1)); %2s average before

       else
           curr_stim_times = all_fish{group(fish), 1}.timeStampStim(find(~contains(all_fish{group(fish), 1}.stimInfo, "VIB_0")));
           curr_stim_ons = curr_stim_times(stim);
           

           
           t = all_fish{group(fish),1}.t;
    %            x = all_fish{fish,1}.x;
           y = all_fish{group(fish),1}.y;
           
           stim_new_time = find(t>floor(curr_stim_ons));
           onset_period = find(t>floor(curr_stim_ons)-duration); 
           marker_stim = stim_new_time(1) - onset_period(1); 
 %                collected_yposition_upper(:, fish, stim) = -y(t>floor(curr_stim_ons)-duration & t<floor(curr_stim_ons)+duration);
           temp_y = -y(t>floor(curr_stim_ons)-duration & t<floor(curr_stim_ons)+duration)+26; 
           collected_yposition(:, fish, stim) = temp_y(1:len_y_trace);
           change_yposition(:, fish, stim) = abs(temp_y(1:len_y_trace)) -abs(nanmean(temp_y(marker_stim-(baseline_sec*15):marker_stim),1)); %2s average before 2s times framerate
       end
    end
   
end
%
% now plotting 

stim_types = unique(all_fish{group(1), 1}.stimInfo); 
stim_types = stim_types(2:end); % we are taking out the VIB 0 one 
stim_types_indice = all_fish{group(1), 1}.stimInfo(find(~contains(all_fish{group(1), 1}.stimInfo, "VIB_0"))); % this is now the stimuli times for the specific stimulus 


% figure('units','centimeters','Position',[2 2 18 24])
figure('units','centimeters','Position',[2 2 12 10])
x_values = 1:len_y_trace; 
labels = all_fish{group(1), 1}.stimInfo(find(~contains(all_fish{group(1), 1}.stimInfo, "VIB_0")));%unique(all_fish{1, 1}.stimInfo);
type_labels = unique(labels);

for stim = 1:length(stim_types)
    current_stim_ind = find(stim_types_indice == stim_types(stim)); 
%     subplot(round(length(stim_types)/2),2,stim)
    subplot(1,length(stim_types),stim)
    H=shadedErrorBar(x_values, squeeze(nanmean(nanmean(collected_yposition(:,:,current_stim_ind),3),2)),squeeze(nanstd(nanmean(collected_yposition(:,:,current_stim_ind),3),0,2)/sqrt(size(group,1))), 'lineProps','k');
%     current_stim_index = find(all_fish{1,1}.t >= stim_times(stim))
    H.mainLine.LineWidth = 1; 
    hold on
    xline(marker_stim, 'Color', 'r' )
    for fish = 1:size(collected_yposition(:,:,stim),2)
        plot(nanmean(collected_yposition(:,fish,current_stim_ind),3))
    end
    hold off
    ylim([-140 0])
    title([name, type_labels(stim)])
end
saveas(gcf, fullfile(folder_path_save, figures_subfolder, name, [name, '_avg_over_stim_fish_yposition.png'])); 
close; 

% now the change
figure('units','centimeters','Position',[2 2 20 8])
change_yposition = -change_yposition; 
for stim = 1:length(stim_types)
    current_stim_ind = find(stim_types_indice == stim_types(stim)); 
%     subplot(round(length(stim_types)/2),2,stim)
    subplot(1,length(stim_types),stim)
    H=shadedErrorBar(x_values, squeeze(nanmean(nanmean(change_yposition(:,:,current_stim_ind),3),2)),squeeze(nanstd(nanmean(change_yposition(:,:,current_stim_ind),3),0,2)/sqrt(size(group,1))), 'lineProps','k');
%     current_stim_index = find(all_fish{1,1}.t >= stim_times(stim))
    H.mainLine.LineWidth = 1; 
    hold on
    xline(marker_stim, 'Color', 'r' )
    for fish = 1:size(change_yposition(:,:,stim),2)
        plot(nanmean(change_yposition(:,fish,current_stim_ind),3))
    end
    hold off
%     ylim([-140 0])
    title([name, type_labels(stim)])
    ylabel('Change of Y Pos')
end
saveas(gcf, fullfile(folder_path_save, figures_subfolder, name, [name, '_change_over_stim_fish_yposition.png'])); 
close; 







if single_ploting
    for fish = 1:size(group,1)
        
            % now plotting the average position for each fish individually but
            % averaged over each stimulus 
            % figure('units','centimeters','Position',[2 2 18 24])
            figure('units','centimeters','Position',[2 2 12 10])
%             figure('units','centimeters','Position',[2 2 20 8])
            x_values = 1:len_y_trace; 
            for stim = 1:length(stim_types) % looping over the types 
                current_stim_ind = find(stim_types_indice == stim_types(stim)); 
            %     subplot(round(length(stim_types)/2),2,stim)
                subplot(1,length(stim_types),stim)
                if ismember(fish,lowerRow)
                    H=shadedErrorBar(x_values, squeeze(nanmean(collected_yposition(:,fish,current_stim_ind),3)),squeeze(nanstd(collected_yposition(:,fish,current_stim_ind),0,3)/sqrt(length(current_stim_ind))), 'lineProps','k');

                %     H=shadedErrorBar(x_values, squeeze(nanmean(collected_yposition_lower(:,:,stim),2)),squeeze(nanstd(collected_yposition_lower(:,:,stim),0,2)/sqrt(size(lowerRow,2))), 'lineProps','k');
                %     current_stim_index = find(all_fish{1,1}.t >= stim_times(stim))
                    hold on
                    xline(marker_stim, 'Color', 'r' )
                    for little_stim = 1:length(current_stim_ind) % looping over the little stim
                        plot(collected_yposition(:,fish,current_stim_ind(little_stim)))
                    end
                    hold off
    %                 ylim([-280 -140])
                    ylim([-140 0])
                    title([name,' Fish ', int2str(fish), ' Lower ', type_labels(stim)])
                else 
                    H=shadedErrorBar(x_values, squeeze(nanmean(collected_yposition(:,fish,current_stim_ind),3)),squeeze(nanstd(collected_yposition(:,fish,current_stim_ind),0,3)/sqrt(length(current_stim_ind))), 'lineProps','k');

                %     H=shadedErrorBar(x_values, squeeze(nanmean(collected_yposition_lower(:,:,stim),2)),squeeze(nanstd(collected_yposition_lower(:,:,stim),0,2)/sqrt(size(lowerRow,2))), 'lineProps','k');
                %     current_stim_index = find(all_fish{1,1}.t >= stim_times(stim))
                    hold on
                    xline(marker_stim, 'Color', 'r' )
                    for little_stim = 1:length(current_stim_ind) % looping over the little stim
                        plot(collected_yposition(:,fish,current_stim_ind(little_stim)))
                    end
                    hold off
                    ylim([-140 0])
                    title([name,' Fish ', int2str(fish), ' Upper ', type_labels(stim)])
                end
            end
            saveas(gcf, fullfile(folder_path_save, figures_subfolder, name, [name, ' Fish_',num2str(all_fish{group(fish), 1}.realNum), '_average_ypos_per_stim.png'])); 
            close; 
        
    end
end
end