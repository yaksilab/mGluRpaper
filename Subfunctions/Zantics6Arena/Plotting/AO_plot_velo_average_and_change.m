function [collected_velo,change_velo, collected_stim] = AO_plot_velo_average_and_change(all_fish, group,name, duration, timebin, folder_path_save, figures_subfolder,single_ploting, baseline_sec)
%AO_plot_velo_average_and_change - Plotting and collecting the trialwise
%velocity and the change in velo
%   Author: Anna Maria Ostenrath
%
%   Syntax:
%       [collected_velo,change_velo, collected_stim] = AO_plot_velo_average_and_change(all_fish, group,name, duration, timebin, folder_path_save, figures_subfolder,single_ploting)
%    
%   Inputs:
%       all_fish - cell array with all the fish data
%       group - list of indices for the fish in this group
%       name - string of the name of the group
%       duration - pre and post duration of stimulus you want to look at
%       timebin - timebin (double) for the velocity data (e.g. 0.5)
%       folder_path_save, figures_subfolder - The save path and folder to
%       save as png
%       single_ploting - logical 1 or 0 if you want to plot the fish
%       individually
%
%   Outputs:
%       collected_velo - array for each fish with average velo per trial in the group with dim time x
%       fish x stim
%       change_velo - array for each fish with change in velo per trial in the group with dim time x
%       fish x stim
%       collected_stim - stimulus time points for each trial (example of the first fish in group)
%
%   Notes: 
%       The baseline is hard coded to 2 s before stim onset (baseline_sec)
%       Line 2 of example
%       Line 3 of example
%
%   Other m-files required: none
%   Subfunctions: none
%   MAT-files required: none
%
%   See also: AO_average_velo_trace_felxibel (looking at full trace),  OTHER_FUNCTION_NAME2
%   Author: Anna Maria Ostenrath 
%   Date : May 2022	

try
    mkdir(fullfile(folder_path_save, figures_subfolder, name))
catch
end

stim_times = all_fish{group(1), 1}.timeStampStim(find(~contains(all_fish{group(1), 1}.stimInfo, "VIB_0")));

n_bins = floor(max(all_fish{group(1),1}.t)/timebin);
new_time = linspace(0, max(all_fish{group(1),1}.t), n_bins);
if timebin < 1
    velo_string = num2str(timebin); 
    velo_string(2) = '_';
    fish_speed = all_fish{group(1),1}.(['binnedVel_' velo_string]);
else
    fish_speed = all_fish{group(1),1}.(['binnedVel_' num2str(timebin)]);
end
len_velo_trace = length(fish_speed(new_time>floor(stim_times(1))-duration & new_time<floor(stim_times(1))+duration));  % double check why you floor it
collected_velo = nan([len_velo_trace, size(group,1), length(stim_times)]);
collected_stim = [];     

change_velo = nan(size(collected_velo)); 
% baseline_sec = 2; % the baseline duration in seconds
for stim = 1: length(stim_times)
 %  
   for fish = 1:size(group,1)

        if timebin < 1
            velo_string = num2str(timebin); 
            velo_string(2) = '_';
            fish_speed = all_fish{group(fish),1}.(['binnedVel_' velo_string]);
        else
            fish_speed = all_fish{group(fish),1}.(['binnedVel_' num2str(timebin)]);
        end
        n_bins = floor(max(all_fish{group(fish),1}.t)/timebin);
        new_time = linspace(0, max(all_fish{group(fish),1}.t), n_bins); 

       curr_stim_times = all_fish{group(fish), 1}.timeStampStim(find(~contains(all_fish{group(fish), 1}.stimInfo, "VIB_0")));
       curr_stim_ons = curr_stim_times(stim); 

       temp_speed = fish_speed(new_time>floor(curr_stim_ons)-duration & new_time<floor(curr_stim_ons)+duration); 
       if size(temp_speed,2) < len_velo_trace
           temp_speed = fish_speed(new_time>floor(curr_stim_ons)-duration & new_time<floor(curr_stim_ons)+duration+1);
       end
       
       stim_new_time = find(new_time>floor(curr_stim_ons));
       onset_period = find(new_time>floor(curr_stim_ons)-duration); 
       marker_stim = stim_new_time(1) - onset_period(1); 
      if group(fish) == group(1)
           collected_stim = [collected_stim, marker_stim];  
       end

       collected_velo(:, fish, stim) = temp_speed(1:len_velo_trace);
       
       change_velo(:, fish, stim) = temp_speed(1:len_velo_trace)-nanmean(temp_speed(marker_stim-round(baseline_sec/timebin):marker_stim),2); 
   end
   
end



%%
%% PLotting       
% figure('units','centimeters','Position',[2 2 18 24])
% x_values = 1:len_velo_trace; 
labels = all_fish{group(1), 1}.stimInfo(find(~contains(all_fish{group(1), 1}.stimInfo, "VIB_0")));%unique(all_fish{1, 1}.stimInfo);
% 
% for stim = 1:length(stim_times)
%     subplot(round(length(stim_times)/2),2,stim)
%     
%     H=shadedErrorBar(x_values, squeeze(nanmean(collected_velo(:,:,stim),2)),squeeze(nanstd(collected_velo(:,:,stim),0,2)/sqrt(size(group,1))), 'lineProps','k');
% %     current_stim_index = find(all_fish{1,1}.t >= stim_times(stim))
%     hold on
%     xline(collected_stim(stim), 'Color', 'r' )
%     for fish = 1:size(collected_velo(:,:,stim),2)
%         plot(collected_velo(:,fish,stim))
%     end
%     hold off
%     ylim([0 10])
%     title([labels(stim)])
%     ylabel([int2str(timebin) ' s binned velo ', name])
% 
% end
stim_types = unique(all_fish{group(1), 1}.stimInfo); 
stim_types = stim_types(2:end); % we are taking out the VIB 0 one 
stim_types_indice = all_fish{group(1), 1}.stimInfo(find(~contains(all_fish{group(1), 1}.stimInfo, "VIB_0"))); % this is now the stimuli times for the specific stimulus 

% figure('units','centimeters','Position',[2 2 18 24])
figure('units','centimeters','Position',[2 2 20 8])
x_values = 1:len_velo_trace; 
type_labels = unique(labels);
%
for stim = 1:length(stim_types)
    current_stim_ind = find(stim_types_indice == stim_types(stim)); 
    subplot(1,length(stim_types),stim)
    H=shadedErrorBar(x_values, squeeze(nanmean(nanmean(collected_velo(:,:,current_stim_ind),3),2)),squeeze(nanstd(nanmean(collected_velo(:,:,current_stim_ind),3),0,2)/sqrt(size(group,1))), 'lineProps','k');
%     current_stim_index = find(all_fish{1,1}.t >= stim_times(stim))
    hold on
    xline(collected_stim(stim), 'Color', 'r' )
%     for stimu = 1:length(current_stim_ind)
%         plot(nanmean(collected_velo(:,:,current_stim_ind(stimu)),2))
%     end
    hold off
    ylim([0 5])
    title([type_labels(stim) , name])
    ylabel([num2str(timebin) ' s binned velo ', name])
    hold off

end
saveas(gcf, fullfile(folder_path_save, figures_subfolder, name, ['Average_velo_per_stim_avg_fish.png'])); 
close; 
% now ploting the change 
figure('units','centimeters','Position',[2 2 20 8])
x_values = 1:len_velo_trace; 
type_labels = unique(labels);
%
for stim = 1:length(stim_types)
    current_stim_ind = find(stim_types_indice == stim_types(stim)); 
    subplot(1,length(stim_types),stim)
    H=shadedErrorBar(x_values, squeeze(nanmean(nanmean(change_velo(:,:,current_stim_ind),3),2)),squeeze(nanstd(nanmean(change_velo(:,:,current_stim_ind),3),0,2)/sqrt(size(group,1))), 'lineProps','k');
%     current_stim_index = find(all_fish{1,1}.t >= stim_times(stim))
    hold on
    xline(collected_stim(stim), 'Color', 'r' )
%     for stimu = 1:length(current_stim_ind)
%         plot(nanmean(collected_velo(:,:,current_stim_ind(stimu)),2))
%     end
    hold off
    ylim([0 3])
    title([type_labels(stim) , name])
    ylabel([num2str(timebin) ' s binned velo change'])
    hold off

end
saveas(gcf, fullfile(folder_path_save, figures_subfolder, name, ['Average_velo_change_per_stim_avg_fish.png'])); 
close;

if single_ploting
    for fish = 1:size(group,1)
        
            figure('units','centimeters','Position',[2 2 20 8])
            x_values = 1:len_velo_trace; 
            type_labels = labels(2:end);
            %
            for stim = 1:length(stim_types)
                current_stim_ind = find(stim_types_indice == stim_types(stim)); 
                subplot(1,length(stim_types),stim)
                H=shadedErrorBar(x_values, squeeze(nanmean(collected_velo(:,fish,current_stim_ind),3)),squeeze(nanstd(collected_velo(:,fish,current_stim_ind),0,3)/sqrt(length(current_stim_ind))), 'lineProps','k');
            %     current_stim_index = find(all_fish{1,1}.t >= stim_times(stim))
                hold on
                xline(collected_stim(stim), 'Color', 'r' )
        %         for stimu = 1:length(current_stim_ind)
        %             plot(collected_velo(:,fish,current_stim_ind(stimu)))
        %         end
                hold off
                ylim([0 5])
                title(['Fish ', int2str(fish),' ',type_labels(stim)])
                ylabel([int2str(timebin) ' s binned velo'])
                hold off

            end
            saveas(gcf, fullfile(folder_path_save, figures_subfolder, name, ['Fish_',int2str(fish), '_average_velo_per_stim.png'])); 
            close; 
        
    end
end
end