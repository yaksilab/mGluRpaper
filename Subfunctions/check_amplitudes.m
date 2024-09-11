function [avg_ampl, collected_ampl] = check_amplitudes(all_dffs_trial, brain_regions_wt, position_wt, ...
    con_names, save_path,name, cmap, thresh_dormedneurons_wt, dorsomed, con_trials, stim_period)

% changed this to not be flexible for brain region anymore
avg_ampl = nan([size(con_trials,1), size(all_dffs_trial,2),1]);
collected_ampl = cell(size(con_trials,1),1); 
for fish = 1:size(all_dffs_trial,2)

    brainnumber = 11; 
    
    current_dff = all_dffs_trial{1,fish};
    % current_dff = all_dffs_wt{1,fish};
    current_hab = current_dff(:,:,find(brain_regions_wt{1,fish} == brainnumber)); 
    % current_hab = current_dff(find(brain_regions_wt{1,fish} == brainnumber),:); 
    
    hab_positions = position_wt{1,fish}(find(brain_regions_wt{1,fish} == brainnumber),:);
    if dorsomed 
        current_hab = current_hab(:,:,find(thresh_dormedneurons_wt{1,fish} ==1));
        hab_positions = hab_positions(find(thresh_dormedneurons_wt{1,fish} ==1),:);
        % hab_resp = hab_resp(find(thresh_dormedneurons_wt{1,fish} ==1),:);
    end
    % hab_resp = resp_list_wt_dff_sign{1,fish}(find(brain_regions_wt{1,fish} == brainnumber),:);
    for con = 1:size(con_trials,1)
        amplitude = nan(size(current_hab,1),1);
        for neuron = 1:size(current_hab,3)
            amplitude(neuron) = squeeze(mean(mean(current_hab(stim_period, con_trials(con,:), neuron),2),1)); 
        end
        collected_ampl{con,1} = [collected_ampl{con,1} ; amplitude];
        avg_ampl(con, fish) = nanmean(amplitude); 
    end
end