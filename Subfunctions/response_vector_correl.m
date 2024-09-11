function [resp_corr_val] = response_vector_correl(all_dfftrialw_wt_FL, brain_regions_wt, position_wt, stim_period, con_trials, resp_list_wt_dff_std, dorsomed, thresh_dormedneurons_wt)


resp_corr_val = nan(size(all_dfftrialw_wt_FL,2),3); 
no_con = size(con_trials,1); 
for fish = 1:size(all_dfftrialw_wt_FL,2)

    brainnumber = 11; 
    current_dff = all_dfftrialw_wt_FL{1,fish};
    % current_dff = all_dffs_wt{1,fish};
    current_hab = current_dff(:,:,find(brain_regions_wt{1,fish} == brainnumber)); 
    % current_hab = current_dff(find(brain_regions_wt{1,fish} == brainnumber),:); 
    hab_positions = position_wt{1,fish}(find(brain_regions_wt{1,fish} == brainnumber),:);
    try
        hab_resp = resp_list_wt_dff_std{1,1}{1,fish}(find(brain_regions_wt{1,fish} == brainnumber),:);
    catch
        hab_resp = resp_list_wt_dff_std{1,1}{1,fish};
    end
    if dorsomed 
        current_hab = current_hab(:,:,find(thresh_dormedneurons_wt{1,fish} ==1));
        hab_positions = hab_positions(find(thresh_dormedneurons_wt{1,fish} ==1),:);
        hab_resp = hab_resp(find(thresh_dormedneurons_wt{1,fish} ==1),:);
    end
    % now I make a little response vector for the different conditions

    resp_vector = nan(size(current_hab,3), no_con);

    for con = 1:no_con
        % resp_vector(:, con) = mean(mean(current_hab(stim_period, con_trials(con,:),:),2),1);
        resp_vector(:, con) = hab_resp(:,con); 
    end
    
    % let us correlated the vector 

    corr_val_LT = corrcoef(resp_vector(:,1),resp_vector(:,2));
    corr_val_LLT = corrcoef(resp_vector(:,1),resp_vector(:,3));
    corr_val_TLT = corrcoef(resp_vector(:,2),resp_vector(:,3));

    resp_corr_val(fish,1) = corr_val_LT(1,2); 
    resp_corr_val(fish,2) = corr_val_LLT(1,2); 
    resp_corr_val(fish,3) = corr_val_TLT(1,2); 
end