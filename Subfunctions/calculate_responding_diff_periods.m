function [resp_list_period, resp_con_list] = calculate_responding_diff_periods(all_dfftrialw_wt_FL, brain_regions_wt, diff_stim_period, con_trials, base_period, no_con, brainnumber, std_factor)
%calculate_responding_diff_periods - Caluclates responding cells with std
%for different periods
%   Author: Anna Maria Ostenrath
%
%   Syntax:
%      [resp_list_period, resp_con_list] = calculate_responding_diff_periods(all_dfftrialw_wt_FL, brain_regions_wt, diff_stim_period, con_trials, base_period, no_con, brainnumber, std_factor)
%       output = function(input1, input2, input3)
%
%   Description:
%       function() - description
%    
%   Inputs:
%       all_dfftrialw_wt_FL - cell array with the dff data (organised: time x trial x neurons) for each fish 
%       brain_regions_wt - cell array with the brain region index for each
%       fish
%       diff_stim_period - cell array with frame indices with the differen
%       periods that should be looked at 
%       con_trials - indicators which trials belong to which condition
%       base_period - duration of the baseline (frames)
%       no_con - number of different conditions (eg. 3)
%       brainnumber - number of the brain region, if 0 then this will look
%       at all the cells 
%       std_factor - how many times the std the positive cells should be
%       over, the standard is 2... inhibited cells are hardcoded at 1
%
%   Outputs:
%       resp_period_list - cell array for each period with each fish with a list of
%       responses for each condition, 1: positive, 0: no response, -1: neg
%       resp_con_list - cell array for each condition with each fish with a list of
%       responses for each condition, 1: positive, 0: no response, -1: neg
%
%
%   Other m-files required: none
%   Subfunctions: none
%   MAT-files required: none
%
%   See also: calculate_resp_cells_AO,  OTHER_FUNCTION_NAME2
%   Author: Anna Maria Ostenrath 
%   Date : November 2024	


resp_list_period = cell(size(diff_stim_period,1),1);
resp_con_list = {};
for fish = 1:size(all_dfftrialw_wt_FL,2)
    
    % cutting the cells down to the brain regions
    current_dff = all_dfftrialw_wt_FL{1,fish};
    % if brainnumber == 0
    %     current_hab = current_dff;
    % else
        current_hab = current_dff(:,:,find(brain_regions_wt{1,fish} == brainnumber)); 
    % end

    % making the empty arrays to fill with responses

    for stim = 1:size(diff_stim_period,1)
        resp_list_period{stim,1}{1,fish} = zeros([size(current_hab,3),no_con]);
    end
    resp_list_con = cell(no_con,1); 
    for con = 1:no_con
        resp_list_con{con,1} = zeros([size(current_hab,3),size(diff_stim_period,1)]);
    end
    for neuron = 1:size(current_hab,3) % looping over neurons

        for stim = 1:size(diff_stim_period,1) % looping over response windows
            stim_period = diff_stim_period(stim,:); 
            for con = 1:no_con % looping over conditions 
                pre_con1 = mean(mean(current_hab(base_period,con_trials(con, :),neuron),2),1); 
                pre_con1_std = std(mean(current_hab(base_period,con_trials(con, :),neuron),2),0,1); 
                post_con1 = mean(mean(current_hab(stim_period,con_trials(con, :),neuron),2),1); 

                % now we compare the mean plus x times std with the post avg to
                % determine pos or neg responding cells 
                if abs(post_con1) > abs(pre_con1)+std_factor*pre_con1_std
                   if post_con1 > pre_con1
                       resp_list_period{stim,1}{1,fish}(neuron, con) = 1; 
                       resp_list_con{con,1}(neuron,stim) = 1; 
                   elseif post_con1 < pre_con1
                       resp_list_period{stim,1}{1,fish}(neuron,con) = -1; 
                       resp_list_con{con,1}(neuron,stim) = -1;
                   end
               % for neg responding cells we also check with only 1 std
                elseif abs(post_con1) > abs(pre_con1)+1*pre_con1_std
                    if post_con1 < pre_con1
                       resp_list_period{stim,1}{1,fish}(neuron,con) = -1; 
                       resp_list_con{con,1}(neuron,stim) = -1;
                   end
                end
            end
        end
    end
%     resp_period_list{fish} = resp_list_period;
    resp_con_list{fish} = resp_list_con;

end


end