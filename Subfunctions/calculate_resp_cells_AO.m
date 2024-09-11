function [pos_resp_list, neg_resp_list, perc_resp_cells, perc_neg_resp_cells, info_list] = calculate_resp_cells_AO(trialwise_zscore, method, baseline_dur, base_period_start, base_period_end, stim_period_end, no_con, con_trials, alpha_factor)
%calculate_resp_cells_AO - Calculates the responding cells for variable
%amount of conditions
%   Author: Anna Maria Ostenrath
%   based on having different conditions of stimuli/ different stim types
%   based on dFF or zscore in the format time x trial x cells
%
%   Syntax:
%       [pos_resp_list, neg_resp_list, perc_resp_cells, perc_neg_resp_cells, info_list] = calculate_resp_cells_AO(trialwise_zscore, method, baseline_dur, base_period_end, stim_period_end, no_con, con_trials)
%       [pos_resp_list, neg_resp_list, perc_resp_cells, perc_neg_resp_cells, info_list] = calculate_resp_cells_AO(trialwise_zscore, method, baseline_dur, base_period_end, stim_period_end, no_con, con_trials, alpha_factor)
%       
%
%   Description:
%       calculate_resp_cells_AO() - calculates the responding cells for
%       each condition with two alternative methods 1: signrank and 2: mean
%       plus factor x std (factor and alpha are fixed to 2 and 0.05 if you don't give anything else). For
%       negative respondig cells will always be a 1*std added
%    
%   Inputs:
%       trialwise_zscore - time x trial x cells of dFF or zscore data
%       method - method different methods of calculating resp cells (1 sign rank, 2 std method)
%       baseline_dur - baseline time in frames before the stim hits (depends on how you created the dFF/ zscore matrix or where your stim hits)
%       base_period_start - start of the baseline period in frames that you want to compare with
%       base_period_end - end of the baseline period in frames that you want to compare with
%       stim_period_end - end of the stimulus period in frames (last fram
%       of period of interest, always starts with baseline_dur plus 1)
%       no_con - number of conditions of the experiment
%       con_trials - variable with dim no_con x trials in condition to know
%       which trial belongs to which condition (e.g. 3x6 matrix for 3 con
%       with 6 stim each, number of trials must be qual in each condition) 
%       alpha_factor - optinal argument to change the alpha or the
%       std_factor depending on what you want 
%    
%
%   Outputs:
%       pos_resp_list - cell array of the cell index that positively responds to the
%       conditions 
%       neg_resp_list - cell array of the cell index that negatively responds to the
%       conditions
%       perc_resp_cells - cell array of the calculate percentage for each condition
%       positive
%       perc_neg_resp_cells - cell array of the calculate percentage for each condition
%       ngeative
%       info_list - is a matrix with cell x no_con where pos responding
%       cells are marked with 1, negative with -1 and non-responding with 0
%
%   Examples: 
%       Line 1 of example
%       Line 2 of example
%       Line 3 of example
%
%   Notes: 
%       Improvement: make it variable to add an alpha or std_factor!
%       Line 2 of example
%       Line 3 of example
%
%   Other m-files required: none
%   Subfunctions: none
%   MAT-files required: none
%
%   See also: responding_and_plotting_mglur (old function with some
%   plotting included),  calculate_resp_cells_3_con_AO (optimized for 3 conditions only)
%   Author: Anna Maria Ostenrath 
%   Date : March 2022	

% conditions
% con1_trials = 1:6; 
% con2_trials = 7:12; 
% con3_trials = 13:18; 

pos_list = zeros([size(trialwise_zscore,3), no_con]); 
neg_list = zeros([size(trialwise_zscore,3), no_con]); 

% here we define some parameters
base_period = base_period_start:base_period_end; % what frames do we average over for baseline period%baseline_period   new_base; %
stim_period = baseline_dur+1: stim_period_end; % period of interest
if ~exist('alpha_factor','var')
    alpha = 0.05; % alpha valus for sign rank
    std_factor = 2; % factor of std for mean plus std calculations
    disp('Standart alpha/std factor used.')
else 
    if method == 1
        alpha = alpha_factor; 
        std_factor = 2; % will be unused
        disp(['Method 1 Signrank alpha chosen ' num2str(alpha_factor)])
    elseif method ==2
        std_factor = alpha_factor; 
        alpha = 0.05; % will be unused
        disp(['Method 2 Std factor chosen ' num2str(alpha_factor)])
    end
  
end

for cell = 1:size(trialwise_zscore,3) % now we loop over each cells
      
    if method == 1
        % so here we now use signrank between base_period and stim_period to
        % find pos and neg resp cells
        for con = 1:no_con % loop over each condition
            pre_con1 = mean(trialwise_zscore(base_period,con_trials(con, :),cell),1); 
            post_con1 = mean(trialwise_zscore(stim_period,con_trials(con, :),cell),1); 
            % left tailed for positive and right tailed for neg responding 
            [p_con1_pos,h_con1_pos,stats_con1_pos] = signrank(pre_con1,post_con1,'tail','left','alpha',alpha);
            [p_con1_neg,h_con1_neg,stats_con1_neg] = signrank(pre_con1,post_con1,'tail','right','alpha',alpha);


            % storing the result in a list for each cell... that will later be
            % combined 
            pos_list(cell,con) = h_con1_pos; 
            neg_list(cell,con) = h_con1_neg; 
        end
        
    elseif method == 2
        % so here we now use std of base_period and stim_period to
        % find pos and neg resp cells
        for con = 1:no_con % loop over each condition

            pre_con1 = mean(mean(trialwise_zscore(base_period,con_trials(con, :),cell),2),1); 
            pre_con1_std = std(mean(trialwise_zscore(base_period,con_trials(con, :),cell),2),0,1); 
            post_con1 = mean(mean(trialwise_zscore(stim_period,con_trials(con, :),cell),2),1); 

            % now we compare the mean plus x times std with the post avg to
            % determine pos or neg responding cells 
            if abs(post_con1) > abs(pre_con1)+std_factor*pre_con1_std
               if post_con1 > pre_con1
                   pos_list(cell, con) = 1; 
               elseif post_con1 < pre_con1
                   neg_list(cell,con) = 1; 
               end
           % for neg responding cells we also check with only 1 std
            elseif abs(post_con1) > abs(pre_con1)+1*pre_con1_std
                if post_con1 < pre_con1
                   neg_list(cell,con) = 1; 
               end
            end
        end
    end
end

info_list = zeros(size(pos_list)); % now we combine the pos and neg in one list together

% now we just save the cell indices in case we want to use them for
% somethig
pos_resp_list = {}; 
neg_resp_list = {}; 

% and here we calculate the percentage
perc_resp_cells = {}; 
perc_neg_resp_cells = {}; 
 
for con = 1:no_con
    info_list(find(pos_list(:,con)==1),con) = 1; 
    info_list(find(neg_list(:,con)==1),con) = -1; 
    
    pos_resp_list{con} = find(pos_list(:,con) == 1);
    neg_resp_list{con} = find(neg_list(:,con) == 1);
    
    perc_resp_cells{con} = length(find(pos_list(:,con) == 1))/length(pos_list(:,con)); 
    perc_neg_resp_cells{con} = length(find(neg_list(:,con) == 1))/length(neg_list(:,con)); 
end




end