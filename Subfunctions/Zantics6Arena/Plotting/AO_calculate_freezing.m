function AO_calculate_freezing(all_fish, groups, n_groups, group_name, stim_times, save_path, start_points,end_points, single_plot, save_name, plotting)
%name - One line description of what the function or script performs (H1 line)
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
%   Other m-files required: me_freezing_fromFP
%   Subfunctions: none
%   MAT-files required: none
%
%   See also: OTHER_FUNCTION_NAME1,  OTHER_FUNCTION_NAME2
%   Author: Anna Maria Ostenrath 
%   Date : June 2023	


%% creating the array to store the individual freezing times in there... 
freezing_perc_list = nan(size(all_fish,1),1); 
t = all_fish{groups{1}(1),1}.t;
start_time = find(t >=start_point);
start_time = start_time(1);
end_time = find(t >=end_point);
end_time = end_time(1);
% now looping over all fish for the current time period
for fish = 1:size(all_fish,1)
    for group = 1:n_groups 
       if ismember(fish,groups{group})
        D = all_fish{fish, 1}.dt; %all_fish{20, 1}.calcBinnedDistance_1; %all_fish{20, 1}.distance % 
        X = all_fish{fish, 1}.x; 
        Y = all_fish{fish, 1}.y; 
        
        [freezing_percentage,freezing_time]=me_freezing_fromFP(X(start_time:end_time),Y(start_time:end_time),D(start_time:end_time),[2 5]); 
       freezing_perc_list(fish) = freezing_percentage; 
       end
       
    end
end

freeze_perc = {};
for group = 1:n_groups
    freeze_perc{group} = freezing_perc_list(groups{group});
    disp(mean(freezing_perc_list(groups{group})))
end
figure
plots = [];
for group = 1:n_groups
    hold on
    s = scatter(ones(size(groups{group},1),1)*group, freeze_perc{group}, 50, 'filled', 'MarkerEdgeColor', cmap(group,:), 'MarkerFaceColor', cmap(group,:));
    plots = [plots s];
end

xlim([0 4])
SEM1 = std(freeze_perc{1}(:,:), 0 ,1)/sqrt(size(freeze_perc{1},1));
SEM2 = std(freeze_perc{2}(:,:), 0 ,1)/sqrt(size(freeze_perc{2},1));
SEM3 = std(freeze_perc{3}(:,:), 0 ,1)/sqrt(size(freeze_perc{3},1));
er = errorbar([1 2 3],[nanmean(freeze_perc{1}), nanmean(freeze_perc{2}), nanmean(freeze_perc{3}(:,:))],[SEM1, SEM2, SEM3]);    
er.Color = [0 0 0];                            
er.LineStyle = 'none'; 
legend(plots, group_name)
title(['Freezing perc'])
xticks([1 2 3])
xticklabels(group_name)

%% Now we do it for after each stimulus 
% n_trials = size(stim_times,1); 
% afterstim = stim_times(3) - stim_times(2); 
% freeze_trial = {};
% for group = 1:n_groups
%     perc_per_trial = nan(size(groups{group},1),n_trials);
%     for fish = 1:size(groups{group},1)
%         for stim = 1:size(stim_times,1)
%             cu_fish = groups{group}(fish); 
%             D = all_fish{cu_fish, 1}.dt; %all_fish{20, 1}.calcBinnedDistance_1; %all_fish{20, 1}.distance % 
%             X = all_fish{cu_fish, 1}.x; 
%             Y = all_fish{cu_fish, 1}.y;
% 
%             start_time = stim_times(stim);
%             end_time = stim_times(stim)+afterstim-5;
%             try
%                 [freezing_percentage,freezing_time]=me_freezing_fromFP(X(start_time:end_time),Y(start_time:end_time),D(start_time:end_time),[2 2]); 
%                 perc_per_trial(fish, stim) = freezing_percentage;
%             catch
%                 perc_per_trial(fish, stim) = 0;
%             end
%         end
% 
%     end
%     freeze_trial{group} = perc_per_trial; 
% end

end