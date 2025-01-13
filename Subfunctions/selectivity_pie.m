function [fish_perc, fish_perc_simple] = selectivity_pie(selective_neurons_hb_wt, name, save_path, cmap1)
%selectivity_pie - Calculate the percentage of selectivity types
%   Author: Anna Maria Ostenrath
%   Syntax:
%       [fish_perc, fish_perc_simple] = selectivity_pie(selective_neurons_hb_wt, name, save_path, cmap1)
%
%   Description:
%       function() - description
%    
%   Inputs:
%       selective_neurons_hb_wt - cell array of all the fish including a
%       cell array of the different types of selectivity (1: only lisght, 2: only vibration, 3: only both)
%       name - string of the group name
%       save_path - string folder path for saving the figure
%       cmap1 - colormap for plotting
%
%   Outputs:
%       fish_perc - array of fish x percentages for each selectivity type
%       fish_perc_simple - array of fish x percentages for each selectivity
%       type where the unimodal (only light and only vibration) are
%       combined
%
%
%   Other m-files required: none
%   Subfunctions: none
%   MAT-files required: none
%
%   See also: OTHER_FUNCTION_NAME1,  OTHER_FUNCTION_NAME2
%   Author: Anna Maria Ostenrath 
%   Date : November 2024	

selective_label = {'Only Light', 'Only Tap', 'Only Both', 'Light and Both', 'Tap and Both' , 'All', 'No Sel'}; 
fish_perc = nan(size(selective_neurons_hb_wt,2), size(selective_neurons_hb_wt{1,1},2));
fish_perc_simple = nan(size(selective_neurons_hb_wt,2), 4);
for fish= 1:size(selective_neurons_hb_wt,2)
    
    
        cur_selectives = selective_neurons_hb_wt{1,fish}();
   
     
    % hab_resp = resp_list_wt_dff_std{1,fish}(find(brain_regions_wt{1,fish} == 11),:)

    total_resp_cells = 0;
    for sel = 1:size(cur_selectives,2)
        
        total_resp_cells = total_resp_cells + size(cur_selectives{sel},1);
        
    end
    perc_response_just_responding_cells = nan(size(cur_selectives,2),1); 
    for sel = 1:size(cur_selectives,2)
        perc_response_just_responding_cells(sel) = size(cur_selectives{sel},1)/total_resp_cells*100;
    end
    % reorganise just for the sake of it beeing pretty

    % figure('units','centimeters','Position',[2 2 20 11])
    % pie(perc_response_just_responding_cells([1 4 2 5 3 6]))
    % lgd = legend(selective_label{[1 4 2 5 3 6]}, 'Location', 'eastoutside');
    % colormap(cmap1([2 5 3 6 4 7],:))
    % title([name, ' Fish ', num2str(fish)])
    % saveas(gcf, fullfile(save_path, 'Indv_Fish', [name, '_Fig20_Fish_' num2str(fish) '_selectivity_spread.png']))
    % saveas(gcf, fullfile(save_path, 'Indv_Fish', [name, '_Fig20_Fish_' num2str(fish)  '_selectivity_spread.svg']))
    % close;
    
    fish_perc(fish,:) = perc_response_just_responding_cells; 

    % now we want the same, but we put 1&2 together, 4&5 and then 3 and
    % then 6
    
    
    perc_response_combined = nan(4,1); 
    changed_sel{1} = [cur_selectives{1}; cur_selectives{2}];
    changed_sel{2} = [cur_selectives{4}; cur_selectives{5}];
    changed_sel{3} = [cur_selectives{3}];
    changed_sel{4} = [cur_selectives{6}];
    perc_response_combined(1) = size(changed_sel{1},1)/total_resp_cells*100;
    perc_response_combined(2) = size(changed_sel{2},1)/total_resp_cells*100;
    perc_response_combined(3) = size(changed_sel{3},1)/total_resp_cells*100;
    perc_response_combined(4) = size(changed_sel{4},1)/total_resp_cells*100;
   
    fish_perc_simple(fish,:) = perc_response_combined;

    
end

err_wt2 = std(fish_perc(:,:), 0 ,1)/sqrt(size(fish_perc,1));
figure('units','centimeters','Position',[2 2 20 11])
hold on
for sel = 1: size(selective_neurons_hb_wt{1,1},2)
    scatter([ones(size(fish_perc,1),1)*sel] , fish_perc(:,sel), 'filled',  'MarkerFaceColor', cmap1(sel+1,:))
end
b4 = bar([1:size(selective_neurons_hb_wt{1,1},2)], mean(fish_perc,1), 0.2,'FaceColor', cmap1(1,:))
set(b4(1), 'facecolor', cmap1(1,:), 'facealpha', 0.3)
% set(b3(2), 'facecolor', map_wt(2,:), 'facealpha', 0.5)
% set(b3(3), 'facecolor', map_wt(3,:), 'facealpha', 0.5)%%b1.CData(1,:) = cmap1(1,:);
er = errorbar([1:size(selective_neurons_hb_wt{1,1},2)],mean(fish_perc,1),-err_wt2,err_wt2);    
er.Color = [0 0 0];                            
er.LineStyle = 'none';  
xticks([1:size(selective_neurons_hb_wt{1,1},2)])
xticklabels(selective_label)
ylabel('% Selectivte')
title([name, ' Sel'])
ylim([0 max(max(fish_perc))+5])
for fish = 1:size(fish_perc,1)
    plot([1:size(selective_neurons_hb_wt{1,1},2)], fish_perc(fish,:), 'Color', 'k')
end
saveas(gcf, fullfile(save_path, [name, '_Fig21_Fish_' num2str(fish) '_selectivity_spread.png']))
saveas(gcf, fullfile(save_path, [name, '_Fig21_Fish_' num2str(fish)  '_selectivity_spread.svg']))
    

end