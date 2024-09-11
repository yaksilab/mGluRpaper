function [perc_resp_wt, perc_resp_het, perc_resp_hom] = plot_perc_resp_cells_AO(resp_list_wt_dff_std, resp_list_het_dff_std, resp_list_hom_dff_std, no_con, map_wt, map_het,map_hom, response_type, con_names)
%plot_perc_resp_cells - One line description of what the function or script performs (H1 line)
%   Author: Anna Maria Ostenrath
%   Optional file header info (to give more details about the function than in the H1 line)
%   Optional file header info (to give more details about the function than in the H1 line)
%
%   Syntax:
%       output = function(input1, input2)
%       output = function(input1, input2, input3)
%
%   Description:
%       plot_perc_resp_cells() - description
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
%   Notes: 
%       Make this variable for groups and conditions
%       Line 2 of example
%       Line 3 of example
%
%   Other m-files required: none
%   Subfunctions: none
%   MAT-files required: none
%
%   See also: OTHER_FUNCTION_NAME1,  OTHER_FUNCTION_NAME2

%   Author: Anna Maria Ostenrath 
%   Date : April 2022	


perc_resp_wt = nan(size(resp_list_wt_dff_std, 2), no_con); 
for fish = 1:size(resp_list_wt_dff_std, 2)
    for con = 1:no_con
        perc_resp_wt(fish,con) = length(find(resp_list_wt_dff_std{fish}(:,con) == response_type))/length(resp_list_wt_dff_std{fish}(:,con)) * 100; 
    end    
end

perc_resp_het = nan(size(resp_list_het_dff_std, 2), no_con); 
for fish = 1:size(resp_list_het_dff_std, 2)
    for con = 1:no_con
        perc_resp_het(fish,con) = length(find(resp_list_het_dff_std{fish}(:,con) == response_type))/length(resp_list_het_dff_std{fish}(:,con)) * 100; 
    end    
end

perc_resp_hom = nan(size(resp_list_hom_dff_std, 2), no_con); 
for fish = 1:size(resp_list_hom_dff_std, 2)
    for con = 1:no_con
        perc_resp_hom(fish,con) = length(find(resp_list_hom_dff_std{fish}(:,con) == response_type))/length(resp_list_hom_dff_std{fish}(:,con)) * 100; 
    end    
end

% low_x = [0.75 1 1.25] ; 
% mod_x = [1.75 2 2.25];
% high_x = [2.75 3 3.25]; 
% 
% err_wt = std(perc_resp_wt(:,:), 0 ,1)/sqrt(size(perc_resp_wt,1));
% err_het = std(perc_resp_het(:,:), 0 ,1)/sqrt(size(perc_resp_het,1));
% err_hom = std(perc_resp_hom(:,:), 0 ,1)/sqrt(size(perc_resp_hom,1));
% 
% % err_wt_mod = std(responding_wt.mod, 0 ,1)/sqrt(length(responding_wt.mod));
% % err_wt_high = std(responding_wt.high, 0 ,1)/sqrt(length(responding_wt.high));
% 
% 
% figure('units','centimeters','Position',[2 2 18 12])
% hold on
% 
% scatter([ones(size(perc_resp_wt,1),1)*low_x(3); ones(size(perc_resp_wt,1),1)*mod_x(3); ones(size(perc_resp_wt,1),1)*high_x(3)] , [perc_resp_wt(:,1); perc_resp_wt(:,2); perc_resp_wt(:,3)], 'filled',  'MarkerFaceColor', map_wt(3,:))
% b3 = bar([low_x(3), mod_x(3), high_x(3)], [mean(perc_resp_wt(:,1)), mean(perc_resp_wt(:,2)), mean(perc_resp_wt(:,3))], 0.2,'FaceColor', map_wt(3,:))
% set(b3(1), 'facecolor', map_wt(1,:), 'facealpha', 0.5)
% % set(b3(2), 'facecolor', map_wt(2,:), 'facealpha', 0.5)
% % set(b3(3), 'facecolor', map_wt(3,:), 'facealpha', 0.5)%%b1.CData(1,:) = cmap1(1,:);
% er = errorbar([low_x(3), mod_x(3), high_x(3)],[mean(perc_resp_wt(:,1)), mean(perc_resp_wt(:,2)), mean(perc_resp_wt(:,3))],-err_wt,err_wt);    
% er.Color = [0 0 0];                            
% er.LineStyle = 'none';  
% 
% scatter([ones(size(perc_resp_het,1),1)*low_x(2); ones(size(perc_resp_het,1),1)*mod_x(2); ones(size(perc_resp_het,1),1)*high_x(2)] , [perc_resp_het(:,1); perc_resp_het(:,2); perc_resp_het(:,3)], 'filled',  'MarkerFaceColor', map_het(3,:))
% b2 = bar([low_x(2), mod_x(2), high_x(2)], [mean(perc_resp_het(:,1)), mean(perc_resp_het(:,2)), mean(perc_resp_het(:,3))], 0.2,'FaceColor', map_wt(3,:))
% set(b2(1), 'facecolor', map_het(1,:), 'facealpha', 0.5)
% % set(b3(2), 'facecolor', map_wt(2,:), 'facealpha', 0.5)
% % set(b3(3), 'facecolor', map_wt(3,:), 'facealpha', 0.5)%%b1.CData(1,:) = cmap1(1,:);
% er = errorbar([low_x(2), mod_x(2), high_x(2)],[mean(perc_resp_het(:,1)), mean(perc_resp_het(:,2)), mean(perc_resp_het(:,3))],-err_het,err_het);    
% er.Color = [0 0 0];                            
% er.LineStyle = 'none';  
% 
% scatter([ones(size(perc_resp_hom,1),1)*low_x(1); ones(size(perc_resp_hom,1),1)*mod_x(1); ones(size(perc_resp_hom,1),1)*high_x(1)] , [perc_resp_hom(:,1); perc_resp_hom(:,2); perc_resp_hom(:,3)], 'filled',  'MarkerFaceColor', map_hom(3,:))
% b1 = bar([low_x(1), mod_x(1), high_x(1)], [mean(perc_resp_hom(:,1)), mean(perc_resp_hom(:,2)), mean(perc_resp_hom(:,3))], 0.2,'FaceColor', map_wt(3,:))
% set(b1(1), 'facecolor', map_hom(1,:), 'facealpha', 0.5)
% % set(b3(2), 'facecolor', map_wt(2,:), 'facealpha', 0.5)
% % set(b3(3), 'facecolor', map_wt(3,:), 'facealpha', 0.5)%%b1.CData(1,:) = cmap1(1,:);
% er = errorbar([low_x(1), mod_x(1), high_x(1)],[mean(perc_resp_hom(:,1)), mean(perc_resp_hom(:,2)), mean(perc_resp_hom(:,3))],-err_hom,err_hom);    
% er.Color = [0 0 0];                            
% er.LineStyle = 'none';  
% 
% xticks([1 2 3])
% xticklabels(con_names)
% legend([b1 b2 b3], {'Hom', 'Het', 'Wt'})
% ylabel('% Resp Cells')
% title('Comparison all cells')
% ylim([0 50])

end