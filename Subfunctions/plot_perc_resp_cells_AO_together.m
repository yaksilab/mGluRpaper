function plot_perc_resp_cells_AO_together(pos_perc_resp_wt, pos_perc_resp_het, pos_perc_resp_hom,neg_perc_resp_wt, neg_perc_resp_het, neg_perc_resp_hom, map_wt, map_het,map_hom, con_names, region_name, save_path, group_names)
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


% low_x = [0.75 1 1.25] ; 
% mod_x = [1.75 2 2.25];
high_x = [2.75 3 3.25]; 

low_x = [ 1.25 1 0.75 ] ; 
mod_x = [2.25 2 1.75];
high_x = [3.25 3 2.75]; 

err_wt = std(pos_perc_resp_wt(:,:), 0 ,1)/sqrt(size(pos_perc_resp_wt,1));
err_het = std(pos_perc_resp_het(:,:), 0 ,1)/sqrt(size(pos_perc_resp_het,1));
err_hom = std(pos_perc_resp_hom(:,:), 0 ,1)/sqrt(size(pos_perc_resp_hom,1));


err_wt2 = std(neg_perc_resp_wt(:,:), 0 ,1)/sqrt(size(neg_perc_resp_wt,1));
err_het2 = std(neg_perc_resp_het(:,:), 0 ,1)/sqrt(size(neg_perc_resp_het,1));
err_hom2 = std(neg_perc_resp_hom(:,:), 0 ,1)/sqrt(size(neg_perc_resp_hom,1));
% err_wt_mod = std(responding_wt.mod, 0 ,1)/sqrt(length(responding_wt.mod));
% err_wt_high = std(responding_wt.high, 0 ,1)/sqrt(length(responding_wt.high));


figure('units','centimeters','Position',[2 2 14 12])
subplot(2,1,1)
hold on

scatter([ones(size(pos_perc_resp_wt,1),1)*low_x(3); ones(size(pos_perc_resp_wt,1),1)*mod_x(3); ones(size(pos_perc_resp_wt,1),1)*high_x(3)] , [pos_perc_resp_wt(:,1); pos_perc_resp_wt(:,2); pos_perc_resp_wt(:,3)], 'filled',  'MarkerFaceColor', map_wt(2,:))
b3 = bar([low_x(3), mod_x(3), high_x(3)], [nanmean(pos_perc_resp_wt(:,1)), nanmean(pos_perc_resp_wt(:,2)), nanmean(pos_perc_resp_wt(:,3))], 0.2,'FaceColor', map_wt(3,:))
set(b3(1), 'facecolor', map_wt(1,:), 'facealpha', 0.5)
% set(b3(2), 'facecolor', map_wt(2,:), 'facealpha', 0.5)
% set(b3(3), 'facecolor', map_wt(3,:), 'facealpha', 0.5)%%b1.CData(1,:) = cmap1(1,:);
er = errorbar([low_x(3), mod_x(3), high_x(3)],[nanmean(pos_perc_resp_wt(:,1)), nanmean(pos_perc_resp_wt(:,2)), nanmean(pos_perc_resp_wt(:,3))],-err_wt,err_wt);    
er.Color = [0 0 0];                            
er.LineStyle = 'none';  

scatter([ones(size(pos_perc_resp_het,1),1)*low_x(2); ones(size(pos_perc_resp_het,1),1)*mod_x(2); ones(size(pos_perc_resp_het,1),1)*high_x(2)] , [pos_perc_resp_het(:,1); pos_perc_resp_het(:,2); pos_perc_resp_het(:,3)], 'filled',  'MarkerFaceColor', map_het(2,:))
b2 = bar([low_x(2), mod_x(2), high_x(2)], [nanmean(pos_perc_resp_het(:,1)), nanmean(pos_perc_resp_het(:,2)), nanmean(pos_perc_resp_het(:,3))], 0.2,'FaceColor', map_wt(3,:))
set(b2(1), 'facecolor', map_het(1,:), 'facealpha', 0.5)
% set(b3(2), 'facecolor', map_wt(2,:), 'facealpha', 0.5)
% set(b3(3), 'facecolor', map_wt(3,:), 'facealpha', 0.5)%%b1.CData(1,:) = cmap1(1,:);
er = errorbar([low_x(2), mod_x(2), high_x(2)],[nanmean(pos_perc_resp_het(:,1)), nanmean(pos_perc_resp_het(:,2)), nanmean(pos_perc_resp_het(:,3))],-err_het,err_het);    
er.Color = [0 0 0];                            
er.LineStyle = 'none';  

scatter([ones(size(pos_perc_resp_hom,1),1)*low_x(1); ones(size(pos_perc_resp_hom,1),1)*mod_x(1); ones(size(pos_perc_resp_hom,1),1)*high_x(1)] , [pos_perc_resp_hom(:,1); pos_perc_resp_hom(:,2); pos_perc_resp_hom(:,3)], 'filled',  'MarkerFaceColor', map_hom(2,:))
b1 = bar([low_x(1), mod_x(1), high_x(1)], [nanmean(pos_perc_resp_hom(:,1)), nanmean(pos_perc_resp_hom(:,2)), nanmean(pos_perc_resp_hom(:,3))], 0.2,'FaceColor', map_wt(3,:))
set(b1(1), 'facecolor', map_hom(1,:), 'facealpha', 0.5)
% set(b3(2), 'facecolor', map_wt(2,:), 'facealpha', 0.5)
% set(b3(3), 'facecolor', map_wt(3,:), 'facealpha', 0.5)%%b1.CData(1,:) = cmap1(1,:);
er = errorbar([low_x(1), mod_x(1), high_x(1)],[nanmean(pos_perc_resp_hom(:,1)), nanmean(pos_perc_resp_hom(:,2)), nanmean(pos_perc_resp_hom(:,3))],-err_hom,err_hom);    
er.Color = [0 0 0];                            
er.LineStyle = 'none';  
xticks([1 2 3])
xticklabels(con_names)
legend([b3 b2 b1], group_names,  'Location', 'Eastoutside')
ylabel('% Excited Cells')
title([region_name, ' Comparison all cells pos'])
% saveas(gcf, fullfile(save_path, [region_name, '_POSperc_resp.png']))

% figure('units','centimeters','Position',[2 2 18 12])
subplot(2,1,2)
hold on
scatter([ones(size(neg_perc_resp_wt,1),1)*low_x(3); ones(size(neg_perc_resp_wt,1),1)*mod_x(3); ones(size(neg_perc_resp_wt,1),1)*high_x(3)] , [neg_perc_resp_wt(:,1); neg_perc_resp_wt(:,2); neg_perc_resp_wt(:,3)], 'filled',  'MarkerFaceColor', map_wt(2,:))
b4 = bar([low_x(3), mod_x(3), high_x(3)], [nanmean(neg_perc_resp_wt(:,1)), nanmean(neg_perc_resp_wt(:,2)), nanmean(neg_perc_resp_wt(:,3))], 0.2,'FaceColor', map_wt(3,:))
set(b4(1), 'facecolor', map_wt(1,:), 'facealpha', 0.5)
% set(b3(2), 'facecolor', map_wt(2,:), 'facealpha', 0.5)
% set(b3(3), 'facecolor', map_wt(3,:), 'facealpha', 0.5)%%b1.CData(1,:) = cmap1(1,:);
er = errorbar([low_x(3), mod_x(3), high_x(3)],[nanmean(neg_perc_resp_wt(:,1)), nanmean(neg_perc_resp_wt(:,2)), nanmean(neg_perc_resp_wt(:,3))],-err_wt2,err_wt2);    
er.Color = [0 0 0];                            
er.LineStyle = 'none';  

scatter([ones(size(neg_perc_resp_het,1),1)*low_x(2); ones(size(neg_perc_resp_het,1),1)*mod_x(2); ones(size(neg_perc_resp_het,1),1)*high_x(2)] , [neg_perc_resp_het(:,1); neg_perc_resp_het(:,2); neg_perc_resp_het(:,3)], 'filled',  'MarkerFaceColor', map_het(2,:))
b5 = bar([low_x(2), mod_x(2), high_x(2)], [nanmean(neg_perc_resp_het(:,1)), nanmean(neg_perc_resp_het(:,2)), nanmean(neg_perc_resp_het(:,3))], 0.2,'FaceColor', map_wt(3,:))
set(b5(1), 'facecolor', map_het(1,:), 'facealpha', 0.5)
% set(b3(2), 'facecolor', map_wt(2,:), 'facealpha', 0.5)
% set(b3(3), 'facecolor', map_wt(3,:), 'facealpha', 0.5)%%b1.CData(1,:) = cmap1(1,:);
er = errorbar([low_x(2), mod_x(2), high_x(2)],[nanmean(neg_perc_resp_het(:,1)), nanmean(neg_perc_resp_het(:,2)), nanmean(neg_perc_resp_het(:,3))],-err_het2,err_het2);    
er.Color = [0 0 0];                            
er.LineStyle = 'none';  

scatter([ones(size(neg_perc_resp_hom,1),1)*low_x(1); ones(size(neg_perc_resp_hom,1),1)*mod_x(1); ones(size(neg_perc_resp_hom,1),1)*high_x(1)] , [neg_perc_resp_hom(:,1); neg_perc_resp_hom(:,2); neg_perc_resp_hom(:,3)], 'filled',  'MarkerFaceColor', map_hom(2,:))
b6 = bar([low_x(1), mod_x(1), high_x(1)], [nanmean(neg_perc_resp_hom(:,1)), nanmean(neg_perc_resp_hom(:,2)), nanmean(neg_perc_resp_hom(:,3))], 0.2,'FaceColor', map_wt(3,:))
set(b6(1), 'facecolor', map_hom(1,:), 'facealpha', 0.5)
% set(b3(2), 'facecolor', map_wt(2,:), 'facealpha', 0.5)
% set(b3(3), 'facecolor', map_wt(3,:), 'facealpha', 0.5)%%b1.CData(1,:) = cmap1(1,:);
er = errorbar([low_x(1), mod_x(1), high_x(1)],[nanmean(neg_perc_resp_hom(:,1)), nanmean(neg_perc_resp_hom(:,2)), nanmean(neg_perc_resp_hom(:,3))],-err_hom2,err_hom2);    
er.Color = [0 0 0];                            
er.LineStyle = 'none';  

xticks([1 2 3])
xticklabels(con_names)
legend([b4 b5 b6], group_names, 'Location', 'Eastoutside')
ylabel('% Inhibited Cells')
title([ region_name, ' Comparison all cells neg'])
saveas(gcf, fullfile(save_path, [region_name, '_NEGperc_resp.png']))
saveas(gcf, fullfile(save_path, [region_name, '_NEGperc_resp.svg']))
% ylim([0 50])

end