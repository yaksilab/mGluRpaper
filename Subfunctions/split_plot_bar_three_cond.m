function [fig1, fig2] = split_plot_bar_three_cond(perc_resp_wt, perc_resp_het, perc_resp_hom, map_wt, map_het, map_hom, con_names, group_names, tit, width, height)

% , save_path, save_name

low_x = [0.8 1.2] ; 
mod_x = [1.8 2.2];
high_x = [2.8 3.2]; 

err_wt = nanstd(perc_resp_wt(:,:), 0 ,1)/sqrt(size(perc_resp_wt,1));
err_het = nanstd(perc_resp_het(:,:), 0 ,1)/sqrt(size(perc_resp_het,1));
err_hom = nanstd(perc_resp_hom(:,:), 0 ,1)/sqrt(size(perc_resp_hom,1));



%%I want to make two out of one figure where I now plot 1 & 2 in one and 2& 3 in one
fig1 = figure('units','centimeters','Position',[2 2 width height])
hold on
scatter([ones(size(perc_resp_wt,1),1)*low_x(1); ones(size(perc_resp_wt,1),1)*mod_x(1); ones(size(perc_resp_wt,1),1)*high_x(1)] , [perc_resp_wt(:,1); perc_resp_wt(:,2); perc_resp_wt(:,3)], 50, 'filled',  'MarkerFaceColor', map_wt(3,:))
b3 = bar([low_x(1), mod_x(1), high_x(1)], [nanmean(perc_resp_wt(:,1)), nanmean(perc_resp_wt(:,2)), nanmean(perc_resp_wt(:,3))], 0.35,'FaceColor', map_wt(3,:))
set(b3(1), 'facecolor', map_wt(1,:), 'facealpha', 0.5)
% set(b3(2), 'facecolor', map_wt(2,:), 'facealpha', 0.5)
% set(b3(3), 'facecolor', map_wt(3,:), 'facealpha', 0.5)%%b1.CData(1,:) = cmap1(1,:);
er = errorbar([low_x(1), mod_x(1), high_x(1)],[nanmean(perc_resp_wt(:,1)), nanmean(perc_resp_wt(:,2)), nanmean(perc_resp_wt(:,3))],-err_wt,err_wt);    
er.Color = [0 0 0];                            
er.LineStyle = 'none';  

scatter([ones(size(perc_resp_het,1),1)*low_x(2); ones(size(perc_resp_het,1),1)*mod_x(2); ones(size(perc_resp_het,1),1)*high_x(2)] , [perc_resp_het(:,1); perc_resp_het(:,2); perc_resp_het(:,3)], 50,'filled',  'MarkerFaceColor', map_het(3,:))
b2 = bar([low_x(2), mod_x(2), high_x(2)], [nanmean(perc_resp_het(:,1)), nanmean(perc_resp_het(:,2)), nanmean(perc_resp_het(:,3))], 0.35,'FaceColor', map_wt(3,:))
set(b2(1), 'facecolor', map_het(1,:), 'facealpha', 0.5)
% set(b3(2), 'facecolor', map_wt(2,:), 'facealpha', 0.5)
% set(b3(3), 'facecolor', map_wt(3,:), 'facealpha', 0.5)%%b1.CData(1,:) = cmap1(1,:);
er = errorbar([low_x(2), mod_x(2), high_x(2)],[nanmean(perc_resp_het(:,1)), nanmean(perc_resp_het(:,2)), nanmean(perc_resp_het(:,3))],-err_het,err_het);    
er.Color = [0 0 0];                            
er.LineStyle = 'none'; 
xticks([1 2 3])
xticklabels(con_names)
legend([b3 b2], group_names{[1 2]})
title(tit)
legend('Location', 'eastoutside')
% ylim([0 65])
% xlim([0.5 2.6])
hold off

fig2 = figure('units','centimeters','Position',[2 2 width height])
hold on
scatter([ones(size(perc_resp_het,1),1)*low_x(1); ones(size(perc_resp_het,1),1)*mod_x(1); ones(size(perc_resp_het,1),1)*high_x(1)] , [perc_resp_het(:,1); perc_resp_het(:,2); perc_resp_het(:,3)], 50,'filled',  'MarkerFaceColor', map_het(3,:))
b2 = bar([low_x(1), mod_x(1), high_x(1)], [nanmean(perc_resp_het(:,1)), nanmean(perc_resp_het(:,2)), nanmean(perc_resp_het(:,3))], 0.35,'FaceColor', map_wt(3,:))
set(b2(1), 'facecolor', map_het(1,:), 'facealpha', 0.5)
% set(b3(2), 'facecolor', map_wt(2,:), 'facealpha', 0.5)
% set(b3(3), 'facecolor', map_wt(3,:), 'facealpha', 0.5)%%b1.CData(1,:) = cmap1(1,:);
er = errorbar([low_x(1), mod_x(1), high_x(1)],[nanmean(perc_resp_het(:,1)), nanmean(perc_resp_het(:,2)), nanmean(perc_resp_het(:,3))],-err_het,err_het);    
er.Color = [0 0 0];                            
er.LineStyle = 'none';  

scatter([ones(size(perc_resp_hom,1),1)*low_x(2); ones(size(perc_resp_hom,1),1)*mod_x(2); ones(size(perc_resp_hom,1),1)*high_x(2)] , [perc_resp_hom(:,1); perc_resp_hom(:,2); perc_resp_hom(:,3)], 50,'filled',  'MarkerFaceColor', map_hom(3,:))
b1 = bar([low_x(2), mod_x(2), high_x(2)], [nanmean(perc_resp_hom(:,1)), nanmean(perc_resp_hom(:,2)), nanmean(perc_resp_hom(:,3))], 0.35,'FaceColor', map_wt(3,:))
set(b1(1), 'facecolor', map_hom(1,:), 'facealpha', 0.5)
% set(b3(2), 'facecolor', map_wt(2,:), 'facealpha', 0.5)
% set(b3(3), 'facecolor', map_wt(3,:), 'facealpha', 0.5)%%b1.CData(1,:) = cmap1(1,:);
er = errorbar([low_x(2), mod_x(2), high_x(2)],[nanmean(perc_resp_hom(:,1)), nanmean(perc_resp_hom(:,2)), nanmean(perc_resp_hom(:,3))],-err_hom,err_hom);    
er.Color = [0 0 0];                            
er.LineStyle = 'none';  
xticks([1 2 3])
xticklabels(con_names)
legend([b2 b1], group_names{[2 3]})
title(tit)
legend('Location', 'eastoutside')
% ylim([0 65])
% xlim([0.5 2.6])
hold off

end