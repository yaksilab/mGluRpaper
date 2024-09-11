function AO_make_bar_scatter_line_plot(data1, data2, SEM_Dm_right_right_hab, SEM_Dm_right_left_hab, cmap1, tit, label1, label2)

scatter(ones(length(data1),1), data1, 'filled', 'MarkerEdgeColor', cmap1(1,:), 'MarkerFaceColor', cmap1(1,:)); 
hold on
scatter(ones(length(data2),1)*2, data2, 'filled', 'MarkerEdgeColor', cmap1(2,:),  'MarkerFaceColor', cmap1(2,:));
% 
b = bar(1, [nanmean(data1)],'stacked')  %,  'FaceColor','flat');
b.FaceColor = cmap1(1,:);
b.FaceAlpha = 0.5;
xlim([0 900]);

b = bar(2, [nanmean(data2)],'stacked')  %,  'FaceColor','flat');
b.FaceColor = cmap1(2,:);
b.FaceAlpha = 0.5;
xlim([0 900]);

er = errorbar([1 2],[nanmean(data1), nanmean(data2)],[SEM_Dm_right_right_hab, SEM_Dm_right_left_hab],[-SEM_Dm_right_right_hab, -SEM_Dm_right_left_hab]);    
er.Color = [0 0 0];                            
er.LineStyle = 'none';  
ylim([0 40]);
xticks([1 2])
xlim([0 3])
xticklabels({label1,label2})
ylabel('% Respoding Cells');

for exp=1:size(data1, 2)
    plot([1 2], [data1(exp), data2(exp)], 'Color', cmap1(1,:))
end
hold off
title(tit)


end