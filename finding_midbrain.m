%% adding the midbrain here... 
new_brain = cell(3,1); 
for group = 1:no_group
    for fish = 1:size(brain_regions{group,1},2)
        cur_pos = positions{group}{1,fish};
        cur_brain = brain_regions{group}{1,fish};
        cut_of_line = mean(cur_pos(find(cur_brain == 11),1))-25; %min(cur_pos(find(cur_brain == 11),1)); % this should be Y

        midbrain = find(cur_pos(:,1) < cut_of_line); 
        midbrain = midbrain(~ismember(midbrain, find((cur_brain == 11))));
        cur_brain(midbrain) = 15;

        figure
        scatter3(cur_pos(:,1), cur_pos(:,2), cur_pos(:,3), 25, cur_brain, 'filled')
        axis equal,
        set (gca,'Zdir','reverse'),  set (gca,'Ydir','reverse'),  set (gca,'Xdir','reverse')
        colormap(cmap)
        colorbar
        
        new_brain{group,1}{fish} = cur_brain;

    end


end
old_brain = brain_regions; 
brain_regions = new_brain;

% %group 3, fish 6 is the example 
% 
% set(gcf, 'InvertHardCopy', 'off');
% set(gcf, 'DefaultFigureRenderer', 'painters');
% set(gcf,'renderer','painters');
% 
% saveas(gcf, fullfile(save_path, 'midbrainexampleCppGfish6_top.png'))
% saveas(gcf, fullfile(save_path, 'midbrainexampleCppGfish6_top.svg'))  
% saveas(gcf, fullfile(save_path, 'midbrainexampleCppGfish6_side.png'))
% saveas(gcf, fullfile(save_path, 'midbrainexampleCppGfish6_side.svg'))  