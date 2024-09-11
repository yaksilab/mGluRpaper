function [all_corr_together, pos_corr_together, neg_corr_together, all_corr_toge_plane, pos_corr_toge_plane, ...
    neg_corr_toge_plane, all_pais, all_pairs_fish] = corr_vs_distance_multi_vs_perplane(all_dffs_wt, brain_regions_wt, position_wt, ...
    all_on, all_off, con_names_plus_ong, save_path,name, cmap,thresh_dormedneurons_wt, dorsomed)

all_corr_together = cell(4,1); % just split into the periods
all_corr_toge_plane = cell(4,4); % split into periods and plane 

pos_corr_together = cell(4,1); % just split into the periods
pos_corr_toge_plane = cell(4,4); % split into periods and plane 

neg_corr_together = cell(4,1); % just split into the periods
neg_corr_toge_plane = cell(4,4); % split into periods and plane 
all_pais = cell(4,1); 
all_pairs_fish = cell(size(all_dffs_wt,2),4); 
for fish = 1:size(all_dffs_wt,2)
    % figure
    % plotsplot = [];
    for period = 1:size(all_on,2)
        
        
        current_dff = all_dffs_wt{1,fish};
        current_hab = current_dff(find(brain_regions_wt{1,fish} == 11),:); 
        
        hab_positions = position_wt{1,fish}(find(brain_regions_wt{1,fish} == 11),:);
         
        % for plane = 1:4 % hardcoding here that it will be only four planes..
        %         % here i can do the actual correlations now 
        %             current_hb_cells = find(hab_positions(:,5) == plane);

        if dorsomed
            current_hab = current_hab(find(thresh_dormedneurons_wt{1,fish} ==1),:);
            hab_positions = hab_positions(find(thresh_dormedneurons_wt{1,fish} ==1),:);

        end
        hab_period_dff = current_hab(:, all_on(period):all_off(period));
        [all_corr_dist, all_corr_dist_pos, all_corr_dist_neg, all_dist_corr_pairs, steps] = corr_vs_distance_AO(hab_period_dff, hab_positions, 1, 1);

        % subplot(1,2,1)
        % hold on
        % p1 = plot(steps(2:end), mean(all_corr_dist), 'LineWidth', 2)
        % title('Cor vs Dis All')
        % plotsplot = [plotsplot , p1];
        % subplot(1,2,2)
        % hold on
        % plot(steps(2:end), mean(all_corr_dist_pos), 'LineWidth', 2)
        % plot(steps(2:end), mean(all_corr_dist_neg), 'LineWidth', 2)
        % title('Cor vs Dis Split')
        % xlabel('Distance')
        % ylabel('Pearson Corr')
        % legend(plotsplot, con_names_plus_ong)
        % title(['Fish NO ' num2str(fish)])

        all_corr_together{period,1} = cat(1, all_corr_together{period,1}, all_corr_dist); 
        pos_corr_together{period,1} = cat(1, pos_corr_together{period,1}, all_corr_dist_pos);
        neg_corr_together{period,1} = cat(1, neg_corr_together{period,1}, all_corr_dist_neg);
        all_pais{period,1} = cat(1, all_pais{period,1}, all_dist_corr_pairs);
        all_pairs_fish{fish, period} = all_dist_corr_pairs; 
        for plane = 1:4 % hard coded 
            current_hb_cells = find(hab_positions(:,5) == plane);
            if length(current_hb_cells) < 2
                disp(['Not enough cells fish ', num2str(fish), ' and plane ', num2str(plane)])
            else
                [all_corr_dist, all_corr_dist_pos, all_corr_dist_neg, all_dist_corr_pairs, steps] = corr_vs_distance_AO(hab_period_dff(current_hb_cells,:), hab_positions(current_hb_cells,:), 0, 1);
    
                all_corr_toge_plane{period,plane} = cat(1, all_corr_toge_plane{period,plane}, all_corr_dist); 
                pos_corr_toge_plane{period,plane} = cat(1, pos_corr_toge_plane{period,plane}, all_corr_dist_pos);
                neg_corr_toge_plane{period,plane} = cat(1, neg_corr_toge_plane{period,plane}, all_corr_dist_neg);
            end
        end
    end

end
% cmap =[0.85, 0.85, 0.85; cmap]; 
% 
% figure('units','centimeters','Position',[2 2 9 15])% subplot(2,1,1)
% all_lines = []; 
% for con = 1: size(all_on,2)
%     hold on
%     H1=shadedErrorBar(steps(2:end),squeeze(nanmean(pos_corr_together{con,1}(:,1:end),1)),squeeze(nanstd(pos_corr_together{con,1}(:,1:end),0,1)/sqrt(size(pos_corr_together{con,1}(:,1:end),1))), 'lineProps','r');
%     H1.mainLine.Color = cmap(con,: ); %'#A23333';
%     H1.patch.FaceColor = cmap(con,: ); '#A23333'; 
%     H1.mainLine.LineWidth = 2; 
%     all_lines = [all_lines, H1.mainLine]; 
%     H2=shadedErrorBar(steps(2:end),squeeze(nanmean(neg_corr_together{con,1}(:,1:end),1)),squeeze(nanstd(neg_corr_together{con,1}(:,1:end),0,1)/sqrt(size(neg_corr_together{con,1}(:,1:end),1))), 'lineProps','r');
%     H2.mainLine.Color = cmap(con,: ); %'#A23333';
%     H2.patch.FaceColor = cmap(con,: ); '#A23333'; 
%     H2.mainLine.LineWidth = 2; 
% end
% legend(all_lines, con_names_plus_ong)
% yline(0, 'Color', 'k')
% xlabel('Distance')
% ylabel('Corr Value')
% title('Corr vs Dist')
% saveas(gcf, fullfile(save_path, [name, '_corr_vs_distance.png']))
% saveas(gcf, fullfile(save_path, [name, '_corr_vs_distance.svg']))
% 
% figure('units','centimeters','Position',[2 2 20 15])% subplot(2,1,1)
% for plane = 1:4
%     subplot(2,2,plane)
%     all_lines = []; 
%     for con = 1: size(all_on,2)
%         hold on
%         H1=shadedErrorBar(steps(2:end),squeeze(nanmean(pos_corr_toge_plane{con,plane}(:,1:end),1)),squeeze(nanstd(pos_corr_toge_plane{con,plane}(:,1:end),0,1)/sqrt(size(pos_corr_toge_plane{con,plane}(:,1:end),1))), 'lineProps','r');
%         H1.mainLine.Color = cmap(con,: ); %'#A23333';
%         H1.patch.FaceColor = cmap(con,: ); '#A23333'; 
%         H1.mainLine.LineWidth = 2; 
%         all_lines = [all_lines, H1.mainLine]; 
%         H2=shadedErrorBar(steps(2:end),squeeze(nanmean(neg_corr_toge_plane{con,plane}(:,1:end),1)),squeeze(nanstd(neg_corr_toge_plane{con,plane}(:,1:end),0,1)/sqrt(size(neg_corr_toge_plane{con,plane}(:,1:end),1))), 'lineProps','r');
%         H2.mainLine.Color = cmap(con,: ); %'#A23333';
%         H2.patch.FaceColor = cmap(con,: ); '#A23333'; 
%         H2.mainLine.LineWidth = 2; 
%     end
%     legend(all_lines, con_names_plus_ong)
%     yline(0, 'Color', 'k')
%     xlabel('Distance')
%     ylabel('Corr Value')
%     title([ 'Corr vs Dist Plane ', num2str(plane)]);
% end
% saveas(gcf, fullfile(save_path, [name, '_corr_vs_distance_planes.png']))
% saveas(gcf, fullfile(save_path, [name, '_corr_vs_distance_planes.svg']))

end