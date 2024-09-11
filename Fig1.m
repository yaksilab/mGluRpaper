%% Fig1 

%load the dataset 
load("X:\anna\Manuscript\FigData\Fig1\Collected_data_all_fish_fig1.mat"); % this is now the data of all the fish
load("X:\anna\code\Repositories\Anna-Code-Collection\everyday functions\beachVibes.mat")
save_path = 'X:\anna\Manuscript\FigData\Fig1'; %replace with your path

% to make the example figures I do Odor fish:4 Tap/Light Fish: 7
% to collect the odor data 
fishO = 6;
hab_resp = collected_data_fish.resp_odor{fishO,1}; 
pos_cells = find(hab_resp== 1);
neg_cells = find(hab_resp== -1);
hab_pos = collected_data_fish.position_odor{fishO,1};
HB.DFFON = hab_pos(pos_cells,:); 
HB.DFFInh = hab_pos(neg_cells,:);


HB.DFFInh_aver = squeeze(mean(collected_data_fish.dff_odor{fishO,1}(1:100,:,neg_cells),2))'; 
HB.DFFON_aver = squeeze(mean(collected_data_fish.dff_odor{fishO,1}(1:100,:,pos_cells),2))';  

tit = 'Odor'
% to collect the Light data 
fishLT = 7; 10 %10%
hab_resp = collected_data_fish.resp_light{fishLT,1}; 
pos_cells = find(hab_resp== 1);
neg_cells = find(hab_resp== -1);
hab_pos =collected_data_fish.position_light{fishLT,1};
HB.DFFON = hab_pos(pos_cells,:); 
HB.DFFInh = hab_pos(neg_cells,:);

% [pos_resp_list_dff_std, neg_resp_list_dff_std, perc_resp_cells_dff_std, perc_neg_resp_cells_dff_std, info_list_dff_std] = calculate_resp_cells_8trials_AO(collected_data_fish.dff_light{fishLT,1}, 2, baseline_dur, base_period_start, base_period_end, stim_period_end, no_con, con_trials, 2);        

% hab_resp = info_list_dff_std(:,1)
HB.DFFInh_aver = squeeze(mean(collected_data_fish.dff_light{fishLT,1}(1:130,1:8 , neg_cells),2))'; 
HB.DFFON_aver = squeeze(mean(collected_data_fish.dff_light{fishLT,1}(1:130,1:8, pos_cells),2))'; 

tit = 'Light'
% to collect the Tap data 
hab_resp = collected_data_fish.resp_tap{fishLT,1}; 
pos_cells = find(hab_resp== 1);
neg_cells = find(hab_resp== -1);
hab_pos = collected_data_fish.position_tap{fishLT,1};
HB.DFFON = hab_pos(pos_cells,:); 
HB.DFFInh = hab_pos(neg_cells,:);


HB.DFFInh_aver = squeeze(mean(collected_data_fish.dff_tap{fishLT,1}(1:130,9:16 , neg_cells),2))'; 
HB.DFFON_aver = squeeze(mean(collected_data_fish.dff_tap{fishLT,1}(1:130,9:16 , pos_cells),2))'; 
tit = 'Tap'

%% Plotting
figure('units','normalized','outerposition',[0 0 1 1])
subplot('Position',[0.12,0.82,0.18,0.18]);
 scatter3(hab_pos(:,1),hab_pos(:,2),hab_pos(:,3), 3, [0.5 0.5 0.5]);
hold on
scatter3(HB.DFFON(:,1),HB.DFFON(:,2),HB.DFFON(:,3), 32,'r', 'filled');
scatter3(HB.DFFInh(:,1),HB.DFFInh(:,2),HB.DFFInh(:,3), 32,'b', 'filled');


axis equal, view (22,80), grid on ; title (tit)
set (gca,'Zdir','reverse'),  set (gca,'Ydir','reverse'),  set (gca,'Xdir','reverse')
set(gca,'XTick',[]), set(gca,'YTick',[]),box off, axis tight, axis off,

subplot('Position',[0.12,0.62,0.18,0.18]);
 scatter3(hab_pos(:,1),hab_pos(:,2),hab_pos(:,3), 3, [0.5 0.5 0.5]);
hold on
scatter3(HB.DFFON(:,1),HB.DFFON(:,2),HB.DFFON(:,3), 32,'r', 'filled');
scatter3(HB.DFFInh(:,1),HB.DFFInh(:,2),HB.DFFInh(:,3), 32,'b', 'filled');


axis equal, view (-72,13), grid on ; title (tit)
set (gca,'Zdir','reverse'),  set (gca,'Ydir','reverse'),  set (gca,'Xdir','reverse')
set(gca,'XTick',[]), set(gca,'YTick',[]),box off, axis tight, axis off,

subplot('Position',[0.01,0.9,0.1,0.07]); imagesc(HB.DFFInh_aver(:,:)), caxis([-20 40]), colormap(beachVibes), title('sorted inhibited ')
subplot('Position',[0.01,0.7,0.1,0.07]); imagesc(HB.DFFON_aver(:,:)), caxis([-20 40]), colormap(beachVibes), title('sorted excited ')
subplot('Position',[0.01,0.86,0.1,0.04]); plot(mean(HB.DFFInh_aver(:,:)), 'b'), axis tight

subplot('Position',[0.01,0.66,0.1,0.04]); plot(mean(HB.DFFON_aver(:,:)), 'r'), axis tight

set(gcf, 'InvertHardCopy', 'off');
set(gcf, 'DefaultFigureRenderer', 'painters');
set(gcf,'renderer','painters');

%%
figure('units','centimeters','Position',[2 2 20 13])
subplot(1,2,1)
 scatter3(hab_pos(:,1),hab_pos(:,2),hab_pos(:,3), 50, [0.5 0.5 0.5], 'filled', 'MarkerFaceAlpha',.8,'MarkerEdgeAlpha',.8);
hold on
scatter3(HB.DFFON(:,1),HB.DFFON(:,2),HB.DFFON(:,3), 80,'r', 'filled');
scatter3(HB.DFFInh(:,1),HB.DFFInh(:,2),HB.DFFInh(:,3), 80,'b', 'filled');

axis equal, view (22,80), grid on ; title (tit)
set (gca,'Zdir','reverse'),  set (gca,'Ydir','reverse'),  set (gca,'Xdir','reverse')
% set(gca,'XTick',[]), set(gca,'YTick',[]),box off, axis tight, axis off,

subplot(1,2,2)
 scatter3(hab_pos(:,1),hab_pos(:,2),hab_pos(:,3), 30, [0.5 0.5 0.5], 'filled','MarkerFaceAlpha',.8,'MarkerEdgeAlpha',.8);
hold on
scatter3(HB.DFFON(:,1),HB.DFFON(:,2),HB.DFFON(:,3), 80,'r', 'filled');
scatter3(HB.DFFInh(:,1),HB.DFFInh(:,2),HB.DFFInh(:,3), 80,'b', 'filled');

axis equal, view (-72,13), grid on ; title (tit)
set (gca,'Zdir','reverse'),  set (gca,'Ydir','reverse'),  set (gca,'Xdir','reverse')
% set(gca,'XTick',[]), set(gca,'YTick',[]),box off, axis tight, axis off,

 set(gcf, 'InvertHardCopy', 'off');
    set(gcf, 'DefaultFigureRenderer', 'painters');
    set(gcf,'renderer','painters');

    %%
    %% All the save paths 
saveas(gcf, fullfile(save_path, 'Odour_Hb_heatmapsinferno.png'))
saveas(gcf, fullfile(save_path, 'Odour_Hb_heatmapsinferno.svg'))   


saveas(gcf, fullfile(save_path, 'Odour_Hb_heatmapshot.png'))
saveas(gcf, fullfile(save_path, 'Odour_Hb_heatmapshot.svg'))  

saveas(gcf, fullfile(save_path, 'Odour_Hb_3D.png'))
saveas(gcf, fullfile(save_path, 'Odour_Hb_3D.svg'))


saveas(gcf, fullfile(save_path, 'Light_Hb_heatmapsinferno.png'))
saveas(gcf, fullfile(save_path, 'Light_Hb_heatmapsinferno.svg'))   

saveas(gcf, fullfile(save_path, 'Light_Hb_heatmapshot.png'))
saveas(gcf, fullfile(save_path, 'Light_Hb_heatmapshot.svg'))  

saveas(gcf, fullfile(save_path, 'Light_Hb_3D.png'))
saveas(gcf, fullfile(save_path, 'Light_Hb_3D.svg'))  


saveas(gcf, fullfile(save_path, 'Tap_Hb_heatmapsinferno.png'))
saveas(gcf, fullfile(save_path, 'Tap_Hb_heatmapsinferno.svg'))   


saveas(gcf, fullfile(save_path, 'Tap_Hb_heatmapshot.png'))
saveas(gcf, fullfile(save_path, 'Tap_Hb_heatmapshot.svg'))  

saveas(gcf, fullfile(save_path, 'Tap_Hb_3D.png'))
saveas(gcf, fullfile(save_path, 'Tap_Hb_3D.svg'))  

%%
%% Now we also want to have perc responding cells
% responding cells were calculated with calculate_resp_cells_AO
% collected_data_fish.resp_odor
cmap1 = [0.85, 0.85, 0.85; 0.5430 0 0; 0, 0.4470, 0.7410;  0.133, 0.128, 0.177; 1.0000    0.8000    0.7961;  0.6902    0.8784    0.9020; 0.7843    0.6353    0.7843];

[neg_perc_resp_o, neg_perc_resp_l, neg_perc_resp_t] = plot_perc_resp_cells_AO(collected_data_fish.resp_odor', collected_data_fish.resp_light', collected_data_fish.resp_tap', 1, cmap1, cmap1,cmap1, -1, {'Stuff'});
close;
[pos_perc_resp_o, pos_perc_resp_l, pos_perc_resp_t] = plot_perc_resp_cells_AO(collected_data_fish.resp_odor', collected_data_fish.resp_light', collected_data_fish.resp_tap', 1, cmap1, cmap1, cmap1, 1, {'Stuff'});
close;

% Odor
err_od_pos = std(pos_perc_resp_o(:,:), 0 ,1)/sqrt(size(pos_perc_resp_o,1));
err_od_neg = std(neg_perc_resp_o(:,:), 0 ,1)/sqrt(size(neg_perc_resp_o,1));

figure('units','centimeters','Position',[2 2 8 10])
hold on 
AO_make_bar_scatter_line_plot(pos_perc_resp_o', neg_perc_resp_o', err_od_pos, err_od_neg, cmap1, 'Odor', 'Exc', 'Inh')
ylim([0 50])
saveas(gcf, fullfile(save_path, 'Od_perc.png'))
saveas(gcf, fullfile(save_path, 'Od_perc.svg')) 

figure('units','centimeters','Position',[2 2 10 5]) 
donout = donut([100-(mean(pos_perc_resp_o) + mean(neg_perc_resp_o)), mean(pos_perc_resp_o),mean(neg_perc_resp_o)],{'Non','Pos', 'Neg'}, cmap1)
legend('Location', 'Eastoutside')
title('Odor')
saveas(gcf, fullfile(save_path, 'Od_perc_pie.png'))
saveas(gcf, fullfile(save_path, 'Od_perc_pie.svg'))

% stats 
[p, h] = quick_statistic_signrank(pos_perc_resp_o, neg_perc_resp_o)

% Light

err_l_pos = std(pos_perc_resp_l(:,:), 0 ,1)/sqrt(size(pos_perc_resp_l,1));
err_l_neg = std(neg_perc_resp_l(:,:), 0 ,1)/sqrt(size(neg_perc_resp_l,1));

figure('units','centimeters','Position',[2 2 8 10])

hold on 
AO_make_bar_scatter_line_plot(pos_perc_resp_l', neg_perc_resp_l', err_l_pos, err_l_neg, cmap1, 'Light', 'Exc', 'Inh')
ylim([0 50])
saveas(gcf, fullfile(save_path, 'Li_perc.png'))
saveas(gcf, fullfile(save_path, 'Li_perc.svg')) 

figure('units','centimeters','Position',[2 2 10 5]) 
donout = donut([100-(mean(pos_perc_resp_l) + mean(neg_perc_resp_l)), mean(pos_perc_resp_l),mean(neg_perc_resp_l)],{'Non','Pos', 'Neg'}, cmap1)
legend('Location', 'Eastoutside')
title('Light')
saveas(gcf, fullfile(save_path, 'Li_perc_pie.png'))
saveas(gcf, fullfile(save_path, 'Li_perc_pie.svg'))

[p, h] = quick_statistic_signrank(pos_perc_resp_l, neg_perc_resp_l)

% Tap
err_t_pos = std(pos_perc_resp_t(:,:), 0 ,1)/sqrt(size(pos_perc_resp_t,1));
err_t_neg = std(neg_perc_resp_t(:,:), 0 ,1)/sqrt(size(neg_perc_resp_t,1));

figure('units','centimeters','Position',[2 2 8 10])

hold on 
AO_make_bar_scatter_line_plot(pos_perc_resp_t', neg_perc_resp_t', err_t_pos, err_t_neg, cmap1, 'Tap', 'Exc', 'Inh')
ylim([0 50])
saveas(gcf, fullfile(save_path, 'Ta_perc.png'))
saveas(gcf, fullfile(save_path, 'Ta_perc.svg')) 

figure('units','centimeters','Position',[2 2 10 5]) 
donout = donut([100-(mean(pos_perc_resp_t) + mean(neg_perc_resp_t)), mean(pos_perc_resp_t),mean(neg_perc_resp_t)],{'Non','Pos', 'Neg'}, cmap1)
legend('Location', 'Eastoutside')
title('Tap')
saveas(gcf, fullfile(save_path, 'Ta_perc_pie.png'))
saveas(gcf, fullfile(save_path, 'Ta_perc_pie.svg'))
[p, h] = quick_statistic_signrank(pos_perc_resp_t, neg_perc_resp_t)
