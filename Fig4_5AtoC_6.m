%% Fig4 

% Loading the data 
load('X:\anna\Manuscript\Data and matlab\FigData\Fig4\Fig4_cppg_injection_data.mat')

dff_fist = cppg_injection_data.dff_fish;
brain_regions = cppg_injection_data.brain_regions; 
positions = cppg_injection_data.positions;
dff_trialwise = cppg_injection_data.dff_trialw; 
fra_rat_two_pho = cppg_injection_data.frame_rate; 
stim_times = cppg_injection_data.stim_times; 

%% Ploting colors
cmap_wt = ['E4E6EB'; 'B0B3B8'; '18191A']; 
cmap_het = ['00FFFF'; '40E0D0'; '008080'];
cmap_hom = ['FF007F'; 'FF007F'; 'A94064'];
% cmap_hom = ['F89B29'; 'FF0F7B'; 'FF0F7B']; cmap_hom = ['FC8EAC';
% 'DE5D83'; 'FF0F7B']; ['F89B29'; 'FC5552'; 'FF0F7B'];
% cmap_hom = ['F89B29'; 'FF0F7B'; 'FF0F7B'];
RGB = hex2rgb(cmap_wt);
map_bef = RGB/255;

RGB = hex2rgb(cmap_het);
map_con = RGB/255;

RGB = hex2rgb(cmap_hom);
map_dru = RGB/255;
cmap1 = [0.85, 0.85, 0.85; 0.5430 0 0; 0, 0.4470, 0.7410;  0.133, 0.128, 0.177; 1.0000    0.8000    0.7961;  0.6902    0.8784    0.9020; 0.7843    0.6353    0.7843];
load('5_cluster_cmap.mat')
% red_map = [];
% % % cmap2 = [];
% for colorcoding = 1:9
%     c = uisetcolor;
%     red_map = [red_map; c];
% end
load("beachVibes.mat")

% load('X:\anna\code\mGluR\Ongoing_Resp_Comparision\5_cluster_cmap.mat')
all_cmap = {map_bef, map_con, map_dru};
cmap2 = [];
for colorcoding = 1:3
    c = uisetcolor;
    cmap2 = [cmap2; c];
end
%% Normalise poisitons

percentage_inclu = 0.40; 
all_planes_y = 1; 
dorsomed_ind = cell(no_group, 1); 
thresh_neu = cell(no_group, 1); 
norm_lat_med = cell(no_group, 1); 
norm_dor_ven = cell(no_group, 1); 
ven_neu  = cell(no_group, 1); 
for group = 1: no_group
    [dorsomed_ind{group, 1}, thresh_neu{group, 1}, norm_lat_med{group, 1}, norm_dor_ven{group, 1}, ven_neu{group, 1}] = calculate_dynamic_dorsomed_neurons(positions{group,1}, brain_regions{group,1}, percentage_inclu, all_planes_y, group_names{1,group}, save_path);
    

end
close all;
%% Responding cells

% Calculate the responding cells for the different conditions
% make bar plots for exc and inhib 
% show some traces and heatmap 

con_trials = results.conIdx; 
con_names = results.conNames; 
std_factor = 2; 
brainnumber = 11;
no_con = 3; 

baseline_dur = floor(5 * fra_rat_two_pho); 
method = 2; % 1: sign rank, 2: std
base_period_start = 2; %floor(2 * fra_rat_two_pho); 
stim_period_end = baseline_dur + 1 + floor(5 * fra_rat_two_pho); %stim_period_end = baseline_dur + floor(2 * fra_rat_two_pho);
base_period_end = baseline_dur;

time_in_s = 10; 
initial_stim = base_period_end:base_period_end+1+floor(time_in_s*fra_rat_two_pho);
middle_stim = initial_stim(end):initial_stim(end)+1+floor(time_in_s*fra_rat_two_pho); 
late_stim = middle_stim(end):middle_stim(end)+1+floor(time_in_s*fra_rat_two_pho); 
after_stim = late_stim(end):late_stim(end)+1+floor(time_in_s*fra_rat_two_pho);
lat_after_stim = after_stim(end):after_stim(end)+1+floor(time_in_s*fra_rat_two_pho);

diff_stim_period = [initial_stim; middle_stim; late_stim; after_stim; lat_after_stim]; 
% diff_stim_period = [initial_stim; middle_stim; late_stim]; 

base_period = base_period_start:base_period_end;
resp_period_list = cell(no_group,1);
resp_con_list = cell(no_group,1);

for group = 1:no_group
    [resp_period_list{group, 1}, resp_con_list{group, 1}] = calculate_responding_diff_periods(dff_trialwise{group,1}, brain_regions{group,1}, diff_stim_period, con_trials, base_period, no_con, brainnumber, std_factor); %positions{group,1}, brain_regions{group,1}


end
% [resp_period_list_wt, resp_con_list_wt] = calculate_responding_diff_periods(all_dfftrialw_wt_FL, brain_regions_wt, diff_stim_period, con_trials, base_period, no_con, brainnumber, std_factor);
% [resp_period_list_het, resp_con_list_het] = calculate_responding_diff_periods(all_dfftrialw_het_FL, brain_regions_het, diff_stim_period, con_trials, base_period, no_con, brainnumber, std_factor);
% [resp_period_list_hom, resp_con_list_hom] = calculate_responding_diff_periods(all_dfftrialw_hom_FL, brain_regions_hom, diff_stim_period, con_trials, base_period, no_con, brainnumber, std_factor);
% stim_period_name = {'0-5', '5-10', '10-15', '15-20', '20-25'};
stim_period_name = {'0-10', '10-20', '20-30', '30-40', '40-50'};

p_val_pos_period = {}; 
p_val_neg_period = {}; 
p_val_pos_period_dormed = {};
p_val_neg_period_dormed = {};
p_val_pos_period_ven = {};
p_val_neg_period_ven = {};
for stim = 1:size(diff_stim_period,1)
    
    region_list_bef = {}; 
    for fish = 1:size(brain_regions{1,1},2)

        region_list_bef{fish} = resp_period_list{1,1}{stim,1}{1,fish};

    end

    region_list_con = {}; 
    for fish = 1:size(brain_regions{2,1},2)
     region_list_con{fish} = resp_period_list{2,1}{stim,1}{1,fish};;

    end


    region_list_dru = {}; 
    for fish = 1:size(brain_regions{3,1},2)
        region_list_dru{fish} = resp_period_list{3,1}{stim,1}{1,fish};;

    end
    
    [neg_perc_resp_wt, neg_perc_resp_het, neg_perc_resp_hom] = plot_perc_resp_cells_AO(region_list_bef, region_list_con, region_list_dru, no_con, map_bef, map_con,map_dru, -1, con_names);
    % close;
    [pos_perc_resp_wt, pos_perc_resp_het, pos_perc_resp_hom] = plot_perc_resp_cells_AO(region_list_bef, region_list_con, region_list_dru, no_con, map_bef, map_con, map_dru, 1, con_names);
    % close;
    plot_perc_resp_cells_AO_together(pos_perc_resp_wt, pos_perc_resp_het, pos_perc_resp_hom,-neg_perc_resp_wt, -neg_perc_resp_het, -neg_perc_resp_hom,  map_bef, map_con,map_dru,con_names, ['all ' stim_period_name{stim}], save_path, group_names)    
    
    p_val_con_pos = {}; 
    p_val_con_neg = {}; 
    for con = 1:3
         [p1, h] = quick_statistic(pos_perc_resp_wt(:,con), pos_perc_resp_het(:,con))
         [p2, h] = quick_statistic(pos_perc_resp_wt(:,con), pos_perc_resp_hom(:,con))
         [p3, h] = quick_statistic(pos_perc_resp_het(:,con), pos_perc_resp_hom(:,con))
        p_val_con_pos{con} = [p1; p2; p3];
    
         [p1, h] = quick_statistic(neg_perc_resp_wt(:,con), neg_perc_resp_het(:,con))
         [p2, h] = quick_statistic(neg_perc_resp_wt(:,con), neg_perc_resp_hom(:,con))
         [p3, h] = quick_statistic(neg_perc_resp_het(:,con), neg_perc_resp_hom(:,con))
         p_val_con_neg{con} = [p1; p2; p3];
    end


    p_val_pos_period{stim} = {p_val_con_pos}; 
    p_val_neg_period{stim} = {p_val_con_neg}; 
    % 
    saveas(gcf, fullfile(save_path, [stim_period_name{stim}, '_Perc_respondingHb.svg']))
     saveas(gcf, fullfile(save_path, [stim_period_name{stim}, '_Perc_respondingHb.png']))

     % Now I also want a split version of the responding ones 
    if stim == 1
         for con = 1:no_con
            figure('units','centimeters','Position',[2 2 8 10])
        
            plot_bar_three_cond([pos_perc_resp_wt(:,con), neg_perc_resp_wt(:,con), neg_perc_resp_wt(:,con)], [pos_perc_resp_het(:,con), neg_perc_resp_het(:,con), neg_perc_resp_het(:,con)],[pos_perc_resp_hom(:,con), neg_perc_resp_hom(:,con), neg_perc_resp_hom(:,con)], map_bef, map_con, map_dru, {'Exc', 'Inh', 'Neg again'}, group_names)
            title(con_names{con})
            legend('Location', 'eastoutside')
            ylim([0 65])
            xlim([0.5 2.6])
            saveas(gcf, fullfile(save_path, [stim_period_name{stim}, '_', con_names{con},'_Perc_respondingHb.svg']))
            saveas(gcf, fullfile(save_path, [stim_period_name{stim}, '_', con_names{con},'_Perc_respondingHb.png']))

            [fig1, fig2] = split_plot_bar_three_cond([pos_perc_resp_wt(:,con), neg_perc_resp_wt(:,con), neg_perc_resp_wt(:,con)], [pos_perc_resp_het(:,con), neg_perc_resp_het(:,con), neg_perc_resp_het(:,con)],[pos_perc_resp_hom(:,con), neg_perc_resp_hom(:,con), neg_perc_resp_hom(:,con)], map_bef, map_con, map_dru, {'Exc', 'Inh', 'Neg again'}, group_names, con_names{con}, 8, 5)
            figure(fig1)
            ylim([0 65])
            xlim([0.5 2.6])

            figure(fig2)
            ylim([0 65])
            xlim([0.5 2.6])

            saveas(fig1, fullfile(save_path, [stim_period_name{stim}, '_', con_names{con},'_Perc_respondingHbCONCOM.svg']))
            saveas(fig1, fullfile(save_path, [stim_period_name{stim}, '_', con_names{con},'_Perc_respondingHbCONCOM.png']))
            saveas(fig2, fullfile(save_path, [stim_period_name{stim}, '_', con_names{con},'_Perc_respondingHbCONDRU.svg']))
            saveas(fig2, fullfile(save_path, [stim_period_name{stim}, '_', con_names{con},'_Perc_respondingHbCONDRU.png']))
            
         end
    end
end
% p_val_10s.p_val_pos = p_val_pos_period
% p_val_10s.p_val_neg = p_val_neg_period
% save(fullfile(save_path, 'pval_resp_10s.mat'), 'p_val_10s')

%% Uni vs. Mulitmodal cells (Selectivity)

%%%% maybe cut out the non positive responding cells 

% first I could just collect all of the dff values both for the heatmap and
% for the scatter plot
stim_per = 1; 
brainnumber = 11; 
collected_dff = cell(no_group,no_con);
selectivity_simple = cell(no_group, 1);
for group = 1:no_group
    for con = 1:no_con
        for fish = 1:size(dff_trialwise{group,1},2)
            
            curr_dff = dff_trialwise{group,1}{1,fish}(:,:,find(brain_regions{group,1}{1,fish} == brainnumber));
            % curr_dff = dff_trialwise{group,1}{1,fish}(:,:,:);

            curr_resp = resp_period_list{group,1}{stim_per,1}{1,fish}; 

            if fish == 1
                collected_dff{group,con} = squeeze(mean(curr_dff(:,con_trials(con,:),:),2)); 
            else
                 collected_dff{group,con} = cat(2, collected_dff{group,con}, squeeze(mean(curr_dff(:,con_trials(con,:),:),2))); 
            end

        end

    end

end

% Now for selectivity 

sel_code = cell(no_group,1);
sel_code_per_fish = cell(no_group,1);
% colorcode_hb = []; 
sel_neuron_p_gro = {no_group,1}; 
% color_code_per_fish_hb = {}; 
onlydor = 0; 

for group = 1:no_group
    selective_neurons_hb = {};
    color_code_per_fish_hb = {};
    colorcode_hb = [];
    for fish = 1:size(dff_trialwise{group,1},2)

        % now I want to colorcode them to what they respond to?
        % hb_neurons = find(brain_regions_ctrl{1,fish} == 11);
        curr_resp = resp_period_list{group,1}{stim_per,1}{1,fish};
        if onlydor
              % hb_neurons = hb_neurons(find(thresh_dormedneurons_wt{1,fish} ==1)); 
              curr_resp = curr_resp(find(thresh_neu{group,1}{1,fish} ==1),:); 
        end
        
        % now I go through the responding list and find which conditions a
        % neuron responds to
        baby_color = [];
        for neuron = 1:size(curr_resp,1)
            cur_neuron = neuron; 
            current_list = curr_resp(cur_neuron, :);
            % current_list = resp_list_ctrl_dff_std{1,fish}(cur_neuron, :); %* (-1); 
            % light only
            if current_list(1) == 1 && current_list(2) ~= 1 %== 0 %  && current_list(3) == 0 
                colorcode_hb = [colorcode_hb; 1];
                baby_color = [baby_color; 1];
            
            % tap only
            elseif current_list(1)  ~= 1 && current_list(2) == 1 %&& current_list(3)== 0 ~= 1 == 0 
                colorcode_hb = [colorcode_hb; 2];
                baby_color = [baby_color; 2];
            
            % elseif current_list(1) == 0 && current_list(2) == 0 %&& current_list(3) == typ 
            %     colorcode_hb = [colorcode_hb; 0];
            %     baby_color = [baby_color; 0];
            % elseif current_list(1) == 1 && current_list(2) == 0 %&& current_list(3) == typ 
            %     colorcode_hb = [colorcode_hb; 4];
            %     baby_color = [baby_color; 4];
            % elseif current_list(1) == 0 && current_list(2) == 1 %&& current_list(3) == typ 
            %     colorcode_hb = [colorcode_hb; 5];
            %     baby_color = [baby_color; 5];
            
            %both respo
            elseif current_list(1) == 1 && current_list(2) == 1 %&& current_list(3) == typ 
                colorcode_hb = [colorcode_hb; 3];
                baby_color = [baby_color; 3];
            % non resp
            else
                colorcode_hb = [colorcode_hb; 0];
                baby_color = [baby_color; 0];
    %             disp('Here')
            end
        end
        selective_neurons_hb{fish} = {find(baby_color == 1), find(baby_color == 2), find(baby_color == 3), find(baby_color == 4), find(baby_color == 5), find(baby_color == 6)};
        color_code_per_fish_hb{fish} = baby_color; 
    end
    sel_code{group, 1} = colorcode_hb; 
    sel_code_per_fish{group,1} = color_code_per_fish_hb; 
    sel_neuron_p_gro{group,1} = selective_neurons_hb; 
end

% find the positive and negative resp cells and if they are tap or light
% only of respond to both 
fish_perc_sel = cell(no_group,1); 

selective_label = {'Only Light', 'Only Tap', 'Only Both', 'Light and Both', 'Tap and Both' , 'All', 'No Sel'}; 
for group = 1:no_group
    [fish_perc_sel{group,1},  fish_perc_simple_wt] = selectivity_pie(sel_neuron_p_gro{group,1}, group_names{group}, save_path, cmap1);

end

%% Plotting 
end_Per = floor(30*fra_rat_two_pho); 
% first the collected dffs for the two conditions per group
for group = 1:no_group
    figure('units','centimeters','Position',[2 2 30 5])
    for con = 1:no_con
        % first we sort from highest to lowest
        
        [~, sortidx] = sort(mean(collected_dff{group,con}(diff_stim_period(stim_per,:),:),1),2,'descend'); 
        subplot(1,no_con,con)
        imagesc(collected_dff{group,con}(1:end_Per,sortidx)')
        colormap(beachVibes)
        caxis([-20 40])
        title([group_names{1,group} ' ' con_names{1,con}])
        % colorbar

    end
    saveas(gcf, fullfile(save_path, [group_names{1,group}, '_heatmapsallcells.svg']))
    saveas(gcf, fullfile(save_path, [group_names{1,group}, '_heatmapsallcells.png']))


    % scattering the avergae activiy in the stim period
    % figure('units','centimeters','Position',[2 2 20 12])
    figure('units','centimeters','Position',[2 2 10 15])
    scatter(mean(collected_dff{group,1}(diff_stim_period(stim_per,:),:),1), mean(collected_dff{group,2}(diff_stim_period(stim_per,:),:),1), 50,sel_code{group, 1} , 'filled'); 
    
    % scatter(mean(collected_dff{group,1}(diff_stim_period(stim_per,:),find(sel_code{group, 1} ~= 0)),1), mean(collected_dff{group,2}(diff_stim_period(stim_per,:),find(sel_code{group, 1} ~= 0)),1), 25,sel_code{group, 1}(find(sel_code{group, 1} ~= 0)), 'filled'); 

    colormap(cmap1(1:4,:))
    title([group_names{1,group}])
    xlabel(['dFF ', con_names{1,1}])
    ylabel(['dFF ', con_names{1,2}])
    xlim([-20 180]);
    ylim([-20 140])
    % saveas(gcf, fullfile(save_path, [group_names{1,group}, '_scatterofmeandff2justresp.svg']))
    % saveas(gcf, fullfile(save_path, [group_names{1,group}, '_scatterofmeandff2justresp.png']))
    saveas(gcf, fullfile(save_path, [group_names{1,group}, '_scatterofmeandff2alsoneg.svg']))
    saveas(gcf, fullfile(save_path, [group_names{1,group}, '_scatterofmeandff2alsoneg.png']))
    
    % a little donut to show the ratio
    perc_all = [length(find(sel_code{group, 1} == 0))/length(sel_code{group, 1}), length(find(sel_code{group, 1} == 1))/length(sel_code{group, 1}), length(find(sel_code{group, 1} == 2))/length(sel_code{group, 1}), length(find(sel_code{group, 1} == 3))/length(sel_code{group, 1}),]
    figure('units','centimeters','Position',[2 2 10 5]) 
    donout = donut(perc_all,{'Non','Light', 'Tap', 'Both'}, cmap1([1, 5,6,7], :));
    legend('Location', 'Eastoutside')
    title([group_names{1,group} ' Selectivity'])
    saveas(gcf, fullfile(save_path, [group_names{1,group}, '_donutspread_selectivity.svg']))
    saveas(gcf, fullfile(save_path, [group_names{1,group}, '_donutspread_selectivity.png']))

end

% Traces
figure('units','centimeters','Position',[2 2 30 5])
for con = 1:no_con
    subplot(1,no_con,con)
    hold on
    plplpl = [];
    for group = 1:no_group
       H1=shadedErrorBar([1:size(mean(collected_dff{group,con}(1:end_Per,:),2),1)],squeeze(mean(collected_dff{group,con}(1:end_Per,:),2)),squeeze(nanstd(collected_dff{group,con}(1:end_Per,:),0,2)/sqrt(size(collected_dff{group,con}(:,:),2))), 'lineProps','r');
        H1.mainLine.Color = all_cmap{1,group}(2,: ); %'#A23333';
        H1.patch.FaceColor = all_cmap{1,group}(2,: ); '#A23333'; 
        % p1 = plot(mean(collected_dff{group,con}(:,:),2), 'Color',all_cmap{1,group}(2,:))
    
        plplpl = [plplpl, H1.mainLine];
    
    end
    legend(plplpl, group_names)
    title(con_names{1,con})
    ylim([-2 25])
    xlim([0 end_Per])
    xline(stim_period(1), 'Color', 'k')
end
saveas(gcf, fullfile(save_path, ['Traces_groups.svg']))
saveas(gcf, fullfile(save_path, ['Traces_groups.png']))
% Now the selectivty per fish 


figure('units','centimeters','Position',[2 2 14 5])
hold on
plotplot = [];
xvalue = [[0.75, 1.75, 2.75, 3.75]; [1 2 3 4]; [1.25 2.25 3.25 4.25]];

for group = 1:no_group
    
    for sel = 1:3 
        scatter(ones(size(fish_perc_sel{group,1},1),1)*xvalue(group,sel), fish_perc_sel{group,1}(:,sel), 'filled',  'MarkerFaceColor',all_cmap{1,group}(2,:))
    end
    err = std(fish_perc_sel{group,1}(:,:), 0 ,1)/sqrt(size(fish_perc_sel{group,1},1));

    b4 = bar(xvalue(group,1:3), nanmean(fish_perc_sel{group,1}(:,1:3),1), 0.2,'FaceColor', all_cmap{1,group}(2,:))
    set(b4(1), 'facecolor', all_cmap{1,group}(2,:), 'facealpha', 0.5)
    % set(b3(2), 'facecolor', map_wt(2,:), 'facealpha', 0.5)
    % set(b3(3), 'facecolor', map_wt(3,:), 'facealpha', 0.5)%%b1.CData(1,:) = cmap1(1,:);
    er = errorbar(xvalue(group,1:3),nanmean(fish_perc_sel{group,1}(:,1:3),1),-err(1:3),err(1:3));    
    er.Color = [0 0 0];                            
    er.LineStyle = 'none'; 

    plotplot = [plotplot, b4];
end
xticks(xvalue(2,1:3))
xticklabels(selective_label)
legend(plotplot, group_names, 'Location', 'Eastoutside')
title('Selectivity')
ylabel('Perc Neurons')
saveas(gcf, fullfile(save_path, [group_names{1,3}, '_Selectivity.svg']))
saveas(gcf, fullfile(save_path, [group_names{1,3}, '_Selectivity.png']))

p_val_sel = {}; 
for sel = 1: 3
    [p1, h] = quick_statistic(fish_perc_sel{1,1}(:,sel), fish_perc_sel{2,1}(:,sel))
    [p2, h] = quick_statistic(fish_perc_sel{1,1}(:,sel), fish_perc_sel{3,1}(:,sel))
    [p3, h] = quick_statistic(fish_perc_sel{2,1}(:,sel), fish_perc_sel{3,1}(:,sel))
    p_val_sel{sel} = [p1; p2; p3]; 
end
pval_sel.all_p = p_val_sel;
% save(fullfile(save_path, 'pval_sel.mat'), 'pval_sel')


% NOw we combine only light and tap tp unimodal
figure('units','centimeters','Position',[2 2 10 5])
hold on
plotplot = [];
xvalue = [[0.75, 1.75, 2.75, 3.75]; [1 2 3 4]; [1.25 2.25 3.25 4.25]];
combined_sel = cell(no_group,1); 
for group = 1:no_group
    combined_sel{group,1} = [fish_perc_sel{group,1}(:,1)+fish_perc_sel{group,1}(:,2), fish_perc_sel{group,1}(:,3)]
    
    scatter(ones(size(fish_perc_sel{group,1},1),1)*xvalue(group,1), combined_sel{group,1}(:,1), 'filled',  'MarkerFaceColor',all_cmap{1,group}(2,:))
    scatter(ones(size(fish_perc_sel{group,1},1),1)*xvalue(group,2), combined_sel{group,1}(:,2), 'filled',  'MarkerFaceColor',all_cmap{1,group}(2,:))

    err = std(combined_sel{group,1}, 0 ,1)/sqrt(size(combined_sel{group,1},1));

    b4 = bar(xvalue(group,1:2), nanmean(combined_sel{group,1},1), 0.2,'FaceColor', all_cmap{1,group}(2,:))
    set(b4(1), 'facecolor', all_cmap{1,group}(2,:), 'facealpha', 0.5)
    % set(b3(2), 'facecolor', map_wt(2,:), 'facealpha', 0.5)
    % set(b3(3), 'facecolor', map_wt(3,:), 'facealpha', 0.5)%%b1.CData(1,:) = cmap1(1,:);
    er = errorbar(xvalue(group,1:2),nanmean(combined_sel{group,1},1),-err,err);    
    er.Color = [0 0 0];                            
    er.LineStyle = 'none'; 

    plotplot = [plotplot, b4];
end
xticks(xvalue(2,1:2))
xticklabels({'Unimodal', 'Mulitmodal'})
legend(plotplot, group_names, 'Location', 'Eastoutside')
ylabel('Perc Neurons')
title('Selectivity Simple')
xlim([0.5 2.6])
saveas(gcf, fullfile(save_path, [group_names{1,3}, '_Selectivity_simple.svg']))
saveas(gcf, fullfile(save_path, [group_names{1,3}, '_Selectivity_simple.png']))

p_val_sel_simple = {}; 
for sel = 1: 2
    [p1, h] = quick_statistic(combined_sel{1,1}(:,sel), combined_sel{2,1}(:,sel))
    [p2, h] = quick_statistic(combined_sel{1,1}(:,sel), combined_sel{3,1}(:,sel))
    [p3, h] = quick_statistic(combined_sel{2,1}(:,sel), combined_sel{3,1}(:,sel))
    p_val_sel_simple{sel} = [p1; p2; p3]; 
end
pval_sel.all_p_simple = p_val_sel_simple;
save(fullfile(save_path, 'pval_sel.mat'), 'pval_sel')

[fig1, fig2] = split_plot_bar_three_cond([combined_sel{1,1}(:,1), combined_sel{1,1}(:,2), combined_sel{1,1}(:,2)], [combined_sel{2,1}(:,1), combined_sel{2,1}(:,2), combined_sel{2,1}(:,2)],[combined_sel{3,1}(:,1), combined_sel{3,1}(:,2), combined_sel{3,1}(:,2)], map_bef, map_con, map_dru, {'Unimodal', 'Multimodal', 'Neg again'}, group_names, {'Selectivity'}, 8 , 5)
figure(fig1)
ylim([0 100])
xlim([0.5 2.6])

figure(fig2)
ylim([0 100])
xlim([0.5 2.6])

saveas(fig1, fullfile(save_path, [group_names{1,3}, '_Selectivity_simple_CON.svg']))
saveas(fig1, fullfile(save_path, [group_names{1,3}, '_Selectivity_simple_CON.png']))
saveas(fig2, fullfile(save_path, [group_names{1,3}, '_Selectivity_simple_CONDRU.svg']))
saveas(fig2, fullfile(save_path, [group_names{1,3}, '_Selectivity_simple_CONDRU.png']))

%% Trial Trial variability? response vector correlations? similarity? 
% for group = 1:no_group
%     [all_cond_saving] = trial_trial_var_groups(dff_trialwise{group,1}, brain_regions{group,1}, fra_rat_two_pho, con_trials, con_names, save_path, group_names{group}, 'Hb', 11)
% end

% maybe i make a vector over the five sec initial response window? start
% out just with like the mean? 
dorsomed = 0; 
resp_corr_val = cell(3,1);
for group = 1:no_group
    [resp_corr_val{group,1}] = response_vector_correl(dff_trialwise{group,1}, brain_regions{group,1}, positions{group,1}, stim_period, con_trials, resp_period_list{group,1}, dorsomed, thresh_neu{group,1});

end


figure('units','centimeters','Position',[2 2 10 5])
plot_bar_three_cond(resp_corr_val{1,1}, resp_corr_val{2,1},resp_corr_val{3,1}, map_bef, map_con, map_dru, {'L vs T', 'L vs LT', 'T vs LT'}, group_names)
title('Corr Val dff' )
ylabel('Corr Values')
xlim([0.5 1.6])
legend('Location', 'Eastoutside')
if dorsomed
    saveas(gcf, fullfile(save_path, ['Groups_respvec_correl_dorsomed.svg']))
    saveas(gcf, fullfile(save_path, ['Groups_respvec_correl_dorsomed.png']))
else
    saveas(gcf, fullfile(save_path, ['Groups_respvec_correl_dff.svg']))
    saveas(gcf, fullfile(save_path, ['Groups_respvec_correl_dff.png']))
end


p_val_respcorr= {}; 

% ylabel('% Perc of late inhibited cells')
for con = 1:3
     [p1, h] = quick_statistic(resp_corr_val{1,1}(:,con), resp_corr_val{2,1}(:,con))
     [p2, h] = quick_statistic(resp_corr_val{1,1}(:,con), resp_corr_val{3,1}(:,con))
     [p3, h] = quick_statistic(resp_corr_val{2,1}(:,con), resp_corr_val{3,1}(:,con))
    p_val_respcorr{con} = [p1; p2; p3];

end
if dorsomed

   p_val_respveccorr.dorsoemd = p_val_respcorr; 
else
    p_val_respveccorr.all = p_val_respcorr; 
end
save(fullfile(save_path, 'p_val_resp_vector.mat'), 'p_val_respveccorr')

[fig1, fig2] = split_plot_bar_three_cond([resp_corr_val{1,1}(:,1), resp_corr_val{1,1}(:,2), resp_corr_val{1,1}(:,2)], [resp_corr_val{2,1}(:,1), resp_corr_val{2,1}(:,2), resp_corr_val{2,1}(:,2)],[resp_corr_val{3,1}(:,1), resp_corr_val{3,1}(:,2), resp_corr_val{3,1}(:,2)], map_bef, map_con, map_dru, {'L vs T', 'LvsLT', 'T vs LT'}, group_names, {'Resp Corr'}, 6, 5)
figure(fig1)
xlim([0.5 1.6])
% xlim([0.5 2.6])
ylabel('Corr Values')
figure(fig2)
xlim([0.5 1.6])
% xlim([0.5 2.6])
ylabel('Corr Values')
saveas(fig1, fullfile(save_path, [group_names{1,3}, '_RespCorr_CON.svg']))
saveas(fig1, fullfile(save_path, [group_names{1,3}, '_RespCorr_CON.png']))
saveas(fig2, fullfile(save_path, [group_names{1,3}, '_RespCorr_CONDRU.svg']))
saveas(fig2, fullfile(save_path, [group_names{1,3}, '_RespCorr_CONDRU.png']))

%% Amplitudes 
dorsomed = 0;
avg_ampl = cell(3,1); 
collected_ampl = cell(3,1);

for group = 1:no_group

    [avg_ampl{group,1}, collected_ampl{group,1}] = check_amplitudes(dff_trialwise{group,1}, brain_regions{group,1}, positions{group,1}, ...
        con_names, save_path,group_names{group}, all_cmap{1,group}, thresh_neu{group,1}, dorsomed, con_trials, diff_stim_period(1,:));

end


% NOw we combine only light and tap tp unimodal
figure('units','centimeters','Position',[2 2 10 5])
hold on
plotplot = [];
xvalue = [[0.75, 1.75, 2.75, 3.75]; [1 2 3 4]; [1.25 2.25 3.25 4.25]];
 
for group = 1:no_group
    
    
    scatter(ones(size(avg_ampl{group,1},2),1)*xvalue(group,1), avg_ampl{group,1}(1,:), 'filled',  'MarkerFaceColor',all_cmap{1,group}(2,:))
    scatter(ones(size(avg_ampl{group,1},2),1)*xvalue(group,2), avg_ampl{group,1}(2,:), 'filled',  'MarkerFaceColor',all_cmap{1,group}(2,:))

    scatter(ones(size(avg_ampl{group,1},2),1)*xvalue(group,3), avg_ampl{group,1}(3,:), 'filled',  'MarkerFaceColor',all_cmap{1,group}(2,:))
    
    err = std(avg_ampl{group,1}(1:3,:), 0 ,2)/sqrt(size(avg_ampl{group,1}(1:3,:),2));

    b4 = bar(xvalue(group,1:3), nanmean(avg_ampl{group,1}(1:3,:),2), 0.2,'FaceColor', all_cmap{1,group}(2,:))
    set(b4(1), 'facecolor', all_cmap{1,group}(2,:), 'facealpha', 0.5)
    % set(b3(2), 'facecolor', map_wt(2,:), 'facealpha', 0.5)
    % set(b3(3), 'facecolor', map_wt(3,:), 'facealpha', 0.5)%%b1.CData(1,:) = cmap1(1,:);
    er = errorbar(xvalue(group,1:3), nanmean(avg_ampl{group,1}(1:3,:),2),-err,err);    
    er.Color = [0 0 0];                            
    er.LineStyle = 'none'; 

    plotplot = [plotplot, b4];
end

xticks(xvalue(2,1:3))
xticklabels({'Light', 'Tap', 'LightTap'})
legend(plotplot, group_names, 'Location', 'Eastoutside')
ylabel('dFF Amplitude')
title('Amplitude')
% ylim([0 60])
xlim([1.5 2.6])
saveas(gcf, fullfile(save_path, [group_names{1,3}, '_Amplitude.svg']))
saveas(gcf, fullfile(save_path, [group_names{1,3}, '_Amplitude.png']))

p_val_ampl = {}; 
for sel = 1: 3
    [p1, h] = quick_statistic(avg_ampl{1,1}(sel,:), avg_ampl{2,1}(sel,:))
    [p2, h] = quick_statistic(avg_ampl{1,1}(sel,:), avg_ampl{3,1}(sel,:))
    [p3, h] = quick_statistic(avg_ampl{2,1}(sel,:), avg_ampl{3,1}(sel,:))
    p_val_ampl{sel} = [p1; p2; p3]; 
end
save(fullfile(save_path, 'p_val_ampl.mat'), 'p_val_ampl')

[fig1, fig2] = split_plot_bar_three_cond(avg_ampl{1,1}', avg_ampl{2,1}', avg_ampl{3,1}', map_bef, map_con, map_dru, {'L', 'T', 'LT'},group_names, 'Amplitude', 5, 5 )
figure(fig1)
ylim([0 20])
xlim([1.5 2.6])

figure(fig2)
ylim([0 20])
% xlim([0.5 2.6])
xlim([1.5 2.6])
saveas(fig1, fullfile(save_path, ['AmplitudeHbCONCOM_T.svg']))
saveas(fig1, fullfile(save_path, ['AmplitudeHbCONCOM_T.png']))
saveas(fig2, fullfile(save_path, ['AmplitudeHbCONDRU_T.svg']))
saveas(fig2, fullfile(save_path, ['AmplitudeHbCONDRU_T.png']))

figure(fig1)
ylim([0 20])
xlim([0.5 1.6])

figure(fig2)
ylim([0 20])
xlim([0.5 1.6])
% xlim([1.5 2.6])
saveas(fig1, fullfile(save_path, ['AmplitudeHbCONCOM_L.svg']))
saveas(fig1, fullfile(save_path, ['AmplitudeHbCONCOM_L.png']))
saveas(fig2, fullfile(save_path, ['AmplitudeHbCONDRU_L.svg']))
saveas(fig2, fullfile(save_path, ['AmplitudeHbCONDRU_L.png']))

%% Fig5
stim_times = results.stim_triggers; 
all_on = [230, stim_times(1)-100, stim_times(9)-100, stim_times(17)-100];
all_off = [2000, stim_times(8)+100, stim_times(16)+100, stim_times(24)+100];
n_clusters = 5; 
con_names_plus_ong = {'Ong', con_names{1}, con_names{2}, con_names{3}};

all_corr_together = cell(no_group,1);
pos_corr_together = cell(no_group,1);
neg_corr_together = cell(no_group,1);
all_corr_toge_plane = cell(no_group,1);
pos_corr_toge_plane = cell(no_group,1);
neg_corr_toge_plane = cell(no_group,1);
all_pairs = cell(no_group,1);
all_pairs_fish = cell(no_group,1);
all_cmap = {map_bef, map_con, map_dru};
for group = 1:no_group
    [all_corr_together{group,1}, pos_corr_together{group,1}, neg_corr_together{group,1}, all_corr_toge_plane{group,1}, pos_corr_toge_plane{group,1}, ...
    neg_corr_toge_plane{group,1}, all_pairs{group,1}, all_pairs_fish{group,1}] = corr_vs_distance_multi_vs_perplane(dff_fist{group,1}, brain_regions{group,1}, positions{group,1}, ...
    all_on, all_off, con_names_plus_ong, save_path, group_names{group}, all_cmap{1,group}, [], 0);

end
steps = floor(linspace(0,100,26))

figure('units','centimeters','Position',[2 2 30 10])% subplot(2,1,1)
all_lines = []; 
for con = 1: size(all_on,2)
    subplot(1,4,con)
    for group = 1:no_group
    hold on
    H1=shadedErrorBar(steps(2:end),squeeze(nanmean(all_corr_together{group,1}{con,1}(:,1:end),1)),squeeze(nanstd(all_corr_together{group,1}{con,1}(:,1:end),0,1)/sqrt(size(all_corr_together{group,1}{con,1}(:,1:end),1))), 'lineProps','r');
    H1.mainLine.Color = all_cmap{1,group}(2,: ); %'#A23333';
    H1.patch.FaceColor = all_cmap{1,group}(2,: ); '#A23333'; 
    H1.mainLine.LineWidth = 2; 
    all_lines = [all_lines, H1.mainLine]; 

    % H2=shadedErrorBar(steps(2:end),squeeze(nanmean(neg_corr_together{group,1}{con,1}(:,1:end),1)),squeeze(nanstd(neg_corr_together{group,1}{con,1}(:,1:end),0,1)/sqrt(size(neg_corr_together{group,1}{con,1}(:,1:end),1))), 'lineProps','r');
    % H2.mainLine.Color = all_cmap{1,group}(2,: ); %'#A23333';
    % H2.patch.FaceColor = all_cmap{1,group}(2,: ); '#A23333'; 
    % H2.mainLine.LineWidth = 2; 

    end
   
    yline(0, 'Color', 'k')
    xlabel('Distance')
    ylabel('Corr Value')
    title(['Corr vs Dist ', con_names_plus_ong{con}])
    ylim([-0.01 0.5])
end
legend(all_lines, group_names)

saveas(gcf, fullfile(save_path, ['Groups_corr_vs_distance.png']))
saveas(gcf, fullfile(save_path, [ 'Groups_corr_vs_distance.svg']))

%% Avg distance for pos and neg
con = 1
avg_dist_pos = cell(no_group,1);
avg_dist_neg = cell(no_group,1);
for group = 1:no_group
    fish_curr_pos = nan(size(all_pairs_fish{group,1},1),1);
    fish_curr_neg = nan(size(all_pairs_fish{group,1},1),1);
    for fish = 1:size(all_pairs_fish{group,1},1)
        curr_pairs = all_pairs_fish{group,1}{fish,con};
        sign_ones = find(curr_pairs(:,3)<= 0.001);
        curr_pairs_sign = curr_pairs(sign_ones, :); 
        posit = find(curr_pairs_sign(:,2)>= 0);
        negat = find(curr_pairs_sign(:,2) < 0);
        fish_curr_pos(fish) = mean(curr_pairs_sign(posit, 1));
        fish_curr_neg(fish) = mean(curr_pairs_sign(negat, 1));
    end
    avg_dist_pos{group,1} = fish_curr_pos; 
    avg_dist_neg{group,1} = fish_curr_neg; 

end
[p, h] = quick_statistic_signrank(avg_dist_pos{1,1}, avg_dist_neg{1,1})
[p, h] = quick_statistic_signrank(avg_dist_pos{2,1}, avg_dist_neg{2,1})
[p, h] = quick_statistic_signrank(avg_dist_pos{3,1}, avg_dist_neg{3,1})

[fig1, fig2] = split_plot_bar_three_cond([avg_dist_pos{1,1}, avg_dist_neg{1,1},  avg_dist_neg{1,1}], [avg_dist_pos{2,1}, avg_dist_neg{2,1},  avg_dist_neg{2,1}],[avg_dist_pos{3,1}, avg_dist_neg{3,1},  avg_dist_neg{3,1}], map_bef, map_con, map_dru, {'Pos', 'Neg', 'Neg again'}, group_names, {'Avg Dist'}, 8, 5)
figure(fig1)
ylim([0 65])
xlim([0.5 2.6])

figure(fig2)
ylim([0 65])
xlim([0.5 2.6])

% saveas(fig1, fullfile(save_path, [stim_period_name{stim}, '_', con_names{con},'_Perc_respondingHbCONCOM.svg']))
% saveas(fig1, fullfile(save_path, [stim_period_name{stim}, '_', con_names{con},'_Perc_respondingHbCONCOM.png']))
% saveas(fig2, fullfile(save_path, [stim_period_name{stim}, '_', con_names{con},'_Perc_respondingHbCONDRU.svg']))
% saveas(fig2, fullfile(save_path, [stim_period_name{stim}, '_', con_names{con},'_Perc_respondingHbCONDRU.png']))


err_con = nanstd(avg_dist_pos{2,1}, 0 ,1)/sqrt(size(avg_dist_pos{2,1},1));
err_con2 = nanstd(avg_dist_neg{2,1}, 0 ,1)/sqrt(size(avg_dist_neg{2,1},1));

% here i want to plot the avg distance of all positive and all negative
% correlations
figure('units','centimeters','Position',[2 2 6 5])
hold on
scatter([ones(size(avg_dist_pos{2,1},1),1)*1] , avg_dist_pos{2,1}, 50, 'filled',  'MarkerFaceColor', map_con(3,:))
er = errorbar([1], mean(avg_dist_pos{2,1}),-err_con,err_con);    
er.Color = [0 0 0];                            
er.LineStyle = 'none';  
scatter([ones(size(avg_dist_neg{2,1},1),1)*2] , avg_dist_neg{2,1}, 50,'filled',  'MarkerFaceColor', map_con(3,:))
er = errorbar([2],mean(avg_dist_neg{2,1}),-err_con2,err_con2);    
er.Color = [0 0 0];                            
er.LineStyle = 'none'; 
for fish = 1:size(avg_dist_pos{2,1},1)
    plot([1 2], [avg_dist_pos{2,1}(fish),avg_dist_neg{2,1}(fish)], 'Color', 'k', 'LineStyle', '--' )

end
xticks([1 2])
xticklabels({'Pos', 'Neg'})
xlim([0 3])
hold off
title('Avg dist p val 0.001')
ylabel('Avg Distance ')
% legend([b3 b2], group_names{[1 2]})
% legend('Location', 'eastoutside')
% ylim([0 65])
saveas(gcf, fullfile(save_path, ['posandnefavgdistp_val0.001.png']))
saveas(gcf, fullfile(save_path, [ 'posandnefavgdistp_val0.001.svg']))
[p_con, h] = quick_statistic_signrank(avg_dist_pos{2,1}, avg_dist_neg{2,1})
save(fullfile(save_path, 'p_val_avgdist.mat'), 'p_con')
%% Anova for corr vs dist

test_corr = {all_corr_together{2,1}{1,1}'; all_corr_together{3,1}{1,1}'}; %all_corr_together{1,1}{1,1}

p_fac = cal_p_fac(steps(2:end), test_corr) % is sign, the first one is if the distance has an impact, the 2nd one is if the groups has a sign impact.

% now judt do this for all of them con
test_corr = {all_corr_together{1,1}{2,1}'; all_corr_together{2,1}{2,1}'; all_corr_together{3,1}{2,1}'};

p_fac = cal_p_fac(steps(2:end), test_corr) % is sign, the first one is if the distance has an impact, the 2nd one is if the groups has a sign impact.

test_corr = {all_corr_together{1,1}{3,1}'; all_corr_together{2,1}{3,1}'; all_corr_together{3,1}{3,1}'};

p_fac = cal_p_fac(steps(2:end), test_corr) % is sign, the first one is if the distance has an impact, the 2nd one is if the groups has a sign impact.

test_corr = {all_corr_together{1,1}{4,1}'; all_corr_together{2,1}{4,1}'; all_corr_together{3,1}{4,1}'};

p_fac = cal_p_fac(steps(2:end), test_corr) % is sign, the first one is if the distance has an impact, the 2nd one is if the groups has a sign impact.

%% Fig6
%% Interaction Analysis

general_ratio = cell(3,1);
general_ratio_resp = cell(3,1);
collected_dff = cell(3,1);
int_code_per_fish = cell(3,1);
interaction_code_coll = cell(3,1);
all_fractions = cell(3,1);
all_fract_fish = cell(3,1);
fraction_resp = cell(3,1);
dorsomed = 0;
plotting = 0; 
method = 4; % Method 1: peak and 2: mean 3: msc 4: int index
stim_period = diff_stim_period(1,:); 
all_fraction_enh = cell(3,1);
all_fraction_supr = cell(3,1);
for group = 1:no_group 

    [general_ratio{group,1}, general_ratio_resp{group,1}, collected_dff{group,1}, int_code_per_fish{group,1}, interaction_code_coll{group,1}, all_fractions{group,1}, all_fract_fish{group,1},fraction_resp{group,1}, all_fraction_enh{group,1}, all_fraction_supr{group,1}] = AO_interaction_analysis(dff_trialwise{group,1}, brain_regions{group,1}, positions{group,1}, resp_period_list{group,1}, no_con, ...
       con_trials, cmap2, save_path, group_names{group}, stim_period, dorsomed, plotting, thresh_neu{group,1},  method);

end


figure('units','centimeters','Position',[2 2 18 6])
plot_bar_three_cond(general_ratio{1,1}*100, general_ratio{2,1}*100, general_ratio{3,1}*100, map_bef, map_con, map_dru, {'Enhanced', 'Non', 'Supression'},group_names)
if method == 1
    title('Hab Cells Interactions PEAK All Cells 10s')
    % legend('Location', 'Eastoutside')
    ylabel('% Hb Cells')
    saveas(gcf, fullfile(save_path, ['Fig24_mixinteractio_withpeak_10s.png']))
    saveas(gcf, fullfile(save_path, ['Fig24_mixinteractio_withpeak_10s.svg']))
elseif method == 2
    title('Hb interaction with mean 10 s')
    ylabel('% Hb Cells')
    saveas(gcf, fullfile(save_path, ['Fig24_mixinteractio_withmean_only_10s.png']))
    saveas(gcf, fullfile(save_path, ['Fig24_mixinteractio_withmean_only_10s.svg']))
elseif method == 4
    title('Hb interaction with interact index 10 s')
    ylabel('% Hb Cells')
    saveas(gcf, fullfile(save_path, ['Fig24_mixinteractio_withIntIdx_10s.png']))
    saveas(gcf, fullfile(save_path, ['Fig24_mixinteractio_withIntIdx_10s.svg']))
end

[fig1, fig2] = split_plot_bar_three_cond(general_ratio{1,1}*100, general_ratio{2,1}*100, general_ratio{3,1}*100, map_bef, map_con, map_dru, {'Enhanced', 'Non', 'Supression'},group_names, 'Hb interaction with interact index 10 s', 10, 5)
figure(fig1)
ylim([0 70])
% xlim([0.5 2.6])

figure(fig2)
ylim([0 70])
% xlim([0.5 2.6])

saveas(fig1, fullfile(save_path, ['Fig24_mixinteractio_withIntIdx_10sHbCONCOM.svg']))
saveas(fig1, fullfile(save_path, ['Fig24_mixinteractio_withIntIdx_10sHbCONCOM.png']))
saveas(fig2, fullfile(save_path, ['Fig24_mixinteractio_withIntIdx_10sHbCONDRU.svg']))
saveas(fig2, fullfile(save_path, ['Fig24_mixinteractio_withIntIdx_10sHbCONDRU.png']))

p_values_all_cells1 = cell(3,1); %nan([3,3,3]); % 3 comparisions Wt-het, wt-hom, het-hom and 3 con: enh, no resp, suppr, 3 directions of rank sum
for clu = 1:3
    
    [p1, h] = quick_statistic(squeeze(general_ratio{1,1}(:,clu)), squeeze(general_ratio{2,1}(:,clu)))
    [p2, h] = quick_statistic(squeeze(general_ratio{1,1}(:,clu)), squeeze(general_ratio{3,1}(:,clu)))
    [p3, h] = quick_statistic(squeeze(general_ratio{2,1}(:,clu)), squeeze(general_ratio{3,1}(:,clu)))

     
    p_values_all_cells1{clu} = [p1;p2;p3]; 
end
if method == 1
    pval_interact.peak_p = p_values_all_cells1;

elseif method == 2
    pval_interact.mean_p = p_values_all_cells1;

end

sub_val = [1 2 3; 4 5 6; 7 8 9];
figure
for group = 1:3
    subplot(3,3,sub_val(group,1))
    hold on
    for con = 1:3
        plot(squeeze(mean(mean(collected_dff{group,1}(:,con_trials(con,:),find(interaction_code_coll{group,1} == 1)),2),3)))
    end
    legend({'L' , 'T', 'Lt'})
    title('Enh')
    ylabel(group_names{1,group})

    subplot(3,3,sub_val(group,2))
    hold on
    for con = 1:3
        plot(squeeze(mean(mean(collected_dff{group,1}(:,con_trials(con,:),find(interaction_code_coll{group,1} == 0)),2),3)))
    end
    legend({'L' , 'T', 'Lt'})
    title('Non')
    subplot(3,3,sub_val(group,3))
    hold on
    for con = 1:3
        plot(squeeze(mean(mean(collected_dff{group,1}(:,con_trials(con,:),find(interaction_code_coll{group,1} == -1)),2),3)))
    end
    legend({'L' , 'T', 'Lt'})
    title('Supr')
end


figure('units','centimeters','Position',[2 2 12 12])
hold on
h2 = cdfplot(all_fractions{2,1})
h2.Color = map_con(2,:); 
h2.LineWidth = 2;
h3 = cdfplot(all_fractions{1,1})
h3.Color = map_bef(2,:); 
h3.LineWidth = 2;
h1 = cdfplot(all_fractions{3,1})
h1.Color = map_dru(2,:); 
h1.LineWidth = 2; 
legend([h1 h2 h3], {'Dru', 'Con', 'Bef'}, 'Location', 'Southeast')
% title('Magnitude of Enhancement')
xlabel('Int Index')
title('II = (LT - max(L,T) /max(L,t)  * 100')
ylim([-0.1 1.1])
xlim([-300 300])
saveas(gcf, fullfile(save_path, ['IntIdx.png']))
saveas(gcf, fullfile(save_path, ['IntIdx.svg']))

[h1, p1] = kstest2(all_fractions{1,1},all_fractions{2,1})
[h2, p2] = kstest2(all_fractions{1,1},all_fractions{3,1})
[h3, p3] = kstest2(all_fractions{2,1},all_fractions{3,1})