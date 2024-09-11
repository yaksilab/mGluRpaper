%% Fig 3

%% CPPG 

% Loading the data 
load('X:\anna\Manuscript\FigData\Fig3\Fig3_cppg_2p_exvivo.mat')
save_path = 'X:\anna\Manuscript\FigData\Fig3';
% Now we can then load this way: 
dff_fish = cppg_2p_exvivo.dff_traces ; 
positions = cppg_2p_exvivo.position  ; 
brain_regions = cppg_2p_exvivo.brain_regions ;
drug_onsets = cppg_2p_exvivo.drug_onsets ;
baseline_time = cppg_2p_exvivo.baseline_time ;

[reacting_cells_per_fish, reacting_rand_per_fish] = AO_compare_region_to_random(dff_fish, positions, brain_regions, drug_onsets, baseline_time);
uni_brain_regions = [1, 2, 3, 4, 5, 6, 7, 9, 11]; 
p_valu_reg_to_ran = nan(size(reacting_cells_per_fish,2), 3); 
for reg = 1:size(reacting_cells_per_fish,2)
    cur_reg = organisation(reg)
    [p, h] = quick_statistic_signrank(reacting_cells_per_fish(:,cur_reg), reacting_rand_per_fish(:,cur_reg))
    p_valu_reg_to_ran(reg,:) = p;

end
[p, h] = quick_statistic_signrank(reacting_cells_per_fish(:,7), reacting_cells_per_fish(:,9)) % testting OB vs Hb

% PLotting the perc responding cells
cmap = [0.85, 0.85, 0.85; ... %grey for 0    
 
      0.4660, 0.6740, 0.1880; ...  %  Greenfor 1
      0, 0.4470, 0.7410; ... %  Blue for 2
  0.3010, 0.7450, 0.9330; ...       % Light Blue for 3
  0.8500, 0.3250, 0.0980; ... % orange for 4
  1, 0, 1; ... % magenta for 5
  0, 0.6, 0.3;... %dark green 6
  0.9290, 0.6940, 0.1250; ...%yellow for 7
  0.5430 0 0; ... % dark red 8
  0.8,0.3,0.3;... %indian red 9
  0.98, 0.122, 0.157; ... %sky blue 10
  0.94, 0.60, 0.108; ... %purple 11
  0.133, 0.128, 0.177; ... %lavendel 12
  0.225, 0. , 0.102 ;...
  0.204, 0.255, 0.255; ...
  0.153, 0.153, 0. ...
  ]  ;
names = categorical({'Dm','Dl', 'Dlv', 'Dc', 'Dd', 'Dmp', 'Dp', 'OB', 'Vd', 'Vv', 'Hab', 'dHab', 'otherHab?', 'latHab', 'Midbrain' });
names = reordercats(names,{'Dm','Dl', 'Dlv', 'Dc', 'Dd', 'Dmp', 'Dp', 'OB', 'Vd', 'Vv', 'Hab', 'dHab', 'otherHab?', 'latHab', 'Midbrain' });
regions_int = [1 2 3 4 5 6 7 8 9 11];  %regions_int = [1 2 3 4 5 6 7 8 9 10 11]; 
% I used different identifiers for the brain regions: 
% % 1: Dm, 2:Dl, 3:Dlv, 4:Dc, 5:Dd, 6:DMp, 7:Dp, 8:OB, 9:Vd, 10 Vv, 11:Hab, 12:dHab, 13:vHab, 14:lat Hab, 15:Midbrain

regions_int = [8 2 3 4 1 5 6 7 9 11];  %regions_int = [1 2 3 4 5 6 7 8 9 10 11]; 

uni_brain_regions = [1, 2, 3, 4, 5, 6, 7, 9, 11]; 
brain_names = {'Dl', 'Da', 'Dd', 'Dm', 'Dp', 'Vd', 'OB', 'Dc', 'Hab'}; 

organisation = [7 2 1 8 4 3 5 6 9]
new_names = {};
for na = 1:9
    cur_reg = organisation(na);
    new_names{na} = brain_names{1, cur_reg}

end

regions_int = organisation; 
figure('units','centimeters','Position',[2 2 12 10])
% first we plot the random line 
R = plot(nanmean(reacting_rand_per_fish(:,regions_int),1)*100, 'Color', 'k')
ylim([0 100])
xlim([0 10])
hold on 
% plot(nanmean(reacting_cells_per_fish(:, regions_int),1)*100, 'Color', cmap(7,:))
H = shadedErrorBar([1:9],nanmean(reacting_cells_per_fish(:, regions_int),1)*100,squeeze(nanstd(reacting_cells_per_fish(:,regions_int),1)*100/sqrt(size(reacting_cells_per_fish,1))), 'lineProps','k')
set(H.mainLine, 'Color', cmap(7,:))

for region = 1:length(regions_int)
    
    x = ones(1,size(reacting_cells_per_fish,1))*region;
    s = scatter(x, reacting_cells_per_fish(:, regions_int(region))*100, 50, 'filled');
    s.LineWidth = 0.1;

end
for fish = 1:size(reacting_cells_per_fish,1)
    plot([1:9], reacting_cells_per_fish(fish, regions_int)*100, 'Color', cmap(1,:))
end
hold off
xticks([1:10])
xticklabels({char(brain_names(regions_int))})
title('CPPG % Responding cells')
ylabel('% Resp Cells per region')
legend([R H.mainLine], {'Shuffled', 'Cells in Region'}, 'Location','northwest')

saveas(gcf, fullfile(save_path, ['Affected_cells_CPPG.png']))
saveas(gcf, fullfile(save_path, ['Affected_cells_CPPG.svg']))

%% L-AP4
% Loading 
load('X:\anna\Manuscript\FigData\Fig3\Fig3_lap4_2p_exvivo.mat')
save_path = 'X:\anna\Manuscript\Data and matlab\FigData\Fig3';
% Now we can then load this way: 
fish_dff = lap4_2p_exvivo.dff_traces ; 
fish_pos = lap4_2p_exvivo.position  ; 
fish_brain_reg = lap4_2p_exvivo.brain_regions ;
fish_red_cells = lap4_2p_exvivo.redcells;
% baseline_time = lap4_2p_exvivo.baseline_time ;
baseline_per = lap4_2p_exvivo.baseline_drug_per(1,:);
drug_per= lap4_2p_exvivo.baseline_drug_per(2,:);

% Affected cells per region

baseline_per = [floor(600*fra_rat_two_pho)-400:floor(600*fra_rat_two_pho)-50]; % minute 8-10; 
drug_per = [4500:4500+length(baseline_per)-1]; 
figure('units','centimeters','Position',[2 2 32 14])
ind = [1 3; 2 4];
brain_regions = [1, 2, 3, 4, 5, 6, 7, 9, 11]; 
brain_names = {'Dl', 'Da', 'Dd', 'Dm', 'Dp', 'Vd', 'OB', 'Dc', 'Hab'}; 
affected_per_region = nan([6, length(brain_regions)]); 
affected_random_per_reg = nan([6, length(brain_regions)]); 
resp_idx_region = cell(length(brain_regions),6);

new_names = {};
for na = 1:9
    cur_reg = organisation(na);
    new_names{na} = brain_names{1, cur_reg}

end

for fish = 1:6    
    for cu_region = 1:size(brain_regions,2)
        region = brain_regions(cu_region);
        current_dff = fish_dff{1,fish};
        current_hab = current_dff(find(fish_brain_reg{1,fish} == region),:); 
        hab_positions = fish_pos{1, fish}(find(fish_brain_reg{1,fish} == region),:);; %position_ctrl{1,fish}(find(brain_regions_ctrl{1,fish} == 11),:);
        hab_period = current_hab ;  %current_hab(:,con_period{con});


        affected_neurons = nan([size(current_hab,1),1]); 
        for neuron = 1:size(current_hab,1)
            pre = mean(current_hab(neuron,baseline_per),2); 
            pre_std = std(current_hab(neuron,baseline_per),0,2); 
            drug = mean(current_hab(neuron,drug_per),2);

            if drug < pre-2*pre_std
                affected_neurons(neuron) = 1; 
            else
                affected_neurons(neuron) = 0;
            end
        end
        resp_idx_region{cu_region,fish} = affected_neurons; 
        subplot(2,3, fish)
        % now lets do ratios
        hold on
        b2 = bar([0.75]+(cu_region-1), [length(find(affected_neurons == 1))/length(affected_neurons)])
        b2.FaceColor = cmap(3,:);
        b2.FaceAlpha = 0.5;
        b3 = bar([1.25]+(cu_region-1), [length(find(affected_neurons == 0))/length(affected_neurons)])
        b3.FaceColor = red_map(1,:);
        b3.FaceAlpha = 0.5;
        legend([b2 b3], {'Affected', 'Not'})
        title(names{fish})
        ylabel('Ratio Affected')
        ylim([0 1])
        affected_per_region(fish, cu_region) = length(find(affected_neurons == 1))/length(affected_neurons); 
     xticks([1:9])
    xticklabels(brain_names)


        % now adding the random part 
        rand_cycles = 100;
        cells_region_excluded = current_dff; %current_dff(find(fish_brain_reg{1,fish}~=region),:); %
        no_cells_of_region = size(current_hab,1); 
        perc_of_rand_pos = []; 
        % perc_of_rand_neg = [];
        disp('Cycle')
        for cycle=1:rand_cycles
            %first I need to make the rand ind
            cycle_ind = randi([1 size(cells_region_excluded,1)],1,no_cells_of_region);
            rand_base_std = std(cells_region_excluded(cycle_ind, baseline_per),0, 2);
            rand_base = mean(cells_region_excluded(cycle_ind, baseline_per),2);
            rand_drug = mean(cells_region_excluded(cycle_ind,drug_per),2);
            positive_reacts = [];
            % negative_racts = [];
            %now I am looking which cells are po/neg repsonding 
            for neuron=1:length(rand_drug)
                if rand_base(neuron)-rand_base_std(neuron)*2 > rand_drug(neuron)
                    positive_reacts = [positive_reacts, neuron];
                % elseif -rand_base_std(cell)*3 > rand_drug(cell)
                %     negative_racts = [negative_racts, cell]; 
                end 
            end
            %now I calculate the % for the rand cycle we are in
            ran_reacting_pos = length(positive_reacts);
            percentage_ind_region_pos_ran = ran_reacting_pos/length(rand_drug);
            perc_of_rand_pos = [perc_of_rand_pos, percentage_ind_region_pos_ran]; 
            % 
            % ran_reacting_neg = length(negative_racts);
            % percentage_ind_region_neg_ran = ran_reacting_neg/length(rand_drug);
            % perc_of_rand_neg = [perc_of_rand_neg, percentage_ind_region_neg_ran];
        end
        % calculate the avergae of that 100 randoms... 
        affected_random_per_reg(fish, cu_region) = mean(perc_of_rand_pos);
        % avg_ran_perc_neg = mean(perc_of_rand_neg); 
    end
   
    
     
end
% saveas(gcf, fullfile(save_path, ['Fig8_otherregions.png']))
% saveas(gcf, fullfile(save_path, ['Fig8_otherregions.svg']))

organisation = [7 2 1 8 4 3 5 6 9]
%  new_names{na} = brain_names{1, cur_reg}
% over all the fish 
figure('units','centimeters','Position',[2 2 12 10])
% plot(nanmean(affected_per_region,1))
hold on
H1=shadedErrorBar([1:9],squeeze(nanmean(affected_per_region(:,organisation),1))*100,squeeze(nanstd(affected_per_region(:,organisation),1)*100/sqrt(size(affected_per_region,1))), 'lineProps','k');
H1.mainLine.Color = red_map(1,: ); %'#A23333';
H1.patch.FaceColor = red_map(1,: ); %'#A23333'; 
H1.mainLine.LineWidth = 2;
hold on
plot(nanmean(affected_random_per_reg(:,organisation),1)*100, 'Color', 'k')
xticks([1:length(brain_regions)]); 
xticklabels(new_names)
xlim([0 10])

p_valu_reg_to_ran = nan(size(brain_regions,2), 3); 
for reg = 1:size(brain_regions,2)
    cur_reg = organisation(reg) 
    scatter(ones(size(affected_per_region,1),1)*reg, affected_per_region(:,cur_reg)*100, 50, 'filled')
    [p, h] = quick_statistic_signrank(affected_per_region(:,cur_reg), affected_random_per_reg(:,cur_reg))
    p_valu_reg_to_ran(reg,:) = p; 
end
[p, h] = quick_statistic_signrank(affected_per_region(:,7), affected_per_region(:,9))
for fish = 1:size(affected_per_region,1)
    plot([1:9], affected_per_region(fish,organisation)*100, 'Color', [.7 .7 .7], 'LineStyle', '--')
end

title('LAP4 % Responding cells')
ylabel('% Resp Cells per region')
saveas(gcf, fullfile(save_path, ['Affected_cells_LAP4.png']))
saveas(gcf, fullfile(save_path, ['Affected_cells_LAP4.svg']))

% Comparing dao positive to rest of Hb
baseline_per = [floor(600*fra_rat_two_pho)-400:floor(600*fra_rat_two_pho)-50]; % minute 8-10; 
drug_per = [4500:4500+length(baseline_per)]; 
figure('units','centimeters','Position',[2 2 32 14])
ind = [1 3; 2 4];
collected_ratio = [];
collected_red_Rat = [];
for fish = 1:6    
    n_clusters = 5; 
    current_dff = fish_dff{1,fish};
    current_hab = current_dff(find(fish_brain_reg{1,fish} == 11),:); 
    hab_positions = fish_pos{1, fish}(find(fish_brain_reg{1,fish} == 11),:);; %position_ctrl{1,fish}(find(brain_regions_ctrl{1,fish} == 11),:);
    hab_period = current_hab ;  %current_hab(:,con_period{con});
    red_list = fish_red_cells{1,fish};
    hab_red = red_list(find(fish_brain_reg{1,fish} == 11),:);
    
    affected_neurons = nan([size(current_hab,1),1]); 
    for neuron = 1:size(current_hab,1)
        pre = mean(current_hab(neuron,baseline_per),2); 
        pre_std = std(current_hab(neuron,baseline_per),0,2); 
        drug = mean(current_hab(neuron,drug_per),2);
        
        if drug < pre-2*pre_std
            affected_neurons(neuron) = 1; 
        else
            affected_neurons(neuron) = 0;
        end
    end
    figure('units','centimeters','Position',[2 2 32 14])
    red_cell_list = affected_neurons(find(hab_red == 1)); 
    nonred_list = affected_neurons(find(hab_red == 0)); 
    % now lets do ratios
    b2 = bar([0.75 1.75], [length(find(red_cell_list == 1))/length(red_cell_list) length(find(red_cell_list == 0))/length(red_cell_list)])
    b2.FaceColor = red_map(2,:);
    b2.FaceAlpha = 0.5;
    hold on 
     b3 = bar([1.25 2.25], [length(find(nonred_list == 1))/length(nonred_list) length(find(nonred_list == 0))/length(nonred_list)])
    b3.FaceColor = red_map(1,:);
    b3.FaceAlpha = 0.5;
    xticks([1 2])
    xticklabels({'Affected', 'Not'})
%     legend([b2 b3])
    title(names{fish})
    ylabel('Ratio Affected')
    ylim([0 1])
%     saveas(gcf, fullfile(save_path, [num2str(fish), '_Fig4_Hbratioresp.png']))
%     saveas(gcf, fullfile(save_path, [num2str(fish), '_Fig4_Hbratioresp.svg']))
    close; 
    collected_ratio = [collected_ratio; [length(find(nonred_list == 1))/length(nonred_list) length(find(nonred_list == 0))/length(nonred_list)]];
    collected_red_Rat = [collected_red_Rat; [length(find(red_cell_list == 1))/length(red_cell_list) length(find(red_cell_list == 0))/length(red_cell_list)]];
    
    if fish == 1
        affected_dff = current_hab(find(nonred_list == 1),:);
    else
        affected_dff = cat(1,affected_dff, current_hab(find(nonred_list == 1),:));
    end
%     collected_ratio = [collected_ratio; [length(find(nonred_list == 1))/length(nonred_list) length(find(nonred_list == 0))/length(nonred_list)]];
%     collected_red_Rat = [collected_red_Rat; [length(find(red_cell_list == 1))/length(red_cell_list) length(find(red_cell_list == 0))/length(red_cell_list)]];
    
end
% saveas(gcf, fullfile(save_path, ['Fig4_Hbratioresp.png']))
% saveas(gcf, fullfile(save_path, ['Fig4_Hbratioresp.svg']))
figure('units','centimeters','Position',[2 2 6 10])
scatter(ones(size(collected_red_Rat,1)), collected_red_Rat(:,1)*100, 50, red_map(2,:), 'filled');
hold on
% scatter(ones(size(collected_red_Rat,1))*1.25, collected_red_Rat(:,2), 50, red_map(2,:), 'filled'); 

scatter(ones(size(collected_ratio,1))*2, collected_ratio(:,1)*100, 50, red_map(1,:), 'filled'); 
% scatter(ones(size(collected_ratio,1))*2.25, collected_ratio(:,2), 50, red_map(1,:), 'filled'); 
for fish = 1:size(collected_red_Rat,1)
    plot([1 2], [collected_red_Rat(fish,1)*100, collected_ratio(fish,1)*100], 'Color', [.7 .7 .7], 'LineStyle', '--')
end
xlim([0 3])
xticks([1 2])
xticklabels({'Red', 'Non Red'})
ylabel('Ratio of affected cells')
title('Treshold 2')
b = bar([1 2], [mean(collected_red_Rat(:,1))*100, mean(collected_ratio(:,1))*100] )
b.FaceColor = red_map(1,:);
b.FaceAlpha = 0.5;
SEM1 = std(collected_red_Rat(:,1), 0 ,1)*100/sqrt(size(collected_red_Rat(:,1),1));
SEM2 = std(collected_ratio(:,1), 0 ,1)*100/sqrt(size(collected_red_Rat(:,1),1));
er = errorbar([1 2],[nanmean(collected_red_Rat(:,1))*100, nanmean(collected_ratio(:,1))*100],[SEM1, SEM2]);    
er.Color = [0 0 0];                            
er.LineStyle = 'none'; 

% saveas(gcf, fullfile(save_path, ['Fig4_Hbratioresp2.png']))
% saveas(gcf, fullfile(save_path, ['Fig4_Hbratioresp2.svg']))
[p, h] = quick_statistic(collected_red_Rat(:,1), collected_ratio(:,1))
saveas(gcf, fullfile(save_path, ['DAOcells_LAP4.png']))
saveas(gcf, fullfile(save_path, ['DAOcells_LAP4.svg']))
