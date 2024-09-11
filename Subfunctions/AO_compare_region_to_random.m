%% statitics for CPPG data

%%
function [reacting_cells_per_fish, reacting_rand_per_fish] = AO_compare_region_to_random(dff_fish, positions, brain_regions, drug_onsets, baseline_time)

uni_brain_regions = [1, 2, 3, 4, 5, 6, 7, 9, 11]; 
reacting_cells_per_fish = zeros(size(dff_fish,2), length(uni_brain_regions));
reacting_rand_per_fish = zeros(size(dff_fish,2), length(uni_brain_regions));

for fish = 1:size(dff_fish,2)
   
    
    dff_all_cells = dff_fish{1,fish};
    
    % now instead of z score do the 2 std from the baseline 
    baseline_start = baseline_time{1,fish}(1); 
    baseline_stop = baseline_time{1,fish}(2); 
    drug_start = 1900; %drug_onsets{1,fish}(1); 1900:2500
    drug_stop = 2500;% drug_onsets{1,fish}(2); 
    % results.NEWindex(find(results.NEWindex==12)) = 11; 
    % results.NEWindex(find(results.NEWindex==14)) = 11; 
    
    % list for storing all the pos/neg percentage of reacting cells for each
    % region
    all_positive_perc_cells_reacting = [];
    all_negative_perc_cells_reacting = [];
    
    % the same, but for the random cells 
    all_positive_ran_perc = [];
    all_negative_ran_perc = [];
    
    results.redoneIndex = brain_regions{1,fish}; 
    for cu_region=1:length(uni_brain_regions)
        region = uni_brain_regions(cu_region);
    % the region of the loop
        baseline = mean(dff_all_cells(find(results.redoneIndex==region),baseline_start:baseline_stop),2);
        std_baseline = std(dff_all_cells(find(results.redoneIndex==region),baseline_start:baseline_stop), 0 ,2);
        drug = mean(dff_all_cells(find(results.redoneIndex==region),drug_start:drug_stop),2);
    
    % here I am trying to see if the cell is increasing or decreasing in
    % activity 
        perc_cells_reacting_pos = []; 
        perc_cells_reacting_neg = []; 
    % 
    %     for cell=1:length(drug)
    %         if std_baseline(cell)*3 < drug(cell)
    %             perc_cells_reacting_pos = [perc_cells_reacting_pos, cell];
    %         elseif -std_baseline(cell)*3 > drug(cell)
    %             perc_cells_reacting_neg = [perc_cells_reacting_neg, cell]; 
    % %         else
    % %             continue
    %         end 
    %     end
        for cell=1:length(drug)
            if baseline(cell)+std_baseline(cell)*2 < drug(cell)
                perc_cells_reacting_pos = [perc_cells_reacting_pos, cell];
            elseif baseline(cell)-std_baseline(cell)*2 > drug(cell)
                perc_cells_reacting_neg = [perc_cells_reacting_neg, cell]; 
    %         else
    %             continue
            end 
        end
    
       % now I am calculating the percentage of cells that are reacting either
       % with increase (positive) or decrease (neg)
        cells_reacting_pos = length(perc_cells_reacting_pos);
        percentage_ind_region_pos = cells_reacting_pos/length(drug);
        all_positive_perc_cells_reacting = [all_positive_perc_cells_reacting, percentage_ind_region_pos]; 
        
        cells_reacting_neg = length(perc_cells_reacting_neg);
        percentage_ind_region_neg = cells_reacting_neg/length(drug);
        all_negative_perc_cells_reacting = [all_negative_perc_cells_reacting, percentage_ind_region_neg];
        % 1: Dm, 2:Dl, 3:Dlv, 4:Dc, 5:Dd, 6:DMp, 7:Dp, 8:OB, 9:Vd, 10 Vv, 11:Hab, 12:dHab, 13:vHab, 14:lat Hab, 15:Midbrain
        %for the random ones... I have to draw it like 100 times and then take
        %the average of that I guess...
        rand_cycles = 100;
        cells_region_excluded = dff_all_cells; %dff_all_cells(find(results.redoneIndex~=region),:); %
        no_cells_of_region = length(drug); 
        perc_of_rand_pos = []; 
        perc_of_rand_neg = [];
        for cycle=1:rand_cycles
            %first I need to make the rand ind
            cycle_ind = randi([1 size(cells_region_excluded,1)],1,no_cells_of_region);
            rand_base = mean(cells_region_excluded(cycle_ind,baseline_start:baseline_stop),2);
            rand_base_std = std(cells_region_excluded(cycle_ind, baseline_start:baseline_stop),0, 2);
            rand_drug = mean(cells_region_excluded(cycle_ind,drug_start:drug_stop),2);
            positive_reacts = [];
            negative_racts = [];
            %now I am looking which cells are po/neg repsonding 
            for cell=1:length(rand_drug)
                if rand_base(cell) + rand_base_std(cell)*2 < rand_drug(cell)
                    positive_reacts = [positive_reacts, cell];
                elseif rand_base(cell)-rand_base_std(cell)*2 > rand_drug(cell)
                    negative_racts = [negative_racts, cell]; 
                end 
            end
            %now I calculate the % for the rand cycle we are in
            ran_reacting_pos = length(positive_reacts);
            percentage_ind_region_pos_ran = ran_reacting_pos/length(rand_drug);
            perc_of_rand_pos = [perc_of_rand_pos, percentage_ind_region_pos_ran]; 
    
            ran_reacting_neg = length(negative_racts);
            percentage_ind_region_neg_ran = ran_reacting_neg/length(rand_drug);
            perc_of_rand_neg = [perc_of_rand_neg, percentage_ind_region_neg_ran];
        end
        % calculate the avergae of that 100 randoms... 
        avg_ran_perc_pos = mean(perc_of_rand_pos);
        avg_ran_perc_neg = mean(perc_of_rand_neg); 
        
        % and add that to our overall list woop
        all_positive_ran_perc = [all_positive_ran_perc, avg_ran_perc_pos ];
        all_negative_ran_perc = [all_negative_ran_perc, avg_ran_perc_neg ];
    end
    reacting_cells_per_fish(fish, :) = all_positive_perc_cells_reacting; 
    reacting_rand_per_fish(fish, :) = all_positive_ran_perc; 
end
        
end   

