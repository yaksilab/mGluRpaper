function [general_ratio, general_ratio_resp, collected_dff, int_code_per_fish, interaction_code_coll, all_fractions, all_fract_fish, fraction_resp, all_fraction_enh, all_fraction_supr] = AO_interaction_analysis(all_dfftrialw_wt,brain_regions_wt, position_wt, resp_list_wt_dff_std, no_con, ...
    con_trials, cmap2, save_path, groupname, stim_period, dorsomed, plotting, thresh_dormedneurons_wt,  method)

% In this function I will take a look at the interation between the two
% stimulus conditions using either method 1 : peak detection or method 2:
% mean over a period 

% some output variables 

int_code_per_fish = {};
general_ratio = nan([size(all_dfftrialw_wt,2), 3]);
general_ratio_resp = nan([size(all_dfftrialw_wt,2), 3]);
collected_dff = nan([size(all_dfftrialw_wt{1,1},1),size(all_dfftrialw_wt{1,1},2),1]);
all_fractions = []; 
all_fract_fish= {};
fraction_resp = [];
all_fraction_enh = [];
all_fraction_supr = [];

for fish = 1:size(all_dfftrialw_wt,2)

    brainnumber = 11; 

    current_dff = all_dfftrialw_wt{1,fish};
    % current_dff = all_dffs_wt{1,fish};
    current_hab = current_dff(:,:,find(brain_regions_wt{1,fish} == brainnumber));
    % current_hab = current_dff;%(:,:,find(brain_regions_wt{1,fish} == brainnumber)); 

    % current_hab = current_dff(find(brain_regions_wt{1,fish} == brainnumber),:); 

    hab_positions = position_wt{1,fish}(find(brain_regions_wt{1,fish} == brainnumber),:);
    try
        hab_resp = resp_list_wt_dff_std{1,1}{1,fish}(find(brain_regions_wt{1,fish} == brainnumber),:);
    catch
        hab_resp = resp_list_wt_dff_std{1,1}{1,fish};

    end
    % now I want to only focus on the cells on the top two planes 
    if dorsomed
        dorsal_hb_cells = find(thresh_dormedneurons_wt{1,fish} ==1); %find(hab_positions(:,5)== 1); %| hab_positions(:,5)== 2 
        current_hab = current_hab(:,:,dorsal_hb_cells);
        hab_positions = hab_positions(dorsal_hb_cells,:); 
        hab_resp = hab_resp(dorsal_hb_cells,:);
    end
    
    % here is now where we splt in the two methods 

    if method == 1

        cell_avg_per_con = nan([size(hab_positions,1), 3]); 
        loc_of_peak = nan([size(hab_positions,1), 3]); 
        special_neurons = []; 
        
        for neuron = 1:size(current_hab,3)
            for con = 1:no_con
    %             cell_avg_per_con(neuron,con) = mean(mean(current_hab(stim_period, con_trials(con,:), neuron),1),2);
                try
                    % find peaks in a bigger window and then make sure they are
                    % in the stimulus window
                    new_stim_p = stim_period; %[stim_period, 34:40];
                    % [temp_avg, temp_loc] = findpeaks(abs(mean(current_hab(new_stim_p, con_trials(con,:), neuron),2)));
                     [temp_avg, temp_loc] = findpeaks((mean(current_hab(new_stim_p, con_trials(con,:), neuron),2)));
                    
                    real_temp_av = mean(current_hab(new_stim_p(temp_loc), con_trials(con,:), neuron),2); 
                    [max_val, loc_max] = max(temp_avg);
                    if ~any(diff(sign(real_temp_av(real_temp_av~=0))))
                        if temp_loc(loc_max) <= length(stim_period)
                            real_max_avg =  real_temp_av(loc_max);
                            cell_avg_per_con(neuron,con) = real_max_avg; 
                            loc_of_peak(neuron,con) = stim_period(temp_loc(loc_max)); 
                        else 
                            new_loc = find(temp_loc < length(stim_period)); 
                            [max_val, loc_max] = max(temp_avg(new_loc));
                            real_max_avg =  real_temp_av(loc_max);
                            cell_avg_per_con(neuron,con) = real_max_avg; 
                            loc_of_peak(neuron,con) = stim_period(temp_loc(loc_max)); 
                        end
                    else
                        if temp_loc(loc_max) <= length(stim_period)
                            real_max_avg =  real_temp_av(loc_max);
                            cell_avg_per_con(neuron,con) = real_max_avg; 
                            loc_of_peak(neuron,con) = stim_period(temp_loc(loc_max)); 
                        else 
                            new_loc = find(temp_loc < length(stim_period)); 
                            [max_val, loc_max] = max(temp_avg(new_loc));
                            real_max_avg =  real_temp_av(loc_max);
                            cell_avg_per_con(neuron,con) = real_max_avg; 
                            loc_of_peak(neuron,con) = stim_period(temp_loc(loc_max)); 
                        end
                        % cur_signs =sign(real_temp_av(real_temp_av~=0)); 
                        % if cur_signs(1) < 0 & mean(current_hab(stim_period(1)-1, con_trials(3,:), neuron),2) > real_temp_av(1) & con == 3 & hab_resp(neuron,3) == 1
                        %     special_neurons = [special_neurons; neuron];
                        %     disp('im special')
                        % end
                    end
    
    
                catch
                    disp(['no peak elongatng the stim period neuron: ', num2str(neuron), ' and con ', num2str(con)])
                    [temp_avg, temp_loc] = findpeaks(mean(current_hab([stim_period, 34:40], con_trials(con,:), neuron),2));
                    [max_val, loc_max] = max(temp_avg);
                    try
                    cell_avg_per_con(neuron,con) = max_val; 
                    new_stim_p = [stim_period, 34:40]; 
                    loc_of_peak(neuron,con) = new_stim_p(temp_loc(loc_max)); 
                    catch
                        disp('No peak maybe neg cell')
                    end
    
                end
                
    %             cell_avg_per_con(neuron,con) = max_val; 
    %             loc_of_peak(neuron,con) = stim_period(temp_loc(loc_max)); 
            end       
        end

    % now we see if the cells are neutral, synergic or supressed
        interact_code = nan(size(current_hab,3),1);
        fraction_of_cell = nan(size(current_hab,3),1);
        fraction_enh = [];
        fraction_supr = [];
        for neuron = 1:size(current_hab,3)
            r_a_r_b = cell_avg_per_con(neuron,1) + cell_avg_per_con(neuron,2); 
            r_ab = cell_avg_per_con(neuron,3); 
  
            fraction_of_cell(neuron) = r_a_r_b-r_ab;  % switched out ratior_ab/r_a_r_b; 
    %         end
            % if cell_avg_per_con(neuron,1) > 0 && cell_avg_per_con(neuron,2) > 0
            if r_a_r_b < r_ab && max([cell_avg_per_con(neuron,1),cell_avg_per_con(neuron,2)])< r_ab %&&  ~any(ismember(sign(cell_avg_per_con(neuron,[1 2])), -1)) %~all(ismember(sign(cell_avg_per_con(neuron,[1 2])), -1))%
                interact_code(neuron) = 1; 
                fraction_enh = [fraction_enh; r_ab- r_a_r_b;];
    %         elseif r_a_r_b < r_ab && sign(cell_avg_per_con(neuron,1)) ~= sign(cell_avg_per_con(neuron,2)) && max([cell_avg_per_con(neuron,1),cell_avg_per_con(neuron,2)]) < r_ab
    %             interaction_code(neuron) = 6;
    %         elseif r_a_r_b < r_ab && all(ismember(sign(cell_avg_per_con(neuron,[1 2])), -1))
    %             interaction_code(neuron) = 8;
            elseif  max([cell_avg_per_con(neuron,1),cell_avg_per_con(neuron,2)])> r_ab %&& ~all(ismember(sign(cell_avg_per_con(neuron,[1 2])), -1))
                interact_code(neuron) = -1;
                fraction_supr = [fraction_supr; max([cell_avg_per_con(neuron,1),cell_avg_per_con(neuron,2)])-r_ab];
            else 
                interact_code(neuron) = 0;
            end
            % else
            % 
            %     interaction_code(neuron) = 0;
            % end
        end

        all_fraction_enh = [all_fraction_enh; fraction_enh];
        all_fraction_supr = [all_fraction_supr; fraction_supr];
        int_code_per_fish{fish,1} = interact_code;
        avg_per_fish{fish} = cell_avg_per_con;
        collected_dff = cat(3,collected_dff, current_hab);
        all_fractions = [all_fractions; fraction_of_cell]; % (find(interaction_code == 1))
        all_fract_fish{fish,1} = fraction_of_cell; %(find(interaction_code == 1)); 
    
        if fish ==1
            interaction_code_coll = interact_code;
            
        else 
            interaction_code_coll = [interaction_code_coll ;interact_code]; 
            
            
        end
        general_ratio(fish,:) = [length(find(interact_code == 1))/length(interact_code); length(find(interact_code == 0))/length(interact_code); length(find(interact_code == -1))/length(interact_code)]; 

        general_ratio_resp(fish,:) = [length(find(interact_code(find(hab_resp(:,3)==1)) == 1))/length(interact_code(find(hab_resp(:,3)==1))); length(find(interact_code(find(hab_resp(:,3)==1)) == 0))/length(interact_code(find(hab_resp(:,3)==1))); length(find(interact_code(find(hab_resp(:,3)==1)) == -1))/length(interact_code(find(hab_resp(:,3)==1)))]; 


    elseif method == 2

        % so first lets check if there is suppression of enhancement
        interact_code = nan(size(current_hab,3),1); %nan(size(last_con_resp,1),1);
        fraction_of_cell = nan(size(current_hab,3),1);
        fraction_enh = [];
        fraction_supr = [];
        for neuron = 1:size(current_hab,3)%size(last_con_resp,1)
            % lightap_avg = mean(mean(current_hab(stim_period, con_trials(3,:), last_con_resp(neuron)),2),1);
            % ligh_avg = mean(mean(current_hab(stim_period, con_trials(1,:), last_con_resp(neuron)),2),1);
            % tap_avg = mean(mean(current_hab(stim_period, con_trials(2,:), last_con_resp(neuron)),2),1);
    
            lightap_avg = mean(mean(current_hab(stim_period, con_trials(3,:), neuron),2),1);
            ligh_avg = mean(mean(current_hab(stim_period, con_trials(1,:), neuron),2),1);
            tap_avg = mean(mean(current_hab(stim_period, con_trials(2,:), neuron),2),1);
    
            r_a_r_b = ligh_avg + tap_avg; 
            r_ab = lightap_avg; 
    
            % if ~any(ismember(sign(cell_avg_per_con(neuron,[1 2])), -1))
                fraction_of_cell(neuron) = r_a_r_b-r_ab;  % switched out ratior_ab/r_a_r_b; 
            % end
    
                 % if ligh_avg < 0 && tap_avg < 0 && lightap_avg < 0
            if r_a_r_b < r_ab && max([ligh_avg, tap_avg]) < r_ab %&& ligh_avg > 0 && tap_avg > 0  here I am checking if the sum is smaller than and also the individual ones are smaller just in case there is a neg
                interact_code(neuron) = 1;
                fraction_enh = [fraction_enh; r_a_r_b-r_ab;];
    
            elseif  r_a_r_b > r_ab  %max([ligh_avg, tap_avg]) > r_ab % here I now check for suppression
                interact_code(neuron) = -1;
                % fraction_supr = [fraction_supr; max([ligh_avg, tap_avg])-r_ab];
                fraction_supr = [fraction_supr; r_a_r_b-r_ab];

            else 
                 interact_code(neuron) = 0;
            end
                 % else
                 %     interact_code(neuron) = 0;
                 % end
             % fraction_supr = [fraction_supr; max([ligh_avg, tap_avg])-r_ab];
             % fraction_enh = [fraction_enh; r_a_r_b-r_ab;];
        end
        all_fraction_enh = [all_fraction_enh; fraction_enh];
        all_fraction_supr = [all_fraction_supr; fraction_supr];
        all_fractions = [all_fractions; fraction_of_cell]; % (find(interaction_code == 1))
        all_fract_fish{fish,1} = fraction_of_cell; %(find(interaction_code == 1));
        fraction_resp = [fraction_resp; fraction_of_cell(find(hab_resp(:,3)==1))];
        general_ratio(fish,:) = [length(find(interact_code == 1))/length(interact_code); length(find(interact_code == 0))/length(interact_code); length(find(interact_code == -1))/length(interact_code)]; 
        general_ratio_resp(fish,:) = [length(find(interact_code(find(hab_resp(:,3)==1)) == 1))/length(interact_code(find(hab_resp(:,3)==1))); length(find(interact_code(find(hab_resp(:,3)==1)) == 0))/length(interact_code(find(hab_resp(:,3)==1))); length(find(interact_code(find(hab_resp(:,3)==1)) == -1))/length(interact_code(find(hab_resp(:,3)==1)))]; 

        int_code_per_fish{fish,1} = interact_code;
        % interact_code_fish{fish,1} = interact_code; 
        collected_dff = cat(3,collected_dff, current_hab);
        if fish ==1
            interaction_code_coll = interact_code;
            
        else 
            interaction_code_coll = [interaction_code_coll ;interact_code]; 
            
            
        end

    elseif method == 3
        % using the mean statistical contrast from the paper "Identifying
        % and Quantifying Multisensory Integration: A Tutorial Review"
        % (Stevenson et al 2014) and these Perrault et al. 2003, 2005; Stanford et al. 2005
        interact_code = nan(size(current_hab,3),1); %nan(size(last_con_resp,1),1);
        fraction_of_cell = nan(size(current_hab,3),1);
        for neuron = 1:size(current_hab,3)%size(last_con_resp,1)
            % lightap_avg = mean(mean(current_hab(stim_period, con_trials(3,:), last_con_resp(neuron)),2),1);
            % ligh_avg = mean(mean(current_hab(stim_period, con_trials(1,:), last_con_resp(neuron)),2),1);
            % tap_avg = mean(mean(current_hab(stim_period, con_trials(2,:), last_con_resp(neuron)),2),1);
    
            % lightap_avg = mean(mean(current_hab(stim_period, con_trials(3,:), neuron),2),1);
            % ligh_avg = mean(mean(current_hab(stim_period, con_trials(1,:), neuron),2),1);
            % tap_avg = mean(mean(current_hab(stim_period, con_trials(2,:), neuron),2),1);
            sum_overalltrials = 0;
            for trial = 1:8 %no trials 
                curr_trial_val = mean(current_hab(stim_period, con_trials(3,trial), neuron),1) - (mean(current_hab(stim_period, con_trials(1,trial), neuron),1) +mean(current_hab(stim_period, con_trials(2,trial), neuron),1));
                sum_overalltrials = sum_overalltrials + curr_trial_val;
    
            end
            msc = sum_overalltrials/8;

            %subadditive (msc < 0), additive (msc = 0), and superadditive (msc > 0) 
            if msc < 0
                interact_code(neuron) = -1; 
            elseif msc > 0 
                interact_code(neuron) = 1;
            elseif msc == 0
                interact_code(neuron) = 0;

            end
            fraction_of_cell(neuron) = msc;
        end
        all_fractions = [all_fractions; fraction_of_cell]; % (find(interaction_code == 1))
        all_fract_fish{fish,1} = fraction_of_cell; %(find(interaction_code == 1));

        general_ratio(fish,:) = [length(find(interact_code == 1))/length(interact_code); length(find(interact_code == 0))/length(interact_code); length(find(interact_code == -1))/length(interact_code)]; 
        general_ratio_resp(fish,:) = [length(find(interact_code(find(hab_resp(:,3)==1)) == 1))/length(interact_code(find(hab_resp(:,3)==1))); length(find(interact_code(find(hab_resp(:,3)==1)) == 0))/length(interact_code(find(hab_resp(:,3)==1))); length(find(interact_code(find(hab_resp(:,3)==1)) == -1))/length(interact_code(find(hab_resp(:,3)==1)))]; 

        int_code_per_fish{fish,1} = interact_code;
        % interact_code_fish{fish,1} = interact_code; 
        collected_dff = cat(3,collected_dff, current_hab);
        if fish ==1
            interaction_code_coll = interact_code;
            
        else 
            interaction_code_coll = [interaction_code_coll ;interact_code]; 
            
            
        end
    elseif method == 4
        % lets do the interactive index instead basically incompassing the
        % ratio  (Meredith and Stein 1983, 1986b). and that review

        interact_code = nan(size(current_hab,3),1); %nan(size(last_con_resp,1),1);
        fraction_of_cell = nan(size(current_hab,3),1);
        for neuron = 1:size(current_hab,3)%size(last_con_resp,1)
            % lightap_avg = mean(mean(current_hab(stim_period, con_trials(3,:), last_con_resp(neuron)),2),1);
            % ligh_avg = mean(mean(current_hab(stim_period, con_trials(1,:), last_con_resp(neuron)),2),1);
            % tap_avg = mean(mean(current_hab(stim_period, con_trials(2,:), last_con_resp(neuron)),2),1);
    
            lightap_avg = mean(mean(current_hab(stim_period, con_trials(3,:), neuron),2),1);
            ligh_avg = mean(mean(current_hab(stim_period, con_trials(1,:), neuron),2),1);
            tap_avg = mean(mean(current_hab(stim_period, con_trials(2,:), neuron),2),1);
            % if tap_avg >0 && ligh_avg >0 && lightap_avg >0
            int_idx = (lightap_avg - max([ligh_avg, tap_avg]))*100 / max([ligh_avg, tap_avg]);
            % int_idx = int_idx * 100; 
            % if abs(int_idx) > 1000
            %     fraction_of_cell(neuron) = NaN; %int_idx;%
            %      interact_code(neuron) = NaN;
            %     disp('Now')
            % else

            %subadditive (int_idx < 0), additive (int_idx  betwee 0 and 100), and superadditive (msc > 100) 
                if int_idx< 0 %-100 < int_idx &&  
                    interact_code(neuron) = -1; 
                elseif int_idx > 100
                    interact_code(neuron) = 1;
                else
                    interact_code(neuron) = 0;
    
                end
                 fraction_of_cell(neuron) = int_idx;
                % else 
                %     interact_code(neuron) = 0;
                %     fraction_of_cell(neuron) = nan;
                % end
               
            end
        % end
        all_fractions = [all_fractions; fraction_of_cell]; % (find(interaction_code == 1))
        all_fract_fish{fish,1} = fraction_of_cell; %(find(interaction_code == 1));
        % general_ratio(fish,:) = [length(find(interact_code == 1))/length(interact_code); length(find(interact_code == 0))/length(interact_code); length(find(interact_code == -1))/length(interact_code)]; 
        general_ratio(fish,:) = [length(find(interact_code == 1))/length(find(~isnan(interact_code))); length(find(interact_code == 0))/length(find(~isnan(interact_code))); length(find(interact_code == -1))/length(find(~isnan(interact_code)))]; 
        general_ratio_resp(fish,:) = [length(find(interact_code(find(hab_resp(:,3)==1)) == 1))/length(interact_code(find(hab_resp(:,3)==1))); length(find(interact_code(find(hab_resp(:,3)==1)) == 0))/length(interact_code(find(hab_resp(:,3)==1))); length(find(interact_code(find(hab_resp(:,3)==1)) == -1))/length(interact_code(find(hab_resp(:,3)==1)))]; 

        int_code_per_fish{fish,1} = interact_code;
        % interact_code_fish{fish,1} = interact_code; 
        collected_dff = cat(3,collected_dff, current_hab);

        if fish ==1
            interaction_code_coll = interact_code;
            
        else 
            interaction_code_coll = [interaction_code_coll ;interact_code]; 
            
            
        end


    end
end



if plotting 
    figure('units','centimeters','Position',[2 2 10 12])
    
    scatter3(hab_positions(:,1), hab_positions(:,2), hab_positions(:,3), 25, interact_code, 'filled')
    colormap(cmap2)
    colorbar
    set(gca, 'Zdir', 'reverse'), set(gca, 'Ydir', 'reverse'), set(gca, 'Xdir', 'reverse')
    view (22,80) 
    axis equal
    title([groupname,'_Fish', num2str(fish) ' All cells'])

end
% end
collected_dff = collected_dff(:,:,2:end); 
end