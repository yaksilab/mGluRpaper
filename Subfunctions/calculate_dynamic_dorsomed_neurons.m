function [dorsomed_index_fish, thresh_dormedneurons, norm_lat_med, norm_dor_ven, ven_neurons, dv_list] = calculate_dynamic_dorsomed_neurons(position_wt, brain_regions_wt, percentage_inclu, all_planes_y, groupname, save_path)
%calculate_dynamic_dorsomed_neurons - Normalise med-lat and dorso-ven axis and find the most dorsomed neurons (H1 line)
%   Author: Anna Maria Ostenrath
%
%   Syntax:
%       [dorsomed_index_fish, thresh_dormedneurons, norm_lat_med, norm_dor_ven, ven_neurons] = calculate_dynamic_dorsomed_neurons(position_wt, brain_regions_wt, percentage_inclu, all_planes_y, groupname, save_path)
%
%   Description:
%       function() - description
%    
%   Inputs:
%       position_wt - cell array with the positions for each
%       fish (1:x, 2:y, 3:z, 4:neuron index, 5:plane index)
%       brain_regions_wt - cell array with the brain region index for each
%       fish
%       percentage_inclu - double of how much to include (e.g 0.4 is 40%)
%       groupname - string with name for each group
%       save_path - string folder path for saving the figure
%
%   Outputs:
%       dorsomed_index_fish - cell array with a list for each fish
%       where each cell is ranked on the dorsomed scale
%       thresh_dormedneurons - cell array with a list for each fish
%       which cells as dorsomed and which not
%       norm_lat_med - cell array normalised lateral-medial position for each fish
%       norm_dor_ven - cell array normalised dorso-ventral position for each fish
%       ven_neurons - cell array with a list for each fish
%       dv_list - cell array with a list for each fish over all the cells
%       indicating the Hb cells that are dorsal or ventral
%       which cells as ventral and which not
%
%
%   Other m-files required: normalise_position_hb
%   Subfunctions: none
%   MAT-files required: none
%
%   See also: normalise_position_hb,  OTHER_FUNCTION_NAME2
%   Author: Anna Maria Ostenrath 
%   Date : November 2024

dorsomed_index_fish = {}; 
thresh_dormedneurons = {}; 
ven_neurons = {}; 
norm_lat_med = {}; 
norm_dor_ven = {}; 

brainnumber = 11;  % hardcoded for Hb
dv_list = {}; % this is a list for the dorsal-ventral neurons.
for fish = 1:size(position_wt,2) 
    % percentage_inclu = 0.20; 
    [norm_latmed, norm_dorven] = normalise_position_hb(position_wt{1,fish}(find(brain_regions_wt{1,fish} == brainnumber),:), all_planes_y); 
    close;
    hb_position = position_wt{1,fish}(find(brain_regions_wt{1,fish} == 11),:); 
    % now we need to find those cells that are dorsomedial ... 
    % medial are the cells close to 0... so it would be the lowest absolute
    % value cells 
    % dorsal would also be the lowest... 
    % can i add them up 
    rank_value = sum([abs(norm_latmed), norm_dorven], 2);
    [~,dorsomed_ind] = sort(rank_value, 'ascend'); 
    treshold = floor(length(dorsomed_ind) * percentage_inclu) ; 
    listcol = zeros(size(hb_position,1),1);
    listcol(dorsomed_ind(1:treshold)) = 1; 
    figure('units','centimeters','Position',[2 2 30 10])
    scatter3(hb_position(:,1), hb_position(:,2), hb_position(:,3), 80, listcol, 'filled')
    % colormap(cmap2)
    % colorbar
    set(gca, 'Zdir', 'reverse'), set(gca, 'Ydir', 'reverse'), set(gca, 'Xdir', 'reverse')
    view (80,7) ;%view (22,80) 
    axis equal
    title([groupname, '_Fish_', num2str(fish)])
    % saveas(gcf, fullfile(save_path, [groupname, '_Fish_', num2str(fish), '_dorsomedcells.png']))
    % saveas(gcf, fullfile(save_path, [groupname, '_Fish_', num2str(fish), '_dorsomedcells.svg']))
    % % close; 
    
    dorsomed_index_fish{fish} = dorsomed_ind; 
    thresh_dormedneurons{fish} = listcol; 

    norm_lat_med{fish} = norm_latmed; 
    norm_dor_ven{fish} = norm_dorven; 

    % lets make the same for ventral cells 
    unique_planes = unique(hb_position(:,5)); 

    listvent = zeros(size(hb_position,1),1);

    ven_planes = unique_planes(find(unique_planes > 2)); 
    for plane = 1:size(ven_planes)
        listvent(find(hb_position(:,5) == ven_planes(plane))) = 1; 

    end
    listvent(find(listcol) == 1) = 0; 
    ven_neurons{fish} = listvent;

    % using the normalised dorsomed axis I can then make the new list to
    % distingish between dor and ven

    thresh_dv = 40; %in micron 
    dv_list_fish = zeros(size(norm_dorven,1),1);
    for neuron = 1:size(norm_dorven,1)
        if norm_dorven(neuron) < thresh_dv
            dv_list_fish(neuron) = 1;

        else
            dv_list_fish(neuron) = 2;
        end

    end
    adapter_dv = zeros(size(position_wt{1,fish},1),1); 
    adapter_dv(find(brain_regions_wt{1,fish} == brainnumber)) = dv_list_fish;  
    dv_list{fish} = adapter_dv;
    % figure('units','centimeters','Position',[2 2 30 10])
    % scatter3(position_wt{1,fish}(:,1), position_wt{1,fish}(:,2), position_wt{1,fish}(:,3), 80, adapter_dv, 'filled')
    % % colormap(cmap2)
    % % colorbar
    % set(gca, 'Zdir', 'reverse'), set(gca, 'Ydir', 'reverse'), set(gca, 'Xdir', 'reverse')
    % view (80,7) ;%view (22,80) 
    % axis equal

    % figure('units','centimeters','Position',[2 2 30 10])
    % scatter3(hb_position(:,1), hb_position(:,2), hb_position(:,3), 80, dv_list_fish, 'filled')
    % % colormap(cmap2)
    % % colorbar
    % set(gca, 'Zdir', 'reverse'), set(gca, 'Ydir', 'reverse'), set(gca, 'Xdir', 'reverse')
    % view (80,7) ;%view (22,80) 
    % axis equal
    % title([groupname, '_Fish_', num2str(fish)])

   %  %%
   %  figure('units','centimeters','Position',[2 2 25 15])
   %  subplot(1,2,1)
   %  scatter3(hb_position(:,1), hb_position(:,2), hb_position(:,3), 80, dv_list_fish, 'filled')
   %  % colormap(cmap2)
   %  axis equal, view (65,90), grid on ; title (tit) % view (22,80)
   %  set (gca,'Zdir','reverse'),  set (gca,'Ydir','reverse'),  set (gca,'Xdir','reverse')
   %  % set(gca,'XTick',[]), set(gca,'YTick',[]),box off, axis tight, axis off,
   % 
   %  subplot(1,2,2)
   % scatter3(hb_position(:,1), hb_position(:,2), hb_position(:,3), 80, dv_list_fish, 'filled')
   %  % colormap(cmap2)
   % 
   %  axis equal, view (-267,8), grid on ; title (tit) % view (-72,13)
   %  set (gca,'Zdir','reverse'),  set (gca,'Ydir','reverse'),  set (gca,'Xdir','reverse')
   %  % set(gca,'XTick',[]), set(gca,'YTick',[]),box off, axis tight, axis off,
   % 
   %   set(gcf, 'InvertHardCopy', 'off');
   %      set(gcf, 'DefaultFigureRenderer', 'painters');
   %      set(gcf,'renderer','painters');
   %      saveas(gcf, fullfile(save_path, [tit, '_3D.png']))
   %  saveas(gcf, fullfile(save_path, [tit, '_3D.svg']))  


end
end