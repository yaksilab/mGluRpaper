function AO_ploting_noveltank_heatmaps_over_fish(all_fish, group, name, folder_path_save, figures_subfolder, single_plotting)
%AO_ploting_noveltank_heatmaps_over_fish - Plots the NTT heatmap early and
%late
%   Author: Anna Maria Ostenrath
%   Optional file header info (to give more details about the function than in the H1 line)
%   Optional file header info (to give more details about the function than in the H1 line)
%
%   Syntax:
%       AO_ploting_noveltank_heatmaps_over_fish(all_fish, group, name, folder_path_save, figures_subfolder, single_plotting)
%       output = function(input1, input2, input3)
%
%   Description:
%       AO_ploting_noveltank_heatmaps_over_fish() - description
%    
%   Inputs:
%       all_fish - cell array with all the fish data
%       group - list for which fish are part of the group
%       name - name of the group
%       folder_path_save, figures_subfolder - part of the saving pathway
%       single_plotting - 1 or 0 if you want to plot the individual fish or
%       not
%
%
%   Outputs:
%       output1 - Description
%       output2 - Description
%
%   Examples: 
%       Line 1 of example
%       Line 2 of example
%       Line 3 of example
%
%   Other m-files required: none
%   Subfunctions: none
%   MAT-files required: none
%
%   See also: OTHER_FUNCTION_NAME1,  OTHER_FUNCTION_NAME2
%   Author: Anna Maria Ostenrath 
%   Date : Nov 2022	

% First we make a subfolder in case it does not exist yet
try
    mkdir(fullfile(folder_path_save, figures_subfolder, name))
catch
end
% here we build the variables to put all the PDFs in
avg_pdf = nan([size(all_fish{group(1), 1}.novelPDF,1), size(all_fish{group(1), 1}.novelPDF,2), 2, size(group,1)]); % 2 for before and after
for fish = 1: size(group,1)  
    avg_pdf(:,:,1,fish) = all_fish{group(fish), 1}.novelPDF(:,:,1); %early
    avg_pdf(:,:,2,fish) = all_fish{group(fish), 1}.novelPDF(:,:,2); %late
end

% plotting part 
if single_plotting 
    for fish = 1:size(group,1)  
        figure('units','centimeters','Position',[2 2 20 8])
        subplot(1,2,1)
        imagesc(avg_pdf(:,:,1,fish))
        title(['Avg Fish ', num2str(fish), ' Early ', name])
        colormap jet %hot
        caxis([0 0.05])
        colorbar

        subplot(1,2,2)
        imagesc(avg_pdf(:,:,2,fish))
        title(['Avg Fish ', num2str(fish), ' Late ', name])
        colormap jet % hot
        caxis([0 0.05])
        colorbar

        
        saveas(gcf, fullfile(folder_path_save, figures_subfolder, name, [name,'_', num2str(all_fish{group(fish), 1}.realNum),'_noveltank_heatmaps_per_fish.png'])); 
        close; 
    end
end


% now the averages
figure('units','centimeters','Position',[2 2 30 10])
subplot(1,2, 1)
imagesc(nanmean(avg_pdf(:,:,1,:),4))
title(['Avg over Fish Early ', name]) 
colormap jet %hot
caxis([0 0.05])
colorbar


subplot(1,2, 2)
imagesc(nanmean(avg_pdf(:,:,2,:),4))
title(['Avg over fish Late ', name])
colormap jet % hot
caxis([0 0.05])
colorbar   
saveas(gcf, fullfile(folder_path_save, figures_subfolder, name, [name, '_noveltank_heatmaps_avg_over_fish.png'])); 
% close; 
% 

end