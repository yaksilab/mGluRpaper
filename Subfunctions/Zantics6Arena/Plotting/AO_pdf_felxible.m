function [group_pdf] = AO_pdf_felxible(all_fish, group, group_name, folder_path_save, figures_subfolder, start_points,end_points, single_plot, save_name, plotting, seconds)
%AO_pdf_felxible - Create and plot PDF for specific time intervals 
%   Author: Anna Maria Ostenrath
%
%   Syntax:
%       [group_pdf] = AO_pdf_felxible(all_fish, group, names, folder_path_save, figures_subfolder, start_points,end_points, single_plot, save_name, plotting)
%       output = function(input1, input2, input3)
%
%   Description:
%       function() - description
%    
%   Inputs:
%       all_fish - cell array with all the fish data
%       group - List of the fish belonging to the group so you can select
%       them from the all_fish struct
%       group_name - Name of the current group eg. 'Wt'
%       folder_path_save &  figures_subfolder- saving paths
%       start_points - list of the start points of the time intervals
%       end_points - list of the end points of the time intervals
%       start_points(1) should correspond to end_points(1)
%       single_plot - 1 or 0 for plotting the individual fish each time
%       point seperately
%       save_name - cell array of names of the time points eg. {'Early
%       NTT', 'Late NTT'}; 
%       plotting - 1 or 0 if you want to plot the average of the indiv time
%       points
%       seconds - 1 or 0 if the time points are in frames or in seconds (1:
%       the timepoint is in seconds, 0: the time point is in frames)
%
%
%   Outputs:
%       group_pdf - matrix of the pdfs dim: xdim(30) x ydim(28) x no of timepoints x no of
%       fish
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
%   See also: AO_aj_getPDF,  AO_aj_pdfPosition
%   Author: Anna Maria Ostenrath 
%   Date : Jan 2023	


% to make the pdf we first have to know to which arena the fish belong
% we have that in the arena variable 

no_time_point = length(start_points); 
group_pdf = nan(30, 28, no_time_point, size(group,1)); % right now this is hard coded for the mesh

% here we make sure we seperate the fish into the right arenas
upperRow = [1 2 3];
lowerRow = [4 5 6];
firstCol = [1 4]; 
secondCol = [2 5];
thirdCol = [3 6];

for fish = 1:size(group,1)
    current_fish = group(fish); 
    
    for time = 1 : no_time_point
        %first we now find the current time interval
        if seconds
            start_time = find(all_fish{current_fish, 1}.t >= start_points(time));
            start_time = start_time(1);
            end_time = find(all_fish{current_fish, 1}.t >= end_points(time));
            end_time = end_time(1);
        else
            start_time = start_points(time); 
            end_time = end_points(time);
        end
        % finding the arena number of our current fish
        arena_num = all_fish{current_fish, 1}.fishNum;
        % extracting the X and Y values according to the current time
        % window
        X = all_fish{current_fish, 1}.x(start_time:end_time); 
        Y = -all_fish{current_fish, 1}.y(start_time:end_time); 
        
        % now we make the grid 
         
        if ismember(arena_num,secondCol)
            x= 130:120/30:240;
        elseif ismember(arena_num,thirdCol)
            x= 250:120/30:360;
        else 
            x= 10:120/30:120; % axis x, which you want to see
        end 
        
        if ismember(arena_num,lowerRow)
            y = -280:120/25:-140; % axis y, which you want to see 168
        else 
            y = -140:120/25: 0; % axis y, which you want to see 26?
        end 
      
                
        pdf = hist3([X ,Y],{x y}); % standard hist3 (calculated for yours axis)
        pdf_normalize = ((pdf./ length(X))); % normalization means devide it by length of the vector
        pdf_normalize =(pdf_normalize)';% fliping the data around
        pdf_normalize = flip(pdf_normalize,1);% fliping the data around
        % figure,surf(Xtemp,Ytemp,pdf_normalize)
        pdf=pdf_normalize;
        group_pdf(:,:,time,fish) = pdf; 

    end
    if single_plot 
        figure('units','centimeters','Position',[2 2 20 8])
        for time = 1 : no_time_point
            subplot(1, no_time_point, time)
            imagesc(group_pdf(:,:,time,fish))
            caxis([0 0.1])
            title([group_name ' Fish' num2str(fish) '  ' save_name{time}])

        end
        saveas(gcf, fullfile(folder_path_save, figures_subfolder, [group_name, '_fish',num2str(fish), '_',save_name{time},'_heatmap_flexible.png']));
        close; 
    end
end

if plotting 
    figure('units','centimeters','Position',[2 2 20 6])
    for time = 1 : no_time_point
        subplot(1, no_time_point, time)
        imagesc(mean(group_pdf(:,:,time,:),4))
        caxis([0 0.05])
        title([group_name ' ' save_name{time}])
        colorbar
        colormap jet
    end
    saveas(gcf, fullfile(folder_path_save, figures_subfolder, [group_name,'_heatmap_flexible.png']));
    
end

end