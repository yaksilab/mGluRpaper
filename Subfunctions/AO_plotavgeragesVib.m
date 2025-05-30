function AO_plotavgeragesVib(collect_dist_vib, vib_stimuli_trigger, metadata, short_cmap, tit)
%AO_plotavgeragesVib - Plot average heatmap and traces over the entire vib
%period
%   Author: Anna Maria Ostenrath
%   Optional file header info (to give more details about the function than in the H1 line)
%   Optional file header info (to give more details about the function than in the H1 line)
%
%   Syntax:
%       AO_plotavgeragesVib(collect_dist_vib, vib_stimuli_trigger, metadata, short_cmap, tit)
%       output = function(input1, input2, input3)
%
%   Description:
%       AO_plotavgeragesVib() - Plot average heatmap and traces over the entire vib
%period
%    
%   Inputs:
%       collect_dist_vib - cell array with data per fish per group
%       vib_stimuli_trigger - list of the stimulus onsets
%       metadata - struct with extra infor like savepath and groupnames
%       short_cmap - colors for plotting (rgb)
%       tit - string for the title
%       input3 - Description
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
%   Date : September 2024	

no_group = size(collect_dist_vib,1)

figure('units','pixel','Position',[100 100 1200 1000])
%for i=1 %for WT group only
for group=1:no_group  %for all groups
        subplot(no_group,1,group), imagesc(collect_dist_vib{group,1}') %select fish from certain group - encoded in T.Groupcounter
        colormap (flipud (hot))
        colorbar
        title ([metadata.GroupName(group,:), ' Vib']) %change according to group name
        ylabel('Fish number')
        xlabel('time (s)')
        box ('off')
        set(gca,'TickDir','out')
        caxis([0 10])
        xline(vib_stimuli_trigger, '--r', 'LineWidth', 2)

end
sgtitle('Vibration')
saveas(gcf, fullfile(metadata.data_save, [tit, '_','Avg_heatmapVib.png']))
saveas(gcf, fullfile(metadata.data_save, [tit, '_','Avg_heatmapVib.svg']))
saveas(gcf, fullfile(metadata.data_save, [tit, '_','Avg_heatmapVib.fig']))

% Now I want the traces per con
no_con = 4;
vibr_conTrials = [[1:10]; [11:20]; [21:30]; [31:40]];
vibr_conNames = {'Vib', 'Light Bef', 'VibLight', 'Light Aft'};
figure
for con = 1:no_con
    subplot(no_con, 1, con)
    hold on
    plplpl = [];
    for group = 1:no_group
        curr_period = [vib_stimuli_trigger(vibr_conTrials(con,1))-59:vib_stimuli_trigger(vibr_conTrials(con,10))+60]; 
        H1 = shadedErrorBar(curr_period, squeeze(mean(collect_dist_vib{group,1}(curr_period,:), 2)), std(collect_dist_vib{group,1}(curr_period,:),[],2)/sqrt(size(squeeze(collect_dist_vib{group,1}(curr_period,:)),2)))
        H1.mainLine.LineWidth = 1;
        H1.mainLine.Color = short_cmap(group,:);
        H1.patch.FaceColor = short_cmap(group,:); 
        H1.patch.EdgeColor= short_cmap(group,:);

        plplpl = [plplpl, H1.mainLine]; 

    end
    title(vibr_conNames{1,con})
    xline(vib_stimuli_trigger(vibr_conTrials(con,:)), 'Color', 'k', 'LineStyle', '--', 'LineWidth', 2)
    
    legend(plplpl, metadata.GroupName)
    ylabel('distance (mm)')
    xlabel('time (s)')
    
end
set(gcf,'units','centimeters','Position',[2 2 30 15])
saveas(gcf, fullfile(metadata.data_save, [tit, '_','Avg_tracesVib.png']))
saveas(gcf, fullfile(metadata.data_save, [tit, '_','Avg_tracesVib.svg']))
saveas(gcf, fullfile(metadata.data_save, [tit, '_','Avg_tracesVib.fig']))

end