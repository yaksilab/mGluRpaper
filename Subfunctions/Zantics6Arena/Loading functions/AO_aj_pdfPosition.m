
function [all_fish] = AO_aj_pdfPosition(all_fish,nROIs, duration) %,idxStimulus,folder_path_save,fish_name,date_exp,figures_subfolder,group1,group2,name_group1,name_group2)
%AO_aj_pdfPosition - Calculates various PDf for different time points
%   Optional file header info (to give more details about the function than in the H1 line)
%   Optional file header info (to give more details about the function than in the H1 line)
%
%   Syntax:
%       output = function(input1, input2)
%       output = function(input1, input2, input3)
%
%   Description:
%       function() - description
%    
%   Inputs:
%       input1 - Description
%       input2 - Description
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
%   See also: OTHER_FUNCTION_NAME1,  OTHER_FUNCTION_NAME2
%   Author: Anna Maria Ostenrath and Ahmed Jamali 
%   Date : March 2022 (updated Nov 2022)	

%% for the first hour spontanious activity
for ROI=1:nROIs
    
    t = all_fish{ROI,1}.t;
    x = all_fish{ROI,1}.x;
    y = all_fish{ROI,1}.y;
    
    if sum(diff(all_fish{ROI, 1}.y))==0 && sum(diff(all_fish{ROI, 1}.x))==0
        pdf1 = ('No movement - no fish?');
        all_fish{ROI,1}.positionPD1 = pdf1;
    elseif isnan(mean(all_fish{ROI, 1}.y,'omitnan'))...
            && isnan(mean(all_fish{ROI, 1}.x,'omitnan'))
        pdf1 = ('Not detected - no fish?');
        all_fish{ROI,1}.positionPD1 = pdf1;
    else
        pdf = AO_aj_getPDF(...
            (x(t>0 & t<=floor(max(t)))),... %-min(all_fish{ROI,1}.x)
            (-y(t>0 & t<=floor(max(t)))),...% -min(all_fish{ROI,1}.y)
            110,...
            ROI,...
            nROIs);
        all_fish{ROI,1}.positionPD_overall = pdf;
    end
end
%% For novel tank
%first 2 min and last 2 min? 
time_duration = 1*60; % 2 min with 15 frame rate
for ROI=1:nROIs
    
    t = all_fish{ROI,1}.t;
    x = all_fish{ROI,1}.x;
    y = all_fish{ROI,1}.y;
    stim_times = all_fish{ROI, 1}.timeStampStim(find(~contains(all_fish{ROI, 1}.stimInfo, "VIB_0"))); 
    
    if sum(diff(all_fish{ROI, 1}.y))==0 && sum(diff(all_fish{ROI, 1}.x))==0
        pdf1 = ('No movement - no fish?');
        all_fish{ROI,1}.positionPD1 = pdf1;
    elseif isnan(mean(all_fish{ROI, 1}.y,'omitnan'))...
            && isnan(mean(all_fish{ROI, 1}.x,'omitnan'))
        pdf1 = ('Not detected - no fish?');
        all_fish{ROI,1}.positionPD1 = pdf1;
    else
        novel_pdf = nan([size(all_fish{ROI,1}.positionPD_overall),2]);
%         novel_pdf = nan([size(all_fish{ROI,1}.positionPD_overall),2]);
        
        curr_stim_ons = stim_times(1); 
        pdf_early = AO_aj_getPDF(...
            (x(t>20& t<20+floor(time_duration))),... %-min(all_fish{ROI,1}.x)
            (-y(t>20& t<20+floor(time_duration))),...% -min(all_fish{ROI,1}.y)
            110,...
            ROI,...
            nROIs);
        close;
%         pdf_late = aj_getPDF_20211206(...
%             (x(t>floor(curr_stim_ons-time_duration) & t<=floor(curr_stim_ons))),... %-min(all_fish{ROI,1}.x)
%             (-y(t>floor(curr_stim_ons-time_duration) & t<=floor(curr_stim_ons))),...% -min(all_fish{ROI,1}.y)
%             110,...
%             ROI,...
%             nROIs);
            pdf_late = AO_aj_getPDF(...
        (x(t>1200-time_duration & t<=1200)),... %-min(all_fish{ROI,1}.x)
        (-y(t>1200-time_duration & t<=1200)),...% -min(all_fish{ROI,1}.y)
        110,...
        ROI,...
        nROIs);
        close; 
        novel_pdf(:,:,1) = pdf_early;
        novel_pdf(:,:,2) = pdf_late;
        
        all_fish{ROI,1}.novelPDF = novel_pdf;
    end
end

%% for each stimulus
for ROI=1:nROIs
    
    t = all_fish{ROI,1}.t;
    x = all_fish{ROI,1}.x;
    y = all_fish{ROI,1}.y;
    stim_times = all_fish{ROI, 1}.timeStampStim(find(~contains(all_fish{ROI, 1}.stimInfo, "VIB_0"))); 
    
    if sum(diff(all_fish{ROI, 1}.y))==0 && sum(diff(all_fish{ROI, 1}.x))==0
        pdf1 = ('No movement - no fish?');
        all_fish{ROI,1}.positionPD1 = pdf1;
    elseif isnan(mean(all_fish{ROI, 1}.y,'omitnan'))...
            && isnan(mean(all_fish{ROI, 1}.x,'omitnan'))
        pdf1 = ('Not detected - no fish?');
        all_fish{ROI,1}.positionPD1 = pdf1;
    else
        stim_pdf = nan([size(all_fish{ROI,1}.positionPD_overall),2, length(stim_times)]);
        for stim = 1:length(stim_times)
            curr_stim_ons = stim_times(stim); 
            pdf_before = AO_aj_getPDF(...
                (x(t>floor(curr_stim_ons)-duration & t<floor(curr_stim_ons))),... %-min(all_fish{ROI,1}.x)
                (-y(t>floor(curr_stim_ons)-duration & t<floor(curr_stim_ons))),...% -min(all_fish{ROI,1}.y)
                110,...
                ROI,...
                nROIs);
            close;
            pdf_after = AO_aj_getPDF(...
                (x(t>floor(curr_stim_ons) & t<=floor(curr_stim_ons)+duration)),... %-min(all_fish{ROI,1}.x)
                (-y(t>floor(curr_stim_ons) & t<=floor(curr_stim_ons)+duration)),...% -min(all_fish{ROI,1}.y)
                110,...
                ROI,...
                nROIs);
            close; 
            stim_pdf(:,:,1,stim) = pdf_before;
            stim_pdf(:,:,2,stim) = pdf_after;
        end
        all_fish{ROI,1}.stimPDF = stim_pdf;
    end
end

%% avg over the stim time compared to base
try
for ROI=1:nROIs
    
    if all_fish{ROI, 1}.selected
        t = all_fish{ROI,1}.t;
        x = all_fish{ROI,1}.x;
        y = all_fish{ROI,1}.y;
        stim_times = all_fish{ROI, 1}.timeStampStim(find(~contains(all_fish{ROI, 1}.stimInfo, "VIB_0"))); 
        stim_types = unique(all_fish{ROI, 1}.stimInfo); 
        stim_types = stim_types(2:end); % we are taking out the VIB 0 one 
        stim_types_indice = all_fish{ROI, 1}.stimInfo(find(~contains(all_fish{ROI, 1}.stimInfo, "VIB_0"))); % this is now the stimuli times for the specific stimulus 

        if sum(diff(all_fish{ROI, 1}.y))==0 && sum(diff(all_fish{ROI, 1}.x))==0
            pdf1 = ('No movement - no fish?');
            all_fish{ROI,1}.positionPD1 = pdf1;
        elseif isnan(mean(all_fish{ROI, 1}.y,'omitnan'))...
                && isnan(mean(all_fish{ROI, 1}.x,'omitnan'))
            pdf1 = ('Not detected - no fish?');
            all_fish{ROI,1}.positionPD1 = pdf1;
        else
            stim_pdf = nan([size(all_fish{ROI,1}.positionPD_overall),2, length(stim_times)]);
            for stim = 1:length(stim_types)
                curr_stim_type = find(stim_types_indice == stim_types(stim)); 

                curr_stim_ons = stim_times(curr_stim_type);
                bef_duration = 4*60; % 4 min 
                % now we average over the 4min before the onset
                pdf_before = aj_getPDF_20211206(...
                    (x(t>floor(curr_stim_ons(1))-bef_duration & t<floor(curr_stim_ons(1)))),... %-min(all_fish{ROI,1}.x)
                    (-y(t>floor(curr_stim_ons(1))-bef_duration & t<floor(curr_stim_ons(1)))),...% -min(all_fish{ROI,1}.y)
                    110,...
                    ROI,...
                    nROIs);
    %             close;
                after_duration = 11*60; %because we look at each stim plus 1min extra 
                pdf_after = aj_getPDF_20211206(...
                    (x(t>floor(curr_stim_ons(1)) & t<=floor(curr_stim_ons(1))+after_duration)),... %-min(all_fish{ROI,1}.x)
                    (-y(t>floor(curr_stim_ons(1)) & t<=floor(curr_stim_ons(1))+after_duration)),...% -min(all_fish{ROI,1}.y)
                    110,...
                    ROI,...
                    nROIs);
    %             close; 
                stim_pdf(:,:,1,stim) = pdf_before;
                stim_pdf(:,:,2,stim) = pdf_after;
            end
            all_fish{ROI,1}.avgstimPDF = stim_pdf;
        end
    else 
        continue
    end
end

catch
end
%% For the stimulus
% for i=1:length(idxStimulus)
%     for ROI=1:nROIs
%         t = all_fish{ROI,1}.t;
%         timeStim = all_fish{ROI,1}.stimulusTime;
%         
%         if sum(diff(all_fish{ROI, 1}.y))==0 &&...
%                 sum(diff(all_fish{ROI, 1}.x))==0
%             pdf1 = ('No movement - no fish?');
%             all_fish{ROI,1}.positionPD1 = pdf1;
%         elseif isnan(mean(all_fish{ROI, 1}.y,'omitnan'))...
%                 && isnan(mean(all_fish{ROI, 1}.x,'omitnan'))
%             pdf1 = ('Not detected - no fish?');
%             all_fish{ROI,1}.positionPD1 = pdf1;
%         else
%             pdf = aj_getPDF_20211206(...
%                 (all_fish{ROI,1}.x(t>=idxStimulus(ROI,i)*timeStim &...
%                 t<idxStimulus(ROI,i)*timeStim)),... %-min(all_fish{ROI,1}.x)
%                 (all_fish{ROI,1}.y(t>=idxStimulus(ROI,i)*timeStim &...
%                 t<idxStimulus(ROI,i)*timeStim)),...% -min(all_fish{ROI,1}.y)
%                 110,...
%                 ROI,...
%                 nROIs);
%             all_fish{ROI,1}.positionPD_overall = pdf;
%         end
%     end
%     
% end
% 



