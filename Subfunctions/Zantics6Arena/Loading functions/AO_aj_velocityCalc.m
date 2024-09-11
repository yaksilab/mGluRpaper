function [all_fish]=AO_aj_velocityCalc(all_fish,nROIs,velBin,folder_path_save,fish_name,conversion_factor)
%AO_aj_velocityCalc - Calculate the velocity for specific time bin
%
%   Syntax:
%       output = function(input1, input2)
%       output = function(input1, input2, input3)
%
%   Description:
%       AO_aj_velocityCalc() - description
%    
%   Inputs:
%       all_fish - Cell array with struct for each fish
%       nROIs - Total number of fish
%       velBin - Time bin for velocity
%
%   Outputs:
%       all_fish - Update Cell array with struct for each fish
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
%   Author: Anna Maria Ostenrath and Ahmed Jamali 
%   Date : March 2022 (updated Nov 2022)

for ROI=1:nROIs %parfor
    %% Load variables
    t           = all_fish{ROI,1}.t;
    dt          = all_fish{ROI,1}.dt;
    s           = all_fish{ROI,1}.distance;
    n_bins      = floor(max(all_fish{ROI,1}.t)/velBin);
    
    %% Calculate speed over time and delta-time
    V1 = nan(1,length(dt));
    for j=1:length(dt)-1
        V1(j)= (s(j+1))/dt(j);
    end
    
    %% Calculate velocity per second
    
    bV_temp=nan(1,n_bins);
%     disp(ROI)
    for i=1:n_bins
%         disp(i)
        try
            bV_temp(i)= sum(s(t>(i-1)*velBin & t<=i*velBin))/...
                        sum(dt(t>(i-1)*velBin & t<=i*velBin));
        catch
            disp('Here is weird thing with dt and length of t')
            disp(ROI)
            disp(i)
            bV_temp(i) = nan; 
        end
    end
    %% Save data inn cell array
    all_fish{ROI,1}.speed_over_time                   =   V1;
    if velBin < 1
        velo_string = num2str(velBin); 
        velo_string(2) = '_';
        all_fish{ROI,1}.(['binnedVel_' velo_string])  =   bV_temp;
        all_fish{ROI,1}.(['speed_over_time_' velo_string])  =   V1;
    else
        all_fish{ROI,1}.(['binnedVel_' num2str(velBin)])  =   bV_temp;
        all_fish{ROI,1}.(['speed_over_time_' num2str(velBin)]) =   V1;
    end
end

% save([folder_path_save,filesep,fish_name,date_exp '_data.mat'],'all_fish','-v7.3');