function [all_fish]=AO_aj_distCalculation(all_fish,nROIs,timebin,folder_path_save,fish_name,date_exp,conversion_factor);
%AO_aj_distCalculation - Calculates x,ydistance
%
%   Syntax:
%       [all_fish]=AO_aj_distCalculation(all_fish,nROIs,timebin,folder_path_save,fish_name,date_exp,conversion_factor);
%       output = function(input1, input2, input3)
%
%   Description:
%       function() - description
%    
%   Inputs:
%       conversion_factor - Pixel conversion factor
%       all_fish - Cell array with struct for each fish
%       nROIs - Total number of fish
%       timebin - Timebin for distance
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

for ROI=1:nROIs 
    t           = all_fish{ROI,1}.t;
    dT          = all_fish{ROI,1}.dt;
    dX          = diff(all_fish{ROI,1}.x);
    dY          = diff(all_fish{ROI,1}.y);
    n_bins2     = floor(max(all_fish{ROI,1}.t)/timebin);
    S2          = all_fish{ROI,1}.distance;
    
    S           = nan(1,length(dT));
    
    for i=1:length(dT)-1
    S (i)   = sqrt(dX(i)^2+dY(i)^2)*conversion_factor;
    end 
    
     

    
    %% Calculate binned distances
        s_temp      = nan(1,n_bins2);
        s_temp2     = nan(1,n_bins2);
        for i=1:n_bins2
            try
                s_temp(i)= sum(S ...
                        (t>(i-1)*timebin & t<=i*timebin));  
            catch 
                disp('Something weird with index again')
                disp(i)
                disp(ROI)
                s_temp(i) = nan; 
            end
            s_temp2(i)= sum(S2 ...
                    (t>(i-1)*timebin & t<=i*timebin));  
        end 
    
    %% Save data inn cell array
    if timebin < 1
        velo_string = num2str(timebin); 
        velo_string(2) = '_';
        all_fish{ROI,1}.(['calculatedRawDistance_' velo_string]) = S       ;
        all_fish{ROI,1}.(['calcBinnedDistance_' velo_string])    = s_temp  ;
        all_fish{ROI,1}.(['bD_' velo_string])                    = s_temp2 ;
    else
        all_fish{ROI,1}.(['calculatedRawDistance_' num2str(timebin)]) = S       ;
        all_fish{ROI,1}.(['calcBinnedDistance_' num2str(timebin)])    = s_temp  ;
        all_fish{ROI,1}.(['bD_' num2str(timebin)])                    = s_temp2 ;
    end
%     all_fish{ROI,1}.(['calculatedRawDistance_' num2str(timebin)]) = S       ;
%     all_fish{ROI,1}.(['calcBinnedDistance_' num2str(timebin)])    = s_temp  ;
%     all_fish{ROI,1}.(['bD_' num2str(timebin)])                    = s_temp2 ;
    
end
%save([folder_path_save,filesep,fish_name,date_exp '_data.mat'],'all_fish','-v7.3');
