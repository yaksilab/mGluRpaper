function [norm_latmed, norm_dorven] = normalise_position_hb(position_data, all_planes_y)
%normalise_position_hb - Normalise the lat-med and dor-ven axis in hb 
%   Author: Anna Maria Ostenrath
%
%   Syntax:
%       [norm_latmed, norm_dorven] = normalise_position_hb(position_data, all_planes_y)
%
%   Description:
%       function() - description
%    
%   Inputs:
%       position_data - position data in the format: neurons x 5 (x, y, z, neuronidx, planeidx)
%       all_planes_y - 1 or 0 if you should use all the planes to caluclate
%       the midline(1) or only plane1&2(0)
%
%   Outputs:
%       norm_latmed - normlaised y position
%       norm_dorven - normalised z position
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
%   Date : Oct 2023



% we will normalise the lat med axis by taking a midline through hb and
% then calculating the distance to that midline. that means that values
% close to 0 (both pos and ned) are medial and close to like 100 are very
% lateral in additiion we have positive and neg values to see what is left
% and what is right but that then depends on the exp 
if all_planes_y
    middle_point = mean(position_data(:,2)); % this should be yposition

else 
     middle_point = mean(position_data(find(position_data(:,5) == 1) | position_data(:,5) == 2),2);
end


% what if I first split them into two and then look for the min and max
% values of each 
[pos_ind,~,~,~] = kmeans(position_data(:,2),2,'Replicates', 50,'Distance','sqeuclidean');
mean_hb1 = mean(position_data(find(pos_ind == 1),2)); 
mean_hb2 = mean(position_data(find(pos_ind == 2),2)); 

middle = mean([mean_hb1, mean_hb2, middle_point]); 

figure('units','centimeters','Position',[2 2 30 10])
scatter3(position_data(:,1), position_data(:,2), position_data(:,3), 50, 'filled')
% colormap(cmap2)
% colorbar
set(gca, 'Zdir', 'reverse'), set(gca, 'Ydir', 'reverse'), set(gca, 'Xdir', 'reverse')
view (80,7) ;%view (22,80) 
axis equal
title('All cells')
hold on 
yline(middle_point)
yline(middle, 'Color', 'b')
% what if I first split them into two and then look for the min and max
% values of each 



norm_latmed = position_data(:,2)-middle; 

% we will normalise the dorso_ventral axis just by taking the highest value
% of all the cells and calculating the distance to that 
top_value = min(position_data(:,3));
norm_dorven = position_data(:,3) - top_value;


end