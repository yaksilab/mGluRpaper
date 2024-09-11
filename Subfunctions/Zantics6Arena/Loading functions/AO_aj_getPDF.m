function[pdf]=AO_aj_getPDF(X,Y,xdimmax,ROI,nROIs)
%AO_aj_getPDF - Calculates PDF for certain interval
%   Author: Anna Maria Ostenrath
%   Optional file header info (to give more details about the function than in the H1 line)
%   Optional file header info (to give more details about the function than in the H1 line)
%
%   Syntax:
%       output = function(input1, input2)
%       output = function(input1, input2, input3)
%
%   Description:
%       AO_aj_getPDF() - description
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
%   Author: Anna Maria Ostenrath and Ahmed Jamali 
%   Date : updated Nov 2022	

numExp = nROIs/6;
for i=1:numExp
lowerRow (i*3-(3-1):i*3) = [4 5 6]+6*(i-1); % make an array for the lower row of ROIs 
end

for i=1:numExp
secondCol (i*2-(2-1):i*2) = [2 5]+6*(i-1); % make an array for the lower row of ROIs 
end

for i=1:numExp
thirdCol (i*2-(2-1):i*2) = [3 6]+6*(i-1); % make an array for the lower row of ROIs 
end 

% if ismember(ROI,secondCol)
%     x= 145:170/25:315;
% elseif ismember(ROI,thirdCol)
%     x= 359:170/25:529;
% else 
%     x= 0:170/25: 170; % axis x, which you want to see
% end 

% the 25 and 37.5 are for the "mesh grid" step size you can adapt that if
% you want more details etc... 
% changed to 15 and 27.5
% if ismember(ROI,secondCol)
%     x= 130:120/15:240;
% elseif ismember(ROI,thirdCol)
%     x= 250:120/15:360;
% else 
%     x= 10:120/15: 120; % axis x, which you want to see
% end 
if ismember(ROI,secondCol)
    x= 130:120/30:240;
elseif ismember(ROI,thirdCol)
    x= 250:120/30:360;
else 
    x= 10:120/30: 120; % axis x, which you want to see
end 
% x= 10:120/30: 120;
% y = -140:120/25: 20;
% if ismember(ROI,lowerRow)
%     y = -310 :170/27.5: -140; % axis y, which you want to see
% else 
%     y = -71 :170/27.5: 100; % axis y, which you want to see
% end 
% if ismember(ROI,lowerRow)
%     y = -280 :120/17.5: -140; % axis y, which you want to see
% else 
%     y = -140:120/17.5: 0; % axis y, which you want to see
% end 
if ismember(ROI,lowerRow)
    y = -280 :120/25: -140; % axis y, which you want to see
else 
    y = -140:120/25: 0; % axis y, which you want to see
end 
% figure, plot (X,Y);

[Xtemp,Ytemp] = meshgrid(x,y); % important for "surf" - makes defined grid

pdf = hist3([X ,Y],{x y}); % standard hist3 (calculated for yours axis)
pdf_normalize = ((pdf./ length(X))); % normalization means devide it by length of the vector
pdf_normalize =(pdf_normalize)';% data_x (or data_y)
pdf_normalize = flip(pdf_normalize,1);
% figure,surf(Xtemp,Ytemp,pdf_normalize)
pdf=pdf_normalize;

end