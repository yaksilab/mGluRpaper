function p_fac = cal_p_fac(dis_blo, pcc_sam_blo_ani)
%CAL_P_FAC - One line description of what the function or script performs (H1 line)
%   Optional file header info (to give more details about the function than in the H1 line)
%   Optional file header info (to give more details about the function than in the H1 line)
%   Optional file header info (to give more details about the function than in the H1 line)
%
%   Syntax:
%       output = CAL_P_FAC(input1, input2)
%       output = CAL_P_FAC(input1, input2, input3)
%
%   Description:
%       CAL_P_FAC() - description
%    
%   Inputs:
%       dis_blo - vector, each element is the distance for one block
%       pcc_sam_blo_ani - cell array, each cell contains the corr matrix
%       for a sample/group where rows are blocks, and colums are animals
%
%   Outputs:
%       p_fac - first element: p-value for the effect of distance on
%       correlation, second element: p-value for the effect of sample/group
%       on correlation
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

%   Author: Aytac Kadir Mutlu
%   Address: Olav Kyrres gate 9, 7030 Trondheim, Norway
%   email: kadir.a.mutlu@ntnu.no
%   Website: https://www.ntnu.edu/kavli
%   Date: 09-Apr-2024; Last revision: 09-Apr-2024
%
%   Copyright (c) 2024, Aytac Kadir Mutlu
%   All rights reserved.
n_blo = length(dis_blo);
n_sam = length(pcc_sam_blo_ani);
n_ani_sam = nan(n_sam, 1);
for i = 1:n_sam
    n_ani_sam(i) = size(pcc_sam_blo_ani{i}, 2);
end

% abl_blo_ani = [zeros(n_blo, n_ani_sam(1)) ones(n_blo, n_ani_sam(2))];
% cor_coe_blo_ani = [pcc_sam_blo_ani{1} pcc_sam_blo_ani{2}];
if size(n_ani_sam,1) > 2
    cat_blo_ani = categorical([1*ones(n_blo, n_ani_sam(1)) 2*ones(n_blo, n_ani_sam(2)) 3*ones(n_blo, ...
    n_ani_sam(3))]);
    cor_coe_blo_ani = [pcc_sam_blo_ani{1} pcc_sam_blo_ani{2} pcc_sam_blo_ani{3}];
else
    cat_blo_ani = categorical([1*ones(n_blo, n_ani_sam(1)) 2*ones(n_blo, n_ani_sam(2))]);
    cor_coe_blo_ani = [pcc_sam_blo_ani{1} pcc_sam_blo_ani{2} ];
end

fac_num = {repmat(dis_blo, sum(n_ani_sam), 1); cat_blo_ani(:)};
fac_num{1} = fac_num{1}(:);
p_fac = anovan(cor_coe_blo_ani(:), fac_num);
end
