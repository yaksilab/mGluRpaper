function acc_epo_xre = cal_acc_epo_xre(act_fra_cel, log_tot_cel_reg, n_reg_tot, n_epo, fra_epo_sam)
%CAL_ACC_EPO_XRE - One line description of what the function or script performs (H1 line)
%   Optional file header info (to give more details about the function than in the H1 line)
%   Optional file header info (to give more details about the function than in the H1 line)
%   Optional file header info (to give more details about the function than in the H1 line)
%
%   Syntax:
%       output = CAL_ACC_EPO_XRE(input1, input2)
%       output = CAL_ACC_EPO_XRE(input1, input2, input3)
%
%   Description:
%       CAL_ACC_EPO_XRE() - description
%    
%   Inputs:
%       log_tot_cel_reg: logical array indicating if a cell is in a region
%       n_reg_tot: number of regions
%       fra_epo_sam: cell array containing the frame numbers per epoch
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

%   Author: Aytac Kadir Mutlu
%   Address: Olav Kyrres gate 9, 7030 Trondheim, Norway
%   email: kadir.a.mutlu@ntnu.no
%   Website: https://www.ntnu.edu/kavli
%   Date: 30-Oct-2023; Last revision: 30-Oct-2023
%
%   Copyright (c) 2023, Aytac Kadir Mutlu
%   All rights reserved.
cor_coe_epo_cel_cel = cell(n_epo, 1);
for j = 1:n_epo
    cor_coe_epo_cel_cel{j} = corr(act_fra_cel(fra_epo_sam{j}, :));
end
%%%%%%%%%%%%% xre
n_cro_tot = nchoosek(n_reg_tot, 2);
cor_coe_epo_xre_cel_cel = cell(n_epo, n_cro_tot);
cor_coe_epo_xre_pai = cell(n_epo, n_cro_tot);
acc_epo_xre_pai = cell(n_epo, n_cro_tot);
acc_epo_xre = nan(n_epo, n_cro_tot);
reg_cro_ind = nchoosek(1:n_reg_tot, 2);
for j = 1:n_epo
    for i = 1:n_cro_tot
        cor_coe_epo_xre_cel_cel{j, i} = cor_coe_epo_cel_cel{j}(log_tot_cel_reg...
            (:, reg_cro_ind(i, 1)), log_tot_cel_reg(:, reg_cro_ind(i, 2)));
        cor_coe_epo_xre_pai{j, i} = ext_upp(cor_coe_epo_xre_cel_cel{j, i});
        acc_epo_xre_pai{j, i} = abs(cor_coe_epo_xre_pai{j, i});
        acc_epo_xre(j, i) = mean(acc_epo_xre_pai{j, i});
    end
end
end

function ele_row = ext_upp(ele_row_col)
%EXT_UPP - One line description of what the function or script performs (H1 line)
%   Optional file header info (to give more details about the function than in the H1 line)
%   Optional file header info (to give more details about the function than in the H1 line)
%   Optional file header info (to give more details about the function than in the H1 line)
%
%   Syntax:
%       output = EXT_UPP(input1, input2)
%       output = EXT_UPP(input1, input2, input3)
%
%   Description:
%       EXT_UPP() - description
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

%   Author: Kadir Mutlu
%   Address: Olav Kyrres gate 9, 7030 Trondheim, Norway
%   email: kadir.a.mutlu@ntnu.no
%   Website: https://www.ntnu.edu/kavli
%   Date: 03-Sep-2021; Last revision: 03-Sep-2021
%
%   Copyright (c) 2021, Kadir Mutlu
%   All rights reserved.
[n_row, n_col] = size(ele_row_col);
tru_row_col = true(n_row, n_col);
log_row_col = triu(tru_row_col, 1);
ele_row = ele_row_col(log_row_col);
end
