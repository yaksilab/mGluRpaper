function [hax_sub, pos_axe, hei_axe, dis_asp_rat, las_pix, rat_wid] = tight_subplot_gen(n_row, n_col, ...
    gap_ver_row, gap_hor_row_col, mar_bot, mar_top, mar_lef, mar_rig, asp_rat_axe, fig_wid, fig_hei, ...
    wid, sca_axe, n_pix_ext)
% tight_subplot creates "subplot" axes with adjustable gaps and margins
%
% [ha, pos] = tight_subplot(Nh, Nw, gap, marg_h, marg_w)
%
%   in:  Nh      number of axes in hight (vertical direction)
%        Nw      number of axes in width (horizontaldirection)
%        gap     gaps between the axes in normalized units (0...1)
%                   or [gap_h gap_w] for different gaps in height and width 
%        marg_h  margins in height in normalized units (0...1)
%                   or [lower upper] for different lower and upper margins 
%        marg_w  margins in width in normalized units (0...1)
%                   or [left right] for different left and right margins 
%
%  out:  ha     array of handles of the axes objects
%                   starting from upper left corner, going row-wise as in
%                   subplot
%        pos    positions of the axes objects; cell array
%
%  Example: ha = tight_subplot(3,2,[.01 .03],[.1 .01],[.01 .01])
%           for ii = 1:6; axes(ha(ii)); plot(randn(10,ii)); end
%           set(ha(1:4),'XTickLabel',''); set(ha,'YTickLabel','')

% Pekka Kumpulainen 21.5.2012   @tut.fi
% Tampere University of Technology / Automation Science and Engineering

% n_row_col, *n_col_row* ... either !

mar_pol = [mar_bot mar_top];
mar_sid = [mar_lef mar_rig];
dis_asp_rat = fig_wid/fig_hei;
n_plo = n_row*n_col;
if wid
    wid_lef = 1 - sum(mar_sid) - sum(gap_hor_row_col(1, 1:end - 1));
    
    wid_ori_axe = (wid_lef/n_col)*ones(n_plo, 1);
    
    wid_axe = wid_ori_axe.*sca_axe;

    hei_ori_axe = wid_axe*dis_asp_rat; % THIS IS THE CORRECT ONE

    %%%!!!! ONLY FOR RASTER MAP, COMMENT OTHERWISE - nei, nei, nei
%     if sca_axe(2) == 0.5
%         hei_ori_axe = ([1; 3; 1; 3].*wid_axe)*dis_asp_rat;
%     else
%         hei_ori_axe = wid_axe*dis_asp_rat;
%     end

    hei_axe = hei_ori_axe./asp_rat_axe; % THIS IS THE CORRECT ONE
    % sum(hei_axe) must be < 1, otherwise overflow 

    % ONLY FOR EQUAL HEIGHT, COMMENT OTHERWISE !!!
    %%%hei_axe = max(hei_axe)*ones(1, n_plo);
    %%%%%hei_axe = min(hei_axe)*ones(1, n_plo);
    %%hei_axe = hei_axe(2)*ones(1, n_plo);
    %hei_axe(1) = hei_axe(2);
    %hei_axe(3) = hei_axe(4);
else
    hei_ori_axe = ((1 - sum(mar_pol) - sum(gap_ver_row(1:end - 1)))/n_row)*ones(n_plo, 1);

    hei_axe = hei_ori_axe.*sca_axe;

    wid_axe = (hei_axe.*asp_rat_axe)/dis_asp_rat;
end

%y_pos = 1 - mar_top - max(hei_axe, [], 'all');
hei_max_row = max(hei_axe, [], 2); % MAKE HEI_ROW_COL !!!!
y_pos = 1 - mar_top - hei_max_row(1);

ii = 0;
for ih = 1:n_row
    px = mar_sid(1);
    for ix = 1:n_col
        ii = ii + 1;% rows, then columns
        hax_sub(ii) = axes('Units', 'normalized', 'Position', [px y_pos wid_axe(ii) hei_axe(ii)], ...
            'XTickLabel', '', 'YTickLabel', '');
        px = px + wid_axe(ii) + gap_hor_row_col(ih, ix);
    end
    %y_pos = y_pos - hei_ori_axe(ii) - gap_ver_row(ih);
    
    if ii < n_plo
        y_pos = y_pos - hei_axe(ii + 1) - gap_ver_row(ih);
    end
    

end
pos_axe = get(hax_sub, 'Position');
hax_sub = hax_sub(:);
if ~iscell(pos_axe)
    pos_axe = {pos_axe};
end
[las_pix, rat_wid] = ext_las_pix(pos_axe, fig_wid, n_pix_ext);
