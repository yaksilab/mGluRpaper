function h_fig = plo_dia_acc(dif_acc_epo_xre_sam, p_acc_epo_xre_xsa)
n_row.fig = 1;
n_col.fig = 1;
gap_ver_row = 0.010*ones(1, n_row.fig);
gap_hor_ave = 0.015;
gap_hor_col = gap_hor_ave*ones(1, n_col.fig);
gap_hor_row_col = repmat(gap_hor_col, n_row.fig, 1);
mar_bot = 0.000;
mar_top = 0.000;
mar_lef = 0.000;
mar_rig = 0.000;
gre = [0.5 0.5 0.5];
n_reg_tot = 12;
n_cro_tot = nchoosek(n_reg_tot, 2);
reg_cro_ind = nchoosek(1:n_reg_tot, 2);
n_row.cal = 1080;
n_col.cal = 789;
asp_rat_cax = n_col.cal/n_row.cal;

asp_rat_axe = asp_rat_cax*ones(1, n_row.fig*n_col.fig);
[h_fig, wid_mon, hei_mon] = fig;
wid = false;
sca_axe = ones(1, n_row.fig*n_col.fig);
n_pix_ext = 7;
[hax_sub, pos_axe, axh, dis_asp_rat, las_pix, rat_wid] = tight_subplot_gen(n_row.fig, n_col.fig, ...
    gap_ver_row, gap_hor_row_col, ...
    mar_bot, mar_top, mar_lef, mar_rig, asp_rat_axe, wid_mon, hei_mon, wid, sca_axe, n_pix_ext);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
one_six = n_col.cal/16;%Dd
eix = n_col.cal/8;
qua = n_col.cal/4;
hal = n_col.cal/2;
fiv_six = (5/16)*n_col.cal;
thr_eix = 3*eix;% dm, dmp
pix_num_hor_lef_reg = [fiv_six; eix; thr_eix; qua; one_six; thr_eix];
pix_num_hor_rig_reg = abs(pix_num_hor_lef_reg - hal) + hal;
pix_num_hor_reg = [pix_num_hor_lef_reg; pix_num_hor_rig_reg] + 7; %%%%!!!!!!!!!!!!!!!!!!!!!!!

adj = 0;
pix_num_hor_reg(5) = pix_num_hor_reg(5) - adj;
pix_num_hor_reg(11) = pix_num_hor_reg(11) + adj;

pix_num_ver_lef_reg = [n_row.cal/16; 0.75*n_row.cal; 0.5*n_row.cal; 0.625*n_row.cal; ...
    0.375*n_row.cal; 0.250*n_row.cal];
pix_num_ver_reg = [pix_num_ver_lef_reg; pix_num_ver_lef_reg] + 70; %%%%
mar_siz = 50;
max_abs_dif_acc_epo = cellfun(@(x) max(abs(x), [], 'all'), dif_acc_epo_xre_sam);
labels = {'Ha' 'Dl' 'Dm' 'Dc' 'Dd' 'Dmp' 'Ha' 'Dl' 'Dm' 'Dc' 'Dd' 'Dmp'};
fon_siz = 24;
sig_lev = 0.05;

j = 1;
scatter(hax_sub(j), pix_num_hor_reg, pix_num_ver_reg, 'Marker', 'o', ...
    'MarkerEdgeColor', gre, 'MarkerFaceColor', gre, 'SizeData', mar_siz);

for i = 1:n_cro_tot
    hold(hax_sub(j), 'on')
    dif = dif_acc_epo_xre_sam{epo}(i);
    reg_one = reg_cro_ind(i, 1);
    reg_two = reg_cro_ind(i, 2);
    if p_acc_epo_xre_xsa{epo, i} <= sig_lev
        if dif > 0
            col = "#00008B";
        else
            col = "#8B0000";
        end
        plo_xre(hax_sub(j), dif, pix_num_hor_reg, pix_num_ver_reg, ...
            reg_one, reg_two, max_abs_dif_acc_epo(epo), col)
    else
        if dif > 0
            col = "#ADD8E6";
        else
            col = "#FF7F7F";
        end
        plo_xre(hax_sub(j), dif, pix_num_hor_reg, pix_num_ver_reg, ...
            reg_one, reg_two, max_abs_dif_acc_epo(epo), col)
    end
    hax_sub(j).Visible = 'off';
end
hax_sub(j).Children = circshift(hax_sub(j).Children, 1);
axes(hax_sub(j))
labelpoints(pix_num_hor_reg, pix_num_ver_reg, labels, 'FontSize', fon_siz);
% should be installed through web?
        
%         scatter(hax_sub(j), x, y, 10, 'k', 'filled')
%         hax_sub(j).XLim(2) = max(x);
%         hax_sub(j).YLim(2) = max(y);

%h_fig = opt_h_fig(h_fig, las_pix, rat_wid);
end

function plo_xre(hax, dif, pix_num_hor_reg, pix_num_ver_reg, reg_one, reg_two, max_abs_dif, col)
lin_wid.two = 1;
abs_dif = abs(dif);
plot(hax, ...
    [pix_num_hor_reg(reg_one); pix_num_hor_reg(reg_two)], ...
    [pix_num_ver_reg(reg_one); pix_num_ver_reg(reg_two)], ...
    'Color', col, 'LineWidth', lin_wid.two + (abs_dif/max_abs_dif)*7)
end
