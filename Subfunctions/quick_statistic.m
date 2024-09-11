function [p, h] = quick_statistic(data1, data2)

[p_r,h_r,stats_r] = ranksum(data1,data2,'tail','right');
[p_b,h_b,stats_b] = ranksum(data1,data2,'tail','both');
[p_l,h_l,stats_l] = ranksum(data1,data2,'tail','left');

p = [p_l p_b p_r];
h = [h_l h_b h_r];

end