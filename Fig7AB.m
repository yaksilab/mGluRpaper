%% Fig 7AB

% Novel Tank
% loading the data with the 1s binns
load('W:\Anna\Ricarda''s behaviour data in Matrix form\variables_1sec_bin.mat')

no_groups = 3; 
group_names = {'Wt', 'Het', 'Hom'}; 

%% Y position

% create the data matrices (base one fabrizios code=
for gr= 1 : no_groups
    y_temp=[];y_temp=y_position_mGluR6a_NTT{gr};
    for fish = 1:size(y_temp,2);
        for minute= 1: size(y_temp,1);

            y_min(gr,minute,fish)=mean(y_temp{minute,fish});
        end
    end
end

timpoint_start = 1;
timpoint_end = 120;

timpoint_start = 780;
timpoint_end = 900;

y_timebin_group1= y_min(1,timpoint_start:timpoint_end, 1:size(y_position_mGluR6a_NTT{1,1},2)); %(group, number of timebin, fish)
y_timebin_group2= y_min(2,timpoint_start:timpoint_end, 1:size(y_position_mGluR6a_NTT{1,2},2));
y_timebin_group3= y_min(3,timpoint_start:timpoint_end, 1:size(y_position_mGluR6a_NTT{1,3},2));

% flip it
y_timebin_group1=squeeze(y_timebin_group1);
y_timebin_group2=squeeze(y_timebin_group2);
y_timebin_group3=squeeze(y_timebin_group3);

%calculate mean and SEM 
ymean1=squeeze(mean(y_timebin_group1));
ymean2=squeeze(mean(y_timebin_group2));
ymean3=squeeze(mean(y_timebin_group3));

ysem1=squeeze(std(y_timebin_group1)/sqrt(length(y_timebin_group1)));
ysem2=squeeze(std(y_timebin_group2)/sqrt(length(y_timebin_group2)));
ysem3=squeeze(std(y_timebin_group3)/sqrt(length(y_timebin_group3)));

% here I could then replace it with the beswarm stuff 
figure()
scatter(ones(1,size(y_timebin_group1,2)),mean(y_timebin_group1),'filled','MarkerFaceAlpha',1/2,'MarkerFaceColor','k');hold on 
scatter(2.*ones(1,size(y_timebin_group2,2)),mean(y_timebin_group2),'filled','MarkerFaceAlpha',1/2,'MarkerFaceColor','r');hold on 
scatter(3.*ones(1,size(y_timebin_group3,2)),mean(y_timebin_group3),'filled','MarkerFaceAlpha',1/2,'MarkerFaceColor','b');hold on  

%Add mean & SEM + slope line
%y1=[ymean1 ymean2];
y1=[mean(ymean1) mean(ymean2) mean(ymean3)];
%y2=[ysem1 ysem2];
y2=[mean(ysem1) mean(ysem2) mean(ysem3)];
hold on;
errorbar(y1,y2,'linewidth',1,'color','k');

%plot design properties 
% ylim([-120 0])
%xlim([0 3])
xlim([0 4])
%xticks([1 2])
xticks([1 2 3])
xtickangle(45)
xticklabels(group_names)
ylabel('Distance from bottom of tank [mm]') 
% yticklabels({'0','20','40','60','80','100','120'})
%sigstar({[1,2]},[nan]);
% sigstar({[1,2] [2,3] [1,3]},[nan,nan,nan]);

[p1, h] = quick_statistic(mean(y_timebin_group1), mean(y_timebin_group2))
[p2, h] = quick_statistic(mean(y_timebin_group1), mean(y_timebin_group3))
[p3, h] = quick_statistic(mean(y_timebin_group2), mean(y_timebin_group3))
% firsth time point
x_spots = [[1.5 2 2.5]; [3.5 4 4.5]];
group_oder1 = [ones(size(y_timebin_group1,2),1)*1.5; ones(size(y_timebin_group2,2),1)*2;ones(size(y_timebin_group3,2),1)*2.5];
combined_data1 = [mean(y_timebin_group1)'; mean(y_timebin_group2)'; mean(y_timebin_group3)'];
grouop_avg1 = y1; 
sems1 = y2;
% second time point
group_oder2 = [ones(size(y_timebin_group1,2),1)*3.5; ones(size(y_timebin_group2,2),1)*4;ones(size(y_timebin_group3,2),1)*4.5]; 
combined_data2 = [mean(y_timebin_group1)'; mean(y_timebin_group2)'; mean(y_timebin_group3)'];
grouop_avg2 = y1; 
sems2 = y2;
group_oder = [group_oder1; group_oder2]; 
combined_data = [combined_data1; combined_data2]; 
grouop_avg = [grouop_avg1, grouop_avg2]; 
sems = [sems1, sems2]; 
figure
x = beeswarm(group_oder,combined_data)
title([' Early vs late '])
hold on
er = errorbar([x_spots(1,:), x_spots(2,:)],grouop_avg, sems)
er.Color = [0 0 0];                            
er.LineStyle = 'none';
xticks([2 4]);
xticklabels({'Early', 'Late'})
saveas(gcf, fullfile('W:\Anna\Ricarda''s behaviour data in Matrix form\', ['beeee.png']))
saveas(gcf, fullfile('W:\Anna\Ricarda''s behaviour data in Matrix form\',['beeee.svg']))
saveas(gcf, fullfile('W:\Anna\Ricarda''s behaviour data in Matrix form\ ',['beeee.fig']))
%% Y traces
timpoint_start = 1;
timpoint_end = 900;
y_timebin_group1= squeeze(y_min(1,timpoint_start:timpoint_end, 1:size(y_position_mGluR6a_NTT{1,1},2))); %(group, number of timebin, fish)
y_timebin_group2= squeeze(y_min(2,timpoint_start:timpoint_end, 1:size(y_position_mGluR6a_NTT{1,2},2)));
y_timebin_group3= squeeze(y_min(3,timpoint_start:timpoint_end, 1:size(y_position_mGluR6a_NTT{1,3},2)));


%calculate mean and SEM 
ymean1=squeeze(mean(y_timebin_group1,2));
ymean2=squeeze(mean(y_timebin_group2,2));
ymean3=squeeze(mean(y_timebin_group3,2));

ysem1=squeeze(std(y_timebin_group1,0,2)/sqrt(size(y_timebin_group1,2)));
ysem2=squeeze(std(y_timebin_group2,0,2)/sqrt(size(y_timebin_group2,2)));
ysem3=squeeze(std(y_timebin_group3,0,2)/sqrt(size(y_timebin_group3,2)));

time = [timpoint_start:timpoint_end];
figure('units','pixel','Position',[100 100 1500 600])
hold on

col=char('b','c','m');

H1 = shadedErrorBar(time, ymean1, ysem1,'lineProps',col(1)), hold on
H1.mainLine.LineWidth = 2;

H2 = shadedErrorBar(time, ymean2, ysem2,'lineProps',col(2)), hold on
H2.mainLine.LineWidth = 2;

H3 = shadedErrorBar(time, ymean3, ysem3,'lineProps',col(3)), hold on
H3.mainLine.LineWidth = 2;
ylabel('distance from bottom (mm)')
xlabel('time (s)')
set(gca,'TickDir','out')
% xline(0, '--k')
% xline(300, '--r')

  
legend([H1.mainLine, H2.mainLine, H3.mainLine], group_names);

%% Heatmaps 
timpoint_start = 1;
timpoint_end = 120;
PDF = heatmaps_mGluR6a_NTT; 
time_range = [timpoint_start:timpoint_end];
avg_PDF_gr = cell(1,no_groups); 
for gr= 1 : no_groups
    
    for fish = 1:size(PDF{1,gr},2)
        coll_pdf = nan(length(timpoint_start:timpoint_end), size(PDF{1,gr}{1,fish},1),size(PDF{1,gr}{1,fish},2));
        for timep = 1:length(timpoint_start:timpoint_end)
            curr_pdf = PDF{1,gr}{time_range(timep),fish};
            coll_pdf(timep, :,:) = curr_pdf; 
        end
        if fish == 1
            avg_fish = squeeze(mean(coll_pdf,1));
        else
            avg_fish = cat(3,avg_fish, squeeze(mean(coll_pdf,1)));
        end
    end
    avg_PDF_gr{1, gr} = avg_fish;
end

scaling=[0 0.005];
figure
for gr= 1 : no_groups
    subplot(1,3,gr)
    imagesc(squeeze(mean(avg_PDF_gr{1,gr},3)))
    colormap jet 
    %colorbar
    caxis (scaling)
    axis equal
end
saveas(gcf, fullfile('W:\Anna\Ricarda''s behaviour data in Matrix form\', ['heatm.png']))
saveas(gcf, fullfile('W:\Anna\Ricarda''s behaviour data in Matrix form\',['heatm.svg']))
saveas(gcf, fullfile('W:\Anna\Ricarda''s behaviour data in Matrix form\ ',['beeee.fig']))

for i=1:size(pdf_mean,1)
subplot(z,x,i+(x*(y-1)))
imagesc(pdf_mean{i,1})
colormap jet 
%colorbar
caxis (scaling)
sum(sum(pdf_mean{i,1}))
axis off
end