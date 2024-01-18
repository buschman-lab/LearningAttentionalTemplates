%Plot Figure 6

load('Figure_6A')
load('Figure_6C')
load('Figure_6D')

addpath(genpath('Violinplot-Matlab-master'))

%% Figure 6A-B

figure
subplot(3,2,1)
n_color=1;
hold on
shadedErrorBar(Figure_6A.window_start_list(4:end),Figure_6A.LIP.mean_r_chosen_contra_boot(4:end),[Figure_6A.LIP.r_chosen_contra_boot_95(4:end),Figure_6A.LIP.r_chosen_contra_boot_5(4:end)]',{'-','color',Figure_6A.color_for_ROI(n_color,:),'LineWidth',2},2)
shadedErrorBar(Figure_6A.window_start_list(4:end),Figure_6A.LIP.mean_r_chosen_ipsi_boot(4:end),[Figure_6A.LIP.r_chosen_ipsi_boot_95(4:end),Figure_6A.LIP.r_chosen_ipsi_boot_5(4:end)]',{'--','color',Figure_6A.color_for_ROI(n_color,:),'LineWidth',2},2)
plot_significance_level(Figure_6A.window_start_list,Figure_6A.LIP.p_chosen_contra_boot,0.8,Figure_6A.color_for_ROI(n_color,:),[0.01,0.05/(length(Figure_6A.window_start_list)-3), 0.01/(length(Figure_6A.window_start_list)-3)])
plot_significance_level(Figure_6A.window_start_list,Figure_6A.LIP.p_chosen_ipsi_boot,0.75,Figure_6A.color_for_ROI(n_color,:),[0.01,0.05/(length(Figure_6A.window_start_list)-3), 0.01/(length(Figure_6A.window_start_list)-3)])
plot_significance_level(Figure_6A.window_start_list,Figure_6A.LIP.p_ch_contra_ipsi,0.85,'k',[0.01,0.05/(length(Figure_6A.window_start_list)-3), 0.01/(length(Figure_6A.window_start_list)-3)])

yl.LabelVerticalAlignment = 'middle';
yl.LabelHorizontalAlignment = 'right';
yl.FontSize=12;
text(-575,0.85,'Contra ≠ Ipsi')
text(-575,0.8,'Contra')
text(-575,0.75,'Ipsi') 
ylim([-0.2 1])
box off
xlabel('Time to target')
xlim([-600 Figure_6A.window_start_list(end)+150])
xticks([-600:200:Figure_6A.window_start_list(end)+150])
yline(0,'--')
ax=gca;
ax.XAxis.FontSize=12;
ax.YAxis.FontSize=12;

subplot(3,2,3)
n_color=2;
hold on
shadedErrorBar(Figure_6A.window_start_list(4:end),Figure_6A.FEF.mean_r_chosen_contra_boot(4:end),[Figure_6A.FEF.r_chosen_contra_boot_95(4:end),Figure_6A.FEF.r_chosen_contra_boot_5(4:end)]',{'-','color',Figure_6A.color_for_ROI(n_color,:),'LineWidth',2},2)
shadedErrorBar(Figure_6A.window_start_list(4:end),Figure_6A.FEF.mean_r_chosen_ipsi_boot(4:end),[Figure_6A.FEF.r_chosen_ipsi_boot_95(4:end),Figure_6A.FEF.r_chosen_ipsi_boot_5(4:end)]',{'--','color',Figure_6A.color_for_ROI(n_color,:),'LineWidth',2},2)
plot_significance_level(Figure_6A.window_start_list,Figure_6A.FEF.p_chosen_contra_boot,0.8,Figure_6A.color_for_ROI(n_color,:),[0.01,0.05/(length(Figure_6A.window_start_list)-3), 0.01/(length(Figure_6A.window_start_list)-3)])
plot_significance_level(Figure_6A.window_start_list,Figure_6A.FEF.p_chosen_ipsi_boot,0.75,Figure_6A.color_for_ROI(n_color,:),[0.01,0.05/(length(Figure_6A.window_start_list)-3), 0.01/(length(Figure_6A.window_start_list)-3)])
plot_significance_level(Figure_6A.window_start_list,Figure_6A.FEF.p_ch_contra_ipsi,0.85,'k',[0.01,0.05/(length(Figure_6A.window_start_list)-3), 0.01/(length(Figure_6A.window_start_list)-3)])

yl.LabelVerticalAlignment = 'middle';
yl.LabelHorizontalAlignment = 'right';
yl.FontSize=12;
text(-575,0.85,'Contra ≠ Ipsi')
text(-575,0.8,'Contra')
text(-575,0.75,'Ipsi') 
ylim([-0.2 1])
box off
xlabel('Time to target')
xlim([-600 Figure_6A.window_start_list(end)+150])
xticks([-600:200:Figure_6A.window_start_list(end)+150])
yline(0,'--')
ax=gca;
ax.XAxis.FontSize=12;
ax.YAxis.FontSize=12;

subplot(3,2,5)
n_color=3;
hold on
shadedErrorBar(Figure_6A.window_start_list(4:end),Figure_6A.PFC.mean_r_chosen_contra_boot(4:end),[Figure_6A.PFC.r_chosen_contra_boot_95(4:end),Figure_6A.PFC.r_chosen_contra_boot_5(4:end)]',{'-','color',Figure_6A.color_for_ROI(n_color,:),'LineWidth',2},2)
shadedErrorBar(Figure_6A.window_start_list(4:end),Figure_6A.PFC.mean_r_chosen_ipsi_boot(4:end),[Figure_6A.PFC.r_chosen_ipsi_boot_95(4:end),Figure_6A.PFC.r_chosen_ipsi_boot_5(4:end)]',{'--','color',Figure_6A.color_for_ROI(n_color,:),'LineWidth',2},2)
plot_significance_level(Figure_6A.window_start_list,Figure_6A.PFC.p_chosen_contra_boot,0.8,Figure_6A.color_for_ROI(n_color,:),[0.01,0.05/(length(Figure_6A.window_start_list)-3), 0.01/(length(Figure_6A.window_start_list)-3)])
plot_significance_level(Figure_6A.window_start_list,Figure_6A.PFC.p_chosen_ipsi_boot,0.75,Figure_6A.color_for_ROI(n_color,:),[0.01,0.05/(length(Figure_6A.window_start_list)-3), 0.01/(length(Figure_6A.window_start_list)-3)])
plot_significance_level(Figure_6A.window_start_list,Figure_6A.PFC.p_ch_contra_ipsi,0.85,'k',[0.01,0.05/(length(Figure_6A.window_start_list)-3), 0.01/(length(Figure_6A.window_start_list)-3)])

yl.LabelVerticalAlignment = 'middle';
yl.LabelHorizontalAlignment = 'right';
yl.FontSize=12;
text(-575,0.85,'Contra ≠ Ipsi')
text(-575,0.8,'Contra')
text(-575,0.75,'Ipsi') 
box off
xlabel('Time to target')
xlim([-600 Figure_6A.window_start_list(end)+150])
xticks([-600:200:Figure_6A.window_start_list(end)+150])
ax=gca;
ax.XAxis.FontSize=12;
ax.YAxis.FontSize=12;
ylim([-0.2 1])

subplot(3,2,2)
n_color=1;
hold on
shadedErrorBar(Figure_6A.window_start_list(4:end),Figure_6A.LIP.mean_r_unchosen_contra_boot(4:end),[Figure_6A.LIP.r_unchosen_contra_boot_95(4:end),Figure_6A.LIP.r_unchosen_contra_boot_5(4:end)]',{'-','color',Figure_6A.color_for_ROI(n_color,:),'LineWidth',2},2)
shadedErrorBar(Figure_6A.window_start_list(4:end),Figure_6A.LIP.mean_r_unchosen_ipsi_boot(4:end),[Figure_6A.LIP.r_unchosen_ipsi_boot_95(4:end),Figure_6A.LIP.r_unchosen_ipsi_boot_5(4:end)]',{'--','color',Figure_6A.color_for_ROI(n_color,:),'LineWidth',2},2)
plot_significance_level(Figure_6A.window_start_list,Figure_6A.LIP.p_unchosen_contra_boot,0.8,Figure_6A.color_for_ROI(n_color,:),[0.01,0.05/(length(Figure_6A.window_start_list)-3), 0.01/(length(Figure_6A.window_start_list)-3)])
plot_significance_level(Figure_6A.window_start_list,Figure_6A.LIP.p_unchosen_ipsi_boot,0.75,Figure_6A.color_for_ROI(n_color,:),[0.01,0.05/(length(Figure_6A.window_start_list)-3), 0.01/(length(Figure_6A.window_start_list)-3)])
plot_significance_level(Figure_6A.window_start_list,Figure_6A.LIP.p_unch_contra_ipsi,0.85,'k',[0.01,0.05/(length(Figure_6A.window_start_list)-3), 0.01/(length(Figure_6A.window_start_list)-3)])

yl.LabelVerticalAlignment = 'middle';
yl.LabelHorizontalAlignment = 'right';
yl.FontSize=12;
text(-575,0.85,'Contra ≠ Ipsi')
text(-575,0.8,'Contra')
text(-575,0.75,'Ipsi') 
ylim([-0.2 1])
box off
xlabel('Time to target')
xlim([-600 Figure_6A.window_start_list(end)+150])
xticks([-600:200:Figure_6A.window_start_list(end)+150])
ax=gca;
ax.XAxis.FontSize=12;
ax.YAxis.FontSize=12;

subplot(3,2,4)
n_color=2;
hold on
shadedErrorBar(Figure_6A.window_start_list(4:end),Figure_6A.FEF.mean_r_unchosen_contra_boot(4:end),[Figure_6A.FEF.r_unchosen_contra_boot_95(4:end),Figure_6A.FEF.r_unchosen_contra_boot_5(4:end)]',{'-','color',Figure_6A.color_for_ROI(n_color,:),'LineWidth',2},2)
shadedErrorBar(Figure_6A.window_start_list(4:end),Figure_6A.FEF.mean_r_unchosen_ipsi_boot(4:end),[Figure_6A.FEF.r_unchosen_ipsi_boot_95(4:end),Figure_6A.FEF.r_unchosen_ipsi_boot_5(4:end)]',{'--','color',Figure_6A.color_for_ROI(n_color,:),'LineWidth',2},2)
plot_significance_level(Figure_6A.window_start_list,Figure_6A.FEF.p_unchosen_contra_boot,0.8,Figure_6A.color_for_ROI(n_color,:),[0.01,0.05/(length(Figure_6A.window_start_list)-3), 0.01/(length(Figure_6A.window_start_list)-3)])
plot_significance_level(Figure_6A.window_start_list,Figure_6A.FEF.p_unchosen_ipsi_boot,0.75,Figure_6A.color_for_ROI(n_color,:),[0.01,0.05/(length(Figure_6A.window_start_list)-3), 0.01/(length(Figure_6A.window_start_list)-3)])
plot_significance_level(Figure_6A.window_start_list,Figure_6A.FEF.p_unch_contra_ipsi,0.85,'k',[0.01,0.05/(length(Figure_6A.window_start_list)-3), 0.01/(length(Figure_6A.window_start_list)-3)])

yl.LabelVerticalAlignment = 'middle';
yl.LabelHorizontalAlignment = 'right';
yl.FontSize=12;
text(-575,0.85,'Contra ≠ Ipsi')
text(-575,0.8,'Contra')
text(-575,0.75,'Ipsi') 
ylim([-0.2 1])
box off
xlabel('Time to target')
xlim([-600 Figure_6A.window_start_list(end)+150])
xticks([-600:200:Figure_6A.window_start_list(end)+150])
yline(0,'--')
ax=gca;
ax.XAxis.FontSize=12;
ax.YAxis.FontSize=12;

subplot(3,2,6)
n_color=3;
hold on
shadedErrorBar(Figure_6A.window_start_list(4:end),Figure_6A.PFC.mean_r_unchosen_contra_boot(4:end),[Figure_6A.PFC.r_unchosen_contra_boot_95(4:end),Figure_6A.PFC.r_unchosen_contra_boot_5(4:end)]',{'-','color',Figure_6A.color_for_ROI(n_color,:),'LineWidth',2},2)
shadedErrorBar(Figure_6A.window_start_list(4:end),Figure_6A.PFC.mean_r_unchosen_ipsi_boot(4:end),[Figure_6A.PFC.r_unchosen_ipsi_boot_95(4:end),Figure_6A.PFC.r_unchosen_ipsi_boot_5(4:end)]',{'--','color',Figure_6A.color_for_ROI(n_color,:),'LineWidth',2},2)
plot_significance_level(Figure_6A.window_start_list,Figure_6A.PFC.p_unchosen_contra_boot,0.8,Figure_6A.color_for_ROI(n_color,:),[0.01,0.05/(length(Figure_6A.window_start_list)-3), 0.01/(length(Figure_6A.window_start_list)-3)])
plot_significance_level(Figure_6A.window_start_list,Figure_6A.PFC.p_unchosen_ipsi_boot,0.75,Figure_6A.color_for_ROI(n_color,:),[0.01,0.05/(length(Figure_6A.window_start_list)-3), 0.01/(length(Figure_6A.window_start_list)-3)])
plot_significance_level(Figure_6A.window_start_list,Figure_6A.PFC.p_unch_contra_ipsi,0.85,'k',[0.01,0.05/(length(Figure_6A.window_start_list)-3), 0.01/(length(Figure_6A.window_start_list)-3)])

yl.LabelVerticalAlignment = 'middle';
yl.LabelHorizontalAlignment = 'right';
yl.FontSize=12;
text(-575,0.85,'Contra ≠ Ipsi')
text(-575,0.8,'Contra')
text(-575,0.75,'Ipsi') 
box off
xlabel('Time to target')
xlim([-600 Figure_6A.window_start_list(end)+150])
xticks([-600:200:Figure_6A.window_start_list(end)+150])
yline(0,'--')
ax=gca;
ax.XAxis.FontSize=12;
ax.YAxis.FontSize=12;
ylim([-0.2 1])
ax.YAxis.FontSize=12;

%% Figure 6C

figure
subplot(3,1,1)
vs = violinplot(Figure_6C.LIP.auc_r_chosen',{'1','2','3','4'},'ViolinColor',Figure_6C.color_for_ROI(1,:),'ShowData',false,'ShowMean',true,'ShowMedian',true);
box off
xlabel('Distance')
ylabel('AUC correlation chosen value vectors')
ax=gca;
ax.XAxis.FontSize=12;
ax.YAxis.FontSize=12;
ylim([2 10])
subplot(3,1,2)
vs = violinplot(Figure_6C.FEF.auc_r_chosen',{'1','2','3','4'},'ViolinColor',Figure_6C.color_for_ROI(2,:),'ShowData',false,'ShowMean',true,'ShowMedian',true);
box off
xlabel('Distance')
ylabel('AUC correlation chosen value vectors')
ax=gca;
ax.XAxis.FontSize=12;
ax.YAxis.FontSize=12;
ylim([2 10])

subplot(3,1,3)
vs = violinplot(Figure_6C.PFC.auc_r_chosen',{'1','2','3','4'},'ViolinColor',Figure_6C.color_for_ROI(3,:),'ShowData',false,'ShowMean',true,'ShowMedian',true);
box off
xlabel('Distance')
ylabel('AUC correlation chosen value vectors')
ax=gca;
ax.XAxis.FontSize=12;
ax.YAxis.FontSize=12;
ylim([2 10])

%% Figure 6D

figure
subplot(3,1,1)
n_color=1;
hold on
shadedErrorBar(Figure_6D.window_start_list(4:end),Figure_6D.LIP.mean_r_chosen_contra_chosen_global_boot(4:end),[Figure_6D.LIP.r_chosen_contra_chosen_global_boot_95(4:end) Figure_6D.LIP.r_chosen_contra_chosen_global_boot_5(4:end)]',{'-','color',Figure_6D.color_for_ROI(1,:),'LineWidth',2},2)
shadedErrorBar(Figure_6D.window_start_list(4:end),Figure_6D.FEF.mean_r_chosen_contra_chosen_global_boot(4:end),[Figure_6D.FEF.r_chosen_contra_chosen_global_boot_95(4:end) Figure_6D.FEF.r_chosen_contra_chosen_global_boot_5(4:end)]',{'-','color',Figure_6D.color_for_ROI(2,:),'LineWidth',2},2)
shadedErrorBar(Figure_6D.window_start_list(4:end),Figure_6D.PFC.mean_r_chosen_contra_chosen_global_boot(4:end),[Figure_6D.PFC.r_chosen_contra_chosen_global_boot_95(4:end) Figure_6D.PFC.r_chosen_contra_chosen_global_boot_5(4:end)]',{'-','color',Figure_6D.color_for_ROI(3,:),'LineWidth',2},2)

plot_significance_level(Figure_6D.window_start_list,Figure_6D.LIP.p_chosen_contra_sim,0.9,Figure_6D.color_for_ROI(1,:),[0.01,0.05/(length(Figure_6D.window_start_list)-3), 0.01/(length(Figure_6D.window_start_list)-3)])
plot_significance_level(Figure_6D.window_start_list,Figure_6D.LIP.p_chosen_contra_diss,0.65,Figure_6D.color_for_ROI(1,:),[0.01,0.05/(length(Figure_6D.window_start_list)-3), 0.01/(length(Figure_6D.window_start_list)-3)])

plot_significance_level(Figure_6D.window_start_list,Figure_6D.FEF.p_chosen_contra_sim,0.85,Figure_6D.color_for_ROI(2,:),[0.01,0.05/(length(Figure_6D.window_start_list)-3), 0.01/(length(Figure_6D.window_start_list)-3)])
plot_significance_level(Figure_6D.window_start_list,Figure_6D.FEF.p_chosen_contra_diss,0.6,Figure_6D.color_for_ROI(2,:),[0.01,0.05/(length(Figure_6D.window_start_list)-3), 0.01/(length(Figure_6D.window_start_list)-3)])

plot_significance_level(Figure_6D.window_start_list,Figure_6D.PFC.p_chosen_contra_sim,0.8,Figure_6D.color_for_ROI(3,:),[0.01,0.05/(length(Figure_6D.window_start_list)-3), 0.01/(length(Figure_6D.window_start_list)-3)])
plot_significance_level(Figure_6D.window_start_list,Figure_6D.PFC.p_chosen_contra_diss,0.55,Figure_6D.color_for_ROI(3,:),[0.01,0.05/(length(Figure_6D.window_start_list)-3), 0.01/(length(Figure_6D.window_start_list)-3)])

yl.LabelVerticalAlignment = 'middle';
yl.LabelHorizontalAlignment = 'right';
yl.FontSize=12;
text(-575,0.95,'Alignment to global')
text(-575,0.9,'LIP')
text(-575,0.85,'FEF')
text(-575,0.8,'PFC')
text(-575,0.7,'Local reliability > alignment to global')
text(-575,0.65,'LIP')
text(-575,0.6,'FEF')
text(-575,0.55,'PFC')
box off
xlabel('Time to target')
xlim([-600 Figure_6D.window_start_list(end)+150])
xticks([-600:200:Figure_6D.window_start_list(end)+150])
yline(0,'--')
xline(300,'-','Color',[0.25 0.25 0.75],'LineWidth',1)
hold on
xl=xline(400,'-',{'Reward'},'Color',[0.25 0.25 0.75],'LineWidth',1);
xl.LabelVerticalAlignment = 'bottom';
xl.LabelHorizontalAlignment = 'right';
xl.FontSize=12;
ax=gca;
ax.XAxis.FontSize=12;
ax.YAxis.FontSize=12;
ylabel('Correlation chosen value vectors')
ylim([-0.2 1])

%% Functions

function plot_significance_level(x,p,a,c,thr)
hold on
for b=1:length(thr)
    this_thr=thr(b);
    for i=4:length(p)-1
        if p(i)<=this_thr && p(i+1)<=this_thr
            plot(x(i):x(i+1),a*ones(1,x(i+1)-x(i)+1),'-','Color',c,'LineWidth',b)
        end
        if p(i-1)>this_thr && p(i)<=this_thr && p(i+1)>this_thr
            plot(x(i),a,'.','Color',c,'LineWidth',b)
        end
    end
        if p(end-1)>this_thr && p(end)<=this_thr
            plot(x(end),a,'.','Color',c,'LineWidth',b)
        end
end
end


