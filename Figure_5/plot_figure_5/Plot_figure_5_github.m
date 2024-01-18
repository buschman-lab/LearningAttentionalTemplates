%Plot Figure 5

load('Figure_5B.mat')
load('Figure_5C.mat')
load('Figure_5D.mat')
load('Figure_5E.mat')

N_loc=4;

%% Figure 5B

figure
subplot(3,1,1)
hold on
shadedErrorBar(Figure_5B.window_start_list,mean(mean(Figure_5B.LIP.Same_belief_value(:,:,:),2),3),[prctile(mean(Figure_5B.LIP.Same_belief_value(:,:,:),2),95,3),prctile(mean(Figure_5B.LIP.Same_belief_value(:,:,:),2),5,3)]',{'Color',Figure_5B.color_for_ROI(1,:),'LineWidth',2},2)
shadedErrorBar(Figure_5B.window_start_list,mean(mean(Figure_5B.LIP.Other_belief_value(:,:,:),2),3),[prctile(mean(Figure_5B.LIP.Other_belief_value(:,:,:),2),95,3),prctile(mean(Figure_5B.LIP.Other_belief_value(:,:,:),2),5,3)]',{'--','Color',Figure_5B.color_for_ROI(1,:),'LineWidth',2},2)

plot_significance_level(Figure_5B.window_start_list,Figure_5B.p_LIP_same_other(:),0.9,Figure_5B.color_for_ROI(1,:),[0.01,0.05/length(Figure_5B.window_start_list), 0.01/length(Figure_5B.window_start_list)])
plot_significance_level(Figure_5B.window_start_list,Figure_5B.p_LIP_other(:),0.95,Figure_5B.color_for_ROI(1,:),[0.01,0.05/length(Figure_5B.window_start_list), 0.01/length(Figure_5B.window_start_list)])
plot_significance_level(Figure_5B.window_start_list,Figure_5B.p_LIP_same(:),1,Figure_5B.color_for_ROI(1,:),[0.01,0.05/length(Figure_5B.window_start_list), 0.01/length(Figure_5B.window_start_list)])

yl=yline(0.5,'--',{'Chance level'},'Color',[0 0 0],'LineWidth',1);
text(-575,0.9,'Within - across')
text(-575,0.95,'Across')
text(-575,1,'Within')
box off
xlabel('Time to targets on')
ylabel('Chosen value decoding accuracy')
xlim([-600 Figure_5B.window_start_list(end)+150])
xticks([-600:200:Figure_5B.window_start_list(end)+150])

subplot(3,1,2)
hold on
shadedErrorBar(Figure_5B.window_start_list,mean(mean(Figure_5B.FEF.Same_belief_value(:,:,:),2),3),[prctile(mean(Figure_5B.FEF.Same_belief_value(:,:,:),2),95,3),prctile(mean(Figure_5B.FEF.Same_belief_value(:,:,:),2),5,3)]',{'Color',Figure_5B.color_for_ROI(2,:),'LineWidth',2},2)
shadedErrorBar(Figure_5B.window_start_list,mean(mean(Figure_5B.FEF.Other_belief_value(:,:,:),2),3),[prctile(mean(Figure_5B.FEF.Other_belief_value(:,:,:),2),95,3),prctile(mean(Figure_5B.FEF.Other_belief_value(:,:,:),2),5,3)]',{'--','Color',Figure_5B.color_for_ROI(2,:),'LineWidth',2},2)

plot_significance_level(Figure_5B.window_start_list,Figure_5B.p_FEF_same_other(:),0.9,Figure_5B.color_for_ROI(2,:),[0.01,0.05/length(Figure_5B.window_start_list), 0.01/length(Figure_5B.window_start_list)])
plot_significance_level(Figure_5B.window_start_list,Figure_5B.p_FEF_other(:),0.95,Figure_5B.color_for_ROI(2,:),[0.01,0.05/length(Figure_5B.window_start_list), 0.01/length(Figure_5B.window_start_list)])
plot_significance_level(Figure_5B.window_start_list,Figure_5B.p_FEF_same(:),1,Figure_5B.color_for_ROI(2,:),[0.01,0.05/length(Figure_5B.window_start_list), 0.01/length(Figure_5B.window_start_list)])

yl=yline(0.5,'--',{'Chance level'},'Color',[0 0 0],'LineWidth',1);
xline(-300,'-','Color',[0.5 0.5 0.5],'LineWidth',1)
text(-575,0.9,'Within - across')
text(-575,0.95,'Across')
text(-575,1,'Within')
box off
xlabel('Time to targets on')
ylabel('Chosen value decoding accuracy')
xlim([-600 Figure_5B.window_start_list(end)+150])
xticks([-600:200:Figure_5B.window_start_list(end)+150])


subplot(3,1,3)
hold on
shadedErrorBar(Figure_5B.window_start_list,mean(mean(Figure_5B.PFC.Same_belief_value(:,:,:),2),3),[prctile(mean(Figure_5B.PFC.Same_belief_value(:,:,:),2),95,3),prctile(mean(Figure_5B.PFC.Same_belief_value(:,:,:),2),5,3)]',{'Color',Figure_5B.color_for_ROI(3,:),'LineWidth',2},2)
shadedErrorBar(Figure_5B.window_start_list,mean(mean(Figure_5B.PFC.Other_belief_value(:,:,:),2),3),[prctile(mean(Figure_5B.PFC.Other_belief_value(:,:,:),2),95,3),prctile(mean(Figure_5B.PFC.Other_belief_value(:,:,:),2),5,3)]',{'--','Color',Figure_5B.color_for_ROI(3,:),'LineWidth',2},2)

plot_significance_level(Figure_5B.window_start_list,Figure_5B.p_PFC_same_other(:),0.9,Figure_5B.color_for_ROI(3,:),[0.01,0.05/length(Figure_5B.window_start_list), 0.01/length(Figure_5B.window_start_list)])
plot_significance_level(Figure_5B.window_start_list,Figure_5B.p_PFC_other(:),0.95,Figure_5B.color_for_ROI(3,:),[0.01,0.05/length(Figure_5B.window_start_list), 0.01/length(Figure_5B.window_start_list)])
plot_significance_level(Figure_5B.window_start_list,Figure_5B.p_PFC_same(:),1,Figure_5B.color_for_ROI(3,:),[0.01,0.05/length(Figure_5B.window_start_list), 0.01/length(Figure_5B.window_start_list)])

yl=yline(0.5,'--',{'Chance level'},'Color',[0 0 0],'LineWidth',1);
text(-575,0.9,'Within - across')
text(-575,0.95,'Across')
text(-575,1,'Within')
box off
xlabel('Time to targets on')
ylabel('Chosen value decoding accuracy')
xlim([-600 Figure_5B.window_start_list(end)+150])
xticks([-600:200:Figure_5B.window_start_list(end)+150])

%% Figure 5C

figure
subplot(3,1,1)
n_color=1;
hold on
shadedErrorBar(Figure_5C.window_start_list(4:end),Figure_5C.LIP.mean_r_reward_global_boot(4:end),[Figure_5C.LIP.r_reward_global_boot_95(4:end), Figure_5C.LIP.r_reward_global_boot_5(4:end)]',{'-','color',[ 0 0.25 0],'LineWidth',2},2)
shadedErrorBar(Figure_5C.window_start_list(4:end),Figure_5C.LIP.mean_r_chosen_global_boot(4:end),[Figure_5C.LIP.r_rchosen_global_boot_95(4:end), Figure_5C.LIP.r_chosen_global_boot_5(4:end)]',{'-','color',Figure_5C.color_for_ROI(n_color,:),'LineWidth',2},2)

plot_significance_level(Figure_5C.window_start_list,Figure_5C.LIP.p_chosen_global_boot,0.8,Figure_5C.color_for_ROI(n_color,:),[0.01,0.05/(length(Figure_5C.window_start_list)-3), 0.01/(length(Figure_5C.window_start_list)-3)])
plot_significance_level(Figure_5C.window_start_list,Figure_5C.LIP.p_reward_global_boot,0.85,[ 0 0.25 0],[0.01,0.05/(length(Figure_5C.window_start_list)-3), 0.01/(length(Figure_5C.window_start_list)-3)])
plot_significance_level(Figure_5C.window_start_list,Figure_5C.LIP.p_ch_rew,0.9,'k',[0.01,0.05/(length(Figure_5C.window_start_list)-3), 0.01/(length(Figure_5C.window_start_list)-3)])

yl.LabelVerticalAlignment = 'middle';
yl.LabelHorizontalAlignment = 'right';
yl.FontSize=12;

text(-575,0.8,'Chosen')
text(-575,0.85,'Reward')
text(-575,0.9,'Reward ≠ Chosen')

box off
xlabel('Time to target')
xlim([-600 Figure_5C.window_start_list(end)+150])
xticks([-600:200:Figure_5C.window_start_list(end)+150])
yline(0,'--')
ax=gca;
ax.XAxis.FontSize=12;
ax.YAxis.FontSize=12;
ylim([-0.2 1])

subplot(3,1,2)
n_color=2;
hold on
shadedErrorBar(Figure_5C.window_start_list(4:end),Figure_5C.FEF.mean_r_reward_global_boot(4:end),[Figure_5C.FEF.r_reward_global_boot_95(4:end), Figure_5C.FEF.r_reward_global_boot_5(4:end)]',{'-','color',[ 0 0.25 0],'LineWidth',2},2)
shadedErrorBar(Figure_5C.window_start_list(4:end),Figure_5C.FEF.mean_r_chosen_global_boot(4:end),[Figure_5C.FEF.r_rchosen_global_boot_95(4:end), Figure_5C.FEF.r_chosen_global_boot_5(4:end)]',{'-','color',Figure_5C.color_for_ROI(n_color,:),'LineWidth',2},2)

plot_significance_level(Figure_5C.window_start_list,Figure_5C.FEF.p_chosen_global_boot,0.8,Figure_5C.color_for_ROI(n_color,:),[0.01,0.05/(length(Figure_5C.window_start_list)-3), 0.01/(length(Figure_5C.window_start_list)-3)])
plot_significance_level(Figure_5C.window_start_list,Figure_5C.FEF.p_reward_global_boot,0.85,[ 0 0.25 0],[0.01,0.05/(length(Figure_5C.window_start_list)-3), 0.01/(length(Figure_5C.window_start_list)-3)])
plot_significance_level(Figure_5C.window_start_list,Figure_5C.FEF.p_ch_rew,0.9,'k',[0.01,0.05/(length(Figure_5C.window_start_list)-3), 0.01/(length(Figure_5C.window_start_list)-3)])

yl.LabelVerticalAlignment = 'middle';
yl.LabelHorizontalAlignment = 'right';
yl.FontSize=12;

text(-575,0.8,'Chosen')
text(-575,0.85,'Reward')
text(-575,0.9,'Reward ≠ Chosen')

box off
xlabel('Time to target')
xlim([-600 Figure_5C.window_start_list(end)+150])
xticks([-600:200:Figure_5C.window_start_list(end)+150])
yline(0,'--')
ax=gca;
ax.XAxis.FontSize=12;
ax.YAxis.FontSize=12;
ylim([-0.2 1])

subplot(3,1,3)
n_color=3;
hold on
shadedErrorBar(Figure_5C.window_start_list(4:end),Figure_5C.PFC.mean_r_reward_global_boot(4:end),[Figure_5C.PFC.r_reward_global_boot_95(4:end), Figure_5C.PFC.r_reward_global_boot_5(4:end)]',{'-','color',[ 0 0.25 0],'LineWidth',2},2)
shadedErrorBar(Figure_5C.window_start_list(4:end),Figure_5C.PFC.mean_r_chosen_global_boot(4:end),[Figure_5C.PFC.r_rchosen_global_boot_95(4:end), Figure_5C.PFC.r_chosen_global_boot_5(4:end)]',{'-','color',Figure_5C.color_for_ROI(n_color,:),'LineWidth',2},2)

plot_significance_level(Figure_5C.window_start_list,Figure_5C.PFC.p_chosen_global_boot,0.8,Figure_5C.color_for_ROI(n_color,:),[0.01,0.05/(length(Figure_5C.window_start_list)-3), 0.01/(length(Figure_5C.window_start_list)-3)])
plot_significance_level(Figure_5C.window_start_list,Figure_5C.PFC.p_reward_global_boot,0.85,[ 0 0.25 0],[0.01,0.05/(length(Figure_5C.window_start_list)-3), 0.01/(length(Figure_5C.window_start_list)-3)])
plot_significance_level(Figure_5C.window_start_list,Figure_5C.PFC.p_ch_rew,0.9,'k',[0.01,0.05/(length(Figure_5C.window_start_list)-3), 0.01/(length(Figure_5C.window_start_list)-3)])

yl.LabelVerticalAlignment = 'middle';
yl.LabelHorizontalAlignment = 'right';
yl.FontSize=12;

text(-575,0.8,'Chosen')
text(-575,0.85,'Reward')
text(-575,0.9,'Reward ≠ Chosen')

box off
xlabel('Time to target')
xlim([-600 Figure_5C.window_start_list(end)+150])
xticks([-600:200:Figure_5C.window_start_list(end)+150])
yline(0,'--')
ax=gca;
ax.XAxis.FontSize=12;
ax.YAxis.FontSize=12;
ylim([-0.2 1])



%% Figure 5D

figure

subplot(3,1,1)
hold on
shadedErrorBar(Figure_5D.window_start_list,nanmean(nanmean(Figure_5D.LIP.Same_belief_chosen_color(:,:,:),2),3),[prctile(nanmean(Figure_5D.LIP.Same_belief_chosen_color(:,:,:),2),95,3),prctile(nanmean(Figure_5D.LIP.Same_belief_chosen_color(:,:,:),2),5,3)]',{'Color',Figure_5D.color_for_ROI(1,:),'LineWidth',2},2)
shadedErrorBar(Figure_5D.window_start_list,nanmean(nanmean(Figure_5D.LIP.Other_belief_chosen_color(:,:,:),2),3),[prctile(nanmean(Figure_5D.LIP.Other_belief_chosen_color(:,:,:),2),95,3),prctile(nanmean(Figure_5D.LIP.Other_belief_chosen_color(:,:,:),2),5,3)]',{'--','Color',Figure_5D.color_for_ROI(1,:),'LineWidth',2},2)

plot_significance_level(Figure_5D.window_start_list,Figure_5D.p_LIP_same_other(:),0.85,Figure_5D.color_for_ROI(1,:),[0.05, 0.01, 0.05/(length(Figure_5D.window_start_list))])
plot_significance_level(Figure_5D.window_start_list,Figure_5D.p_LIP_other(:),0.9,Figure_5D.color_for_ROI(1,:),[0.05, 0.01, 0.05/(length(Figure_5D.window_start_list))])
plot_significance_level(Figure_5D.window_start_list,Figure_5D.p_LIP_same(:),0.95,Figure_5D.color_for_ROI(1,:),[0.05, 0.01, 0.05/(length(Figure_5D.window_start_list))])

yl=yline(1/2,'--',{'Chance level'},'Color',[0 0 0],'LineWidth',1);
yl.LabelVerticalAlignment = 'bottom';
yl.LabelHorizontalAlignment = 'left';
yl.FontSize=12;
ylim([0.1 1])
text(-575,0.85,'Within - Across')
text(-575,0.9,'Across')
text(-575,0.95,'Withins')
box off
xlabel('Time to targets on')
ylabel('Chosen color decoding accuracy')
xlim([-600 Figure_5D.window_start_list(end)+150])
xticks([-600:200:Figure_5D.window_start_list(end)+150])

subplot(3,1,2)
hold on
shadedErrorBar(Figure_5D.window_start_list,nanmean(nanmean(Figure_5D.FEF.Same_belief_chosen_color(:,:,:),2),3),[prctile(nanmean(Figure_5D.FEF.Same_belief_chosen_color(:,:,:),2),95,3),prctile(nanmean(Figure_5D.FEF.Same_belief_chosen_color(:,:,:),2),5,3)]',{'Color',Figure_5D.color_for_ROI(2,:),'LineWidth',2},2)
shadedErrorBar(Figure_5D.window_start_list,nanmean(nanmean(Figure_5D.FEF.Other_belief_chosen_color(:,:,:),2),3),[prctile(nanmean(Figure_5D.FEF.Other_belief_chosen_color(:,:,:),2),95,3),prctile(nanmean(Figure_5D.FEF.Other_belief_chosen_color(:,:,:),2),5,3)]',{'--','Color',Figure_5D.color_for_ROI(2,:),'LineWidth',2},2)

plot_significance_level(Figure_5D.window_start_list,Figure_5D.p_FEF_same_other(:),0.85,Figure_5D.color_for_ROI(2,:),[0.05, 0.01, 0.05/(length(Figure_5D.window_start_list))])
plot_significance_level(Figure_5D.window_start_list,Figure_5D.p_FEF_other(:),0.9,Figure_5D.color_for_ROI(2,:),[0.05, 0.01, 0.05/(length(Figure_5D.window_start_list))])
plot_significance_level(Figure_5D.window_start_list,Figure_5D.p_FEF_same(:),0.95,Figure_5D.color_for_ROI(2,:),[0.05, 0.01, 0.05/(length(Figure_5D.window_start_list))])

yl=yline(1/2,'--',{'Chance level'},'Color',[0 0 0],'LineWidth',1);
yl.LabelVerticalAlignment = 'bottom';
yl.LabelHorizontalAlignment = 'left';
yl.FontSize=12;
ylim([0.1 1])
text(-575,0.85,'Within - Across')
text(-575,0.9,'Across')
text(-575,0.95,'Withins')
box off
xlabel('Time to targets on')
ylabel('Chosen color decoding accuracy')
xlim([-600 Figure_5D.window_start_list(end)+150])
xticks([-600:200:Figure_5D.window_start_list(end)+150])

subplot(3,1,3)
hold on
shadedErrorBar(Figure_5D.window_start_list,nanmean(nanmean(Figure_5D.PFC.Same_belief_chosen_color(:,:,:),2),3),[prctile(nanmean(Figure_5D.PFC.Same_belief_chosen_color(:,:,:),2),95,3),prctile(nanmean(Figure_5D.PFC.Same_belief_chosen_color(:,:,:),2),5,3)]',{'Color',Figure_5D.color_for_ROI(3,:),'LineWidth',2},2)
shadedErrorBar(Figure_5D.window_start_list,nanmean(nanmean(Figure_5D.PFC.Other_belief_chosen_color(:,:,:),2),3),[prctile(nanmean(Figure_5D.PFC.Other_belief_chosen_color(:,:,:),2),95,3),prctile(nanmean(Figure_5D.PFC.Other_belief_chosen_color(:,:,:),2),5,3)]',{'--','Color',Figure_5D.color_for_ROI(3,:),'LineWidth',2},2)

plot_significance_level(Figure_5D.window_start_list,Figure_5D.p_PFC_same_other(:),0.85,Figure_5D.color_for_ROI(3,:),[0.05, 0.01, 0.05/(length(Figure_5D.window_start_list))])
plot_significance_level(Figure_5D.window_start_list,Figure_5D.p_PFC_other(:),0.9,Figure_5D.color_for_ROI(3,:),[0.05, 0.01, 0.05/(length(Figure_5D.window_start_list))])
plot_significance_level(Figure_5D.window_start_list,Figure_5D.p_PFC_same(:),0.95,Figure_5D.color_for_ROI(3,:),[0.05, 0.01, 0.05/(length(Figure_5D.window_start_list))])

yl=yline(1/2,'--',{'Chance level'},'Color',[0 0 0],'LineWidth',1);
yl.LabelVerticalAlignment = 'bottom';
yl.LabelHorizontalAlignment = 'left';
yl.FontSize=12;
ylim([0.1 1])
text(-575,0.85,'Within - Across')
text(-575,0.9,'Across')
text(-575,0.95,'Withins')
box off
xlabel('Time to targets on')
ylabel('Chosen color decoding accuracy')
xlim([-600 Figure_5D.window_start_list(end)+150])
xticks([-600:200:Figure_5D.window_start_list(end)+150])

%% Figure 5E


figure

for loc=2 %Figure 5E %1:N_loc (see all locs: fig S5E)
    subplot(3,4,loc)
    hold on
    
    shadedErrorBar(Figure_5E.window_start_list,mean(mean(Figure_5E.LIP.Same_belief_choice(:,:,:,loc),2),3),[prctile(mean(Figure_5E.LIP.Same_belief_choice(:,:,:,loc),2),95,3),prctile(mean(Figure_5E.LIP.Same_belief_choice(:,:,:,loc),2),5,3)]',{'Color',Figure_5E.color_for_ROI(1,:),'LineWidth',2},2)
    shadedErrorBar(Figure_5E.window_start_list,mean(mean(Figure_5E.LIP.Other_belief_choice(:,:,:,loc),2),3),[prctile(mean(Figure_5E.LIP.Other_belief_choice(:,:,:,loc),2),95,3),prctile(mean(Figure_5E.LIP.Other_belief_choice(:,:,:,loc),2),5,3)]',{'--','Color',Figure_5E.color_for_ROI(1,:),'LineWidth',2},2)
    
    plot_significance_level(Figure_5E.window_start_list,Figure_5E.p_LIP_same_other(:,loc),1,Figure_5E.color_for_ROI(1,:),[0.01,0.05/(length(Figure_5E.window_start_list))/4, 0.01/(length(Figure_5E.window_start_list))/4])
    plot_significance_level(Figure_5E.window_start_list,Figure_5E.p_LIP_other(:,loc),1.05,Figure_5E.color_for_ROI(1,:),[0.01,0.05/(length(Figure_5E.window_start_list))/4, 0.01/(length(Figure_5E.window_start_list))/4])
    plot_significance_level(Figure_5E.window_start_list,Figure_5E.p_LIP_same(:,loc),1.1,Figure_5E.color_for_ROI(1,:),[0.01,0.05/(length(Figure_5E.window_start_list))/4, 0.01/(length(Figure_5E.window_start_list))/4])
    
    yl=yline(0.5,'--',{'Chance level'},'Color',[0 0 0],'LineWidth',1);
    yl.LabelVerticalAlignment = 'middle';
    yl.LabelHorizontalAlignment = 'right';
    yl.FontSize=12;
    ylim([0.35 1.15])
    text(-575,1,'Within - across')
    text(-575,1.05,'Across')
    text(-575,1.1,'Within')
    box off
    xlabel('Time to response')
    ylabel('Choice decoding accuracy')
    xlim([-600 Figure_5E.window_start_list(end)+150])
    xticks([-600:200:Figure_5E.window_start_list(end)+150])
    subplot(3,4,loc+4)
    hold on
    
    shadedErrorBar(Figure_5E.window_start_list,mean(mean(Figure_5E.FEF.Same_belief_choice(:,:,:,loc),2),3),[prctile(mean(Figure_5E.FEF.Same_belief_choice(:,:,:,loc),2),95,3),prctile(mean(Figure_5E.FEF.Same_belief_choice(:,:,:,loc),2),5,3)]',{'Color',Figure_5E.color_for_ROI(2,:),'LineWidth',2},2)
    shadedErrorBar(Figure_5E.window_start_list,mean(mean(Figure_5E.FEF.Other_belief_choice(:,:,:,loc),2),3),[prctile(mean(Figure_5E.FEF.Other_belief_choice(:,:,:,loc),2),95,3),prctile(mean(Figure_5E.FEF.Other_belief_choice(:,:,:,loc),2),5,3)]',{'--','Color',Figure_5E.color_for_ROI(2,:),'LineWidth',2},2)
    
    plot_significance_level(Figure_5E.window_start_list,Figure_5E.p_FEF_same_other(:,loc),1,Figure_5E.color_for_ROI(2,:),[0.01,0.05/(length(Figure_5E.window_start_list))/4, 0.01/(length(Figure_5E.window_start_list))/4])
    plot_significance_level(Figure_5E.window_start_list,Figure_5E.p_FEF_other(:,loc),1.05,Figure_5E.color_for_ROI(2,:),[0.01,0.05/(length(Figure_5E.window_start_list))/4, 0.01/(length(Figure_5E.window_start_list))/4])
    plot_significance_level(Figure_5E.window_start_list,Figure_5E.p_FEF_same(:,loc),1.1,Figure_5E.color_for_ROI(2,:),[0.01,0.05/(length(Figure_5E.window_start_list))/4, 0.01/(length(Figure_5E.window_start_list))/4])
    
    yl=yline(0.5,'--',{'Chance level'},'Color',[0 0 0],'LineWidth',1);
    yl.LabelVerticalAlignment = 'middle';
    yl.LabelHorizontalAlignment = 'right';
    yl.FontSize=12;
    ylim([0.35 1.15])
    text(-575,1,'Within - across')
    text(-575,1.05,'Across')
    text(-575,1.1,'Within')
    box off
    xlabel('Time to response')
    ylabel('Choice decoding accuracy')
    xlim([-600 Figure_5E.window_start_list(end)+150])
    xticks([-600:200:Figure_5E.window_start_list(end)+150])
    subplot(3,4,loc+8)
    hold on
    
    shadedErrorBar(Figure_5E.window_start_list,mean(mean(Figure_5E.PFC.Same_belief_choice(:,:,:,loc),2),3),[prctile(mean(Figure_5E.PFC.Same_belief_choice(:,:,:,loc),2),95,3),prctile(mean(Figure_5E.PFC.Same_belief_choice(:,:,:,loc),2),5,3)]',{'Color',Figure_5E.color_for_ROI(3,:),'LineWidth',2},2)
    shadedErrorBar(Figure_5E.window_start_list,mean(mean(Figure_5E.PFC.Other_belief_choice(:,:,:,loc),2),3),[prctile(mean(Figure_5E.PFC.Other_belief_choice(:,:,:,loc),2),95,3),prctile(mean(Figure_5E.PFC.Other_belief_choice(:,:,:,loc),2),5,3)]',{'--','Color',Figure_5E.color_for_ROI(3,:),'LineWidth',2},2)
    
    plot_significance_level(Figure_5E.window_start_list,Figure_5E.p_PFC_same_other(:,loc),1,Figure_5E.color_for_ROI(3,:),[0.01,0.05/(length(Figure_5E.window_start_list))/4, 0.01/(length(Figure_5E.window_start_list))/4])
    plot_significance_level(Figure_5E.window_start_list,Figure_5E.p_PFC_other(:,loc),1.05,Figure_5E.color_for_ROI(3,:),[0.01,0.05/(length(Figure_5E.window_start_list))/4, 0.01/(length(Figure_5E.window_start_list))/4])
    plot_significance_level(Figure_5E.window_start_list,Figure_5E.p_PFC_same(:,loc),1.1,Figure_5E.color_for_ROI(3,:),[0.01,0.05/(length(Figure_5E.window_start_list))/4, 0.01/(length(Figure_5E.window_start_list))/4])
    
    yl=yline(0.5,'--',{'Chance level'},'Color',[0 0 0],'LineWidth',1);
    yl.LabelVerticalAlignment = 'middle';
    yl.LabelHorizontalAlignment = 'right';
    yl.FontSize=12;
    ylim([0.35 1.15])
    text(-575,1,'Within - across')
    text(-575,1.05,'Across')
    text(-575,1.1,'Within')
    box off
    xlabel('Time to response')
    ylabel('Choice decoding accuracy')
    xlim([-600 Figure_5E.window_start_list(end)+150])
    xticks([-600:200:Figure_5E.window_start_list(end)+150])
    
end


%% Functions

function plot_significance_level(x,p,a,c,thr)
hold on
for b=1:length(thr)
    this_thr=thr(b);
    for i=1:length(p)-1
        if p(i)<=this_thr && p(i+1)<=this_thr
            plot(x(i):x(i+1),a*ones(1,x(i+1)-x(i)+1),'-','Color',c,'LineWidth',b)
        end
        if i>1
            if p(i-1)>=this_thr && p(i)<=this_thr && p(i+1)>=this_thr
                plot(x(i),a,'.','Color',c,'LineWidth',b)
            end
        end
    end
    if p(end-1)>=this_thr && p(end)<=this_thr
        plot(x(end),a,'.','Color',c,'LineWidth',b)
    end
end


end
