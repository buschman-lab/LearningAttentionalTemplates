%Plot Figure 4

load('Figure_4.mat')

addpath(genpath('Violinplot-Matlab-master'))

%% Figure 4A

figure
subplot(2,2,1)
vs = violinplot([mean_update_model_pos_sess', mean_update_pos_sess'],{'Behavioral','Neural'},'ViolinColor',[0 0 0],'ShowData',true,'ShowMean',true,'ShowMedian',true);
hold on
if ttest(mean_update_model_pos_sess,0,'Tail','right')
    plot(1,0.5,'k*')
end
if ttest(mean_update_pos_sess,0,'Tail','right')
    plot(2,0.5,'k*')
end
yline(0,'--')
box off
ylabel(sprintf('Update toward the chosen color'))
title('+RPE')

subplot(2,2,3)

vs = violinplot(mean_update_model_neg_sess',{'Behavioral'},'ViolinColor',[0 0 0],'ShowData',true,'ShowMean',true,'ShowMedian',true);
hold on
if ttest(mean_update_model_neg_sess,0,'Tail','right')
    plot(1,0.06,'k*')
end
yline(0,'--')
box off
ylabel(sprintf('Update away from the chosen color'))
title('-RPE')

subplot(2,2,4)
vs = violinplot(mean_update_neg_sess',{'Neural'},'ViolinColor',[0 0 0],'ShowData',true,'ShowMean',true,'ShowMedian',true);
if ttest(mean_update_neg_sess,0,'Tail','right')
    plot(1,0.4,'k*')
end
yline(0,'--')
box off
ylabel(sprintf('Update away from  chosen color'))
title('-RPE')

%% Figure 4B

figure
hold on
bar(1,nanmean(All_update_pos(All_true_update_pos<0.1)),'Facecolor',[0.5 0.5 0.5])
bar(2,nanmean(All_update_pos(All_true_update_pos>=0.1 & All_true_update_pos<0.2)),'Facecolor',[0.5 0.5 0.5])
bar(3,nanmean(All_update_pos(All_true_update_pos>=0.2)),'Facecolor',[0.5 0.5 0.5])
errorbar(1,nanmean(All_update_pos(All_true_update_pos<=0.1)),nanstd(All_update_pos(All_true_update_pos<0.1))/sqrt(sum(All_true_update_pos<=0.1)-1),'k','LineWidth',1)
errorbar(2,nanmean(All_update_pos(All_true_update_pos>=0.1 & All_true_update_pos<=0.2)),nanstd(All_update_pos(All_true_update_pos>=0.1 & All_true_update_pos<=0.2))/sqrt(sum(All_true_update_pos>=0.1 & All_true_update_pos<=0.2)-1),'k','LineWidth',1)
errorbar(3,nanmean(All_update_pos(All_true_update_pos>=0.2)),nanstd(All_update_pos(All_true_update_pos>=0.2))/sqrt(sum(All_true_update_pos>=0.2)-1),'k','LineWidth',1)
box off
ylabel('Mean behavioral update toward chosen color')
xlabel('Neural update toward chosen color')
xticks([1 2 3])
xticklabels({'≤0.1','[0.1 0.2]','≥0.2'})

figure
vs = violinplot([z_dist_angular_change_pos_sess'],{'Neural ~ beahvioral'},'ViolinColor',[0 0 0],'ShowData',true,'ShowMean',true,'ShowMedian',true);
if ttest(z_dist_angular_change_pos_sess,0,'Tail','left')
    plot(1,1.5,'k*')
end
yline(0,'--')
box off
ylabel(sprintf('Z-score distance between updates'))

%% Figure 4C

All_prc=prctile(All_prev_rpe_r(All_prev_rpe_r>0),0:10:100);

for i=1:length(All_prc)-4
    mean_All_update(i)=nanmean(abs(All_angular_change_r(All_prev_rpe_r>=All_prc(i) & All_prev_rpe_r<=All_prc(i+4) & All_prev_rpe_r>0)));
    sem_All_update(i)=nanstd(abs(All_angular_change_r(All_prev_rpe_r>=All_prc(i) & All_prev_rpe_r<=All_prc(i+4) & All_prev_rpe_r>0)))/sqrt(sum(All_prev_rpe_r>=All_prc(i) & All_prev_rpe_r<=All_prc(i+4) & All_prev_rpe_r>0)-1);
    
    
    mean_All_true_update(i)=nanmean(abs(All_true_angular_change_r(All_prev_rpe_r>=All_prc(i) & All_prev_rpe_r<=All_prc(i+4) & All_prev_rpe_r>0)));
    sem_All_true_update(i)=nanstd(abs(All_true_angular_change_r(All_prev_rpe_r>=All_prc(i) & All_prev_rpe_r<=All_prc(i+4) & All_prev_rpe_r>0)))/sqrt(sum(All_prev_rpe_r>=All_prc(i) & All_prev_rpe_r<=All_prc(i+4) & All_prev_rpe_r>0)-1);
    
end


figure
subplot(2,2,2)
shadedErrorBar(All_prc(3:end-2),mean_All_update,sem_All_update,{'-','color',[0 0 0],'LineWidth',2},2)
xlabel('RPE')
ylabel(sprintf('Neural update (RPE > 0)'))
box off
title('All')

subplot(2,2,1)
shadedErrorBar(All_prc(3:end-2),mean_All_true_update,sem_All_true_update,{'-','color',[0 0 0],'LineWidth',2},2)
xlabel('RPE')
ylabel(sprintf('Behavioral update (RPE > 0)'))
box off
title('All')

All_angle=prctile(abs(All_dist_prev_decoded_angle_cc_r(All_prev_rpe_r>0)),0:10:100);
All_angle2=prctile(abs(All_dist_prev_true_peak_cc_r(All_prev_rpe_r>0)),0:10:100);

for i=1:length(All_angle)-4
    mean_All_update_angle(i)=nanmean(abs(All_angular_change_r(abs(All_dist_prev_decoded_angle_cc_r)>=All_angle(i) & abs(All_dist_prev_decoded_angle_cc_r)<=All_angle(i+4) & All_prev_rpe_r>0)));
    sem_All_update_angle(i)=nanstd(abs(All_angular_change_r(abs(All_dist_prev_decoded_angle_cc_r)>=All_angle(i) & abs(All_dist_prev_decoded_angle_cc_r)<=All_angle(i+4) & All_prev_rpe_r>0)))/sqrt(sum(abs(All_dist_prev_decoded_angle_cc_r)>=All_angle(i) & abs(All_dist_prev_decoded_angle_cc_r)<=All_angle(i+4) & All_prev_rpe_r>0)-1);
    
    mean_All_true_update_angle(i)=nanmean(abs(All_true_angular_change_r(abs(All_dist_prev_true_peak_cc_r)>=All_angle2(i) & abs(All_dist_prev_true_peak_cc_r)<=All_angle2(i+4) & All_prev_rpe_r>0)));
    sem_All_true_update_angle(i)=nanstd(abs(All_true_angular_change_r(abs(All_dist_prev_true_peak_cc_r)>=All_angle2(i) & abs(All_dist_prev_true_peak_cc_r)<=All_angle2(i+4) & All_prev_rpe_r>0)))/sqrt(sum(abs(All_dist_prev_true_peak_cc_r)>=All_angle2(i) & abs(All_dist_prev_true_peak_cc_r)<=All_angle2(i+4) & All_prev_rpe_r>0)-1);
    
end

subplot(2,2,4)
shadedErrorBar(All_angle(3:end-2),mean_All_update_angle,sem_All_update_angle,{'-','color',[0 0 0],'LineWidth',2},2)
xlabel(sprintf('Neural |d(ET(t-1), CC(t-1))|'))
ylabel(sprintf('Neural update magnitude'))
box off
xticks([0 pi/4 pi/2 3*pi/4])
xticklabels({'0','π/4','π/2','3π/4'})
xlim([0 pi-pi/8])

subplot(2,2,3)
shadedErrorBar(All_angle(3:end-2),mean_All_true_update_angle,sem_All_true_update_angle,{'-','color',[0 0 0],'LineWidth',2},2)
xlabel(sprintf('Model |d(ET(t-1), CC(t-1))|'))
ylabel(sprintf('Behavioral update magnitude'))
box off
xticks([0 pi/4 pi/2 3*pi/4])
xticklabels({'0','π/4','π/2','3π/4'})
xlim([0 pi-pi/8])

figure

subplot(2,1,1)
vs = violinplot([pearson_corr_rpe_model', pearson_corr_rpe_data'],{'Behavioral','Neural'},'ViolinColor',[0 0 0],'ShowData',true,'ShowMean',true,'ShowMedian',true);
if ttest(pearson_corr_rpe_model,0,'Tail','right')
    plot(1,0.8,'k*')
end
if ttest(pearson_corr_rpe_data,0,'Tail','right')
    plot(2,0.8,'k*')
end
yline(0,'--')
box off
ylabel(sprintf('Pearson correlation with |RPE|'))


subplot(2,1,2)
vs = violinplot([pearson_corr_angle_model', pearson_corr_angle_data'],{'Behavioral','Neural'},'ViolinColor',[0 0 0],'ShowData',true,'ShowMean',true,'ShowMedian',true);
hold on
if ttest(pearson_corr_angle_model,0,'Tail','right')
    plot(1,0.8,'k*')
end
if ttest(pearson_corr_angle_data,0,'Tail','right')
    plot(2,0.8,'k*')
end
yline(0,'--')
box off
ylabel(sprintf('Pearson correlation with |d(ET(t-1),CC(t-1))|'))


