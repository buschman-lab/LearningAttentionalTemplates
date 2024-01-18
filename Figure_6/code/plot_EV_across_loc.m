clear all

window_size=200;
for i=1:25
    initial_window=-400;
    event_list{i}='target';
    window_start_list(i)=initial_window+(i-1)*50;
end

event='target';

fsroot='/Volumes/buschman';
task='Learning_Attentional_Templates';
subtask='exploreexploit/Restricted_FEF/Reset_RW_model';

subsubtask_save=sprintf('Classifier_%s_%d_to_%s_%d',event,window_start_list(1),event,window_start_list(end));
dirstem2 = fullfile('Projects',task,'Analysis','Electrophy_analysis',subtask,subsubtask_save);
load_path=fullfile(fsroot,dirstem2);


count_ROI=ones(length(window_start_list),1);

arrayID_list=[3,8,12,14,18,20,22,26,28,30,32,34,36,38,40,42,44];

color_for_ROI=[0 0.447 0.741; 0.85 0.325 0.098; 0.929 0.694 0.125]; %LIP FEF PFC


color_each_loc(1,:,1)=color_for_ROI(1,:)+0.25;
color_each_loc(1,:,2)=color_for_ROI(1,:)+0.19;
color_each_loc(1,:,3)=color_for_ROI(1,:)+0.06;
color_each_loc(1,:,4)=color_for_ROI(1,:);

color_each_loc(2,:,1)=color_for_ROI(2,:)+0.15;
color_each_loc(2,:,2)=color_for_ROI(2,:)+0.09;
color_each_loc(2,:,3)=color_for_ROI(2,:)-0.03;
color_each_loc(2,:,4)=color_for_ROI(2,:)-0.09;

color_each_loc(3,:,1)=color_for_ROI(3,:)+0.07;
color_each_loc(3,:,2)=color_for_ROI(3,:)+0.01;
color_each_loc(3,:,3)=color_for_ROI(3,:)-0.06;
color_each_loc(3,:,4)=color_for_ROI(3,:)-0.12;

%%

%LIP
load_name=sprintf('glm_explained_variance_split_%s','LIP');
load(fullfile(load_path,load_name),'contra*','p*')

LIP.contra_z_R_square_across_loc=contra_z_R_square_across_loc;
LIP.contra_contra_z_R_square_across_loc=contra_contra_z_R_square_across_loc;
LIP.contra_ipsi_z_R_square_across_loc=contra_ipsi_z_R_square_across_loc;
LIP.contra_ipsi_diag_z_R_square_across_loc=contra_ipsi_diag_z_R_square_across_loc;

LIP.p_contra_z_R_square_across_loc=p_contra_z_R_square_across_loc;
LIP.p_contra_contra_z_R_square_across_loc=p_contra_contra_z_R_square_across_loc;
LIP.p_contra_ipsi_z_R_square_across_loc=p_contra_ipsi_z_R_square_across_loc;
LIP.p_contra_ipsi_diag_z_R_square_across_loc=p_contra_ipsi_diag_z_R_square_across_loc;
  
LIP.p_contra_to_contra_contra=p_contra_to_contra_contra;
LIP.p_contra_to_contra_ipsi=p_contra_to_contra_ipsi;
LIP.p_contra_to_contra_ipsi_diag=p_contra_to_contra_ipsi_diag;


clear contra* p*

%FEF
load_name=sprintf('glm_explained_variance_split_%s','FEF');
load(fullfile(load_path,load_name),'contra*','p*')

FEF.contra_z_R_square_across_loc=contra_z_R_square_across_loc;
FEF.contra_contra_z_R_square_across_loc=contra_contra_z_R_square_across_loc;
FEF.contra_ipsi_z_R_square_across_loc=contra_ipsi_z_R_square_across_loc;
FEF.contra_ipsi_diag_z_R_square_across_loc=contra_ipsi_diag_z_R_square_across_loc;

FEF.p_contra_z_R_square_across_loc=p_contra_z_R_square_across_loc;
FEF.p_contra_contra_z_R_square_across_loc=p_contra_contra_z_R_square_across_loc;
FEF.p_contra_ipsi_z_R_square_across_loc=p_contra_ipsi_z_R_square_across_loc;
FEF.p_contra_ipsi_diag_z_R_square_across_loc=p_contra_ipsi_diag_z_R_square_across_loc;
  
FEF.p_contra_to_contra_contra=p_contra_to_contra_contra;
FEF.p_contra_to_contra_ipsi=p_contra_to_contra_ipsi;
FEF.p_contra_to_contra_ipsi_diag=p_contra_to_contra_ipsi_diag;


clear contra* p*

%PFC
load_name=sprintf('glm_explained_variance_split_%s','PFC');
load(fullfile(load_path,load_name),'contra*','p*')

PFC.contra_z_R_square_across_loc=contra_z_R_square_across_loc;
PFC.contra_contra_z_R_square_across_loc=contra_contra_z_R_square_across_loc;
PFC.contra_ipsi_z_R_square_across_loc=contra_ipsi_z_R_square_across_loc;
PFC.contra_ipsi_diag_z_R_square_across_loc=contra_ipsi_diag_z_R_square_across_loc;

PFC.p_contra_z_R_square_across_loc=p_contra_z_R_square_across_loc;
PFC.p_contra_contra_z_R_square_across_loc=p_contra_contra_z_R_square_across_loc;
PFC.p_contra_ipsi_z_R_square_across_loc=p_contra_ipsi_z_R_square_across_loc;
PFC.p_contra_ipsi_diag_z_R_square_across_loc=p_contra_ipsi_diag_z_R_square_across_loc;
  
PFC.p_contra_to_contra_contra=p_contra_to_contra_contra;
PFC.p_contra_to_contra_ipsi=p_contra_to_contra_ipsi;
PFC.p_contra_to_contra_ipsi_diag=p_contra_to_contra_ipsi_diag;


clear contra* p*


%%

  LIP.auc_r_chosen(1,:)=sum(LIP.contra_z_R_square_across_loc(:,4:end),2);
  LIP.auc_r_chosen(2,:)=sum(LIP.contra_contra_z_R_square_across_loc(:,4:end),2);
  LIP.auc_r_chosen(3,:)=sum(LIP.contra_ipsi_z_R_square_across_loc(:,4:end),2);
  LIP.auc_r_chosen(4,:)=sum(LIP.contra_ipsi_diag_z_R_square_across_loc(:,4:end),2);
  
  FEF.auc_r_chosen(1,:)=sum(FEF.contra_z_R_square_across_loc(:,4:end),2);
  FEF.auc_r_chosen(2,:)=sum(FEF.contra_contra_z_R_square_across_loc(:,4:end),2);
  FEF.auc_r_chosen(3,:)=sum(FEF.contra_ipsi_z_R_square_across_loc(:,4:end),2);
  FEF.auc_r_chosen(4,:)=sum(FEF.contra_ipsi_diag_z_R_square_across_loc(:,4:end),2);

  PFC.auc_r_chosen(1,:)=sum(PFC.contra_z_R_square_across_loc(:,4:end),2);
  PFC.auc_r_chosen(2,:)=sum(PFC.contra_contra_z_R_square_across_loc(:,4:end),2);
  PFC.auc_r_chosen(3,:)=sum(PFC.contra_ipsi_z_R_square_across_loc(:,4:end),2);
  PFC.auc_r_chosen(4,:)=sum(PFC.contra_ipsi_diag_z_R_square_across_loc(:,4:end),2);
  
%%
p_LIP_auc_across_loc=NaN(4,4);
p_FEF_auc_across_loc=NaN(4,4);
p_PFC_auc_across_loc=NaN(4,4);

for i=1:3
    for j=i+1:4
        [~, p_LIP_auc_across_loc(i,j)]=ttest(LIP.auc_r_chosen(i,:)-LIP.auc_r_chosen(j,:),0);
        [~, p_FEF_auc_across_loc(i,j)]=ttest(FEF.auc_r_chosen(i,:)-FEF.auc_r_chosen(j,:),0);
        [~, p_PFC_auc_across_loc(i,j)]=ttest(PFC.auc_r_chosen(i,:)-PFC.auc_r_chosen(j,:),0);
    end
end


figure
subplot(3,1,1)
shadedErrorBar(1:4,nanmean(LIP.auc_r_chosen,2),nanstd(LIP.auc_r_chosen,0,2)./sqrt(sum(~isnan(LIP.auc_r_chosen),2)-1),{'-','color',color_for_ROI(1,:),'LineWidth',2},2)
box off
xlabel('Distance')
ylabel('AUC correlation chosen value vectors')
ax=gca;
ax.XAxis.FontSize=12;
ax.YAxis.FontSize=12;

subplot(3,1,2)
shadedErrorBar(1:4,nanmean(FEF.auc_r_chosen,2),nanstd(FEF.auc_r_chosen,0,2)./sqrt(sum(~isnan(FEF.auc_r_chosen),2)-1),{'-','color',color_for_ROI(2,:),'LineWidth',2},2)
box off
xlabel('Distance')
ylabel('AUC correlation chosen value vectors')
ax=gca;
ax.XAxis.FontSize=12;
ax.YAxis.FontSize=12;

subplot(3,1,3)
shadedErrorBar(1:4,nanmean(PFC.auc_r_chosen,2),nanstd(PFC.auc_r_chosen,0,2)./sqrt(sum(~isnan(PFC.auc_r_chosen),2)-1),{'-','color',color_for_ROI(3,:),'LineWidth',2},2)
box off
xlabel('Distance')
ylabel('AUC correlation chosen value vectors')
ax=gca;
ax.XAxis.FontSize=12;
ax.YAxis.FontSize=12;


%%

figure
subplot(1,4,1)
hold on

shadedErrorBar(window_start_list(4:end),nanmean(LIP.contra_z_R_square_across_loc(:,4:end),1),nanstd(LIP.contra_z_R_square_across_loc(:,4:end),0,1)./sqrt(sum(~isnan(LIP.contra_z_R_square_across_loc(:,4:end)),1)-1),{'-','color',color_for_ROI(1,:),'LineWidth',2},2)
shadedErrorBar(window_start_list(4:end),nanmean(FEF.contra_z_R_square_across_loc(:,4:end),1),nanstd(FEF.contra_z_R_square_across_loc(:,4:end),0,1)./sqrt(sum(~isnan(FEF.contra_z_R_square_across_loc(:,4:end)),1)-1),{'-','color',color_for_ROI(2,:),'LineWidth',2},2)
shadedErrorBar(window_start_list(4:end),nanmean(PFC.contra_z_R_square_across_loc(:,4:end),1),nanstd(PFC.contra_z_R_square_across_loc(:,4:end),0,1)./sqrt(sum(~isnan(PFC.contra_z_R_square_across_loc(:,4:end)),1)-1),{'-','color',color_for_ROI(3,:),'LineWidth',2},2)

plot_significance_level(window_start_list,LIP.p_contra_z_R_square_across_loc,1.05,color_for_ROI(1,:),[0.01,0.05/(length(window_start_list)-4), 0.01/(length(window_start_list)-4)])
plot_significance_level(window_start_list,FEF.p_contra_z_R_square_across_loc,1.1,color_for_ROI(2,:),[0.01,0.05/(length(window_start_list)-4), 0.01/(length(window_start_list)-4)])
plot_significance_level(window_start_list,PFC.p_contra_z_R_square_across_loc,1.15,color_for_ROI(3,:),[0.01,0.05/(length(window_start_list)-4), 0.01/(length(window_start_list)-4)])

subplot(1,4,2)
hold on
shadedErrorBar(window_start_list(4:end),nanmean(LIP.contra_contra_z_R_square_across_loc(:,4:end),1),nanstd(LIP.contra_contra_z_R_square_across_loc(:,4:end),0,1)./sqrt(sum(~isnan(LIP.contra_contra_z_R_square_across_loc(:,4:end)),1)-1),{'-','color',color_for_ROI(1,:),'LineWidth',2},2)
shadedErrorBar(window_start_list(4:end),nanmean(FEF.contra_contra_z_R_square_across_loc(:,4:end),1),nanstd(FEF.contra_contra_z_R_square_across_loc(:,4:end),0,1)./sqrt(sum(~isnan(FEF.contra_contra_z_R_square_across_loc(:,4:end)),1)-1),{'-','color',color_for_ROI(2,:),'LineWidth',2},2)
shadedErrorBar(window_start_list(4:end),nanmean(PFC.contra_contra_z_R_square_across_loc(:,4:end),1),nanstd(PFC.contra_contra_z_R_square_across_loc(:,4:end),0,1)./sqrt(sum(~isnan(PFC.contra_contra_z_R_square_across_loc(:,4:end)),1)-1),{'-','color',color_for_ROI(3,:),'LineWidth',2},2)

plot_significance_level(window_start_list,LIP.p_contra_contra_z_R_square_across_loc,1.05,color_for_ROI(1,:),[0.01,0.05/(length(window_start_list)-4), 0.01/(length(window_start_list)-4)])
plot_significance_level(window_start_list,FEF.p_contra_contra_z_R_square_across_loc,1.1,color_for_ROI(2,:),[0.01,0.05/(length(window_start_list)-4), 0.01/(length(window_start_list)-4)])
plot_significance_level(window_start_list,PFC.p_contra_contra_z_R_square_across_loc,1.15,color_for_ROI(3,:),[0.01,0.05/(length(window_start_list)-4), 0.01/(length(window_start_list)-4)])

plot_significance_level(window_start_list,LIP.p_contra_to_contra_contra,0.85,color_for_ROI(1,:),[0.01,0.05/(length(window_start_list)-4), 0.01/(length(window_start_list)-4)])
plot_significance_level(window_start_list,FEF.p_contra_to_contra_contra,0.9,color_for_ROI(2,:),[0.01,0.05/(length(window_start_list)-4), 0.01/(length(window_start_list)-4)])
plot_significance_level(window_start_list,PFC.p_contra_to_contra_contra,0.95,color_for_ROI(3,:),[0.01,0.05/(length(window_start_list)-4), 0.01/(length(window_start_list)-4)])


subplot(1,4,3)
hold on
shadedErrorBar(window_start_list(4:end),nanmean(LIP.contra_ipsi_z_R_square_across_loc(:,4:end),1),nanstd(LIP.contra_ipsi_z_R_square_across_loc(:,4:end),0,1)./sqrt(sum(~isnan(LIP.contra_ipsi_z_R_square_across_loc(:,4:end)),1)-1),{'-','color',color_for_ROI(1,:),'LineWidth',2},2)
shadedErrorBar(window_start_list(4:end),nanmean(FEF.contra_ipsi_z_R_square_across_loc(:,4:end),1),nanstd(FEF.contra_ipsi_z_R_square_across_loc(:,4:end),0,1)./sqrt(sum(~isnan(FEF.contra_ipsi_z_R_square_across_loc(:,4:end)),1)-1),{'-','color',color_for_ROI(2,:),'LineWidth',2},2)
shadedErrorBar(window_start_list(4:end),nanmean(PFC.contra_ipsi_z_R_square_across_loc(:,4:end),1),nanstd(PFC.contra_ipsi_z_R_square_across_loc(:,4:end),0,1)./sqrt(sum(~isnan(PFC.contra_ipsi_z_R_square_across_loc(:,4:end)),1)-1),{'-','color',color_for_ROI(3,:),'LineWidth',2},2)

plot_significance_level(window_start_list,LIP.p_contra_ipsi_z_R_square_across_loc,1.05,color_for_ROI(1,:),[0.01,0.05/(length(window_start_list)-4), 0.01/(length(window_start_list)-4)])
plot_significance_level(window_start_list,FEF.p_contra_ipsi_z_R_square_across_loc,1.1,color_for_ROI(2,:),[0.01,0.05/(length(window_start_list)-4), 0.01/(length(window_start_list)-4)])
plot_significance_level(window_start_list,PFC.p_contra_ipsi_z_R_square_across_loc,1.15,color_for_ROI(3,:),[0.01,0.05/(length(window_start_list)-4), 0.01/(length(window_start_list)-4)])

plot_significance_level(window_start_list,LIP.p_contra_to_contra_ipsi,0.85,color_for_ROI(1,:),[0.01,0.05/(length(window_start_list)-4), 0.01/(length(window_start_list)-4)])
plot_significance_level(window_start_list,FEF.p_contra_to_contra_ipsi,0.9,color_for_ROI(2,:),[0.01,0.05/(length(window_start_list)-4), 0.01/(length(window_start_list)-4)])
plot_significance_level(window_start_list,PFC.p_contra_to_contra_ipsi,0.95,color_for_ROI(3,:),[0.01,0.05/(length(window_start_list)-4), 0.01/(length(window_start_list)-4)])

subplot(1,4,4)
hold on
shadedErrorBar(window_start_list(4:end),nanmean(LIP.contra_ipsi_diag_z_R_square_across_loc(:,4:end),1),nanstd(LIP.contra_ipsi_diag_z_R_square_across_loc(:,4:end),0,1)./sqrt(sum(~isnan(LIP.contra_ipsi_diag_z_R_square_across_loc(:,4:end)),1)-1),{'-','color',color_for_ROI(1,:),'LineWidth',2},2)
shadedErrorBar(window_start_list(4:end),nanmean(FEF.contra_ipsi_diag_z_R_square_across_loc(:,4:end),1),nanstd(FEF.contra_ipsi_diag_z_R_square_across_loc(:,4:end),0,1)./sqrt(sum(~isnan(FEF.contra_ipsi_diag_z_R_square_across_loc(:,4:end)),1)-1),{'-','color',color_for_ROI(2,:),'LineWidth',2},2)
shadedErrorBar(window_start_list(4:end),nanmean(PFC.contra_ipsi_diag_z_R_square_across_loc(:,4:end),1),nanstd(PFC.contra_ipsi_diag_z_R_square_across_loc(:,4:end),0,1)./sqrt(sum(~isnan(PFC.contra_ipsi_diag_z_R_square_across_loc(:,4:end)),1)-1),{'-','color',color_for_ROI(3,:),'LineWidth',2},2)

plot_significance_level(window_start_list,LIP.p_contra_ipsi_diag_z_R_square_across_loc,1.05,color_for_ROI(1,:),[0.01,0.05/(length(window_start_list)-4), 0.01/(length(window_start_list)-4)])
plot_significance_level(window_start_list,FEF.p_contra_ipsi_diag_z_R_square_across_loc,1.1,color_for_ROI(2,:),[0.01,0.05/(length(window_start_list)-4), 0.01/(length(window_start_list)-4)])
plot_significance_level(window_start_list,PFC.p_contra_ipsi_diag_z_R_square_across_loc,1.15,color_for_ROI(3,:),[0.01,0.05/(length(window_start_list)-4), 0.01/(length(window_start_list)-4)])

plot_significance_level(window_start_list,LIP.p_contra_to_contra_ipsi_diag,0.85,color_for_ROI(1,:),[0.01,0.05/(length(window_start_list)-4), 0.01/(length(window_start_list)-4)])
plot_significance_level(window_start_list,FEF.p_contra_to_contra_ipsi_diag,0.9,color_for_ROI(2,:),[0.01,0.05/(length(window_start_list)-4), 0.01/(length(window_start_list)-4)])
plot_significance_level(window_start_list,PFC.p_contra_to_contra_ipsi_diag,0.95,color_for_ROI(3,:),[0.01,0.05/(length(window_start_list)-4), 0.01/(length(window_start_list)-4)])


%%
figure
subplot(1,4,1)
hold on

shadedErrorBar(window_start_list(4:end),nanmean(LIP.contra_z_R_square_across_loc(:,4:end),1),nanstd(LIP.contra_z_R_square_across_loc(:,4:end),0,1)./sqrt(sum(~isnan(LIP.contra_z_R_square_across_loc(:,4:end)),1)-1),{'-','color',color_each_loc(1,:,1),'LineWidth',2},2)
shadedErrorBar(window_start_list(4:end),nanmean(LIP.contra_contra_z_R_square_across_loc(:,4:end),1),nanstd(LIP.contra_contra_z_R_square_across_loc(:,4:end),0,1)./sqrt(sum(~isnan(LIP.contra_contra_z_R_square_across_loc(:,4:end)),1)-1),{'-','color',color_each_loc(1,:,2),'LineWidth',2},2)
shadedErrorBar(window_start_list(4:end),nanmean(LIP.contra_ipsi_z_R_square_across_loc(:,4:end),1),nanstd(LIP.contra_ipsi_z_R_square_across_loc(:,4:end),0,1)./sqrt(sum(~isnan(LIP.contra_ipsi_z_R_square_across_loc(:,4:end)),1)-1),{'-','color',color_each_loc(1,:,3),'LineWidth',2},2)
shadedErrorBar(window_start_list(4:end),nanmean(LIP.contra_ipsi_diag_z_R_square_across_loc(:,4:end),1),nanstd(LIP.contra_ipsi_diag_z_R_square_across_loc(:,4:end),0,1)./sqrt(sum(~isnan(LIP.contra_ipsi_diag_z_R_square_across_loc(:,4:end)),1)-1),{'-','color',color_each_loc(1,:,4),'LineWidth',2},2)

plot_significance_level(window_start_list,LIP.p_contra_z_R_square_across_loc,1.05,color_each_loc(1,:,1),[0.01,0.05/(length(window_start_list)-4), 0.01/(length(window_start_list)-4)])
plot_significance_level(window_start_list,LIP.p_contra_contra_z_R_square_across_loc,1,color_each_loc(1,:,2),[0.01,0.05/(length(window_start_list)-4), 0.01/(length(window_start_list)-4)])
plot_significance_level(window_start_list,LIP.p_contra_ipsi_z_R_square_across_loc,0.95,color_each_loc(1,:,3),[0.01,0.05/(length(window_start_list)-4), 0.01/(length(window_start_list)-4)])
plot_significance_level(window_start_list,LIP.p_contra_ipsi_diag_z_R_square_across_loc,0.9,color_each_loc(1,:,4),[0.01,0.05/(length(window_start_list)-4), 0.01/(length(window_start_list)-4)])

text(-575,1.05,'0')
text(-575,1.0,'1')
text(-575,0.95,'2')
text(-575,0.9,'3')

yl.LabelVerticalAlignment = 'middle';
yl.LabelHorizontalAlignment = 'right';
yl.FontSize=12;
box off
xlabel('Time to target')
xlim([-600 window_start_list(end)+150])
xticks([-600:200:window_start_list(end)+150])
yline(0,'--')
ax.YAxis.FontSize=12;
ylabel('z-scored EV chosen value')
ylim([-0.2 1.1])
ax.XAxis.FontSize=12;
ax.YAxis.FontSize=12;


subplot(1,4,2)
hold on
shadedErrorBar(window_start_list(4:end),nanmean(FEF.contra_z_R_square_across_loc(:,4:end),1),nanstd(FEF.contra_z_R_square_across_loc(:,4:end),0,1)./sqrt(sum(~isnan(FEF.contra_z_R_square_across_loc(:,4:end)),1)-1),{'-','color',color_each_loc(2,:,1),'LineWidth',2},2)
shadedErrorBar(window_start_list(4:end),nanmean(FEF.contra_contra_z_R_square_across_loc(:,4:end),1),nanstd(FEF.contra_contra_z_R_square_across_loc(:,4:end),0,1)./sqrt(sum(~isnan(FEF.contra_contra_z_R_square_across_loc(:,4:end)),1)-1),{'-','color',color_each_loc(2,:,2),'LineWidth',2},2)
shadedErrorBar(window_start_list(4:end),nanmean(FEF.contra_ipsi_z_R_square_across_loc(:,4:end),1),nanstd(FEF.contra_ipsi_z_R_square_across_loc(:,4:end),0,1)./sqrt(sum(~isnan(FEF.contra_ipsi_z_R_square_across_loc(:,4:end)),1)-1),{'-','color',color_each_loc(2,:,3),'LineWidth',2},2)
shadedErrorBar(window_start_list(4:end),nanmean(FEF.contra_ipsi_diag_z_R_square_across_loc(:,4:end),1),nanstd(FEF.contra_ipsi_diag_z_R_square_across_loc(:,4:end),0,1)./sqrt(sum(~isnan(FEF.contra_ipsi_diag_z_R_square_across_loc(:,4:end)),1)-1),{'-','color',color_each_loc(2,:,4),'LineWidth',2},2)

plot_significance_level(window_start_list,FEF.p_contra_z_R_square_across_loc,1.05,color_each_loc(2,:,1),[0.01,0.05/(length(window_start_list)-4), 0.01/(length(window_start_list)-4)])
plot_significance_level(window_start_list,FEF.p_contra_contra_z_R_square_across_loc,1,color_each_loc(2,:,2),[0.01,0.05/(length(window_start_list)-4), 0.01/(length(window_start_list)-4)])
plot_significance_level(window_start_list,FEF.p_contra_ipsi_z_R_square_across_loc,0.95,color_each_loc(2,:,3),[0.01,0.05/(length(window_start_list)-4), 0.01/(length(window_start_list)-4)])
plot_significance_level(window_start_list,FEF.p_contra_ipsi_diag_z_R_square_across_loc,0.9,color_each_loc(2,:,4),[0.01,0.05/(length(window_start_list)-4), 0.01/(length(window_start_list)-4)])

text(-575,1.05,'0')
text(-575,1.0,'1')
text(-575,0.95,'2')
text(-575,0.9,'3')

yl.LabelVerticalAlignment = 'middle';
yl.LabelHorizontalAlignment = 'right';
yl.FontSize=12;
box off
xlabel('Time to target')
xlim([-600 window_start_list(end)+150])
xticks([-600:200:window_start_list(end)+150])
yline(0,'--')
ax.YAxis.FontSize=12;
ylabel('z-scored EV chosen value')
ylim([-0.2 1.1])
ax.XAxis.FontSize=12;
ax.YAxis.FontSize=12;


subplot(1,4,3)
hold on
shadedErrorBar(window_start_list(4:end),nanmean(PFC.contra_z_R_square_across_loc(:,4:end),1),nanstd(PFC.contra_z_R_square_across_loc(:,4:end),0,1)./sqrt(sum(~isnan(PFC.contra_z_R_square_across_loc(:,4:end)),1)-1),{'-','color',color_each_loc(3,:,1),'LineWidth',2},2)
shadedErrorBar(window_start_list(4:end),nanmean(PFC.contra_contra_z_R_square_across_loc(:,4:end),1),nanstd(PFC.contra_contra_z_R_square_across_loc(:,4:end),0,1)./sqrt(sum(~isnan(PFC.contra_contra_z_R_square_across_loc(:,4:end)),1)-1),{'-','color',color_each_loc(3,:,2),'LineWidth',2},2)
shadedErrorBar(window_start_list(4:end),nanmean(PFC.contra_ipsi_z_R_square_across_loc(:,4:end),1),nanstd(PFC.contra_ipsi_z_R_square_across_loc(:,4:end),0,1)./sqrt(sum(~isnan(PFC.contra_ipsi_z_R_square_across_loc(:,4:end)),1)-1),{'-','color',color_each_loc(3,:,3),'LineWidth',2},2)
shadedErrorBar(window_start_list(4:end),nanmean(PFC.contra_ipsi_diag_z_R_square_across_loc(:,4:end),1),nanstd(PFC.contra_ipsi_diag_z_R_square_across_loc(:,4:end),0,1)./sqrt(sum(~isnan(PFC.contra_ipsi_diag_z_R_square_across_loc(:,4:end)),1)-1),{'-','color',color_each_loc(3,:,4),'LineWidth',2},2)

plot_significance_level(window_start_list,PFC.p_contra_z_R_square_across_loc,1.05,color_each_loc(3,:,1),[0.01,0.05/(length(window_start_list)-4), 0.01/(length(window_start_list)-4)])
plot_significance_level(window_start_list,PFC.p_contra_contra_z_R_square_across_loc,1,color_each_loc(3,:,2),[0.01,0.05/(length(window_start_list)-4), 0.01/(length(window_start_list)-4)])
plot_significance_level(window_start_list,PFC.p_contra_ipsi_z_R_square_across_loc,0.95,color_each_loc(3,:,3),[0.01,0.05/(length(window_start_list)-4), 0.01/(length(window_start_list)-4)])
plot_significance_level(window_start_list,PFC.p_contra_ipsi_diag_z_R_square_across_loc,0.9,color_each_loc(3,:,4),[0.01,0.05/(length(window_start_list)-4), 0.01/(length(window_start_list)-4)])

text(-575,1.05,'0')
text(-575,1.0,'1')
text(-575,0.95,'2')
text(-575,0.9,'3')

yl.LabelVerticalAlignment = 'middle';
yl.LabelHorizontalAlignment = 'right';
yl.FontSize=12;
box off
xlabel('Time to target')
xlim([-600 window_start_list(end)+150])
xticks([-600:200:window_start_list(end)+150])
yline(0,'--')
ax.YAxis.FontSize=12;
ylabel('z-scored EV chosen value')
ylim([-0.2 1.1])
ax.XAxis.FontSize=12;
ax.YAxis.FontSize=12;


subplot(1,4,4)
hold on
shadedErrorBar(1:4,nanmean(LIP.auc_r_chosen,2),nanstd(LIP.auc_r_chosen,0,2)./sqrt(sum(~isnan(LIP.auc_r_chosen),2)-1),{'-','color',color_for_ROI(1,:),'LineWidth',2},2)
shadedErrorBar(1:4,nanmean(FEF.auc_r_chosen,2),nanstd(FEF.auc_r_chosen,0,2)./sqrt(sum(~isnan(FEF.auc_r_chosen),2)-1),{'-','color',color_for_ROI(2,:),'LineWidth',2},2)
shadedErrorBar(1:4,nanmean(PFC.auc_r_chosen,2),nanstd(PFC.auc_r_chosen,0,2)./sqrt(sum(~isnan(PFC.auc_r_chosen),2)-1),{'-','color',color_for_ROI(3,:),'LineWidth',2},2)
box off
xlabel('Distance')
ylabel('AUC correlation chosen value vectors')
ax=gca;
ax.XAxis.FontSize=12;
ax.YAxis.FontSize=12;


%% Functions

function plot_significance_level(x,p,a,c,thr)
hold on
for b=1:length(thr)
    this_thr=thr(b);
    for i=4:length(p)-1
        if p(i)<this_thr && p(i+1)<this_thr
            plot(x(i):x(i+1),a*ones(1,x(i+1)-x(i)+1),'-','Color',c,'LineWidth',b)
        end
        if p(i-1)>=this_thr && p(i)<this_thr && p(i+1)>=this_thr
            plot(x(i),a,'.','Color',c,'LineWidth',b)
        end
    end
end
end



function p = z_test_function_bootstrap_2tail(dist,null)

if sum(~isnan(dist))>2
    
    dist=dist(~isnan(dist));
    m = mean(dist);
    s = std(dist);
    z = (m-null)/s;
    p=2*(1-normcdf(abs(z)));
    
else
    p=NaN;
end

end



