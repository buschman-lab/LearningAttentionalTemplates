
clear all

event='target';

window_size=200;
for i=1:25
    initial_window=-400;
    event_list{i}='target';
    window_start_list(i)=initial_window+(i-1)*50;
end

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


%% load

%LIP
load_name=sprintf('glm_explained_variance_split_boot_pref_%s','LIP');
load(fullfile(load_path,load_name),'r_chosen_pref_boot','p_chosen_pref_boot',...
    'r_unchosen_pref_boot','p_unchosen_ref_boot',...
    'r_chosen_unpref_boot','p_chosen_unpref_boot',...
    'r_unchosen_unpref_boot','p_unchosen_unpref_boot');

LIP.r_chosen_pref_boot=r_chosen_pref_boot;
LIP.p_chosen_pref_boot=p_chosen_pref_boot;
LIP.r_unchosen_pref_boot=r_unchosen_pref_boot;
LIP.p_unchosen_pref_boot=p_unchosen_ref_boot;

LIP.r_chosen_unpref_boot=r_chosen_unpref_boot;
LIP.p_chosen_unpref_boot=p_chosen_unpref_boot;
LIP.r_unchosen_unpref_boot=r_unchosen_unpref_boot;
LIP.p_unchosen_unpref_boot=p_unchosen_unpref_boot;

clear *boot* *diss*

%FEF
load_name=sprintf('glm_explained_variance_split_boot_pref_%s','FEF');
load(fullfile(load_path,load_name),'r_chosen_pref_boot','p_chosen_pref_boot',...
    'r_unchosen_pref_boot','p_unchosen_ref_boot',...
    'r_chosen_unpref_boot','p_chosen_unpref_boot',...
    'r_unchosen_unpref_boot','p_unchosen_unpref_boot');

FEF.r_chosen_pref_boot=r_chosen_pref_boot;
FEF.p_chosen_pref_boot=p_chosen_pref_boot;
FEF.r_unchosen_pref_boot=r_unchosen_pref_boot;
FEF.p_unchosen_pref_boot=p_unchosen_ref_boot;

FEF.r_chosen_unpref_boot=r_chosen_unpref_boot;
FEF.p_chosen_unpref_boot=p_chosen_unpref_boot;
FEF.r_unchosen_unpref_boot=r_unchosen_unpref_boot;
FEF.p_unchosen_unpref_boot=p_unchosen_unpref_boot;

clear *boot* *diss*

%PFC
load_name=sprintf('glm_explained_variance_split_boot_pref_%s','PFC');
load(fullfile(load_path,load_name),'r_chosen_pref_boot','p_chosen_pref_boot',...
    'r_unchosen_pref_boot','p_unchosen_ref_boot',...
    'r_chosen_unpref_boot','p_chosen_unpref_boot',...
    'r_unchosen_unpref_boot','p_unchosen_unpref_boot');

PFC.r_chosen_pref_boot=r_chosen_pref_boot;
PFC.p_chosen_pref_boot=p_chosen_pref_boot;
PFC.r_unchosen_pref_boot=r_unchosen_pref_boot;
PFC.p_unchosen_pref_boot=p_unchosen_ref_boot;

PFC.r_chosen_unpref_boot=r_chosen_unpref_boot;
PFC.p_chosen_unpref_boot=p_chosen_unpref_boot;
PFC.r_unchosen_unpref_boot=r_unchosen_unpref_boot;
PFC.p_unchosen_unpref_boot=p_unchosen_unpref_boot;

clear *boot* *diss*


%%

for t=1:length(window_start_list)
    
    LIP.p_ch_pref_unpref(t)=z_test_function_bootstrap(mean(LIP.r_chosen_pref_boot(t,:,:),2)-mean(LIP.r_chosen_unpref_boot(t,:,:),2),0);
    FEF.p_ch_pref_unpref(t)=z_test_function_bootstrap(mean(FEF.r_chosen_pref_boot(t,:,:),2)-mean(FEF.r_chosen_unpref_boot(t,:,:),2),0);
    PFC.p_ch_pref_unpref(t)=z_test_function_bootstrap(mean(PFC.r_chosen_pref_boot(t,:,:),2)-mean(PFC.r_chosen_unpref_boot(t,:,:),2),0);

    LIP.p_unch_pref_unpref(t)=z_test_function_bootstrap(mean(LIP.r_unchosen_pref_boot(t,:,:),2)-mean(LIP.r_unchosen_unpref_boot(t,:,:),2),0);
    FEF.p_unch_pref_unpref(t)=z_test_function_bootstrap(mean(FEF.r_unchosen_pref_boot(t,:,:),2)-mean(FEF.r_unchosen_unpref_boot(t,:,:),2),0);
    PFC.p_unch_pref_unpref(t)=z_test_function_bootstrap(mean(PFC.r_unchosen_pref_boot(t,:,:),2)-mean(PFC.r_unchosen_unpref_boot(t,:,:),2),0);
    
    p_LIP_FEF_unch(t)=z_test_function_bootstrap_2_sample(mean(LIP.r_unchosen_pref_boot(t,:,:),2),mean(FEF.r_unchosen_pref_boot(t,:,:),2));
    p_LIP_PFC_unch(t)=z_test_function_bootstrap_2_sample(mean(LIP.r_unchosen_pref_boot(t,:,:),2),mean(PFC.r_unchosen_pref_boot(t,:,:),2));
    p_FEF_PFC_unch(t)=z_test_function_bootstrap_2_sample(mean(FEF.r_unchosen_pref_boot(t,:,:),2),mean(PFC.r_unchosen_pref_boot(t,:,:),2));
end

%%
for nb=1:size(LIP.r_unchosen_pref_boot,3)
    LIP.time_to_peak_unch(nb)=window_start_list(find(mean(LIP.r_unchosen_pref_boot(6:end,:,nb),2)==max(mean(LIP.r_unchosen_pref_boot(6:end,:,nb),2)),1,'first'))+5*50;
    FEF.time_to_peak_unch(nb)=window_start_list(find(mean(FEF.r_unchosen_pref_boot(6:end,:,nb),2)==max(mean(FEF.r_unchosen_pref_boot(6:end,:,nb),2)),1,'first'))+5*50;
    PFC.time_to_peak_unch(nb)=window_start_list(find(mean(PFC.r_unchosen_pref_boot(6:end,:,nb),2)==max(mean(PFC.r_unchosen_pref_boot(6:end,:,nb),2)),1,'first'))+5*50;
end

%%

p_LIP_FEF=z_test_function_bootstrap_2_sample(LIP.time_to_peak_unch,FEF.time_to_peak_unch);
p_LIP_PFC=z_test_function_bootstrap_2_sample(LIP.time_to_peak_unch,PFC.time_to_peak_unch);
p_FEF_PFC=z_test_function_bootstrap_2_sample(FEF.time_to_peak_unch,PFC.time_to_peak_unch);


%%

figure
subplot(3,1,1)
n_color=1;
hold on
shadedErrorBar(window_start_list(4:end),mean(mean(LIP.r_chosen_pref_boot(4:end,:,:),2),3),[prctile(mean(LIP.r_chosen_pref_boot(4:end,:,:),2),95,3),prctile(mean(LIP.r_chosen_pref_boot(4:end,:,:),2),5,3)]',{'-','color',color_for_ROI(n_color,:),'LineWidth',2},2)
shadedErrorBar(window_start_list(4:end),mean(mean(LIP.r_unchosen_pref_boot(4:end,:,:),2),3),[prctile(mean(LIP.r_unchosen_pref_boot(4:end,:,:),2),95,3),prctile(mean(LIP.r_unchosen_pref_boot(4:end,:,:),2),5,3)]',{'--','color',color_for_ROI(n_color,:),'LineWidth',2},2)
plot_significance_level(window_start_list,LIP.p_chosen_pref_boot,0.8,color_for_ROI(n_color,:),[0.01,0.05/(length(window_start_list)-3), 0.01/(length(window_start_list)-3)])
plot_significance_level(window_start_list,LIP.p_unchosen_pref_boot,0.75,color_for_ROI(n_color,:),[0.01,0.05/(length(window_start_list)-3), 0.01/(length(window_start_list)-3)])

yl.LabelVerticalAlignment = 'middle';
yl.LabelHorizontalAlignment = 'right';
yl.FontSize=12;

text(-575,0.75,'Unchosen')
text(-575,0.8,'Chosen')
text(-575,0.7,'Corr(ch,unch)')
ylim([-0.2 1])

box off
xlabel('Time to target')
xlim([-600 window_start_list(end)+150])
xticks([-600:200:window_start_list(end)+150])
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

subplot(3,1,2)
n_color=2;
hold on
shadedErrorBar(window_start_list(4:end),mean(mean(FEF.r_chosen_pref_boot(4:end,:,:),2),3),[prctile(mean(FEF.r_chosen_pref_boot(4:end,:,:),2),95,3),prctile(mean(FEF.r_chosen_pref_boot(4:end,:,:),2),5,3)]',{'-','color',color_for_ROI(n_color,:),'LineWidth',2},2)
shadedErrorBar(window_start_list(4:end),mean(mean(FEF.r_unchosen_pref_boot(4:end,:,:),2),3),[prctile(mean(FEF.r_unchosen_pref_boot(4:end,:,:),2),95,3),prctile(mean(FEF.r_unchosen_pref_boot(4:end,:,:),2),5,3)]',{'--','color',color_for_ROI(n_color,:),'LineWidth',2},2)
plot_significance_level(window_start_list,FEF.p_chosen_pref_boot,0.8,color_for_ROI(n_color,:),[0.01,0.05/(length(window_start_list)-3), 0.01/(length(window_start_list)-3)])
plot_significance_level(window_start_list,FEF.p_unchosen_pref_boot,0.75,color_for_ROI(n_color,:),[0.01,0.05/(length(window_start_list)-3), 0.01/(length(window_start_list)-3)])

yl.LabelVerticalAlignment = 'middle';
yl.LabelHorizontalAlignment = 'right';
yl.FontSize=12;
text(-575,0.75,'Unchosen')
text(-575,0.8,'Chosen')
text(-575,0.7,'Corr(ch,unch)')
ylim([-0.2 1])


box off
xlabel('Time to target')
xlim([-600 window_start_list(end)+150])
xticks([-600:200:window_start_list(end)+150])
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

subplot(3,1,3)
n_color=3;
hold on
shadedErrorBar(window_start_list(4:end),mean(mean(PFC.r_chosen_pref_boot(4:end,:,:),2),3),[prctile(mean(PFC.r_chosen_pref_boot(4:end,:,:),2),95,3),prctile(mean(PFC.r_chosen_pref_boot(4:end,:,:),2),5,3)]',{'-','color',color_for_ROI(n_color,:),'LineWidth',2},2)
shadedErrorBar(window_start_list(4:end),mean(mean(PFC.r_unchosen_pref_boot(4:end,:,:),2),3),[prctile(mean(PFC.r_unchosen_pref_boot(4:end,:,:),2),95,3),prctile(mean(PFC.r_unchosen_pref_boot(4:end,:,:),2),5,3)]',{'--','color',color_for_ROI(n_color,:),'LineWidth',2},2)
plot_significance_level(window_start_list,PFC.p_chosen_pref_boot,0.8,color_for_ROI(n_color,:),[0.01,0.05/(length(window_start_list)-3), 0.01/(length(window_start_list)-3)])
plot_significance_level(window_start_list,PFC.p_unchosen_pref_boot,0.75,color_for_ROI(n_color,:),[0.01,0.05/(length(window_start_list)-3), 0.01/(length(window_start_list)-3)])

yl.LabelVerticalAlignment = 'middle';
yl.LabelHorizontalAlignment = 'right';
yl.FontSize=12;
text(-575,0.75,'Unchosen')
text(-575,0.8,'Chosen')
text(-575,0.7,'Corr(ch,unch)')

box off
xlabel('Time to target')
xlim([-600 window_start_list(end)+150])
xticks([-600:200:window_start_list(end)+150])
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
ylim([-0.2 1])


%%

figure
subplot(1,3,2)

hold on
shadedErrorBar(window_start_list(4:end),mean(mean(LIP.r_chosen_pref_boot(4:end,:,:),2),3),[prctile(mean(LIP.r_chosen_pref_boot(4:end,:,:),2),95,3),prctile(mean(LIP.r_chosen_pref_boot(4:end,:,:),2),5,3)]',{'-','color',color_for_ROI(1,:),'LineWidth',2},2)
shadedErrorBar(window_start_list(4:end),mean(mean(FEF.r_chosen_pref_boot(4:end,:,:),2),3),[prctile(mean(FEF.r_chosen_pref_boot(4:end,:,:),2),95,3),prctile(mean(FEF.r_chosen_pref_boot(4:end,:,:),2),5,3)]',{'-','color',color_for_ROI(2,:),'LineWidth',2},2)
shadedErrorBar(window_start_list(4:end),mean(mean(PFC.r_chosen_pref_boot(4:end,:,:),2),3),[prctile(mean(PFC.r_chosen_pref_boot(4:end,:,:),2),95,3),prctile(mean(PFC.r_chosen_pref_boot(4:end,:,:),2),5,3)]',{'-','color',color_for_ROI(3,:),'LineWidth',2},2)

plot_significance_level(window_start_list,LIP.p_chosen_pref_boot,0.8,color_for_ROI(1,:),[0.01,0.05/(length(window_start_list)-3), 0.01/(length(window_start_list)-3)])
plot_significance_level(window_start_list,FEF.p_chosen_pref_boot,0.85,color_for_ROI(2,:),[0.01,0.05/(length(window_start_list)-3), 0.01/(length(window_start_list)-3)])
plot_significance_level(window_start_list,PFC.p_chosen_pref_boot,0.9,color_for_ROI(3,:),[0.01,0.05/(length(window_start_list)-3), 0.01/(length(window_start_list)-3)])


text(-575,0.8,'LIP')
text(-575,0.85,'FEF')
text(-575,0.9,'lPFC')
ylim([-0.2 1])

box off
xlabel('Time to target')
xlim([-600 window_start_list(end)+150])
xticks([-600:200:window_start_list(end)+150])
yline(0,'--')
hold on
ax=gca;
ax.XAxis.FontSize=12;
ax.YAxis.FontSize=12;
ylabel('Local chosen value split-half reliability')

subplot(1,3,3)
hold on
shadedErrorBar(window_start_list(4:end),mean(mean(LIP.r_unchosen_pref_boot(4:end,:,:),2),3),[prctile(mean(LIP.r_unchosen_pref_boot(4:end,:,:),2),95,3),prctile(mean(LIP.r_unchosen_pref_boot(4:end,:,:),2),5,3)]',{'-','color',color_for_ROI(1,:),'LineWidth',2},2)
shadedErrorBar(window_start_list(4:end),mean(mean(FEF.r_unchosen_pref_boot(4:end,:,:),2),3),[prctile(mean(FEF.r_unchosen_pref_boot(4:end,:,:),2),95,3),prctile(mean(FEF.r_unchosen_pref_boot(4:end,:,:),2),5,3)]',{'-','color',color_for_ROI(2,:),'LineWidth',2},2)
shadedErrorBar(window_start_list(4:end),mean(mean(PFC.r_unchosen_pref_boot(4:end,:,:),2),3),[prctile(mean(PFC.r_unchosen_pref_boot(4:end,:,:),2),95,3),prctile(mean(PFC.r_unchosen_pref_boot(4:end,:,:),2),5,3)]',{'-','color',color_for_ROI(3,:),'LineWidth',2},2)

plot_significance_level(window_start_list,LIP.p_unchosen_pref_boot,0.8,color_for_ROI(1,:),[0.01,0.05/(length(window_start_list)-3), 0.01/(length(window_start_list)-3)])
plot_significance_level(window_start_list,FEF.p_unchosen_pref_boot,0.85,color_for_ROI(2,:),[0.01,0.05/(length(window_start_list)-3), 0.01/(length(window_start_list)-3)])
plot_significance_level(window_start_list,PFC.p_unchosen_pref_boot,0.9,color_for_ROI(3,:),[0.01,0.05/(length(window_start_list)-3), 0.01/(length(window_start_list)-3)])
plot_significance_level(window_start_list,p_FEF_PFC_unch,0.75,(color_for_ROI(2,:)+color_for_ROI(3,:))./2,[0.01,0.05/(length(window_start_list)-3), 0.01/(length(window_start_list)-3)])

text(-575,0.8,'LIP')
text(-575,0.85,'FEF')
text(-575,0.9,'lPFC')
text(-575,0.75,'FEF ≠ lPFC')
ylabel('Local unchosen value split-half reliability')

ylim([-0.2 1])


box off
xlabel('Time to target')
xlim([-600 window_start_list(end)+150])
xticks([-600:200:window_start_list(end)+150])
yline(0,'--')
hold on
ax=gca;
ax.XAxis.FontSize=12;
ax.YAxis.FontSize=12;

%%

figure
vs = violinplot([LIP.time_to_peak_unch' FEF.time_to_peak_unch' PFC.time_to_peak_unch'],{'LIP','FEF','lPFC'},'ViolinColor',color_for_ROI,'ShowData',false,'ShowMean',true,'ShowMedian',true);
box off
xlabel('ROI')
ylabel('Unchosen value time to peak')
ax=gca;
ax.XAxis.FontSize=12;
ax.YAxis.FontSize=12;


%%
figure
subplot(3,2,1)
n_color=1;
hold on
shadedErrorBar(window_start_list(4:end),mean(mean(LIP.r_chosen_pref_boot(4:end,:,:),2),3),[prctile(mean(LIP.r_chosen_pref_boot(4:end,:,:),2),95,3),prctile(mean(LIP.r_chosen_pref_boot(4:end,:,:),2),5,3)]',{'-','color',color_for_ROI(n_color,:),'LineWidth',2},2)
shadedErrorBar(window_start_list(4:end),mean(mean(LIP.r_chosen_unpref_boot(4:end,:,:),2),3),[prctile(mean(LIP.r_chosen_unpref_boot(4:end,:,:),2),95,3),prctile(mean(LIP.r_chosen_unpref_boot(4:end,:,:),2),5,3)]',{'--','color',color_for_ROI(n_color,:),'LineWidth',2},2)
plot_significance_level(window_start_list,LIP.p_chosen_pref_boot,0.8,color_for_ROI(n_color,:),[0.01,0.05/(length(window_start_list)-3), 0.01/(length(window_start_list)-3)])
plot_significance_level(window_start_list,LIP.p_chosen_unpref_boot,0.75,color_for_ROI(n_color,:),[0.01,0.05/(length(window_start_list)-3), 0.01/(length(window_start_list)-3)])

yl.LabelVerticalAlignment = 'middle';
yl.LabelHorizontalAlignment = 'right';
yl.FontSize=12;

text(-575,0.85,'pref ≠ unpref')
text(-575,0.8,'pref')
text(-575,0.75,'unpref') 
ylim([-0.2 1])

box off
xlabel('Time to target')
xlim([-600 window_start_list(end)+150])
xticks([-600:200:window_start_list(end)+150])
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

subplot(3,2,3)
n_color=2;
hold on
shadedErrorBar(window_start_list(4:end),mean(mean(FEF.r_chosen_pref_boot(4:end,:,:),2),3),[prctile(mean(FEF.r_chosen_pref_boot(4:end,:,:),2),95,3),prctile(mean(FEF.r_chosen_pref_boot(4:end,:,:),2),5,3)]',{'-','color',color_for_ROI(n_color,:),'LineWidth',2},2)
shadedErrorBar(window_start_list(4:end),mean(mean(FEF.r_chosen_unpref_boot(4:end,:,:),2),3),[prctile(mean(FEF.r_chosen_unpref_boot(4:end,:,:),2),95,3),prctile(mean(FEF.r_chosen_unpref_boot(4:end,:,:),2),5,3)]',{'--','color',color_for_ROI(n_color,:),'LineWidth',2},2)
plot_significance_level(window_start_list,FEF.p_chosen_pref_boot,0.8,color_for_ROI(n_color,:),[0.01,0.05/(length(window_start_list)-3), 0.01/(length(window_start_list)-3)])
plot_significance_level(window_start_list,FEF.p_chosen_unpref_boot,0.75,color_for_ROI(n_color,:),[0.01,0.05/(length(window_start_list)-3), 0.01/(length(window_start_list)-3)])

yl.LabelVerticalAlignment = 'middle';
yl.LabelHorizontalAlignment = 'right';
yl.FontSize=12;
text(-575,0.85,'pref ≠ unpref')
text(-575,0.8,'pref')
text(-575,0.75,'unpref') 
ylim([-0.2 1])
box off
xlabel('Time to target')
xlim([-600 window_start_list(end)+150])
xticks([-600:200:window_start_list(end)+150])
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

subplot(3,2,5)
n_color=3;
hold on
shadedErrorBar(window_start_list(4:end),mean(mean(PFC.r_chosen_pref_boot(4:end,:,:),2),3),[prctile(mean(PFC.r_chosen_pref_boot(4:end,:,:),2),95,3),prctile(mean(PFC.r_chosen_pref_boot(4:end,:,:),2),5,3)]',{'-','color',color_for_ROI(n_color,:),'LineWidth',2},2)
shadedErrorBar(window_start_list(4:end),mean(mean(PFC.r_chosen_unpref_boot(4:end,:,:),2),3),[prctile(mean(PFC.r_chosen_unpref_boot(4:end,:,:),2),95,3),prctile(mean(PFC.r_chosen_unpref_boot(4:end,:,:),2),5,3)]',{'--','color',color_for_ROI(n_color,:),'LineWidth',2},2)
plot_significance_level(window_start_list,PFC.p_chosen_pref_boot,0.8,color_for_ROI(n_color,:),[0.01,0.05/(length(window_start_list)-3), 0.01/(length(window_start_list)-3)])
plot_significance_level(window_start_list,PFC.p_chosen_unpref_boot,0.75,color_for_ROI(n_color,:),[0.01,0.05/(length(window_start_list)-3), 0.01/(length(window_start_list)-3)])

yl.LabelVerticalAlignment = 'middle';
yl.LabelHorizontalAlignment = 'right';
yl.FontSize=12;
text(-575,0.85,'pref ≠ unpref')

text(-575,0.8,'pref')
text(-575,0.75,'unpref') 

box off
xlabel('Time to target')
xlim([-600 window_start_list(end)+150])
xticks([-600:200:window_start_list(end)+150])
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
ylim([-0.2 1])

subplot(3,2,2)
n_color=1;
hold on
shadedErrorBar(window_start_list(4:end),mean(mean(LIP.r_unchosen_pref_boot(4:end,:,:),2),3),[prctile(mean(LIP.r_unchosen_pref_boot(4:end,:,:),2),95,3),prctile(mean(LIP.r_unchosen_pref_boot(4:end,:,:),2),5,3)]',{'-','color',color_for_ROI(n_color,:),'LineWidth',2},2)
shadedErrorBar(window_start_list(4:end),mean(mean(LIP.r_unchosen_unpref_boot(4:end,:,:),2),3),[prctile(mean(LIP.r_unchosen_unpref_boot(4:end,:,:),2),95,3),prctile(mean(LIP.r_unchosen_unpref_boot(4:end,:,:),2),5,3)]',{'--','color',color_for_ROI(n_color,:),'LineWidth',2},2)
plot_significance_level(window_start_list,LIP.p_unchosen_pref_boot,0.8,color_for_ROI(n_color,:),[0.01,0.05/(length(window_start_list)-3), 0.01/(length(window_start_list)-3)])
plot_significance_level(window_start_list,LIP.p_unchosen_unpref_boot,0.75,color_for_ROI(n_color,:),[0.01,0.05/(length(window_start_list)-3), 0.01/(length(window_start_list)-3)])

yl.LabelVerticalAlignment = 'middle';
yl.LabelHorizontalAlignment = 'right';
yl.FontSize=12;
text(-575,0.85,'pref ≠ unpref')

text(-575,0.8,'pref')
text(-575,0.75,'unpref') 
ylim([-0.2 1])

box off
xlabel('Time to target')
xlim([-600 window_start_list(end)+150])
xticks([-600:200:window_start_list(end)+150])
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

subplot(3,2,4)
n_color=2;
hold on
shadedErrorBar(window_start_list(4:end),mean(mean(FEF.r_unchosen_pref_boot(4:end,:,:),2),3),[prctile(mean(FEF.r_unchosen_pref_boot(4:end,:,:),2),95,3),prctile(mean(FEF.r_unchosen_pref_boot(4:end,:,:),2),5,3)]',{'-','color',color_for_ROI(n_color,:),'LineWidth',2},2)
shadedErrorBar(window_start_list(4:end),mean(mean(FEF.r_unchosen_unpref_boot(4:end,:,:),2),3),[prctile(mean(FEF.r_unchosen_unpref_boot(4:end,:,:),2),95,3),prctile(mean(FEF.r_unchosen_unpref_boot(4:end,:,:),2),5,3)]',{'--','color',color_for_ROI(n_color,:),'LineWidth',2},2)
plot_significance_level(window_start_list,FEF.p_unchosen_pref_boot,0.8,color_for_ROI(n_color,:),[0.01,0.05/(length(window_start_list)-3), 0.01/(length(window_start_list)-3)])
plot_significance_level(window_start_list,FEF.p_unchosen_unpref_boot,0.75,color_for_ROI(n_color,:),[0.01,0.05/(length(window_start_list)-3), 0.01/(length(window_start_list)-3)])

yl.LabelVerticalAlignment = 'middle';
yl.LabelHorizontalAlignment = 'right';
yl.FontSize=12;
text(-575,0.85,'pref ≠ unpref')

text(-575,0.8,'pref')
text(-575,0.75,'unpref') 
ylim([-0.2 1])


box off
xlabel('Time to target')
xlim([-600 window_start_list(end)+150])
xticks([-600:200:window_start_list(end)+150])
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

subplot(3,2,6)
n_color=3;
hold on
shadedErrorBar(window_start_list(4:end),mean(mean(PFC.r_unchosen_pref_boot(4:end,:,:),2),3),[prctile(mean(PFC.r_unchosen_pref_boot(4:end,:,:),2),95,3),prctile(mean(PFC.r_unchosen_pref_boot(4:end,:,:),2),5,3)]',{'-','color',color_for_ROI(n_color,:),'LineWidth',2},2)
shadedErrorBar(window_start_list(4:end),mean(mean(PFC.r_unchosen_unpref_boot(4:end,:,:),2),3),[prctile(mean(PFC.r_unchosen_unpref_boot(4:end,:,:),2),95,3),prctile(mean(PFC.r_unchosen_unpref_boot(4:end,:,:),2),5,3)]',{'--','color',color_for_ROI(n_color,:),'LineWidth',2},2)
plot_significance_level(window_start_list,PFC.p_unchosen_pref_boot,0.8,color_for_ROI(n_color,:),[0.01,0.05/(length(window_start_list)-3), 0.01/(length(window_start_list)-3)])
plot_significance_level(window_start_list,PFC.p_unchosen_unpref_boot,0.75,color_for_ROI(n_color,:),[0.01,0.05/(length(window_start_list)-3), 0.01/(length(window_start_list)-3)])

yl.LabelVerticalAlignment = 'middle';
yl.LabelHorizontalAlignment = 'right';
yl.FontSize=12;
text(-575,0.85,'pref ≠ unpref')

text(-575,0.8,'pref')
text(-575,0.75,'unpref') 

box off
xlabel('Time to target')
xlim([-600 window_start_list(end)+150])
xticks([-600:200:window_start_list(end)+150])
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
ylim([-0.2 1])

ax.YAxis.FontSize=12;

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

function p = z_test_function_bootstrap(dist,null)

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

function p = z_test_function_bootstrap_2_sample(dist1,dist2)

if sum(isnan(dist1))==0 && sum(isnan(dist2))==0
    
    m1 = mean(dist1);
    m2 = mean(dist2);
    
    s1 = std(dist1);
    s2 = std(dist2);
    
    z = (m1-m2)/sqrt(s1^2 + s2^2);
    p=2*(1-normcdf(abs(z)));
    
else
    p=NaN;
end

end


