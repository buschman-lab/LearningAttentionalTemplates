
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

%%

color_for_ROI_2(1,:,1)=color_for_ROI(1,:)+0.25;
color_for_ROI_2(1,:,2)=color_for_ROI(1,:);

color_for_ROI_2(2,:,1)=color_for_ROI(2,:)+0.15;
color_for_ROI_2(2,:,2)=color_for_ROI(2,:)-0.09;

color_for_ROI_2(3,:,1)=color_for_ROI(3,:)+0.07;
color_for_ROI_2(3,:,2)=color_for_ROI(3,:)-0.12;



%% load

%LIP

load_name=sprintf('glm_explained_variance_split_boot_%s','LIP');
load(fullfile(load_path,load_name),'*_chosen*')

LIP.r_chosen_contra_boot=r_chosen_contra_boot;
LIP.r_chosen_contra_contra_boot=r_chosen_contra_contra_boot;
LIP.r_chosen_ipsi_boot=r_chosen_ipsi_boot;
LIP.r_chosen_contra_ipsi_boot=r_chosen_contra_ipsi_boot;
LIP.r_chosen_contra_contra_diag_boot=r_chosen_contra_contra_diag_boot;
LIP.p_chosen_contra_boot=p_chosen_contra_boot;

for t=1:length(window_start_list)
    for n=1:size(r_chosen_contra_boot,3)
        if mean(r_chosen_contra_boot(t,:,n),2)>0 %&& mean(r_chosen_ipsi_boot(t,:,n),2)>0 
            LIP.sim_r_contra_ipsi_chosen(t,n)=mean(r_chosen_contra_ipsi_boot(t,:,n),2)./mean(r_chosen_contra_boot(t,:,n),2);
            LIP.diss_r_contra_ipsi_chosen(t,n)=mean(r_chosen_contra_boot(t,:,n),2)-mean(r_chosen_contra_ipsi_boot(t,:,n),2);
        else
            LIP.sim_r_contra_ipsi_chosen(t,n)=NaN;
            LIP.diss_r_contra_ipsi_chosen(t,n)=NaN;
        end
        if mean(r_chosen_contra_boot(t,:,n),2)>0 %&& mean(r_chosen_ipsi_boot(t,:,n),2)>0
            LIP.sim_r_contra_contra_diag_chosen(t,n)=mean(r_chosen_contra_contra_diag_boot(t,:,n),2)./mean(r_chosen_contra_boot(t,:,n),2);
            LIP.diss_r_contra_contra_diag_chosen(t,n)=mean(r_chosen_contra_boot(t,:,n),2)-mean(r_chosen_contra_contra_diag_boot(t,:,n),2);
        else
            LIP.sim_r_contra_contra_diag_chosen(t,n)=NaN;
            LIP.diss_r_contra_contra_diag_chosen(t,n)=NaN;
        end
        if mean(r_chosen_contra_boot(t,:,n),2)>0  
            LIP.sim_r_contra_contra_chosen(t,n)=mean(r_chosen_contra_contra_boot(t,:,n),2)./mean(r_chosen_contra_boot(t,:,n),2);
            LIP.diss_r_contra_contra_chosen(t,n)=mean(r_chosen_contra_boot(t,:,n),2)-mean(r_chosen_contra_contra_boot(t,:,n),2);
        else
            LIP.sim_r_contra_contra_chosen(t,n)=NaN;
            LIP.diss_r_contra_contra_chosen(t,n)=NaN;
        end
    end
    LIP.p_chosen_contra_ipsi_diss(t) = z_test_function_bootstrap(LIP.diss_r_contra_ipsi_chosen(t,:),0);
    LIP.p_chosen_contra_ipsi_sim(t) = z_test_function_bootstrap(LIP.sim_r_contra_ipsi_chosen(t,:),0);
    LIP.p_chosen_contra_contra_diag_diss(t) = z_test_function_bootstrap(LIP.diss_r_contra_contra_diag_chosen(t,:),0);
    LIP.p_chosen_contra_contra_diag_sim(t) = z_test_function_bootstrap(LIP.sim_r_contra_contra_diag_chosen(t,:),0);
    LIP.p_chosen_contra_contra_diss(t) = z_test_function_bootstrap(LIP.diss_r_contra_contra_chosen(t,:),0);
    LIP.p_chosen_contra_contra_sim(t) = z_test_function_bootstrap(LIP.sim_r_contra_contra_chosen(t,:),0);
end

clear *_chosen*

%FEF
load_name=sprintf('glm_explained_variance_split_boot_%s','FEF');
load(fullfile(load_path,load_name),'*_chosen*')

FEF.r_chosen_contra_boot=r_chosen_contra_boot;
FEF.r_chosen_contra_contra_boot=r_chosen_contra_contra_boot;
FEF.r_chosen_ipsi_boot=r_chosen_ipsi_boot;
FEF.r_chosen_contra_ipsi_boot=r_chosen_contra_ipsi_boot;
FEF.r_chosen_contra_contra_diag_boot=r_chosen_contra_contra_diag_boot;
FEF.p_chosen_contra_boot=p_chosen_contra_boot;

for t=1:length(window_start_list)
    for n=1:size(r_chosen_contra_boot,3)
        if mean(r_chosen_contra_boot(t,:,n),2)>0 %&& mean(r_chosen_ipsi_boot(t,:,n),2)>0 
            FEF.sim_r_contra_ipsi_chosen(t,n)=mean(r_chosen_contra_ipsi_boot(t,:,n),2)./mean(r_chosen_contra_boot(t,:,n),2);
            FEF.diss_r_contra_ipsi_chosen(t,n)=mean(r_chosen_contra_boot(t,:,n),2)-mean(r_chosen_contra_ipsi_boot(t,:,n),2);
        else
            FEF.sim_r_contra_ipsi_chosen(t,n)=NaN;
            FEF.diss_r_contra_ipsi_chosen(t,n)=NaN;
        end
        if mean(r_chosen_contra_boot(t,:,n),2)>0 %&& mean(r_chosen_ipsi_boot(t,:,n),2)>0
            FEF.sim_r_contra_contra_diag_chosen(t,n)=mean(r_chosen_contra_contra_diag_boot(t,:,n),2)./mean(r_chosen_contra_boot(t,:,n),2);
            FEF.diss_r_contra_contra_diag_chosen(t,n)=mean(r_chosen_contra_boot(t,:,n),2)-mean(r_chosen_contra_contra_diag_boot(t,:,n),2);
        else
            FEF.sim_r_contra_contra_diag_chosen(t,n)=NaN;
            FEF.diss_r_contra_contra_diag_chosen(t,n)=NaN;
        end
        if mean(r_chosen_contra_boot(t,:,n),2)>0  
            FEF.sim_r_contra_contra_chosen(t,n)=mean(r_chosen_contra_contra_boot(t,:,n),2)./mean(r_chosen_contra_boot(t,:,n),2);
            FEF.diss_r_contra_contra_chosen(t,n)=mean(r_chosen_contra_boot(t,:,n),2)-mean(r_chosen_contra_contra_boot(t,:,n),2);
        else
            FEF.sim_r_contra_contra_chosen(t,n)=NaN;
            FEF.diss_r_contra_contra_chosen(t,n)=NaN;
        end
    end
    FEF.p_chosen_contra_ipsi_diss(t) = z_test_function_bootstrap(FEF.diss_r_contra_ipsi_chosen(t,:),0);
    FEF.p_chosen_contra_ipsi_sim(t) = z_test_function_bootstrap(FEF.sim_r_contra_ipsi_chosen(t,:),0);
    FEF.p_chosen_contra_contra_diag_diss(t) = z_test_function_bootstrap(FEF.diss_r_contra_contra_diag_chosen(t,:),0);
    FEF.p_chosen_contra_contra_diag_sim(t) = z_test_function_bootstrap(FEF.sim_r_contra_contra_diag_chosen(t,:),0);
    FEF.p_chosen_contra_contra_diss(t) = z_test_function_bootstrap(FEF.diss_r_contra_contra_chosen(t,:),0);
    FEF.p_chosen_contra_contra_sim(t) = z_test_function_bootstrap(FEF.sim_r_contra_contra_chosen(t,:),0);
end

clear *_chosen*

%PFC
load_name=sprintf('glm_explained_variance_split_boot_%s','PFC');
load(fullfile(load_path,load_name),'*_chosen*')

PFC.r_chosen_contra_boot=r_chosen_contra_boot;
PFC.r_chosen_contra_contra_boot=r_chosen_contra_contra_boot;
PFC.r_chosen_ipsi_boot=r_chosen_ipsi_boot;
PFC.r_chosen_contra_ipsi_boot=r_chosen_contra_ipsi_boot;
PFC.r_chosen_contra_contra_diag_boot=r_chosen_contra_contra_diag_boot;
PFC.p_chosen_contra_boot=p_chosen_contra_boot;

for t=1:length(window_start_list)
    for n=1:size(r_chosen_contra_boot,3)
        if mean(r_chosen_contra_boot(t,:,n),2)>0 %&& mean(r_chosen_ipsi_boot(t,:,n),2)>0 
            PFC.sim_r_contra_ipsi_chosen(t,n)=mean(r_chosen_contra_ipsi_boot(t,:,n),2)./mean(r_chosen_contra_boot(t,:,n),2);
            PFC.diss_r_contra_ipsi_chosen(t,n)=mean(r_chosen_contra_boot(t,:,n),2)-mean(r_chosen_contra_ipsi_boot(t,:,n),2);
        else
            PFC.sim_r_contra_ipsi_chosen(t,n)=NaN;
            PFC.diss_r_contra_ipsi_chosen(t,n)=NaN;
        end
        if mean(r_chosen_contra_boot(t,:,n),2)>0 %&& mean(r_chosen_ipsi_boot(t,:,n),2)>0
            PFC.sim_r_contra_contra_diag_chosen(t,n)=mean(r_chosen_contra_contra_diag_boot(t,:,n),2)./mean(r_chosen_contra_boot(t,:,n),2);
            PFC.diss_r_contra_contra_diag_chosen(t,n)=mean(r_chosen_contra_boot(t,:,n),2)-mean(r_chosen_contra_contra_diag_boot(t,:,n),2);
        else
            PFC.sim_r_contra_contra_diag_chosen(t,n)=NaN;
            PFC.diss_r_contra_contra_diag_chosen(t,n)=NaN;
        end
        if mean(r_chosen_contra_boot(t,:,n),2)>0  
            PFC.sim_r_contra_contra_chosen(t,n)=mean(r_chosen_contra_contra_boot(t,:,n),2)./mean(r_chosen_contra_boot(t,:,n),2);
            PFC.diss_r_contra_contra_chosen(t,n)=mean(r_chosen_contra_boot(t,:,n),2)-mean(r_chosen_contra_contra_boot(t,:,n),2);
        else
            PFC.sim_r_contra_contra_chosen(t,n)=NaN;
            PFC.diss_r_contra_contra_chosen(t,n)=NaN;
        end
    end
    PFC.p_chosen_contra_ipsi_diss(t) = z_test_function_bootstrap(PFC.diss_r_contra_ipsi_chosen(t,:),0);
    PFC.p_chosen_contra_ipsi_sim(t) = z_test_function_bootstrap(PFC.sim_r_contra_ipsi_chosen(t,:),0);
    PFC.p_chosen_contra_contra_diag_diss(t) = z_test_function_bootstrap(PFC.diss_r_contra_contra_diag_chosen(t,:),0);
    PFC.p_chosen_contra_contra_diag_sim(t) = z_test_function_bootstrap(PFC.sim_r_contra_contra_diag_chosen(t,:),0);
    PFC.p_chosen_contra_contra_diss(t) = z_test_function_bootstrap(PFC.diss_r_contra_contra_chosen(t,:),0);
    PFC.p_chosen_contra_contra_sim(t) = z_test_function_bootstrap(PFC.sim_r_contra_contra_chosen(t,:),0);
end



%%
for nb=1:size(LIP.r_chosen_contra_boot,3)
  LIP.auc_r_chosen(1,nb)=sum(mean(LIP.r_chosen_contra_boot(4:end,:,nb),2));
  LIP.auc_r_chosen(2,nb)=sum(mean(LIP.r_chosen_contra_contra_boot(4:end,:,nb),2));
  LIP.auc_r_chosen(3,nb)=sum(mean(LIP.r_chosen_contra_ipsi_boot(4:end,:,nb),2));
  LIP.auc_r_chosen(4,nb)=sum(mean(LIP.r_chosen_contra_contra_diag_boot(4:end,:,nb),2));
  
  FEF.auc_r_chosen(1,nb)=sum(mean(FEF.r_chosen_contra_boot(4:end,:,nb),2));
  FEF.auc_r_chosen(2,nb)=sum(mean(FEF.r_chosen_contra_contra_boot(4:end,:,nb),2));
  FEF.auc_r_chosen(3,nb)=sum(mean(FEF.r_chosen_contra_ipsi_boot(4:end,:,nb),2));
  FEF.auc_r_chosen(4,nb)=sum(mean(FEF.r_chosen_contra_contra_diag_boot(4:end,:,nb),2));
  
  PFC.auc_r_chosen(1,nb)=sum(mean(PFC.r_chosen_contra_boot(4:end,:,nb),2));
  PFC.auc_r_chosen(2,nb)=sum(mean(PFC.r_chosen_contra_contra_boot(4:end,:,nb),2));
  PFC.auc_r_chosen(3,nb)=sum(mean(PFC.r_chosen_contra_ipsi_boot(4:end,:,nb),2));
  PFC.auc_r_chosen(4,nb)=sum(mean(PFC.r_chosen_contra_contra_diag_boot(4:end,:,nb),2));
  
end

%%
for nb=1:size(LIP.r_chosen_contra_boot,3)
  LIP.auc_early_r_chosen(1,nb)=sum(mean(LIP.r_chosen_contra_boot(4:15,:,nb),2));
  LIP.auc_early_r_chosen(2,nb)=sum(mean(LIP.r_chosen_contra_contra_boot(4:15,:,nb),2));
  LIP.auc_early_r_chosen(3,nb)=sum(mean(LIP.r_chosen_contra_ipsi_boot(4:15,:,nb),2));
  LIP.auc_early_r_chosen(4,nb)=sum(mean(LIP.r_chosen_contra_contra_diag_boot(4:15,:,nb),2));
  
  FEF.auc_early_r_chosen(1,nb)=sum(mean(FEF.r_chosen_contra_boot(4:15,:,nb),2));
  FEF.auc_early_r_chosen(2,nb)=sum(mean(FEF.r_chosen_contra_contra_boot(4:15,:,nb),2));
  FEF.auc_early_r_chosen(3,nb)=sum(mean(FEF.r_chosen_contra_ipsi_boot(4:15,:,nb),2));
  FEF.auc_early_r_chosen(4,nb)=sum(mean(FEF.r_chosen_contra_contra_diag_boot(4:15,:,nb),2));
  
  PFC.auc_early_r_chosen(1,nb)=sum(mean(PFC.r_chosen_contra_boot(4:15,:,nb),2));
  PFC.auc_early_r_chosen(2,nb)=sum(mean(PFC.r_chosen_contra_contra_boot(4:15,:,nb),2));
  PFC.auc_early_r_chosen(3,nb)=sum(mean(PFC.r_chosen_contra_ipsi_boot(4:15,:,nb),2));
  PFC.auc_early_r_chosen(4,nb)=sum(mean(PFC.r_chosen_contra_contra_diag_boot(4:15,:,nb),2));
  
end

%%
p_LIP_auc_across_loc=NaN(4,4);
p_FEF_auc_across_loc=NaN(4,4);
p_PFC_auc_across_loc=NaN(4,4);

for i=1:3
    for j=i+1:4
        p_LIP_auc_across_loc(i,j)=z_test_function_bootstrap(LIP.auc_r_chosen(i,:)-LIP.auc_r_chosen(j,:),0);
        p_FEF_auc_across_loc(i,j)=z_test_function_bootstrap(FEF.auc_r_chosen(i,:)-FEF.auc_r_chosen(j,:),0);
        p_PFC_auc_across_loc(i,j)=z_test_function_bootstrap(PFC.auc_r_chosen(i,:)-PFC.auc_r_chosen(j,:),0);
    end
end

%%

p_LIP_auc_early_across_loc=NaN(4,4);
p_FEF_auc_early_across_loc=NaN(4,4);
p_PFC_auc_early_across_loc=NaN(4,4);


for i=1:3
    for j=i+1:4
        p_LIP_auc_early_across_loc(i,j)=z_test_function_bootstrap(LIP.auc_early_r_chosen(i,:)-LIP.auc_early_r_chosen(j,:),0);
        p_FEF_auc_early_across_loc(i,j)=z_test_function_bootstrap(FEF.auc_early_r_chosen(i,:)-FEF.auc_early_r_chosen(j,:),0);
        p_PFC_auc_early_across_loc(i,j)=z_test_function_bootstrap(PFC.auc_early_r_chosen(i,:)-PFC.auc_early_r_chosen(j,:),0);
    end
end

        
%%

figure
subplot(3,1,1)
vs = violinplot(LIP.auc_r_chosen',{'1','2','3','4'},'ViolinColor',color_for_ROI(1,:),'ShowData',false,'ShowMean',true,'ShowMedian',true);
box off
xlabel('Distance')
ylabel('AUC correlation chosen value vectors')
ax=gca;
ax.XAxis.FontSize=12;
ax.YAxis.FontSize=12;
ylim([2 10])
subplot(3,1,2)
vs = violinplot(FEF.auc_r_chosen',{'1','2','3','4'},'ViolinColor',color_for_ROI(2,:),'ShowData',false,'ShowMean',true,'ShowMedian',true);
box off
xlabel('Distance')
ylabel('AUC correlation chosen value vectors')
ax=gca;
ax.XAxis.FontSize=12;
ax.YAxis.FontSize=12;
ylim([2 10])

subplot(3,1,3)
vs = violinplot(PFC.auc_r_chosen',{'1','2','3','4'},'ViolinColor',color_for_ROI(3,:),'ShowData',false,'ShowMean',true,'ShowMedian',true);
box off
xlabel('Distance')
ylabel('AUC correlation chosen value vectors')
ax=gca;
ax.XAxis.FontSize=12;
ax.YAxis.FontSize=12;
ylim([2 10])

%%
figure
subplot(3,1,1)
vs = violinplot(LIP.auc_early_r_chosen',{'1','2','3','4'},'ViolinColor',color_for_ROI(1,:),'ShowData',false,'ShowMean',true,'ShowMedian',true);
box off
xlabel('Distance')
ylabel('AUC correlation chosen value vectors')
ax=gca;
ax.XAxis.FontSize=12;
ax.YAxis.FontSize=12;
ylim([-1 3.5])

subplot(3,1,2)
vs = violinplot(FEF.auc_early_r_chosen',{'1','2','3','4'},'ViolinColor',color_for_ROI(2,:),'ShowData',false,'ShowMean',true,'ShowMedian',true);
box off
xlabel('Distance')
ylabel('AUC correlation chosen value vectors')
ax=gca;
ax.XAxis.FontSize=12;
ax.YAxis.FontSize=12;
ylim([-1 3.5])

subplot(3,1,3)
vs = violinplot(PFC.auc_early_r_chosen',{'1','2','3','4'},'ViolinColor',color_for_ROI(3,:),'ShowData',false,'ShowMean',true,'ShowMedian',true);
box off
xlabel('Distance')
ylabel('AUC correlation chosen value vectors')
ax=gca;
ax.XAxis.FontSize=12;
ax.YAxis.FontSize=12;
ylim([-1 3.5])


%%
figure
subplot(3,1,1)

vs = violinplot(LIP.auc_early_r_chosen'./mean(LIP.auc_early_r_chosen(1,:)),{'1','2','3','4'},'ViolinColor',color_for_ROI_2(1,:,2),'ShowData',false,'ViolinAlpha',0.5,'HalfViolin','left','BoxColor',[0 0 0],'MedianColor',[0.75 0.75 0.75]);
hold on
vs = violinplot(LIP.auc_r_chosen'./mean(LIP.auc_r_chosen(1,:)),{'1','2','3','4'},'ViolinColor',color_for_ROI_2(1,:,1),'ShowData',false,'ViolinAlpha',0.5,'HalfViolin','right','BoxColor',[0.5 0.5 0.5],'MedianColor',[1 1 1]);

box off
xlabel('Distance')
ylabel('auc_early correlation chosen value vectors')
ax=gca;
ax.XAxis.FontSize=12;
ax.YAxis.FontSize=12;
ylim([-1 2])
subplot(3,1,2)
vs = violinplot(FEF.auc_early_r_chosen'./mean(FEF.auc_early_r_chosen(1,:)),{'1','2','3','4'},'ViolinColor',color_for_ROI_2(2,:,2),'ShowData',false,'ViolinAlpha',0.5,'HalfViolin','left','BoxColor',[0 0 0],'MedianColor',[0.75 0.75 0.75]);
hold on
vs = violinplot(FEF.auc_r_chosen'./mean(FEF.auc_r_chosen(1,:)),{'1','2','3','4'},'ViolinColor',color_for_ROI_2(2,:,1),'ShowData',false,'ViolinAlpha',0.5,'HalfViolin','right','BoxColor',[0.5 0.5 0.5],'MedianColor',[1 1 1]);

box off
xlabel('Distance')
ylabel('auc_early correlation chosen value vectors')
ax=gca;
ax.XAxis.FontSize=12;
ax.YAxis.FontSize=12;
ylim([-1 2])

subplot(3,1,3)
vs = violinplot(PFC.auc_early_r_chosen'./mean(PFC.auc_early_r_chosen(1,:)),{'1','2','3','4'},'ViolinColor',color_for_ROI_2(3,:,2),'ShowData',false,'ViolinAlpha',0.5,'HalfViolin','left','BoxColor',[0 0 0],'MedianColor',[0.75 0.75 0.75]);
hold on
vs = violinplot(PFC.auc_r_chosen'./mean(PFC.auc_r_chosen(1,:)),{'1','2','3','4'},'ViolinColor',color_for_ROI_2(3,:,1),'ShowData',false,'ViolinAlpha',0.5,'HalfViolin','right','BoxColor',[0.5 0.5 0.5],'MedianColor',[1 1 1]);

box off
xlabel('Distance')
ylabel('auc_early correlation chosen value vectors')
ax=gca;
ax.XAxis.FontSize=12;
ax.YAxis.FontSize=12;
ylim([-1 2])

%%

figure
hold on
shadedErrorBar([1:4],mean(LIP.auc_r_chosen,2),[prctile(LIP.auc_r_chosen,95,2) prctile(LIP.auc_r_chosen,5,2)]',{'-','color',color_for_ROI(1,:),'LineWidth',2},2)
shadedErrorBar([1:4],mean(FEF.auc_r_chosen,2),[prctile(FEF.auc_r_chosen,95,2) prctile(FEF.auc_r_chosen,5,2)]',{'-','color',color_for_ROI(2,:),'LineWidth',2},2)
shadedErrorBar([1:4],mean(PFC.auc_r_chosen,2),[prctile(PFC.auc_r_chosen,95,2) prctile(PFC.auc_r_chosen,5,2)]',{'-','color',color_for_ROI(3,:),'LineWidth',2},2)
box off

%%
figure
subplot(1,4,1)
hold on
shadedErrorBar(window_start_list(4:end),mean(mean(LIP.r_chosen_contra_boot(4:end,:,:),2),3),[prctile(mean(LIP.r_chosen_contra_boot(4:end,:,:),2),95,3) prctile(mean(LIP.r_chosen_contra_boot(4:end,:,:),2),5,3)]',{'-','color',color_each_loc(1,:,1),'LineWidth',2},2)
shadedErrorBar(window_start_list(4:end),mean(mean(LIP.r_chosen_contra_contra_boot(4:end,:,:),2),3),[prctile(mean(LIP.r_chosen_contra_contra_boot(4:end,:,:),2),95,3) prctile(mean(LIP.r_chosen_contra_contra_boot(4:end,:,:),2),5,3)]',{'-','color',color_each_loc(1,:,2),'LineWidth',2},2)
shadedErrorBar(window_start_list(4:end),mean(mean(LIP.r_chosen_contra_ipsi_boot(4:end,:,:),2),3),[prctile(mean(LIP.r_chosen_contra_ipsi_boot(4:end,:,:),2),95,3) prctile(mean(LIP.r_chosen_contra_ipsi_boot(4:end,:,:),2),5,3)]',{'-','color',color_each_loc(1,:,3),'LineWidth',2},2)
shadedErrorBar(window_start_list(4:end),mean(mean(LIP.r_chosen_contra_contra_diag_boot(4:end,:,:),2),3),[prctile(mean(LIP.r_chosen_contra_contra_diag_boot(4:end,:,:),2),95,3) prctile(mean(LIP.r_chosen_contra_contra_diag_boot(4:end,:,:),2),5,3)]',{'-','color',color_each_loc(1,:,4),'LineWidth',2},2)

plot_significance_level(window_start_list,LIP.p_chosen_contra_boot,0.95,color_each_loc(1,:,1),[0.01,0.05/(length(window_start_list)-3), 0.01/(length(window_start_list)-3)])
plot_significance_level(window_start_list,LIP.p_chosen_contra_contra_sim,0.9,color_each_loc(1,:,2),[0.01,0.05/(length(window_start_list)-3), 0.01/(length(window_start_list)-3)])
plot_significance_level(window_start_list,LIP.p_chosen_contra_ipsi_sim,0.85,color_each_loc(1,:,3),[0.01,0.05/(length(window_start_list)-3), 0.01/(length(window_start_list)-3)])
plot_significance_level(window_start_list,LIP.p_chosen_contra_contra_diag_sim,0.8,color_each_loc(1,:,4),[0.01,0.05/(length(window_start_list)-3), 0.01/(length(window_start_list)-3)])


yl.LabelVerticalAlignment = 'middle';
yl.LabelHorizontalAlignment = 'right';
yl.FontSize=12;
box off
xlabel('Time to target')
xlim([-600 window_start_list(end)+150])
xticks([-600:200:window_start_list(end)+150])
yline(0,'--')
ax=gca;
ax.XAxis.FontSize=12;
ax.YAxis.FontSize=12;
ylabel('Correlation chosen value vectors')
ylim([-0.2 1])

subplot(1,4,2)
hold on
shadedErrorBar(window_start_list(4:end),mean(mean(FEF.r_chosen_contra_boot(4:end,:,:),2),3),[prctile(mean(FEF.r_chosen_contra_boot(4:end,:,:),2),95,3) prctile(mean(FEF.r_chosen_contra_boot(4:end,:,:),2),5,3)]',{'-','color',color_each_loc(2,:,1),'LineWidth',2},2)
shadedErrorBar(window_start_list(4:end),mean(mean(FEF.r_chosen_contra_contra_boot(4:end,:,:),2),3),[prctile(mean(FEF.r_chosen_contra_contra_boot(4:end,:,:),2),95,3) prctile(mean(FEF.r_chosen_contra_contra_boot(4:end,:,:),2),5,3)]',{'-','color',color_each_loc(2,:,2),'LineWidth',2},2)
shadedErrorBar(window_start_list(4:end),mean(mean(FEF.r_chosen_contra_ipsi_boot(4:end,:,:),2),3),[prctile(mean(FEF.r_chosen_contra_ipsi_boot(4:end,:,:),2),95,3) prctile(mean(FEF.r_chosen_contra_ipsi_boot(4:end,:,:),2),5,3)]',{'-','color',color_each_loc(2,:,3),'LineWidth',2},2)
shadedErrorBar(window_start_list(4:end),mean(mean(FEF.r_chosen_contra_contra_diag_boot(4:end,:,:),2),3),[prctile(mean(FEF.r_chosen_contra_contra_diag_boot(4:end,:,:),2),95,3) prctile(mean(FEF.r_chosen_contra_contra_diag_boot(4:end,:,:),2),5,3)]',{'-','color',color_each_loc(2,:,4),'LineWidth',2},2)

plot_significance_level(window_start_list,FEF.p_chosen_contra_boot,0.95,color_each_loc(2,:,1),[0.01,0.05/(length(window_start_list)-3), 0.01/(length(window_start_list)-3)])
plot_significance_level(window_start_list,FEF.p_chosen_contra_contra_sim,0.9,color_each_loc(2,:,2),[0.01,0.05/(length(window_start_list)-3), 0.01/(length(window_start_list)-3)])
plot_significance_level(window_start_list,FEF.p_chosen_contra_ipsi_sim,0.85,color_each_loc(2,:,3),[0.01,0.05/(length(window_start_list)-3), 0.01/(length(window_start_list)-3)])
plot_significance_level(window_start_list,FEF.p_chosen_contra_contra_diag_sim,0.8,color_each_loc(2,:,4),[0.01,0.05/(length(window_start_list)-3), 0.01/(length(window_start_list)-3)])

yl.LabelVerticalAlignment = 'middle';
yl.LabelHorizontalAlignment = 'right';
yl.FontSize=12;
box off
xlabel('Time to target')
xlim([-600 window_start_list(end)+150])
xticks([-600:200:window_start_list(end)+150])
yline(0,'--')
ax=gca;
ax.XAxis.FontSize=12;
ax.YAxis.FontSize=12;
ylabel('Correlation chosen value vectors')
ylim([-0.2 1])

subplot(1,4,3)
hold on
shadedErrorBar(window_start_list(4:end),mean(mean(PFC.r_chosen_contra_boot(4:end,:,:),2),3),[prctile(mean(PFC.r_chosen_contra_boot(4:end,:,:),2),95,3) prctile(mean(PFC.r_chosen_contra_boot(4:end,:,:),2),5,3)]',{'-','color',color_each_loc(3,:,1),'LineWidth',2},2)
shadedErrorBar(window_start_list(4:end),mean(mean(PFC.r_chosen_contra_contra_boot(4:end,:,:),2),3),[prctile(mean(PFC.r_chosen_contra_contra_boot(4:end,:,:),2),95,3) prctile(mean(PFC.r_chosen_contra_contra_boot(4:end,:,:),2),5,3)]',{'-','color',color_each_loc(3,:,2),'LineWidth',2},2)
shadedErrorBar(window_start_list(4:end),mean(mean(PFC.r_chosen_contra_ipsi_boot(4:end,:,:),2),3),[prctile(mean(PFC.r_chosen_contra_ipsi_boot(4:end,:,:),2),95,3) prctile(mean(PFC.r_chosen_contra_ipsi_boot(4:end,:,:),2),5,3)]',{'-','color',color_each_loc(3,:,3),'LineWidth',2},2)
shadedErrorBar(window_start_list(4:end),mean(mean(PFC.r_chosen_contra_contra_diag_boot(4:end,:,:),2),3),[prctile(mean(PFC.r_chosen_contra_contra_diag_boot(4:end,:,:),2),95,3) prctile(mean(PFC.r_chosen_contra_contra_diag_boot(4:end,:,:),2),5,3)]',{'-','color',color_each_loc(3,:,4),'LineWidth',2},2)

plot_significance_level(window_start_list,PFC.p_chosen_contra_boot,0.95,color_each_loc(3,:,1),[0.01,0.05/(length(window_start_list)-3), 0.01/(length(window_start_list)-3)])
plot_significance_level(window_start_list,PFC.p_chosen_contra_contra_sim,0.9,color_each_loc(3,:,2),[0.01,0.05/(length(window_start_list)-3), 0.01/(length(window_start_list)-3)])
plot_significance_level(window_start_list,PFC.p_chosen_contra_ipsi_sim,0.85,color_each_loc(3,:,3),[0.01,0.05/(length(window_start_list)-3), 0.01/(length(window_start_list)-3)])
plot_significance_level(window_start_list,PFC.p_chosen_contra_contra_diag_sim,0.8,color_each_loc(3,:,4),[0.01,0.05/(length(window_start_list)-3), 0.01/(length(window_start_list)-3)])

yl.LabelVerticalAlignment = 'middle';
yl.LabelHorizontalAlignment = 'right';
yl.FontSize=12;
box off
xlabel('Time to target')
xlim([-600 window_start_list(end)+150])
xticks([-600:200:window_start_list(end)+150])
yline(0,'--')
ax=gca;
ax.XAxis.FontSize=12;
ax.YAxis.FontSize=12;
ylabel('Correlation chosen value vectors')
ylim([-0.2 1])

subplot(1,4,4)
hold on
shadedErrorBar([1:4],mean(LIP.auc_r_chosen,2),[prctile(LIP.auc_r_chosen,95,2) prctile(LIP.auc_r_chosen,5,2)]',{'-','color',color_for_ROI(1,:),'LineWidth',2},2)
shadedErrorBar([1:4],mean(FEF.auc_r_chosen,2),[prctile(FEF.auc_r_chosen,95,2) prctile(FEF.auc_r_chosen,5,2)]',{'-','color',color_for_ROI(2,:),'LineWidth',2},2)
shadedErrorBar([1:4],mean(PFC.auc_r_chosen,2),[prctile(PFC.auc_r_chosen,95,2) prctile(PFC.auc_r_chosen,5,2)]',{'-','color',color_for_ROI(3,:),'LineWidth',2},2)
box off
xlabel('Distance')
ylabel('AUC correlation chosen value vectors')
ax=gca;
ax.XAxis.FontSize=12;
ax.YAxis.FontSize=12;



%%
n_color=1;

figure
subplot(3,3,1)
hold on
shadedErrorBar(window_start_list(4:end),mean(mean(LIP.r_chosen_contra_boot(4:end,:,:),2),3),[prctile(mean(LIP.r_chosen_contra_boot(4:end,:,:),2),95,3) prctile(mean(LIP.r_chosen_contra_boot(4:end,:,:),2),5,3)]',{'-','color',[0.5 0.5 0.5],'LineWidth',2},2)
shadedErrorBar(window_start_list(4:end),mean(mean(LIP.r_chosen_contra_contra_boot(4:end,:,:),2),3),[prctile(mean(LIP.r_chosen_contra_contra_boot(4:end,:,:),2),95,3) prctile(mean(LIP.r_chosen_contra_contra_boot(4:end,:,:),2),5,3)]',{'-','color',color_for_ROI(n_color,:),'LineWidth',2},2)
plot_significance_level(window_start_list,LIP.p_chosen_contra_contra_sim,0.9,color_for_ROI(n_color,:),[0.01,0.05/(length(window_start_list)-3), 0.01/(length(window_start_list)-3)])
plot_significance_level(window_start_list,LIP.p_chosen_contra_contra_diss,0.95,color_for_ROI(n_color,:),[0.01,0.05/(length(window_start_list)-3), 0.01/(length(window_start_list)-3)])
plot_significance_level(window_start_list,p_chosen_contra_boot,0.85,[0.5 0.5 0.5],[0.01,0.05/(length(window_start_list)-3), 0.01/(length(window_start_list)-3)])
text(-575,0.9,'sim')
text(-575,0.95,'diss')
text(-575,0.85,'contra')

yl.LabelVerticalAlignment = 'middle';
yl.LabelHorizontalAlignment = 'right';
yl.FontSize=12;
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
ylabel('Correlation chosen value vectors')
ylim([-0.2 1])

subplot(3,3,2)
hold on
shadedErrorBar(window_start_list(4:end),mean(mean(LIP.r_chosen_contra_boot(4:end,:,:),2),3),[prctile(mean(LIP.r_chosen_contra_boot(4:end,:,:),2),95,3) prctile(mean(LIP.r_chosen_contra_boot(4:end,:,:),2),5,3)]',{'-','color',[0.5 0.5 0.5],'LineWidth',2},2)
shadedErrorBar(window_start_list(4:end),mean(mean(LIP.r_chosen_contra_ipsi_boot(4:end,:,:),2),3),[prctile(mean(LIP.r_chosen_contra_ipsi_boot(4:end,:,:),2),95,3) prctile(mean(LIP.r_chosen_contra_ipsi_boot(4:end,:,:),2),5,3)]',{'-','color',color_for_ROI(n_color,:),'LineWidth',2},2)
shadedErrorBar(window_start_list(4:end),mean(mean(LIP.r_chosen_ipsi_boot(4:end,:,:),2),3),[prctile(mean(LIP.r_chosen_ipsi_boot(4:end,:,:),2),95,3) prctile(mean(LIP.r_chosen_ipsi_boot(4:end,:,:),2),5,3)]',{'--','color',[0.5 0.5 0.5],'LineWidth',2},2)
plot_significance_level(window_start_list,LIP.p_chosen_contra_ipsi_sim,0.9,color_for_ROI(n_color,:),[0.01,0.05/(length(window_start_list)-3), 0.01/(length(window_start_list)-3)])
plot_significance_level(window_start_list,LIP.p_chosen_contra_ipsi_diss,0.95,color_for_ROI(n_color,:),[0.01,0.05/(length(window_start_list)-3), 0.01/(length(window_start_list)-3)])
plot_significance_level(window_start_list,p_chosen_contra_boot,0.85,[0.5 0.5 0.5],[0.01,0.05/(length(window_start_list)-3), 0.01/(length(window_start_list)-3)])
plot_significance_level(window_start_list,p_chosen_ipsi_boot,0.8,[0.5 0.5 0.5],[0.01,0.05/(length(window_start_list)-3), 0.01/(length(window_start_list)-3)])

text(-575,0.9,'sim')
text(-575,0.95,'diss')
text(-575,0.85,'contra')
text(-575,0.8,'ipsi')

text(-575,0.9,'sim')
text(-575,0.95,'diss')
yl.LabelVerticalAlignment = 'middle';
yl.LabelHorizontalAlignment = 'right';
yl.FontSize=12;
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
ylabel('Correlation chosen value vectors')
ylim([-0.2 1])

subplot(3,3,3)

hold on
shadedErrorBar(window_start_list(4:end),mean(mean(LIP.r_chosen_contra_boot(4:end,:,:),2),3),[prctile(mean(LIP.r_chosen_contra_boot(4:end,:,:),2),95,3) prctile(mean(LIP.r_chosen_contra_boot(4:end,:,:),2),5,3)]',{'-','color',[0.5 0.5 0.5],'LineWidth',2},2)
shadedErrorBar(window_start_list(4:end),mean(mean(LIP.r_chosen_contra_contra_diag_boot(4:end,:,:),2),3),[prctile(mean(LIP.r_chosen_contra_contra_diag_boot(4:end,:,:),2),95,3) prctile(mean(LIP.r_chosen_contra_contra_diag_boot(4:end,:,:),2),5,3)]',{'-','color',color_for_ROI(n_color,:),'LineWidth',2},2)
shadedErrorBar(window_start_list(4:end),mean(mean(LIP.r_chosen_ipsi_boot(4:end,:,:),2),3),[prctile(mean(LIP.r_chosen_ipsi_boot(4:end,:,:),2),95,3) prctile(mean(LIP.r_chosen_ipsi_boot(4:end,:,:),2),5,3)]',{'--','color',[0.5 0.5 0.5],'LineWidth',2},2)
plot_significance_level(window_start_list,LIP.p_chosen_contra_contra_diag_sim,0.9,color_for_ROI(n_color,:),[0.01,0.05/(length(window_start_list)-3), 0.01/(length(window_start_list)-3)])
plot_significance_level(window_start_list,LIP.p_chosen_contra_contra_diag_diss,0.95,color_for_ROI(n_color,:),[0.01,0.05/(length(window_start_list)-3), 0.01/(length(window_start_list)-3)])
plot_significance_level(window_start_list,p_chosen_contra_boot,0.85,[0.5 0.5 0.5],[0.01,0.05/(length(window_start_list)-3), 0.01/(length(window_start_list)-3)])
plot_significance_level(window_start_list,p_chosen_ipsi_boot,0.8,[0.5 0.5 0.5],[0.01,0.05/(length(window_start_list)-3), 0.01/(length(window_start_list)-3)])

text(-575,0.9,'sim')
text(-575,0.95,'diss')
text(-575,0.85,'contra')
text(-575,0.8,'ipsi')

yl.LabelVerticalAlignment = 'middle';
yl.LabelHorizontalAlignment = 'right';
yl.FontSize=12;
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
ylabel('Correlation chosen value vectors')
ylim([-0.2 1])

n_color=2;
subplot(3,3,4)
hold on
shadedErrorBar(window_start_list(4:end),mean(mean(FEF.r_chosen_contra_boot(4:end,:,:),2),3),[prctile(mean(FEF.r_chosen_contra_boot(4:end,:,:),2),95,3) prctile(mean(FEF.r_chosen_contra_boot(4:end,:,:),2),5,3)]',{'-','color',[0.5 0.5 0.5],'LineWidth',2},2)
shadedErrorBar(window_start_list(4:end),mean(mean(FEF.r_chosen_contra_contra_boot(4:end,:,:),2),3),[prctile(mean(FEF.r_chosen_contra_contra_boot(4:end,:,:),2),95,3) prctile(mean(FEF.r_chosen_contra_contra_boot(4:end,:,:),2),5,3)]',{'-','color',color_for_ROI(n_color,:),'LineWidth',2},2)
plot_significance_level(window_start_list,FEF.p_chosen_contra_contra_sim,0.9,color_for_ROI(n_color,:),[0.01,0.05/(length(window_start_list)-3), 0.01/(length(window_start_list)-3)])
plot_significance_level(window_start_list,FEF.p_chosen_contra_contra_diss,0.95,color_for_ROI(n_color,:),[0.01,0.05/(length(window_start_list)-3), 0.01/(length(window_start_list)-3)])
plot_significance_level(window_start_list,p_chosen_contra_boot,0.85,[0.5 0.5 0.5],[0.01,0.05/(length(window_start_list)-3), 0.01/(length(window_start_list)-3)])
text(-575,0.9,'sim')
text(-575,0.95,'diss')
text(-575,0.85,'contra')

yl.LabelVerticalAlignment = 'middle';
yl.LabelHorizontalAlignment = 'right';
yl.FontSize=12;
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
ylabel('Correlation chosen value vectors')
ylim([-0.2 1])

subplot(3,3,5)
hold on
shadedErrorBar(window_start_list(4:end),mean(mean(FEF.r_chosen_contra_boot(4:end,:,:),2),3),[prctile(mean(FEF.r_chosen_contra_boot(4:end,:,:),2),95,3) prctile(mean(FEF.r_chosen_contra_boot(4:end,:,:),2),5,3)]',{'-','color',[0.5 0.5 0.5],'LineWidth',2},2)
shadedErrorBar(window_start_list(4:end),mean(mean(FEF.r_chosen_contra_ipsi_boot(4:end,:,:),2),3),[prctile(mean(FEF.r_chosen_contra_ipsi_boot(4:end,:,:),2),95,3) prctile(mean(FEF.r_chosen_contra_ipsi_boot(4:end,:,:),2),5,3)]',{'-','color',color_for_ROI(n_color,:),'LineWidth',2},2)
shadedErrorBar(window_start_list(4:end),mean(mean(FEF.r_chosen_ipsi_boot(4:end,:,:),2),3),[prctile(mean(FEF.r_chosen_ipsi_boot(4:end,:,:),2),95,3) prctile(mean(FEF.r_chosen_ipsi_boot(4:end,:,:),2),5,3)]',{'--','color',[0.5 0.5 0.5],'LineWidth',2},2)
plot_significance_level(window_start_list,FEF.p_chosen_contra_ipsi_sim,0.9,color_for_ROI(n_color,:),[0.01,0.05/(length(window_start_list)-3), 0.01/(length(window_start_list)-3)])
plot_significance_level(window_start_list,FEF.p_chosen_contra_ipsi_diss,0.95,color_for_ROI(n_color,:),[0.01,0.05/(length(window_start_list)-3), 0.01/(length(window_start_list)-3)])
plot_significance_level(window_start_list,p_chosen_contra_boot,0.85,[0.5 0.5 0.5],[0.01,0.05/(length(window_start_list)-3), 0.01/(length(window_start_list)-3)])
plot_significance_level(window_start_list,p_chosen_ipsi_boot,0.8,[0.5 0.5 0.5],[0.01,0.05/(length(window_start_list)-3), 0.01/(length(window_start_list)-3)])

text(-575,0.9,'sim')
text(-575,0.95,'diss')
text(-575,0.85,'contra')
text(-575,0.8,'ipsi')

text(-575,0.9,'sim')
text(-575,0.95,'diss')
yl.LabelVerticalAlignment = 'middle';
yl.LabelHorizontalAlignment = 'right';
yl.FontSize=12;
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
ylabel('Correlation chosen value vectors')
ylim([-0.2 1])

subplot(3,3,6)

hold on
shadedErrorBar(window_start_list(4:end),mean(mean(FEF.r_chosen_contra_boot(4:end,:,:),2),3),[prctile(mean(FEF.r_chosen_contra_boot(4:end,:,:),2),95,3) prctile(mean(FEF.r_chosen_contra_boot(4:end,:,:),2),5,3)]',{'-','color',[0.5 0.5 0.5],'LineWidth',2},2)
shadedErrorBar(window_start_list(4:end),mean(mean(FEF.r_chosen_contra_contra_diag_boot(4:end,:,:),2),3),[prctile(mean(FEF.r_chosen_contra_contra_diag_boot(4:end,:,:),2),95,3) prctile(mean(FEF.r_chosen_contra_contra_diag_boot(4:end,:,:),2),5,3)]',{'-','color',color_for_ROI(n_color,:),'LineWidth',2},2)
shadedErrorBar(window_start_list(4:end),mean(mean(FEF.r_chosen_ipsi_boot(4:end,:,:),2),3),[prctile(mean(FEF.r_chosen_ipsi_boot(4:end,:,:),2),95,3) prctile(mean(FEF.r_chosen_ipsi_boot(4:end,:,:),2),5,3)]',{'--','color',[0.5 0.5 0.5],'LineWidth',2},2)
plot_significance_level(window_start_list,FEF.p_chosen_contra_contra_diag_sim,0.9,color_for_ROI(n_color,:),[0.01,0.05/(length(window_start_list)-3), 0.01/(length(window_start_list)-3)])
plot_significance_level(window_start_list,FEF.p_chosen_contra_contra_diag_diss,0.95,color_for_ROI(n_color,:),[0.01,0.05/(length(window_start_list)-3), 0.01/(length(window_start_list)-3)])
plot_significance_level(window_start_list,p_chosen_contra_boot,0.85,[0.5 0.5 0.5],[0.01,0.05/(length(window_start_list)-3), 0.01/(length(window_start_list)-3)])
plot_significance_level(window_start_list,p_chosen_ipsi_boot,0.8,[0.5 0.5 0.5],[0.01,0.05/(length(window_start_list)-3), 0.01/(length(window_start_list)-3)])

text(-575,0.9,'sim')
text(-575,0.95,'diss')
text(-575,0.85,'contra')
text(-575,0.8,'ipsi')

yl.LabelVerticalAlignment = 'middle';
yl.LabelHorizontalAlignment = 'right';
yl.FontSize=12;
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
ylabel('Correlation chosen value vectors')
ylim([-0.2 1])

n_color=3;
subplot(3,3,7)
hold on
shadedErrorBar(window_start_list(4:end),mean(mean(PFC.r_chosen_contra_boot(4:end,:,:),2),3),[prctile(mean(PFC.r_chosen_contra_boot(4:end,:,:),2),95,3) prctile(mean(PFC.r_chosen_contra_boot(4:end,:,:),2),5,3)]',{'-','color',[0.5 0.5 0.5],'LineWidth',2},2)
shadedErrorBar(window_start_list(4:end),mean(mean(PFC.r_chosen_contra_contra_boot(4:end,:,:),2),3),[prctile(mean(PFC.r_chosen_contra_contra_boot(4:end,:,:),2),95,3) prctile(mean(PFC.r_chosen_contra_contra_boot(4:end,:,:),2),5,3)]',{'-','color',color_for_ROI(n_color,:),'LineWidth',2},2)
plot_significance_level(window_start_list,PFC.p_chosen_contra_contra_sim,0.9,color_for_ROI(n_color,:),[0.01,0.05/(length(window_start_list)-3), 0.01/(length(window_start_list)-3)])
plot_significance_level(window_start_list,PFC.p_chosen_contra_contra_diss,0.95,color_for_ROI(n_color,:),[0.01,0.05/(length(window_start_list)-3), 0.01/(length(window_start_list)-3)])
plot_significance_level(window_start_list,p_chosen_contra_boot,0.85,[0.5 0.5 0.5],[0.01,0.05/(length(window_start_list)-3), 0.01/(length(window_start_list)-3)])
text(-575,0.9,'sim')
text(-575,0.95,'diss')
text(-575,0.85,'contra')

yl.LabelVerticalAlignment = 'middle';
yl.LabelHorizontalAlignment = 'right';
yl.FontSize=12;
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
ylabel('Correlation chosen value vectors')
ylim([-0.2 1])

subplot(3,3,8)
hold on
shadedErrorBar(window_start_list(4:end),mean(mean(PFC.r_chosen_contra_boot(4:end,:,:),2),3),[prctile(mean(PFC.r_chosen_contra_boot(4:end,:,:),2),95,3) prctile(mean(PFC.r_chosen_contra_boot(4:end,:,:),2),5,3)]',{'-','color',[0.5 0.5 0.5],'LineWidth',2},2)
shadedErrorBar(window_start_list(4:end),mean(mean(PFC.r_chosen_contra_ipsi_boot(4:end,:,:),2),3),[prctile(mean(PFC.r_chosen_contra_ipsi_boot(4:end,:,:),2),95,3) prctile(mean(PFC.r_chosen_contra_ipsi_boot(4:end,:,:),2),5,3)]',{'-','color',color_for_ROI(n_color,:),'LineWidth',2},2)
shadedErrorBar(window_start_list(4:end),mean(mean(PFC.r_chosen_ipsi_boot(4:end,:,:),2),3),[prctile(mean(PFC.r_chosen_ipsi_boot(4:end,:,:),2),95,3) prctile(mean(PFC.r_chosen_ipsi_boot(4:end,:,:),2),5,3)]',{'--','color',[0.5 0.5 0.5],'LineWidth',2},2)
plot_significance_level(window_start_list,PFC.p_chosen_contra_ipsi_sim,0.9,color_for_ROI(n_color,:),[0.01,0.05/(length(window_start_list)-3), 0.01/(length(window_start_list)-3)])
plot_significance_level(window_start_list,PFC.p_chosen_contra_ipsi_diss,0.95,color_for_ROI(n_color,:),[0.01,0.05/(length(window_start_list)-3), 0.01/(length(window_start_list)-3)])
plot_significance_level(window_start_list,p_chosen_contra_boot,0.85,[0.5 0.5 0.5],[0.01,0.05/(length(window_start_list)-3), 0.01/(length(window_start_list)-3)])
plot_significance_level(window_start_list,p_chosen_ipsi_boot,0.8,[0.5 0.5 0.5],[0.01,0.05/(length(window_start_list)-3), 0.01/(length(window_start_list)-3)])

text(-575,0.9,'sim')
text(-575,0.95,'diss')
text(-575,0.85,'contra')
text(-575,0.8,'ipsi')

text(-575,0.9,'sim')
text(-575,0.95,'diss')
yl.LabelVerticalAlignment = 'middle';
yl.LabelHorizontalAlignment = 'right';
yl.FontSize=12;
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
ylabel('Correlation chosen value vectors')
ylim([-0.2 1])

subplot(3,3,9)

hold on
shadedErrorBar(window_start_list(4:end),mean(mean(PFC.r_chosen_contra_boot(4:end,:,:),2),3),[prctile(mean(PFC.r_chosen_contra_boot(4:end,:,:),2),95,3) prctile(mean(PFC.r_chosen_contra_boot(4:end,:,:),2),5,3)]',{'-','color',[0.5 0.5 0.5],'LineWidth',2},2)
shadedErrorBar(window_start_list(4:end),mean(mean(PFC.r_chosen_contra_contra_diag_boot(4:end,:,:),2),3),[prctile(mean(PFC.r_chosen_contra_contra_diag_boot(4:end,:,:),2),95,3) prctile(mean(PFC.r_chosen_contra_contra_diag_boot(4:end,:,:),2),5,3)]',{'-','color',color_for_ROI(n_color,:),'LineWidth',2},2)
shadedErrorBar(window_start_list(4:end),mean(mean(PFC.r_chosen_ipsi_boot(4:end,:,:),2),3),[prctile(mean(PFC.r_chosen_ipsi_boot(4:end,:,:),2),95,3) prctile(mean(PFC.r_chosen_ipsi_boot(4:end,:,:),2),5,3)]',{'--','color',[0.5 0.5 0.5],'LineWidth',2},2)
plot_significance_level(window_start_list,PFC.p_chosen_contra_contra_diag_sim,0.9,color_for_ROI(n_color,:),[0.01,0.05/(length(window_start_list)-3), 0.01/(length(window_start_list)-3)])
plot_significance_level(window_start_list,PFC.p_chosen_contra_contra_diag_diss,0.95,color_for_ROI(n_color,:),[0.01,0.05/(length(window_start_list)-3), 0.01/(length(window_start_list)-3)])
plot_significance_level(window_start_list,p_chosen_contra_boot,0.85,[0.5 0.5 0.5],[0.01,0.05/(length(window_start_list)-3), 0.01/(length(window_start_list)-3)])
plot_significance_level(window_start_list,p_chosen_ipsi_boot,0.8,[0.5 0.5 0.5],[0.01,0.05/(length(window_start_list)-3), 0.01/(length(window_start_list)-3)])

text(-575,0.9,'sim')
text(-575,0.95,'diss')
text(-575,0.85,'contra')
text(-575,0.8,'ipsi')

yl.LabelVerticalAlignment = 'middle';
yl.LabelHorizontalAlignment = 'right';
yl.FontSize=12;
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
ylabel('Correlation chosen value vectors')
ylim([-0.2 1])


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

function p = z_test_function_bootstrap(dist,null)

if sum(~isnan(dist))>2
    
    dist=dist(~isnan(dist));
    m = mean(dist);
    s = std(dist);
    z = (m-null)/s;
    p=1-normcdf(z);
    
else
    p=NaN;
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


