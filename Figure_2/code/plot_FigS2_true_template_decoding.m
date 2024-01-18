clear all;
fsroot='/Volumes/buschman';
event='target';

for i=1:7
    initial_window=0;
    event_list{i}='reward_end';
    window_start_list(i)=initial_window+(i-1)*50;
end

for i=1:20
    initial_window=-600;
    event_list{i+7}='target';
    window_start_list(i+7)=initial_window+(i-1)*50;
end

task='Learning_Attentional_Templates';
subtask='exploreexploit/Restricted_FEF/Reset_RW_model';

subsubtask_classifier=sprintf('Classifier_%s_%d_to_%s_%d',event_list{1},window_start_list(1),event_list{end},window_start_list(end));

dirstem = fullfile('Projects',task,'Analysis','Electrophy_analysis',subtask,subsubtask_classifier);
data_path_clasifier = fullfile(fsroot,dirstem);

color_for_ROI=[0 0.447 0.741; 0.85 0.325 0.098; 0.929 0.694 0.125]; %LIP FEF PFC

load('colors')

color_template(1,:)=colors(1,:);
color_template(2,:)=colors(14,:);
color_template(3,:)=colors(27,:);

this_combination=2;

LOAD=1;
if this_combination==1
    save_name=sprintf('Pseudo_pop_true_template_prog_combined_time_results_more_trials_%s_600ms_proj','PFC');
elseif this_combination==2
    save_name=sprintf('Pseudo_pop_true_template_prog_combined_time_results_more_trials_%s_900ms_proj','PFC');
end

%%

if LOAD==1

for n_tt=1:100
    
    ROI='LIP';
    if this_combination==1
        results_name=sprintf('Pseudo_pop_true_template_prog_combined_time_results_more_trials_%s_600ms_%d',ROI,n_tt);
    elseif this_combination==2
        results_name=sprintf('Pseudo_pop_true_template_prog_combined_time_results_more_trials_%s_900ms_%d',ROI,n_tt);
        results_NN_name=sprintf('Pseudo_pop_true_template_prog_combined_time_results_more_trials_%s_900ms_NN_%d',ROI,n_tt);
    end
    
    if exist(sprintf('%s.mat',fullfile(data_path_clasifier,results_name)))==2
        
        load(fullfile(data_path_clasifier,results_name),'Classification_correct','W')
        load(fullfile(data_path_clasifier,results_NN_name),'Classification_net_correct')
        
        LIP.Classification_net_correct(:,:,n_tt)=Classification_net_correct;
        LIP.Classification_correct(:,:,n_tt)=Classification_correct;
        LIP.W(:,:,n_tt)=W;
        
%         LIP.Classifier_projection(:,:,:,n_tt) = get_projection_classifier_combined_time_prog_NN(fsroot, 'LIP', n_tt, this_combination);
        
        clear Classification*_correct W 
        
    else
        LIP.Classification_net_correct(:,:,n_tt)=NaN;
        LIP.Classification_correct(:,:,n_tt)=NaN;
        LIP.W(:,:,n_tt)=NaN;
        sprintf('%s %d missing',ROI,n_tt)
        
    end
    
    %FEF
    ROI='FEF';
    
    if this_combination==1
        results_name=sprintf('Pseudo_pop_true_template_prog_combined_time_results_more_trials_%s_600ms_%d',ROI,n_tt);
    elseif this_combination==2
        results_name=sprintf('Pseudo_pop_true_template_prog_combined_time_results_more_trials_%s_900ms_%d',ROI,n_tt);
        results_NN_name=sprintf('Pseudo_pop_true_template_prog_combined_time_results_more_trials_%s_900ms_NN_%d',ROI,n_tt);
    end
    
    if exist(sprintf('%s.mat',fullfile(data_path_clasifier,results_name)))==2
        
        load(fullfile(data_path_clasifier,results_name),'Classification_correct','W')
                load(fullfile(data_path_clasifier,results_NN_name),'Classification_net_correct')
        
        FEF.Classification_net_correct(:,:,n_tt)=Classification_net_correct;

        FEF.Classification_correct(:,:,n_tt)=Classification_correct;
        FEF.W(:,:,n_tt)=W;
        
%         FEF.Classifier_projection(:,:,:,n_tt) = get_projection_classifier_combined_time_prog_NN(fsroot, 'FEF', n_tt, this_combination);
        
        clear Classification*_correct W
        
    else
        FEF.Classification_net_correct(:,:,n_tt)=NaN;
        FEF.Classification_correct(:,:,n_tt)=NaN;
        FEF.W(:,:,:,n_tt)=NaN;
        
        sprintf('%s %d missing',ROI,n_tt)
        
    end
    
    %PFC
    ROI='PFC';
    if this_combination==1
        results_name=sprintf('Pseudo_pop_true_template_prog_combined_time_results_more_trials_%s_600ms_%d',ROI,n_tt);
    elseif this_combination==2
        results_name=sprintf('Pseudo_pop_true_template_prog_combined_time_results_more_trials_%s_900ms_%d',ROI,n_tt);
        results_NN_name=sprintf('Pseudo_pop_true_template_prog_combined_time_results_more_trials_%s_900ms_NN_%d',ROI,n_tt);
    end
    
    if exist(sprintf('%s.mat',fullfile(data_path_clasifier,results_name)))==2
        
        
        load(fullfile(data_path_clasifier,results_name),'Classification_correct','W')
        load(fullfile(data_path_clasifier,results_NN_name),'Classification_net_correct')
        
        PFC.Classification_net_correct(:,:,n_tt)=Classification_net_correct;
        PFC.Classification_correct(:,:,n_tt)=Classification_correct;
        PFC.W(:,:,n_tt)=W;
        
%         PFC.Classifier_projection(:,:,:,n_tt) = get_projection_classifier_combined_time_prog_NN(fsroot, 'PFC', n_tt, this_combination);
        
        clear Classification*_correct W
        
    else
        PFC.Classification_net_correct(:,:,n_tt)=NaN;
        PFC.Classification_correct(:,:,n_tt)=NaN;
        PFC.W(:,:,:,n_tt)=NaN;
        sprintf('%s %d missing',ROI,n_tt)
        
    end
    
    
end

% LIP.Mean_classifier_projection=mean(LIP.Classifier_projection,1);
% FEF.Mean_classifier_projection=mean(FEF.Classifier_projection,1);
% PFC.Mean_classifier_projection=mean(PFC.Classifier_projection,1);



LIP_mean_classification(:,:)=nanmean(LIP.Classification_correct,2);
FEF_mean_classification(:,:)=nanmean(FEF.Classification_correct,2);
PFC_mean_classification(:,:)=nanmean(PFC.Classification_correct,2);

LIP_mean_classification_net(:,:)=nanmean(LIP.Classification_net_correct,2);
FEF_mean_classification_net(:,:)=nanmean(FEF.Classification_net_correct,2);
PFC_mean_classification_net(:,:)=nanmean(PFC.Classification_net_correct,2);

% save(fullfile(data_path_clasifier,save_name));

else
    
    load(fullfile(data_path_clasifier,save_name));
    
    
end

%%
figure
subplot(2,3,1)
vs = violinplot(LIP_mean_classification',{'1/3','2/3','3/3'},'ViolinColor',color_for_ROI(1,:),'ShowData',false);
hold on
yline(1/3,'--',{'Chance level'},'Color',[0.5 0.5 0.5],'LabelVerticalAlignment','bottom','LabelHorizontalAlignment','center')

box off
title('LIP')
ylabel('Expected template classification accuracy')
xlabel('Progression in block')
ylim([0.25 1])

subplot(2,3,2)
vs = violinplot(FEF_mean_classification',{'1/3','2/3','3/3'},'ViolinColor',color_for_ROI(2,:),'ShowData',false);
hold on
yline(1/3,'--',{'Chance level'},'Color',[0.5 0.5 0.5],'LabelVerticalAlignment','bottom','LabelHorizontalAlignment','center')
box off
title('FEF')
ylabel('Expected template classification accuracy')
xlabel('Progression in block')
ylim([0.25 1])

subplot(2,3,3)
vs = violinplot(PFC_mean_classification',{'1/3','2/3','3/3'},'ViolinColor',color_for_ROI(3,:),'ShowData',false);
hold on
yline(1/3,'--',{'Chance level'},'Color',[0.5 0.5 0.5],'LabelVerticalAlignment','bottom','LabelHorizontalAlignment','center')
box off
title('PFC')
ylabel('Expected template classification accuracy')
xlabel('Progression in block')
ylim([0.25 1])

subplot(2,3,4)
vs = violinplot(LIP_mean_classification_net',{'1/3','2/3','3/3'},'ViolinColor',color_for_ROI(1,:),'ShowData',false);
hold on
yline(1/3,'--',{'Chance level'},'Color',[0.5 0.5 0.5],'LabelVerticalAlignment','bottom','LabelHorizontalAlignment','center')
box off
title('LIP')
ylabel('Expected template classification accuracy NN')
xlabel('Progression in block')
ylim([0.25 1])

subplot(2,3,5)
vs = violinplot(FEF_mean_classification_net',{'1/3','2/3','3/3'},'ViolinColor',color_for_ROI(2,:),'ShowData',false);
hold on
yline(1/3,'--',{'Chance level'},'Color',[0.5 0.5 0.5],'LabelVerticalAlignment','bottom','LabelHorizontalAlignment','center')
box off
title('FEF')
ylabel('Expected template classification accuracy NN')
xlabel('Progression in block')
ylim([0.25 1])

subplot(2,3,6)
vs = violinplot(PFC_mean_classification_net',{'1/3','2/3','3/3'},'ViolinColor',color_for_ROI(3,:),'ShowData',false);
hold on
yline(1/3,'--',{'Chance level'},'Color',[0.5 0.5 0.5],'LabelVerticalAlignment','bottom','LabelHorizontalAlignment','center')
box off
title('PFC')
ylabel('Expected template classification accuracy NN')
xlabel('Progression in block')
ylim([0.25 1])


%%
figure
subplot(3,1,1)
vs = violinplot(LIP_mean_classification_net',{'1/3','2/3','3/3'},'ViolinColor',color_for_ROI(1,:),'ShowData',false);
hold on
yline(1/3,'--',{'Chance level'},'Color',[0.5 0.5 0.5],'LabelVerticalAlignment','middle','LabelHorizontalAlignment','center')
box off
title('LIP','FontSize',12)
ylabel(sprintf('Expected template\nclassification accuracy'),'FontSize',12)
xlabel('Progression in block','FontSize',12)
ylim([0 1])

subplot(3,1,2)
vs = violinplot(FEF_mean_classification_net',{'1/3','2/3','3/3'},'ViolinColor',color_for_ROI(2,:),'ShowData',false);
hold on
yline(1/3,'--',{'Chance level'},'Color',[0.5 0.5 0.5],'LabelVerticalAlignment','middle','LabelHorizontalAlignment','center')
box off
title('FEF','FontSize',12)
ylabel(sprintf('Expected template\nclassification accuracy'),'FontSize',12)
xlabel('Progression in block','FontSize',12)
ylim([0 1])

subplot(3,1,3)
vs = violinplot(PFC_mean_classification_net',{'1/3','2/3','3/3'},'ViolinColor',color_for_ROI(3,:),'ShowData',false);
hold on
yline(1/3,'--',{'Chance level'},'Color',[0.5 0.5 0.5],'LabelVerticalAlignment','middle','LabelHorizontalAlignment','center')
box off
title('PFC','FontSize',12)
ylabel(sprintf('Expected template\nclassification accuracy'),'FontSize',12)
xlabel('Progression in block','FontSize',12)
ylim([0 1])

%%
figure
subplot(1,3,1)
vs = violinplot(LIP_mean_classification_net',{'1/3','2/3','3/3'},'ViolinColor',color_for_ROI(1,:),'ShowData',false,'ViolinAlpha',0.5);
hold on
yline(1/3,'--',{'Chance level'},'Color',[0.5 0.5 0.5],'LabelVerticalAlignment','middle','LabelHorizontalAlignment','center')
box off
title('LIP','FontSize',12)
ylabel(sprintf('Expected template\nclassification accuracy'),'FontSize',12)
xlabel('Progression in block','FontSize',12)
ylim([0 1])

subplot(1,3,2)
vs = violinplot(FEF_mean_classification_net',{'1/3','2/3','3/3'},'ViolinColor',color_for_ROI(2,:),'ShowData',false,'ViolinAlpha',0.5);
hold on
yline(1/3,'--',{'Chance level'},'Color',[0.5 0.5 0.5],'LabelVerticalAlignment','middle','LabelHorizontalAlignment','center')
box off
title('FEF','FontSize',12)
ylabel(sprintf('Expected template\nclassification accuracy'),'FontSize',12)
xlabel('Progression in block','FontSize',12)
ylim([0 1])

subplot(1,3,3)
vs = violinplot(PFC_mean_classification_net',{'1/3','2/3','3/3'},'ViolinColor',color_for_ROI(3,:),'ShowData',false,'ViolinAlpha',0.5);
hold on
yline(1/3,'--',{'Chance level'},'Color',[0.5 0.5 0.5],'LabelVerticalAlignment','middle','LabelHorizontalAlignment','center')
box off
title('PFC','FontSize',12)
ylabel(sprintf('Expected template\nclassification accuracy'),'FontSize',12)
xlabel('Progression in block','FontSize',12)
ylim([0 1])


%%
figure
subplot(3,1,1)
vs = violinplot(LIP_mean_classification_net',{'1/3','2/3','3/3'},'ViolinColor',color_for_ROI(1,:),'ShowData',false,'ViolinAlpha',0.5);
hold on
yline(1/3,'--',{'Chance level'},'Color',[0.5 0.5 0.5],'LabelVerticalAlignment','middle','LabelHorizontalAlignment','center')
box off
% title('LIP','FontSize',12)
ylabel(sprintf('Expected template\nclassification accuracy'),'FontSize',12)
xlabel('Progression in block','FontSize',12)
ylim([0 1])

subplot(3,1,2)
vs = violinplot(FEF_mean_classification_net',{'1/3','2/3','3/3'},'ViolinColor',color_for_ROI(2,:),'ShowData',false,'ViolinAlpha',0.5);
hold on
yline(1/3,'--',{'Chance level'},'Color',[0.5 0.5 0.5],'LabelVerticalAlignment','middle','LabelHorizontalAlignment','center')
box off
% title('FEF','FontSize',12)
ylabel(sprintf('Expected template\nclassification accuracy'),'FontSize',12)
xlabel('Progression in block','FontSize',12)
ylim([0 1])

subplot(3,1,3)
vs = violinplot(PFC_mean_classification_net',{'1/3','2/3','3/3'},'ViolinColor',color_for_ROI(3,:),'ShowData',false,'ViolinAlpha',0.5);
hold on
yline(1/3,'--',{'Chance level'},'Color',[0.5 0.5 0.5],'LabelVerticalAlignment','middle','LabelHorizontalAlignment','center')
box off
% title('PFC','FontSize',12)
ylabel(sprintf('Expected template\nclassification accuracy'),'FontSize',12)
xlabel('Progression in block','FontSize',12)
ylim([0 1])


%%
for k=1:3
    p_LIP(k)=z_test_function_bootstrap(LIP_mean_classification_net(k,:),1/3);
    p_FEF(k)=z_test_function_bootstrap(FEF_mean_classification_net(k,:),1/3);
    p_PFC(k)=z_test_function_bootstrap(PFC_mean_classification_net(k,:),1/3);
end


%%
% p_LIP=z_test_function_bootstrap(LIP_mean_classification_net(3,:)-LIP_mean_classification_net(1,:),0)
% p_FEF=z_test_function_bootstrap(FEF_mean_classification_net(3,:)-FEF_mean_classification_net(1,:),0)
% p_PFC=z_test_function_bootstrap(PFC_mean_classification_net(3,:)-PFC_mean_classification_net(1,:),0)
% 
% %%
% for n=1:100
%    b_LIP(:,n)=glmfit(zscore(LIP_mean_classification_net(:,n)),zscore([1 2 3]));
%    b_FEF(:,n)=glmfit(zscore(FEF_mean_classification_net(:,n)),zscore([1 2 3]));
%    b_PFC(:,n)=glmfit(zscore(PFC_mean_classification_net(:,n)),zscore([1 2 3]));
% end
% 
% %%
% p_lin_LIP=z_test_function_bootstrap(b_LIP(2,:),0)
% p_lin_FEF=z_test_function_bootstrap(b_FEF(2,:),0)
% p_lin_PFC=z_test_function_bootstrap(b_PFC(2,:),0)

% %% 3d ellipses
% 
% Mean_LIP_proj(:,:)=mean(LIP.Mean_classifier_projection,4);
% Std_LIP_proj(:,:)=std(LIP.Mean_classifier_projection,0,4);
% 
% Mean_FEF_proj(:,:)=mean(FEF.Mean_classifier_projection,4);
% Std_FEF_proj(:,:)=std(FEF.Mean_classifier_projection,0,4);
% 
% Mean_PFC_proj(:,:)=mean(PFC.Mean_classifier_projection,4);
% Std_PFC_proj(:,:)=std(PFC.Mean_classifier_projection,0,4);
% 
% load('colors')
% 
% color_template(1,:)=colors(1,:);
% color_template(2,:)=colors(14,:);
% color_template(3,:)=colors(28,:);
% 
% %%
% f=figure
% subplot(3,1,1)
% hold on
% for j=1:3
%     [X, Y, Z]=ellipsoid(Mean_LIP_proj(j,1),Mean_LIP_proj(j,3),Mean_LIP_proj(j,2),Std_LIP_proj(j,1),Std_LIP_proj(j,3),Std_LIP_proj(j,2));%
%     line([0.5 Mean_LIP_proj(j,1)], [0.5, Mean_LIP_proj(j,3)], [0.5 Mean_LIP_proj(j,2)],'Color',color_for_ROI(1,:),'LineWidth',2)
%     plot3(Mean_LIP_proj(j,1),Mean_LIP_proj(j,3),Mean_LIP_proj(j,2),'o','color',color_template(j,:),'MarkerFacecolor',color_template(j,:),'MarkerSize',10)
%     surf(X,Y,Z,'Facecolor',color_template(j,:),'FaceAlpha',0.35,'EdgeColor','none')
% end
% xlabel('Pink template')
% ylabel('Blue template')
% zlabel('Brown template')
% ylim([0 1])
% xlim([0 1])
% zlim([0 1])
% grid on
% xticks(0:0.2:1)
% yticks(0:0.2:1)
% zticks(0:0.2:1)
% view([-0.5 3 1])
% 
% subplot(3,1,2)
% hold on
% for j=1:3
%     [X, Y, Z]=ellipsoid(Mean_FEF_proj(j,1),Mean_FEF_proj(j,3),Mean_FEF_proj(j,2),Std_FEF_proj(j,1),Std_FEF_proj(j,3),Std_FEF_proj(j,2));%
%     line([0.5 Mean_FEF_proj(j,1)], [0.5, Mean_FEF_proj(j,3)], [0.5 Mean_FEF_proj(j,2)],'Color',color_for_ROI(2,:),'LineWidth',2)
%     plot3(Mean_FEF_proj(j,1),Mean_FEF_proj(j,3),Mean_FEF_proj(j,2),'o','color',color_template(j,:),'MarkerFacecolor',color_template(j,:),'MarkerSize',10)
%     surf(X,Y,Z,'Facecolor',color_template(j,:),'FaceAlpha',0.35,'EdgeColor','none')
% end
% xlabel('Pink template')
% ylabel('Blue template')
% zlabel('Brown template')
% ylim([0 1])
% xlim([0 1])
% zlim([0 1])
% grid on
% xticks(0:0.2:1)
% yticks(0:0.2:1)
% zticks(0:0.2:1)
% view([-0.5 3 1])
% 
% subplot(3,1,3)
% hold on
% for j=1:3
%     [X, Y, Z]=ellipsoid(Mean_PFC_proj(j,1),Mean_PFC_proj(j,3),Mean_PFC_proj(j,2),Std_PFC_proj(j,1),Std_PFC_proj(j,3),Std_PFC_proj(j,2));%
%     line([0.5 Mean_PFC_proj(j,1)], [0.5, Mean_PFC_proj(j,3)], [0.5 Mean_PFC_proj(j,2)],'Color',color_for_ROI(3,:),'LineWidth',2)
%     plot3(Mean_PFC_proj(j,1),Mean_PFC_proj(j,3),Mean_PFC_proj(j,2),'o','color',color_template(j,:),'MarkerFacecolor',color_template(j,:),'MarkerSize',10)
%     surf(X,Y,Z,'Facecolor',color_template(j,:),'FaceAlpha',0.35,'EdgeColor','none')
% end
% xticks(0:0.2:1)
% yticks(0:0.2:1)
% zticks(0:0.2:1)
% 
% xlabel('Pink template')
% ylabel('Blue template')
% zlabel('Brown template')
% ylim([0 1])
% xlim([0 1])
% zlim([0 1])
% grid on
% 
% view([-0.5 3 1])
% 
% % saveas(f,'projection_dec','epsc')
% 
% 
% 
%
function p = z_test_function_bootstrap(dist,null)

m = mean(dist);
s = std(dist);
z = (m-null)/s;
p=1-normcdf(z);

end