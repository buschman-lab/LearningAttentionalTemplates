clear all;
fsroot='/Volumes/buschman';
event='target';

task='Learning_Attentional_Templates';
% subtask='exploreexploit/Reset_RW_model';
subtask='exploreexploit/Restricted_FEF/Reset_RW_model';

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


subsubtask_classifier=sprintf('Classifier_%s_%d_to_%s_%d',event_list{1},window_start_list(1),event_list{end},window_start_list(end));

dirstem = fullfile('Projects',task,'Analysis','Electrophy_analysis',subtask,subsubtask_classifier);
data_path_clasifier = fullfile(fsroot,dirstem);

color_for_ROI=[0 0.447 0.741; 0.85 0.325 0.098; 0.929 0.694 0.125]; %LIP FEF PFC

color_for_ROI_2(1,:,1)=color_for_ROI(1,:)+0.25;
color_for_ROI_2(1,:,2)=color_for_ROI(1,:);

color_for_ROI_2(2,:,1)=color_for_ROI(2,:)+0.15;
color_for_ROI_2(2,:,2)=color_for_ROI(2,:)-0.09;

color_for_ROI_2(3,:,1)=color_for_ROI(3,:)+0.07;
color_for_ROI_2(3,:,2)=color_for_ROI(3,:)-0.12;


load('colors')

color_template(1,:)=colors(1,:);
color_template(2,:)=colors(14,:);
color_template(3,:)=colors(27,:);

this_combination=2;

LOAD=1;
% if this_combination==1
%     save_name=sprintf('Pseudo_pop_peak_belief_prog_combined_time_results_4bins_%s_600ms_proj','PFC');
% elseif this_combination==2
%     save_name=sprintf('Pseudo_pop_peak_belief_prog_combined_time_results_4bins_%s_900ms_proj','PFC');
% end

%%

if LOAD==1

for n_tt=1:100
    
    ROI='LIP';
    if this_combination==2
        results_name=sprintf('Pseudo_pop_peak_belief_prog_combined_time_results_4bins_%s_900ms_%d',ROI,n_tt);
        results_NN_name=sprintf('Pseudo_pop_peak_belief_prog_combined_time_results_4bins_%s_900ms_NN_%d',ROI,n_tt);
%         results_NN_each_prog_name=sprintf('Pseudo_pop_peak_belief_each_prog_combined_time_results_4bins_%s_900ms_NN_%d',ROI,n_tt);
%         results_proj_name=sprintf('Pseudo_pop_peak_belief_prog_combined_time_results_4bins_%s_900ms_%d_proj',ROI,n_tt);
    end
    
    if exist(sprintf('%s.mat',fullfile(data_path_clasifier,results_name)))==2
        
        load(fullfile(data_path_clasifier,results_name),'Classification_correct')
        load(fullfile(data_path_clasifier,results_NN_name),'Classification_net_correct', 'Classification_net_correct_across_time')
        
        LIP.Classification_net_correct(:,:,n_tt)=Classification_net_correct;
        LIP.Classification_net_correct_across_time(:,:,:,n_tt)=Classification_net_correct_across_time;
        LIP.Classification_correct(:,:,n_tt)=Classification_correct;
        
        clear Classification*_correct*
        
%         load(fullfile(data_path_clasifier,results_NN_each_prog_name),'Classification_net_correct')
%         
%         LIP.Classification_net_correct_each_prog(:,:,n_tt)=Classification_net_correct;
%         
%         load(fullfile(data_path_clasifier,results_proj_name),'Classifier_proba_val');
%         LIP.Classifier_projection(:,:,:,n_tt) = Classifier_proba_val;
%         
%         clear Classification*_correct Classifier_proba_val
        
    else
        sprintf('%s %d missing',ROI,n_tt)
        
    end
    
    %FEF
    ROI='FEF';
    
    if this_combination==2
        results_name=sprintf('Pseudo_pop_peak_belief_prog_combined_time_results_4bins_%s_900ms_%d',ROI,n_tt);
        results_NN_name=sprintf('Pseudo_pop_peak_belief_prog_combined_time_results_4bins_%s_900ms_NN_%d',ROI,n_tt);
%         results_NN_each_prog_name=sprintf('Pseudo_pop_peak_belief_each_prog_combined_time_results_4bins_%s_900ms_NN_%d',ROI,n_tt);
%         results_proj_name=sprintf('Pseudo_pop_peak_belief_prog_combined_time_results_4bins_%s_900ms_%d_proj',ROI,n_tt);
    end
    
    if exist(sprintf('%s.mat',fullfile(data_path_clasifier,results_name)))==2
        
        load(fullfile(data_path_clasifier,results_name),'Classification_correct')
        load(fullfile(data_path_clasifier,results_NN_name),'Classification_net_correct', 'Classification_net_correct_across_time')
        
        FEF.Classification_net_correct(:,:,n_tt)=Classification_net_correct;
        FEF.Classification_net_correct_across_time(:,:,:,n_tt)=Classification_net_correct_across_time;
        FEF.Classification_correct(:,:,n_tt)=Classification_correct;
        
        clear Classification*_correct*
        
%         load(fullfile(data_path_clasifier,results_NN_each_prog_name),'Classification_net_correct')
%         
%         FEF.Classification_net_correct_each_prog(:,:,n_tt)=Classification_net_correct;
%         
%         load(fullfile(data_path_clasifier,results_proj_name),'Classifier_proba_val');
%         FEF.Classifier_projection(:,:,:,n_tt) = Classifier_proba_val;
%         
%         clear Classification*_correct Classifier_proba_val
        
    else
        sprintf('%s %d missing',ROI,n_tt)
        
    end
    
    %PFC
    ROI='PFC';
    if this_combination==2
        results_name=sprintf('Pseudo_pop_peak_belief_prog_combined_time_results_4bins_%s_900ms_%d',ROI,n_tt);
        results_NN_name=sprintf('Pseudo_pop_peak_belief_prog_combined_time_results_4bins_%s_900ms_NN_%d',ROI,n_tt);
%         results_NN_each_prog_name=sprintf('Pseudo_pop_peak_belief_each_prog_combined_time_results_4bins_%s_900ms_NN_%d',ROI,n_tt);
%         results_proj_name=sprintf('Pseudo_pop_peak_belief_prog_combined_time_results_4bins_%s_900ms_%d_proj',ROI,n_tt);
    end
    
    if exist(sprintf('%s.mat',fullfile(data_path_clasifier,results_name)))==2
        
        load(fullfile(data_path_clasifier,results_name),'Classification_correct')
        load(fullfile(data_path_clasifier,results_NN_name),'Classification_net_correct', 'Classification_net_correct_across_time')
        
        PFC.Classification_net_correct(:,:,n_tt)=Classification_net_correct;
        PFC.Classification_net_correct_across_time(:,:,:,n_tt)=Classification_net_correct_across_time;
        PFC.Classification_correct(:,:,n_tt)=Classification_correct;
        
        clear Classification*_correct*
        
%         load(fullfile(data_path_clasifier,results_NN_each_prog_name),'Classification_net_correct')
%         
%         PFC.Classification_net_correct_each_prog(:,:,n_tt)=Classification_net_correct;
%         
%         load(fullfile(data_path_clasifier,results_proj_name),'Classifier_proba_val');
%         PFC.Classifier_projection(:,:,:,n_tt) = Classifier_proba_val;
%         
%         clear Classification*_correct Classifier_proba_val
        
    else
        sprintf('%s %d missing',ROI,n_tt)
        
    end
    
    
end

% LIP.Mean_classifier_projection=mean(LIP.Classifier_projection,1);
% FEF.Mean_classifier_projection=mean(FEF.Classifier_projection,1);
% PFC.Mean_classifier_projection=mean(PFC.Classifier_projection,1);

LIP_mean_classification(:,:)=mean(LIP.Classification_correct,2);
FEF_mean_classification(:,:)=mean(FEF.Classification_correct,2);
PFC_mean_classification(:,:)=mean(PFC.Classification_correct,2);

LIP_mean_classification_net(:,:)=mean(LIP.Classification_net_correct,2);
FEF_mean_classification_net(:,:)=mean(FEF.Classification_net_correct,2);
PFC_mean_classification_net(:,:)=mean(PFC.Classification_net_correct,2);

for tp=1:length(window_start_list)
    for n_tt=1:100
        LIP_across_time(tp,n_tt)=mean(mean(LIP.Classification_net_correct_across_time(:,:,tp,n_tt),1),2);
        FEF_across_time(tp,n_tt)=mean(mean(FEF.Classification_net_correct_across_time(:,:,tp,n_tt),1),2);
        PFC_across_time(tp,n_tt)=mean(mean(PFC.Classification_net_correct_across_time(:,:,tp,n_tt),1),2);
    end
end

% LIP_mean_classification_net_each_prog(:,:)=mean(LIP.Classification_net_correct_each_prog,2);
% FEF_mean_classification_net_each_prog(:,:)=mean(FEF.Classification_net_correct_each_prog,2);
% PFC_mean_classification_net_each_prog(:,:)=mean(PFC.Classification_net_correct_each_prog,2);

% save(fullfile(data_path_clasifier,save_name));

else
    
    load(fullfile(data_path_clasifier,save_name));
    
    
end

%%
figure
subplot(2,3,1)
vs = violinplot(LIP_mean_classification',{'1/3','2/3','3/3'},'ViolinColor',color_for_ROI(1,:),'ShowData',false);
hold on
yline(1/4,'--',{'Chance level'},'Color',[0.5 0.5 0.5],'LabelVerticalAlignment','bottom','LabelHorizontalAlignment','center')

box off
title('LIP')
ylabel('Expected template classification accuracy')
xlabel('Progression in block')
ylim([0 0.8])

subplot(2,3,2)
vs = violinplot(FEF_mean_classification',{'1/3','2/3','3/3'},'ViolinColor',color_for_ROI(2,:),'ShowData',false);
hold on
yline(1/4,'--',{'Chance level'},'Color',[0.5 0.5 0.5],'LabelVerticalAlignment','bottom','LabelHorizontalAlignment','center')
box off
title('FEF')
ylabel('Expected template classification accuracy')
xlabel('Progression in block')
ylim([0 0.8])

subplot(2,3,3)
vs = violinplot(PFC_mean_classification',{'1/3','2/3','3/3'},'ViolinColor',color_for_ROI(3,:),'ShowData',false);
hold on
yline(1/4,'--',{'Chance level'},'Color',[0.5 0.5 0.5],'LabelVerticalAlignment','bottom','LabelHorizontalAlignment','center')
box off
title('PFC')
ylabel('Expected template classification accuracy')
xlabel('Progression in block')
ylim([0 0.8])

subplot(2,3,4)
vs = violinplot(LIP_mean_classification_net',{'1/3','2/3','3/3'},'ViolinColor',color_for_ROI(1,:),'ShowData',false);
hold on
yline(1/4,'--',{'Chance level'},'Color',[0.5 0.5 0.5],'LabelVerticalAlignment','bottom','LabelHorizontalAlignment','center')
box off
title('LIP')
ylabel('Expected template classification accuracy NN')
xlabel('Progression in block')
ylim([0 0.8])

subplot(2,3,5)
vs = violinplot(FEF_mean_classification_net',{'1/3','2/3','3/3'},'ViolinColor',color_for_ROI(2,:),'ShowData',false);
hold on
yline(1/4,'--',{'Chance level'},'Color',[0.5 0.5 0.5],'LabelVerticalAlignment','bottom','LabelHorizontalAlignment','center')
box off
title('FEF')
ylabel('Expected template classification accuracy NN')
xlabel('Progression in block')
ylim([0 0.8])

subplot(2,3,6)
vs = violinplot(PFC_mean_classification_net',{'1/3','2/3','3/3'},'ViolinColor',color_for_ROI(3,:),'ShowData',false);
hold on
yline(1/4,'--',{'Chance level'},'Color',[0.5 0.5 0.5],'LabelVerticalAlignment','bottom','LabelHorizontalAlignment','center')
box off
title('PFC')
ylabel('Expected template classification accuracy NN')
xlabel('Progression in block')
ylim([0 0.8])


%%
figure
subplot(3,1,1)
vs = violinplot(LIP_mean_classification_net',{'1/3','2/3','3/3'},'ViolinColor',color_for_ROI(1,:),'ShowData',false);
hold on
yline(1/4,'--',{'Chance level'},'Color',[0.5 0.5 0.5],'LabelVerticalAlignment','middle','LabelHorizontalAlignment','center')
box off
title('LIP','FontSize',12)
ylabel(sprintf('Expected template\nclassification accuracy'),'FontSize',12)
xlabel('Progression in block','FontSize',12)
ylim([0 0.8])

subplot(3,1,2)
vs = violinplot(FEF_mean_classification_net',{'1/3','2/3','3/3'},'ViolinColor',color_for_ROI(2,:),'ShowData',false);
hold on
yline(1/4,'--',{'Chance level'},'Color',[0.5 0.5 0.5],'LabelVerticalAlignment','middle','LabelHorizontalAlignment','center')
box off
title('FEF','FontSize',12)
ylabel(sprintf('Expected template\nclassification accuracy'),'FontSize',12)
xlabel('Progression in block','FontSize',12)
ylim([0 0.8])

subplot(3,1,3)
vs = violinplot(PFC_mean_classification_net',{'1/3','2/3','3/3'},'ViolinColor',color_for_ROI(3,:),'ShowData',false);
hold on
yline(1/4,'--',{'Chance level'},'Color',[0.5 0.5 0.5],'LabelVerticalAlignment','middle','LabelHorizontalAlignment','center')
box off
title('PFC','FontSize',12)
ylabel(sprintf('Expected template\nclassification accuracy'),'FontSize',12)
xlabel('Progression in block','FontSize',12)
ylim([0 0.8])

%%
figure
subplot(1,3,1)
vs = violinplot(LIP_mean_classification_net',{'1/3','2/3','3/3'},'ViolinColor',color_for_ROI(1,:),'ShowData',false,'ViolinAlpha',0.5);
hold on
yline(1/4,'--',{'Chance level'},'Color',[0.5 0.5 0.5],'LabelVerticalAlignment','middle','LabelHorizontalAlignment','center')
box off
title('LIP','FontSize',12)
ylabel(sprintf('Expected template\nclassification accuracy'),'FontSize',12)
xlabel('Progression in block','FontSize',12)
ylim([0 1])

subplot(1,3,2)
vs = violinplot(FEF_mean_classification_net',{'1/3','2/3','3/3'},'ViolinColor',color_for_ROI(2,:),'ShowData',false,'ViolinAlpha',0.5);
hold on
yline(1/4,'--',{'Chance level'},'Color',[0.5 0.5 0.5],'LabelVerticalAlignment','middle','LabelHorizontalAlignment','center')
box off
title('FEF','FontSize',12)
ylabel(sprintf('Expected template\nclassification accuracy'),'FontSize',12)
xlabel('Progression in block','FontSize',12)
ylim([0 1])

subplot(1,3,3)
vs = violinplot(PFC_mean_classification_net',{'1/3','2/3','3/3'},'ViolinColor',color_for_ROI(3,:),'ShowData',false,'ViolinAlpha',0.5);
hold on
yline(1/4,'--',{'Chance level'},'Color',[0.5 0.5 0.5],'LabelVerticalAlignment','middle','LabelHorizontalAlignment','center')
box off
title('PFC','FontSize',12)
ylabel(sprintf('Expected template\nclassification accuracy'),'FontSize',12)
xlabel('Progression in block','FontSize',12)
ylim([0 1])


%% Figure 2D

% figure
% subplot(3,1,1)
% vs = violinplot(LIP_mean_classification_net',{'1/3','2/3','3/3'},'ViolinColor',color_for_ROI_2(1,:,2),'ShowData',false,'ViolinAlpha',0.5,'HalfViolin','left','BoxColor',[0 0 0],'MedianColor',[0.75 0.75 0.75]);
% hold on
% vs = violinplot(LIP_mean_classification_net_each_prog',{'1/3','2/3','3/3'},'ViolinColor',color_for_ROI_2(1,:,1),'ShowData',false,'ViolinAlpha',0.5,'HalfViolin','right','BoxColor',[0.5 0.5 0.5],'MedianColor',[1 1 1]);
% 
% yline(1/4,'--',{'Chance level'},'Color',[0.5 0.5 0.5],'LabelVerticalAlignment','middle','LabelHorizontalAlignment','center')
% box off
% % title('LIP','FontSize',12)
% ylabel(sprintf('Expected template\nclassification accuracy'),'FontSize',12)
% xlabel('Progression in block','FontSize',12)
% ylim([0 0.8])
% 
% subplot(3,1,2)
% vs = violinplot(FEF_mean_classification_net',{'1/3','2/3','3/3'},'ViolinColor',color_for_ROI_2(2,:,2),'ShowData',false,'ViolinAlpha',0.5,'HalfViolin','left','BoxColor',[0 0 0],'MedianColor',[0.75 0.75 0.75]);
% hold on
% vs = violinplot(FEF_mean_classification_net_each_prog',{'1/3','2/3','3/3'},'ViolinColor',color_for_ROI_2(2,:,1),'ShowData',false,'ViolinAlpha',0.5,'HalfViolin','right','BoxColor',[0.5 0.5 0.5],'MedianColor',[1 1 1]);
% 
% yline(1/4,'--',{'Chance level'},'Color',[0.5 0.5 0.5],'LabelVerticalAlignment','middle','LabelHorizontalAlignment','center')
% box off
% % title('FEF','FontSize',12)
% ylabel(sprintf('Expected template\nclassification accuracy'),'FontSize',12)
% xlabel('Progression in block','FontSize',12)
% ylim([0 0.8])
% 
% subplot(3,1,3)
% vs = violinplot(PFC_mean_classification_net',{'1/3','2/3','3/3'},'ViolinColor',color_for_ROI_2(3,:,2),'ShowData',false,'ViolinAlpha',0.5,'HalfViolin','left','BoxColor',[0 0 0],'MedianColor',[0.75 0.75 0.75]);
% hold on
% vs = violinplot(PFC_mean_classification_net_each_prog',{'1/3','2/3','3/3'},'ViolinColor',color_for_ROI_2(3,:,1),'ShowData',false,'ViolinAlpha',0.5,'HalfViolin','right','BoxColor',[0.5 0.5 0.5],'MedianColor',[1 1 1]);
% 
% yline(1/4,'--',{'Chance level'},'Color',[0.5 0.5 0.5],'LabelVerticalAlignment','middle','LabelHorizontalAlignment','center')
% box off
% % title('PFC','FontSize',12)
% ylabel(sprintf('Expected template\nclassification accuracy'),'FontSize',12)
% xlabel('Progression in block','FontSize',12)
% ylim([0 0.8])

%%
for k=1:3
    p_LIP(k,1)=z_test_function_bootstrap(LIP_mean_classification_net(k,:),1/4);
%     p_LIP(k,2)=z_test_function_bootstrap(LIP_mean_classification_net_each_prog(k,:),1/3);
    p_FEF(k,1)=z_test_function_bootstrap(FEF_mean_classification_net(k,:),1/4);
%     p_FEF(k,2)=z_test_function_bootstrap(FEF_mean_classification_net_each_prog(k,:),1/3);
    p_PFC(k,1)=z_test_function_bootstrap(PFC_mean_classification_net(k,:),1/4);
%     p_PFC(k,2)=z_test_function_bootstrap(PFC_mean_classification_net_each_prog(k,:),1/3);
end


for tp=1:length(window_start_list)
    p_LIP_across_time(tp)=z_test_function_bootstrap(LIP_across_time(tp,:),1/4);
    p_FEF_across_time(tp)=z_test_function_bootstrap(FEF_across_time(tp,:),1/4);
    p_PFC_across_time(tp)=z_test_function_bootstrap(PFC_across_time(tp,:),1/4);
end



%% 3d ellipses
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
% %% Figure 2C
% 
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
% ylim([0 0.8])
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
% ylim([0 0.8])
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
% ylim([0 0.8])
% xlim([0 1])
% zlim([0 1])
% grid on
% 
% view([-0.5 3 1])
% 
% % saveas(f,'projection_dec','epsc')
% 
% %%
% 
% figure
% subplot(3,1,1)
% hold on
% shadedErrorBar(window_start_list(8:end),mean(LIP_across_time(8:end,:),2)',[prctile(LIP_across_time(8:end,:),95,2),prctile(LIP_across_time(8:end,:),5,2)]',{'color',color_for_ROI(1,:),'LineWidth',2},2)
% plot_significance_level(window_start_list(8:end),p_LIP_across_time(8:end),0.85,color_for_ROI(1,:),[0.01,0.05/(length(window_start_list)-7), 0.01/(length(window_start_list)-7)])
% 
% yl=yline(1/4,'--',{'Chance level'},'Color',[0 0 0],'LineWidth',1);
% yl.LabelVerticalAlignment = 'middle';
% yl.LabelHorizontalAlignment = 'center';
% yl.FontSize=12;
% xline(-300,'-','Color',[0.5 0.5 0.5],'LineWidth',1)
% xl=xline(-400,'-',{'Fixation'},'Color',[0.5 0.5 0.5],'LineWidth',1);
% xl.LabelVerticalAlignment = 'bottom';
% xl.LabelHorizontalAlignment = 'left';
% xl.FontSize=12;
% 
% xline(300,'-','Color',[0.25 0.75 0.25],'LineWidth',1)
% xl=xline(200,'-',{'Response'},'Color',[0.25 0.75 0.25],'LineWidth',1);
% xl.LabelVerticalAlignment = 'bottom';
% xl.LabelHorizontalAlignment = 'left';
% xl.FontSize=12;
% ylabel('Expected template decoding accuracy','FontSize',14)
% xlabel(sprintf('Time to %s',event_list{end}'))
% ax=gca;
% ax.XAxis.FontSize=12;
% ax.YAxis.FontSize=12;
% box off
% ylim([0 0.8])
% xlim([-600 350])
% 
% subplot(3,1,2)
% hold on
% shadedErrorBar(window_start_list(8:end),mean(FEF_across_time(8:end,:),2)',[prctile(FEF_across_time(8:end,:),95,2),prctile(FEF_across_time(8:end,:),5,2)]',{'color',color_for_ROI(2,:),'LineWidth',2},2)
% plot_significance_level(window_start_list(8:end),p_FEF_across_time(8:end),0.85,color_for_ROI(2,:),[0.01,0.05/(length(window_start_list)-7), 0.01/(length(window_start_list)-7)])
% 
% yl=yline(1/4,'--',{'Chance level'},'Color',[0 0 0],'LineWidth',1);
% yl.LabelVerticalAlignment = 'middle';
% yl.LabelHorizontalAlignment = 'center';
% yl.FontSize=12;
% xline(-300,'-','Color',[0.5 0.5 0.5],'LineWidth',1)
% xl=xline(-400,'-',{'Fixation'},'Color',[0.5 0.5 0.5],'LineWidth',1);
% xl.LabelVerticalAlignment = 'bottom';
% xl.LabelHorizontalAlignment = 'left';
% xl.FontSize=12;
% 
% xline(300,'-','Color',[0.25 0.75 0.25],'LineWidth',1)
% xl=xline(200,'-',{'Response'},'Color',[0.25 0.75 0.25],'LineWidth',1);
% xl.LabelVerticalAlignment = 'bottom';
% xl.LabelHorizontalAlignment = 'left';
% xl.FontSize=12;
% ylabel('Expected template decoding accuracy','FontSize',14)
% xlabel(sprintf('Time to %s',event_list{end}'))
% ax=gca;
% ax.XAxis.FontSize=12;
% ax.YAxis.FontSize=12;
% box off
% ylim([0 0.8])
% xlim([-600 350])
% 
% subplot(3,1,3)
% hold on
% shadedErrorBar(window_start_list(8:end),mean(PFC_across_time(8:end,:),2)',[prctile(PFC_across_time(8:end,:),95,2),prctile(PFC_across_time(8:end,:),5,2)]',{'color',color_for_ROI(3,:),'LineWidth',2},2)
% plot_significance_level(window_start_list(8:end),p_PFC_across_time(8:end),0.85,color_for_ROI(3,:),[0.01,0.05/(length(window_start_list)-7), 0.01/(length(window_start_list)-7)])
% 
% yl=yline(1/4,'--',{'Chance level'},'Color',[0 0 0],'LineWidth',1);
% yl.LabelVerticalAlignment = 'middle';
% yl.LabelHorizontalAlignment = 'center';
% yl.FontSize=12;
% xline(-300,'-','Color',[0.5 0.5 0.5],'LineWidth',1)
% xl=xline(-400,'-',{'Fixation'},'Color',[0.5 0.5 0.5],'LineWidth',1);
% xl.LabelVerticalAlignment = 'bottom';
% xl.LabelHorizontalAlignment = 'left';
% xl.FontSize=12;
% 
% xline(300,'-','Color',[0.25 0.75 0.25],'LineWidth',1)
% xl=xline(200,'-',{'Response'},'Color',[0.25 0.75 0.25],'LineWidth',1);
% xl.LabelVerticalAlignment = 'bottom';
% xl.LabelHorizontalAlignment = 'left';
% xl.FontSize=12;
% ylabel('Expected template decoding accuracy','FontSize',14)
% xlabel(sprintf('Time to %s',event_list{end}'))
% ax=gca;
% ax.XAxis.FontSize=12;
% ax.YAxis.FontSize=12;
% box off
% ylim([0 0.8])
% xlim([-600 350])
% 


%%
p = z_test_function_bootstrap(LIP_across_time(20,:)-LIP_across_time(8,:),0);


%% functions



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
end
end


function p = z_test_function_bootstrap(dist,null)

m = mean(dist);
s = std(dist);
z = (m-null)/s;
p=1-normcdf(z);

end




