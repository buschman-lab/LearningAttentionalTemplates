clear all;

fsroot='/Volumes/buschman';


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

N_boot=100;

dirstem = fullfile('Projects',task,'Analysis','Electrophy_analysis',subtask,subsubtask_classifier);
data_path_clasifier = fullfile(fsroot,dirstem);

color_for_ROI=[0 0.447 0.741; 0.85 0.325 0.098; 0.929 0.694 0.125]; %LIP FEF PFC

load('colors')

color_belief(1,:)=colors(1,:);
color_belief(2,:)=colors(14,:);
color_belief(3,:)=colors(27,:);


%%
%LIP
ROI='LIP';

for n_tt=1:N_boot
    
    for this_time=1:length(window_start_list)
        
        results_name=sprintf('Pseudo_pop_peak_belief_prog_across_time_results_more_trials_NN_%s_%d_%d',ROI,this_time,n_tt);
        
        if exist(sprintf('%s.mat',fullfile(data_path_clasifier,results_name)))==2
            
            load(fullfile(data_path_clasifier,results_name),'Classification_net_correct')
            
            LIP.Classification_net_correct(this_time,:,:,:,n_tt)=Classification_net_correct;
            
            clear Classification_net_correct W
            
        else
            sprintf('%s time %d boot %d is missing',ROI,this_time,n_tt)
        end
        
    end
end

%FEF
ROI='FEF';

for n_tt=1:N_boot
    
    for this_time=1:length(window_start_list)
        
        results_name=sprintf('Pseudo_pop_peak_belief_prog_across_time_results_more_trials_NN_%s_%d_%d',ROI,this_time,n_tt);
        
        if exist(sprintf('%s.mat',fullfile(data_path_clasifier,results_name)))==2
            
            load(fullfile(data_path_clasifier,results_name),'Classification_net_correct')
            
            FEF.Classification_net_correct(this_time,:,:,:,n_tt)=Classification_net_correct;
            
            clear Classification_net_correct W
            
        else
            sprintf('%s time %d boot %d is missing',ROI,this_time,n_tt)
        end
        
    end
end

%PFC
ROI='PFC';

for n_tt=1:N_boot
    
    for this_time=1:length(window_start_list)
        
        results_name=sprintf('Pseudo_pop_peak_belief_prog_across_time_results_more_trials_NN_%s_%d_%d',ROI,this_time,n_tt);
        
        if exist(sprintf('%s.mat',fullfile(data_path_clasifier,results_name)))==2
            
            load(fullfile(data_path_clasifier,results_name),'Classification_net_correct')
            
            PFC.Classification_net_correct(this_time,:,:,:,n_tt)=Classification_net_correct;
            
            clear Classification_net_correct W
            
        else
            sprintf('%s time %d boot %d is missing',ROI,this_time,n_tt)
        end
        
    end
end



%%
for i=1:length(window_start_list)
    for k=1:3
        p_LIP(i,k)=z_test_function_bootstrap(mean(LIP.Classification_net_correct(i,k,:,i,:),3),1/3);
        p_FEF(i,k)=z_test_function_bootstrap(mean(FEF.Classification_net_correct(i,k,:,i,:),3),1/3);
        p_PFC(i,k)=z_test_function_bootstrap(mean(PFC.Classification_net_correct(i,k,:,i,:),3),1/3);
    end
end

%%
for i=1:length(window_start_list)
    p_LIP_early_late(i)=z_test_function_bootstrap(mean(LIP.Classification_net_correct(i,3,:,i,:),3)-(mean(LIP.Classification_net_correct(i,1,:,i,:),3)),0);
    p_FEF_early_late(i)=z_test_function_bootstrap(mean(FEF.Classification_net_correct(i,3,:,i,:),3)-(mean(FEF.Classification_net_correct(i,1,:,i,:),3)),0);
    p_PFC_early_late(i)=z_test_function_bootstrap(mean(PFC.Classification_net_correct(i,3,:,i,:),3)-(mean(PFC.Classification_net_correct(i,1,:,i,:),3)),0);
end

%%
for i=1:length(window_start_list)
    for k=1:3
        LIP_across_time(i,k,:)=mean(LIP.Classification_net_correct(i,k,:,i,:),3);
        FEF_across_time(i,k,:)=mean(FEF.Classification_net_correct(i,k,:,i,:),3);
        PFC_across_time(i,k,:)=mean(PFC.Classification_net_correct(i,k,:,i,:),3);
    end
end

for i=1:length(window_start_list)
    for j=1:length(window_start_list)
        LIP_cross_temporal_decoding(i,j,:)=mean(mean(LIP.Classification_net_correct(i,:,:,j,:),2),3);
        FEF_cross_temporal_decoding(i,j,:)=mean(mean(FEF.Classification_net_correct(i,:,:,j,:),2),3);
        PFC_cross_temporal_decoding(i,j,:)=mean(mean(PFC.Classification_net_correct(i,:,:,j,:),2),3);
    end
end


%%

for i=1:length(window_start_list)
    for j=1:length(window_start_list)
        p_LIP_cross_temporal_decoding(i,j)=z_test_function_bootstrap(mean(mean(LIP.Classification_net_correct(i,:,:,j,:),2),3),1/3);
        p_FEF_cross_temporal_decoding(i,j)=z_test_function_bootstrap(mean(mean(FEF.Classification_net_correct(i,:,:,j,:),2),3),1/3);
        p_PFC_cross_temporal_decoding(i,j)=z_test_function_bootstrap(mean(mean(PFC.Classification_net_correct(i,:,:,j,:),2),3),1/3);
    end
end

%%

figure

    subplot(3,1,1)
    hold on
    imagesc(window_start_list(8:end),window_start_list(8:end),mean(LIP_cross_temporal_decoding(8:end,8:end,:),3),[0.33 0.7])
    colorbar
    colormap hot
    contour(window_start_list(8:end),window_start_list(8:end),(p_LIP_cross_temporal_decoding(8:end,8:end)),[0.05 0.01 0.001],'-k','ShowText','on')
    ax=gca;
    ax.XAxis.FontSize=12;
    ax.YAxis.FontSize=12;
    xlabel('Training time to tagets onset')
    ylabel('Testing time to tagets onset')
    
    title('LIP')
    subplot(3,1,2)
    hold on
    imagesc(window_start_list(8:end),window_start_list(8:end),mean(FEF_cross_temporal_decoding(8:end,8:end,:),3),[0.33 0.7])
    colorbar
    colormap hot
    contour(window_start_list(8:end),window_start_list(8:end),(p_FEF_cross_temporal_decoding(8:end,8:end)),[0.05 0.01 0.001],'-k','ShowText','on')
    
    title('FEF')
    ax=gca;
    ax.XAxis.FontSize=12;
    ax.YAxis.FontSize=12;
    xlabel('Training time to tagets onset')
    ylabel('Testing time to tagets onset')
    
    subplot(3,1,3)
    hold on
    imagesc(window_start_list(8:end),window_start_list(8:end),mean(PFC_cross_temporal_decoding(8:end,8:end,:),3),[0.33 0.7])
    colorbar
    colormap hot
    contour(window_start_list(8:end),window_start_list(8:end),(p_PFC_cross_temporal_decoding(8:end,8:end)),[0.05 0.01 0.001],'-k','ShowText','on')
    
    title('PFC')
    ax=gca;
    ax.XAxis.FontSize=12;
    ax.YAxis.FontSize=12;
    xlabel('Training time to tagets onset')
    ylabel('Testing time to tagets onset')

%%
% for i=1:27
%     for n=1:100
%         mdl=fitglm(1:3,LIP_across_time(i,:,n));
%         T_LIP(i,n)=mdl.Coefficients.tStat(2);
%
%         mdl=fitglm(1:3,FEF_across_time(i,:,n));
%         T_FEF(i,n)=mdl.Coefficients.tStat(2);
%
%         mdl=fitglm(1:3,PFC_across_time(i,:,n));
%         T_PFC(i,n)=mdl.Coefficients.tStat(2);
%     end
% end
%
%
% %%
% figure
%
% subplot(1,3,1)
% errorbar(mean(T_LIP(1:7,:),2),std(T_LIP(1:7,:),0,2),'Color',color_for_ROI(1,:),'LineWidth',2)
%


%%

% figure
% subplot(1,3,1)
% for k=1:3
%     errorbar(window_start_list+10*(k-2),mean(LIP_across_time(:,k,:),3),std(LIP_across_time(:,k,:),0,3),'.','Color',color_for_ROI(1,:),'LineWidth',2)
%     hold on
% end
% yline(1/3,'--k')
% ylabel('Accuracy','FontSize',16)
% ax=gca;
% ax.XAxis.FontSize=12;
% ax.YAxis.FontSize=12;
% box off
%
% subplot(1,3,2)
% for k=1:3
%     errorbar(window_start_list+10*(k-2),mean(FEF_across_time(:,k,:),3),std(FEF_across_time(:,k,:),0,3),'.','Color',color_for_ROI(2,:),'LineWidth',2)
%     hold on
% end
% yline(1/3,'--k')
% ylabel('Accuracy','FontSize',16)
% ax=gca;
% ax.XAxis.FontSize=12;
% ax.YAxis.FontSize=12;
% box off
%
% subplot(1,3,3)
% for k=1:3
%     errorbar(window_start_list+10*(k-2),mean(PFC_across_time(:,k,:),3),std(PFC_across_time(:,k,:),0,3),'.','Color',color_for_ROI(3,:),'LineWidth',2)
%     hold on
% end
% yline(1/3,'--k')
%
% ylabel('Accuracy','FontSize',16)
% ax=gca;
% ax.XAxis.FontSize=12;
% ax.YAxis.FontSize=12;
% box off
%
%
%%
for i=1:3
    color_for_ROI_early(i,:)=rgb2hsv(color_for_ROI(i,:));
    color_for_ROI_early(i,:)=hsv2rgb(color_for_ROI_early(i,1),color_for_ROI_early(i,2),1);
    color_for_ROI_late(i,:)=rgb2hsv(color_for_ROI(i,:));
    color_for_ROI_late(i,:)=hsv2rgb(color_for_ROI_late(i,1),color_for_ROI_late(i,2),0.55);
    
end

%%

figure
subplot(3,2,1)
hold on
shadedErrorBar(window_start_list(1:7),mean(LIP_across_time(1:7,1,:),3)',[prctile(LIP_across_time(1:7,1,:),95,3),prctile(LIP_across_time(1:7,1,:),5,3)]',{'color',color_for_ROI_early(1,:),'LineWidth',2},2)
shadedErrorBar(window_start_list(1:7),mean(LIP_across_time(1:7,2,:),3)',[prctile(LIP_across_time(1:7,2,:),95,3),prctile(LIP_across_time(1:7,2,:),5,3)]',{'color',color_for_ROI(1,:),'LineWidth',2},2)
shadedErrorBar(window_start_list(1:7),mean(LIP_across_time(1:7,3,:),3)',[prctile(LIP_across_time(1:7,3,:),95,3),prctile(LIP_across_time(1:7,3,:),5,3)]',{'color',color_for_ROI_late(1,:),'LineWidth',2},2)
plot_significance_level(window_start_list(1:7),p_LIP(1:7,1),0.85,color_for_ROI_early(1,:),[0.01,0.05/(length(window_start_list)), 0.01/(length(window_start_list))])
plot_significance_level(window_start_list(1:7),p_LIP(1:7,2),0.9,color_for_ROI(1,:),[0.01,0.05/(length(window_start_list)), 0.01/(length(window_start_list))])
plot_significance_level(window_start_list(1:7),p_LIP(1:7,3),0.95,color_for_ROI_late(1,:),[0.01,0.05/(length(window_start_list)), 0.01/(length(window_start_list))])

yl=yline(1/3,'--',{'Chance level'},'Color',[0 0 0],'LineWidth',1);
yl.LabelVerticalAlignment = 'middle';
yl.LabelHorizontalAlignment = 'center';
yl.FontSize=12;

ylabel('Expected template decoding accuracy','FontSize',14)
xlabel(sprintf('Time to end of reward'))
ax=gca;
ax.XAxis.FontSize=12;
ax.YAxis.FontSize=12;
box off
%legend({'Early in block','Late in block'},'FontSize',12,'Location','South')
ylim([0 1])

subplot(3,2,2)
hold on
shadedErrorBar(window_start_list(8:end),mean(LIP_across_time(8:end,1,:),3)',[prctile(LIP_across_time(8:end,1,:),95,3),prctile(LIP_across_time(8:end,1,:),5,3)]',{'color',color_for_ROI_early(1,:),'LineWidth',2},2)
shadedErrorBar(window_start_list(8:end),mean(LIP_across_time(8:end,2,:),3)',[prctile(LIP_across_time(8:end,2,:),95,3),prctile(LIP_across_time(8:end,2,:),5,3)]',{'color',color_for_ROI(1,:),'LineWidth',2},2)
shadedErrorBar(window_start_list(8:end),mean(LIP_across_time(8:end,3,:),3)',[prctile(LIP_across_time(8:end,3,:),95,3),prctile(LIP_across_time(8:end,3,:),5,3)]',{'color',color_for_ROI_late(1,:),'LineWidth',2},2)
plot_significance_level(window_start_list(8:end),p_LIP(8:end,1),0.85,color_for_ROI_early(1,:),[0.01,0.05/(length(window_start_list)), 0.01/(length(window_start_list))])
plot_significance_level(window_start_list(8:end),p_LIP(8:end,2),0.9,color_for_ROI(1,:),[0.01,0.05/(length(window_start_list)), 0.01/(length(window_start_list))])
plot_significance_level(window_start_list(8:end),p_LIP(8:end,3),0.95,color_for_ROI_late(1,:),[0.01,0.05/(length(window_start_list)), 0.01/(length(window_start_list))])

yl=yline(1/3,'--',{'Chance level'},'Color',[0 0 0],'LineWidth',1);
yl.LabelVerticalAlignment = 'middle';
yl.LabelHorizontalAlignment = 'center';
yl.FontSize=12;
xline(-300,'-','Color',[0.5 0.5 0.5],'LineWidth',1)
xl=xline(-400,'-',{'Fixation'},'Color',[0.5 0.5 0.5],'LineWidth',1);
xl.LabelVerticalAlignment = 'bottom';
xl.LabelHorizontalAlignment = 'left';
xl.FontSize=12;

xline(300,'-','Color',[0.25 0.75 0.25],'LineWidth',1)
xl=xline(200,'-',{'Response'},'Color',[0.25 0.75 0.25],'LineWidth',1);
xl.LabelVerticalAlignment = 'bottom';
xl.LabelHorizontalAlignment = 'left';
xl.FontSize=12;
ylabel('Expected template decoding accuracy','FontSize',14)
xlabel(sprintf('Time to %s',event_list{end}'))
ax=gca;
ax.XAxis.FontSize=12;
ax.YAxis.FontSize=12;
box off
ylim([0 1])
xlim([-600 350])

subplot(3,2,3)
hold on
shadedErrorBar(window_start_list(1:7),mean(FEF_across_time(1:7,1,:),3)',[prctile(FEF_across_time(1:7,1,:),95,3),prctile(FEF_across_time(1:7,1,:),5,3)]',{'color',color_for_ROI_early(2,:),'LineWidth',2},2)
shadedErrorBar(window_start_list(1:7),mean(FEF_across_time(1:7,2,:),3)',[prctile(FEF_across_time(1:7,2,:),95,3),prctile(FEF_across_time(1:7,2,:),5,3)]',{'color',color_for_ROI(2,:),'LineWidth',2},2)
shadedErrorBar(window_start_list(1:7),mean(FEF_across_time(1:7,3,:),3)',[prctile(FEF_across_time(1:7,3,:),95,3),prctile(FEF_across_time(1:7,3,:),5,3)]',{'color',color_for_ROI_late(2,:),'LineWidth',2},2)
plot_significance_level(window_start_list(1:7),p_FEF(1:7,1),0.85,color_for_ROI_early(2,:),[0.01,0.05/(length(window_start_list)), 0.01/(length(window_start_list))])
plot_significance_level(window_start_list(1:7),p_FEF(1:7,2),0.9,color_for_ROI(2,:),[0.01,0.05/(length(window_start_list)), 0.01/(length(window_start_list))])
plot_significance_level(window_start_list(1:7),p_FEF(1:7,3),0.95,color_for_ROI_late(2,:),[0.01,0.05/(length(window_start_list)), 0.01/(length(window_start_list))])

yl=yline(1/3,'--',{'Chance level'},'Color',[0 0 0],'LineWidth',1);
yl.LabelVerticalAlignment = 'middle';
yl.LabelHorizontalAlignment = 'center';
yl.FontSize=12;

ylabel('Expected template decoding accuracy','FontSize',14)
xlabel(sprintf('Time to end of reward'))
ax=gca;
ax.XAxis.FontSize=12;
ax.YAxis.FontSize=12;
box off
%legend({'Early in block','Late in block'},'FontSize',12,'Location','South')
ylim([0 1])

subplot(3,2,4)
hold on
shadedErrorBar(window_start_list(8:end),mean(FEF_across_time(8:end,1,:),3)',[prctile(FEF_across_time(8:end,1,:),95,3),prctile(FEF_across_time(8:end,1,:),5,3)]',{'color',color_for_ROI_early(2,:),'LineWidth',2},2)
shadedErrorBar(window_start_list(8:end),mean(FEF_across_time(8:end,2,:),3)',[prctile(FEF_across_time(8:end,2,:),95,3),prctile(FEF_across_time(8:end,2,:),5,3)]',{'color',color_for_ROI(2,:),'LineWidth',2},2)
shadedErrorBar(window_start_list(8:end),mean(FEF_across_time(8:end,3,:),3)',[prctile(FEF_across_time(8:end,3,:),95,3),prctile(FEF_across_time(8:end,3,:),5,3)]',{'color',color_for_ROI_late(2,:),'LineWidth',2},2)
plot_significance_level(window_start_list(8:end),p_FEF(8:end,1),0.85,color_for_ROI_early(2,:),[0.01,0.05/(length(window_start_list)), 0.01/(length(window_start_list))])
plot_significance_level(window_start_list(8:end),p_FEF(8:end,2),0.9,color_for_ROI(2,:),[0.01,0.05/(length(window_start_list)), 0.01/(length(window_start_list))])
plot_significance_level(window_start_list(8:end),p_FEF(8:end,3),0.95,color_for_ROI_late(2,:),[0.01,0.05/(length(window_start_list)), 0.01/(length(window_start_list))])
yl=yline(1/3,'--',{'Chance level'},'Color',[0 0 0],'LineWidth',1);
yl.LabelVerticalAlignment = 'middle';
yl.LabelHorizontalAlignment = 'center';
yl.FontSize=12;

xline(-300,'-','Color',[0.5 0.5 0.5],'LineWidth',1)
xl=xline(-400,'-',{'Fixation'},'Color',[0.5 0.5 0.5],'LineWidth',1);
xl.LabelVerticalAlignment = 'bottom';
xl.LabelHorizontalAlignment = 'left';
xl.FontSize=12;

xline(300,'-','Color',[0.25 0.75 0.25],'LineWidth',1)
xl=xline(200,'-',{'Response'},'Color',[0.25 0.75 0.25],'LineWidth',1);
xl.LabelVerticalAlignment = 'bottom';
xl.LabelHorizontalAlignment = 'left';
xl.FontSize=12;
ylabel('Expected template decoding accuracy','FontSize',14)
xlabel(sprintf('Time to %s',event_list{end}'))
ax=gca;
ax.XAxis.FontSize=12;
ax.YAxis.FontSize=12;
box off
% %legend({'Early in block','Late in block'},'FontSize',12,'Location','South')
ylim([0 1])

xlim([-600 350])

subplot(3,2,5)
hold on;
shadedErrorBar(window_start_list(1:7),mean(PFC_across_time(1:7,1,:),3)',[prctile(PFC_across_time(1:7,1,:),95,3),prctile(PFC_across_time(1:7,1,:),5,3)]',{'color',color_for_ROI_early(3,:),'LineWidth',2},2)
shadedErrorBar(window_start_list(1:7),mean(PFC_across_time(1:7,2,:),3)',[prctile(PFC_across_time(1:7,2,:),95,3),prctile(PFC_across_time(1:7,2,:),5,3)]',{'color',color_for_ROI(3,:),'LineWidth',2},2)
shadedErrorBar(window_start_list(1:7),mean(PFC_across_time(1:7,3,:),3)',[prctile(PFC_across_time(1:7,3,:),95,3),prctile(PFC_across_time(1:7,3,:),5,3)]',{'color',color_for_ROI_late(3,:),'LineWidth',2},2)
plot_significance_level(window_start_list(1:7),p_PFC(1:7,1),0.85,color_for_ROI_early(3,:),[0.01,0.05/(length(window_start_list)), 0.01/(length(window_start_list))])
plot_significance_level(window_start_list(1:7),p_PFC(1:7,2),0.9,color_for_ROI(3,:),[0.01,0.05/(length(window_start_list)), 0.01/(length(window_start_list))])
plot_significance_level(window_start_list(1:7),p_PFC(1:7,3),0.95,color_for_ROI_late(3,:),[0.01,0.05/(length(window_start_list)), 0.01/(length(window_start_list))])

yl=yline(1/3,'--',{'Chance level'},'Color',[0 0 0],'LineWidth',1);
yl.LabelVerticalAlignment = 'middle';
yl.LabelHorizontalAlignment = 'center';
yl.FontSize=12;
ylabel('Expected template decoding accuracy','FontSize',14)
xlabel(sprintf('Time to end of reward'))
ax=gca;
ax.XAxis.FontSize=12;
ax.YAxis.FontSize=12;
box off
%legend({'Early in block','Late in block'},'FontSize',12,'Location','South')
ylim([0 1])

subplot(3,2,6)
hold on
shadedErrorBar(window_start_list(8:end),mean(PFC_across_time(8:end,1,:),3)',[prctile(PFC_across_time(8:end,1,:),95,3),prctile(PFC_across_time(8:end,1,:),5,3)]',{'color',color_for_ROI_early(3,:),'LineWidth',2},2)
shadedErrorBar(window_start_list(8:end),mean(PFC_across_time(8:end,2,:),3)',[prctile(PFC_across_time(8:end,2,:),95,3),prctile(PFC_across_time(8:end,2,:),5,3)]',{'color',color_for_ROI(3,:),'LineWidth',2},2)
shadedErrorBar(window_start_list(8:end),mean(PFC_across_time(8:end,3,:),3)',[prctile(PFC_across_time(8:end,3,:),95,3),prctile(PFC_across_time(8:end,3,:),5,3)]',{'color',color_for_ROI_late(3,:),'LineWidth',2},2)
plot_significance_level(window_start_list(8:end),p_PFC(8:end,1),0.85,color_for_ROI_early(3,:),[0.01,0.05/(length(window_start_list)), 0.01/(length(window_start_list))])
plot_significance_level(window_start_list(8:end),p_PFC(8:end,2),0.9,color_for_ROI(3,:),[0.01,0.05/(length(window_start_list)), 0.01/(length(window_start_list))])
plot_significance_level(window_start_list(8:end),p_PFC(8:end,3),0.95,color_for_ROI_late(3,:),[0.01,0.05/(length(window_start_list)), 0.01/(length(window_start_list))])

yl=yline(1/3,'--',{'Chance level'},'Color',[0 0 0],'LineWidth',1);
yl.LabelVerticalAlignment = 'middle';
yl.LabelHorizontalAlignment = 'center';
yl.FontSize=12;

xline(-300,'-','Color',[0.5 0.5 0.5],'LineWidth',1)
xl=xline(-400,'-',{'Fixation'},'Color',[0.5 0.5 0.5],'LineWidth',1);
xl.LabelVerticalAlignment = 'bottom';
xl.LabelHorizontalAlignment = 'left';
xl.FontSize=12;

xline(300,'-','Color',[0.25 0.75 0.25],'LineWidth',1)
xl=xline(200,'-',{'Response'},'Color',[0.25 0.75 0.25],'LineWidth',1);
xl.LabelVerticalAlignment = 'bottom';
xl.LabelHorizontalAlignment = 'left';
xl.FontSize=12;
ylabel('Expected template decoding accuracy','FontSize',14)
xlabel(sprintf('Time to %s',event_list{end}'))
ax=gca;
ax.XAxis.FontSize=12;
ax.YAxis.FontSize=12;
box off
ylim([0 1])
xlim([-600 350])


%%
% figure
% subplot(3,3,1)
% plot(window_start_list(8:end),mean(PFC_across_time(8:end,1,:),3),'Color',color_for_ROI_early(3,:),'LineWidth',2)
% hold on
% plot(window_start_list(8:end),mean(LIP_across_time(8:end,1,:),3),'Color',color_for_ROI_early(1,:),'LineWidth',2)
% hold on
% shadedErrorBar(window_start_list(8:end),mean(PFC_across_time(8:end,1,:),3),std(PFC_across_time(8:end,1,:),0,3),{'color',color_for_ROI_early(3,:),'LineWidth',2},2)
% hold on
% shadedErrorBar(window_start_list(8:end),mean(LIP_across_time(8:end,1,:),3),std(LIP_across_time(8:end,1,:),0,3),{'color',color_for_ROI_early(1,:),'LineWidth',2},2)
% hold on
% % for i=8:length(window_start_list)
% %     if p_LIP_PFC(i,1)<0.05/length(window_start_list(8:end)) %Bonferoni
% %         plot(window_start_list(i),0.8,'*','Color',(color_for_ROI_early(3,:)+color_for_ROI_early(1,:))/2,'LineWidth',2)
% %     end
% % end
% hold on
% yl=yline(1/3,'--',{'Chance level'},'Color',[0 0 0],'LineWidth',1);
%
% yl.LabelVerticalAlignment = 'middle';
% yl.LabelHorizontalAlignment = 'center';
% yl.FontSize=12;
% ylim([0.3 0.8])
% switch event_list{end}
%     case 'target'
%
% hold on
% xline(-300,'-','Color',[0.5 0.5 0.5],'LineWidth',1)
% hold on
% xl=xline(-400,'-',{'Fixation'},'Color',[0.5 0.5 0.5],'LineWidth',1);
% xl.LabelVerticalAlignment = 'bottom';
% xl.LabelHorizontalAlignment = 'left';
% xl.FontSize=12;
%
% hold on
% xline(300,'-','Color',[0.25 0.75 0.25],'LineWidth',1)
% hold on
% xl=xline(200,'-',{'Response'},'Color',[0.25 0.75 0.25],'LineWidth',1);
% xl.LabelVerticalAlignment = 'bottom';
% xl.LabelHorizontalAlignment = 'left';
% xl.FontSize=12;
% end
% ylabel('Expected template decoding accuracy','FontSize',14)
% xlabel(sprintf('Time to %s',event_list{end}'))
% ax=gca;
% ax.XAxis.FontSize=12;
% ax.YAxis.FontSize=12;
% box off
% %legend({'PFC early','LIP early'},'FontSize',16,'Location','South')
% ylim([0.18 0.82])
% xlim([-600 350])
%
% subplot(3,3,3)
% plot(window_start_list(8:end),mean(PFC_across_time(8:end,3,:),3),'Color',color_for_ROI_late(3,:),'LineWidth',2)
% hold on
% plot(window_start_list(8:end),mean(LIP_across_time(8:end,3,:),3),'Color',color_for_ROI_late(1,:),'LineWidth',2)
% hold on
% shadedErrorBar(window_start_list(8:end),mean(PFC_across_time(8:end,3,:),3),std(PFC_across_time(8:end,3,:),0,3),{'color',color_for_ROI_late(3,:),'LineWidth',2},2)
% hold on
% shadedErrorBar(window_start_list(8:end),mean(LIP_across_time(8:end,3,:),3),std(LIP_across_time(8:end,3,:),0,3),{'color',color_for_ROI_late(1,:),'LineWidth',2},2)
% hold on
% % for i=8:length(window_start_list)
% %     if p_LIP_PFC(i,3)<0.05/length(window_start_list(8:end)) %Bonferoni
% %         plot(window_start_list(i),0.8,'*','Color',(color_for_ROI_late(3,:)+color_for_ROI_late(1,:))/2,'LineWidth',2)
% %     end
% % end
% hold on
% yl=yline(1/3,'--',{'Chance level'},'Color',[0 0 0],'LineWidth',1);
% yl.LabelVerticalAlignment = 'middle';
% yl.LabelHorizontalAlignment = 'center';
% yl.FontSize=12;
% switch event_list{end}
%     case 'target'
%
% hold on
% xline(-300,'-','Color',[0.5 0.5 0.5],'LineWidth',1)
% hold on
% xl=xline(-400,'-',{'Fixation'},'Color',[0.5 0.5 0.5],'LineWidth',1);
% xl.LabelVerticalAlignment = 'bottom';
% xl.LabelHorizontalAlignment = 'left';
% xl.FontSize=12;
%
% hold on
% xline(300,'-','Color',[0.25 0.75 0.25],'LineWidth',1)
% hold on
% xl=xline(200,'-',{'Response'},'Color',[0.25 0.75 0.25],'LineWidth',1);
% xl.LabelVerticalAlignment = 'bottom';
% xl.LabelHorizontalAlignment = 'left';
% xl.FontSize=12;
% end
% ylabel('Expected template decoding accuracy','FontSize',14)
% xlabel(sprintf('Time to %s',event_list{end}'))
% ax=gca;
% ax.XAxis.FontSize=12;
% ax.YAxis.FontSize=12;
% box off
% %legend({'PFC late','LIP late'},'FontSize',16,'Location','South')
% ylim([0.18 0.82])
% xlim([-600 350])
%
%
% subplot(3,3,2)
% plot(window_start_list(8:end),mean(PFC_across_time(8:end,2,:),3),'Color',color_for_ROI(3,:),'LineWidth',2)
% hold on
% plot(window_start_list(8:end),mean(LIP_across_time(8:end,2,:),3),'Color',color_for_ROI(1,:),'LineWidth',2)
% hold on
% shadedErrorBar(window_start_list(8:end),mean(PFC_across_time(8:end,2,:),3),std(PFC_across_time(8:end,2,:),0,3),{'color',color_for_ROI(3,:),'LineWidth',2},2)
% hold on
% shadedErrorBar(window_start_list(8:end),mean(LIP_across_time(8:end,2,:),3),std(LIP_across_time(8:end,2,:),0,3),{'color',color_for_ROI(1,:),'LineWidth',2},2)
% hold on
% % for i=8:length(window_start_list)
% %     if p_LIP_PFC(i,2)<0.05/length(window_start_list(8:end)) %Bonferoni
% %         plot(window_start_list(i),0.8,'*','Color',(color_for_ROI(3,:)+color_for_ROI(1,:))/2,'LineWidth',2)
% %     end
% % end
% hold on
% yl=yline(1/3,'--',{'Chance level'},'Color',[0 0 0],'LineWidth',1);
% yl.LabelVerticalAlignment = 'middle';
% yl.LabelHorizontalAlignment = 'center';
% yl.FontSize=12;
% switch event_list{end}
%     case 'target'
%
% hold on
% xline(-300,'-','Color',[0.5 0.5 0.5],'LineWidth',1)
% hold on
% xl=xline(-400,'-',{'Fixation'},'Color',[0.5 0.5 0.5],'LineWidth',1);
% xl.LabelVerticalAlignment = 'bottom';
% xl.LabelHorizontalAlignment = 'left';
% xl.FontSize=12;
%
% hold on
% xline(300,'-','Color',[0.25 0.75 0.25],'LineWidth',1)
% hold on
% xl=xline(200,'-',{'Response'},'Color',[0.25 0.75 0.25],'LineWidth',1);
% xl.LabelVerticalAlignment = 'bottom';
% xl.LabelHorizontalAlignment = 'left';
% xl.FontSize=12;
% end
% ylabel('Expected template decoding accuracy','FontSize',14)
% xlabel(sprintf('Time to %s',event_list{end}'))
% ax=gca;
% ax.XAxis.FontSize=12;
% ax.YAxis.FontSize=12;
% box off
% %legend({'PFC mid','LIP mid'},'FontSize',16,'Location','South')
% ylim([0.18 0.82])
% xlim([-600 350])
%
%
% subplot(3,3,4)
% plot(window_start_list(8:end),mean(FEF_across_time(8:end,1,:),3),'Color',color_for_ROI_early(2,:),'LineWidth',2)
% hold on
% plot(window_start_list(8:end),mean(LIP_across_time(8:end,1,:),3),'Color',color_for_ROI_early(1,:),'LineWidth',2)
% hold on
% shadedErrorBar(window_start_list(8:end),mean(FEF_across_time(8:end,1,:),3),std(FEF_across_time(8:end,1,:),0,3),{'color',color_for_ROI_early(2,:),'LineWidth',2},2)
% hold on
% shadedErrorBar(window_start_list(8:end),mean(LIP_across_time(8:end,1,:),3),std(LIP_across_time(8:end,1,:),0,3),{'color',color_for_ROI_early(1,:),'LineWidth',2},2)
% hold on
% % for i=8:length(window_start_list)
% %     if p_LIP_FEF(i,1)<0.05/length(window_start_list(8:end)) %Bonferoni
% %         plot(window_start_list(i),0.8,'*','Color',(color_for_ROI_early(2,:)+color_for_ROI_early(1,:))/2,'LineWidth',2)
% %     end
% % end
% hold on
% yl=yline(1/3,'--',{'Chance level'},'Color',[0 0 0],'LineWidth',1);
%
% yl.LabelVerticalAlignment = 'middle';
% yl.LabelHorizontalAlignment = 'center';
% yl.FontSize=12;
% ylim([0.3 0.8])
% switch event_list{end}
%     case 'target'
%
% hold on
% xline(-300,'-','Color',[0.5 0.5 0.5],'LineWidth',1)
% hold on
% xl=xline(-400,'-',{'Fixation'},'Color',[0.5 0.5 0.5],'LineWidth',1);
% xl.LabelVerticalAlignment = 'bottom';
% xl.LabelHorizontalAlignment = 'left';
% xl.FontSize=12;
%
% hold on
% xline(300,'-','Color',[0.25 0.75 0.25],'LineWidth',1)
% hold on
% xl=xline(200,'-',{'Response'},'Color',[0.25 0.75 0.25],'LineWidth',1);
% xl.LabelVerticalAlignment = 'bottom';
% xl.LabelHorizontalAlignment = 'left';
% xl.FontSize=12;
% end
% ylabel('Expected template decoding accuracy','FontSize',14)
% xlabel(sprintf('Time to %s',event_list{end}'))
% ax=gca;
% ax.XAxis.FontSize=12;
% ax.YAxis.FontSize=12;
% box off
% %legend({'FEF early','LIP early'},'FontSize',16,'Location','South')
% ylim([0.18 0.82])
% xlim([-600 350])
%
% subplot(3,3,6)
% plot(window_start_list(8:end),mean(FEF_across_time(8:end,3,:),3),'Color',color_for_ROI_late(2,:),'LineWidth',2)
% hold on
% plot(window_start_list(8:end),mean(LIP_across_time(8:end,3,:),3),'Color',color_for_ROI_late(1,:),'LineWidth',2)
% hold on
% shadedErrorBar(window_start_list(8:end),mean(FEF_across_time(8:end,3,:),3),std(FEF_across_time(8:end,3,:),0,3),{'color',color_for_ROI_late(2,:),'LineWidth',2},2)
% hold on
% shadedErrorBar(window_start_list(8:end),mean(LIP_across_time(8:end,3,:),3),std(LIP_across_time(8:end,3,:),0,3),{'color',color_for_ROI_late(1,:),'LineWidth',2},2)
% hold on
% % for i=8:length(window_start_list)
% %     if p_LIP_FEF(i,3)<0.05/length(window_start_list(8:end)) %Bonferoni
% %         plot(window_start_list(i),0.8,'*','Color',(color_for_ROI_late(2,:)+color_for_ROI_late(1,:))/2,'LineWidth',2)
% %     end
% % end
% hold on
% yl=yline(1/3,'--',{'Chance level'},'Color',[0 0 0],'LineWidth',1);
% yl.LabelVerticalAlignment = 'middle';
% yl.LabelHorizontalAlignment = 'center';
% yl.FontSize=12;
% switch event_list{end}
%     case 'target'
%
% hold on
% xline(-300,'-','Color',[0.5 0.5 0.5],'LineWidth',1)
% hold on
% xl=xline(-400,'-',{'Fixation'},'Color',[0.5 0.5 0.5],'LineWidth',1);
% xl.LabelVerticalAlignment = 'bottom';
% xl.LabelHorizontalAlignment = 'left';
% xl.FontSize=12;
%
% hold on
% xline(300,'-','Color',[0.25 0.75 0.25],'LineWidth',1)
% hold on
% xl=xline(200,'-',{'Response'},'Color',[0.25 0.75 0.25],'LineWidth',1);
% xl.LabelVerticalAlignment = 'bottom';
% xl.LabelHorizontalAlignment = 'left';
% xl.FontSize=12;
% end
% ylabel('Expected template decoding accuracy','FontSize',14)
% xlabel(sprintf('Time to %s',event_list{end}'))
% ax=gca;
% ax.XAxis.FontSize=12;
% ax.YAxis.FontSize=12;
% box off
% %legend({'FEF late','LIP late'},'FontSize',16,'Location','South')
% ylim([0.18 0.82])
% xlim([-600 350])
%
%
% subplot(3,3,5)
% plot(window_start_list(8:end),mean(FEF_across_time(8:end,2,:),3),'Color',color_for_ROI(2,:),'LineWidth',2)
% hold on
% plot(window_start_list(8:end),mean(LIP_across_time(8:end,2,:),3),'Color',color_for_ROI(1,:),'LineWidth',2)
% hold on
% shadedErrorBar(window_start_list(8:end),mean(FEF_across_time(8:end,2,:),3),std(FEF_across_time(8:end,2,:),0,3),{'color',color_for_ROI(2,:),'LineWidth',2},2)
% hold on
% shadedErrorBar(window_start_list(8:end),mean(LIP_across_time(8:end,2,:),3),std(LIP_across_time(8:end,2,:),0,3),{'color',color_for_ROI(1,:),'LineWidth',2},2)
% hold on
% % for i=8:length(window_start_list)
% %     if p_LIP_FEF(i,2)<0.05/length(window_start_list(8:end)) %Bonferoni
% %         plot(window_start_list(i),0.8,'*','Color',(color_for_ROI(2,:)+color_for_ROI(1,:))/2,'LineWidth',2)
% %     end
% % end
% hold on
% yl=yline(1/3,'--',{'Chance level'},'Color',[0 0 0],'LineWidth',1);
% yl.LabelVerticalAlignment = 'middle';
% yl.LabelHorizontalAlignment = 'center';
% yl.FontSize=12;
% switch event_list{end}
%     case 'target'
%
% hold on
% xline(-300,'-','Color',[0.5 0.5 0.5],'LineWidth',1)
% hold on
% xl=xline(-400,'-',{'Fixation'},'Color',[0.5 0.5 0.5],'LineWidth',1);
% xl.LabelVerticalAlignment = 'bottom';
% xl.LabelHorizontalAlignment = 'left';
% xl.FontSize=12;
%
% hold on
% xline(300,'-','Color',[0.25 0.75 0.25],'LineWidth',1)
% hold on
% xl=xline(200,'-',{'Response'},'Color',[0.25 0.75 0.25],'LineWidth',1);
% xl.LabelVerticalAlignment = 'bottom';
% xl.LabelHorizontalAlignment = 'left';
% xl.FontSize=12;
% end
% ylabel('Expected template decoding accuracy','FontSize',14)
% xlabel(sprintf('Time to %s',event_list{end}'))
% ax=gca;
% ax.XAxis.FontSize=12;
% ax.YAxis.FontSize=12;
% box off
% %legend({'FEF mid','LIP mid'},'FontSize',16,'Location','South')
% ylim([0.18 0.82])
% xlim([-600 350])
%
%
% subplot(3,3,7)
% plot(window_start_list(8:end),mean(FEF_across_time(8:end,1,:),3),'Color',color_for_ROI_early(2,:),'LineWidth',2)
% hold on
% plot(window_start_list(8:end),mean(PFC_across_time(8:end,1,:),3),'Color',color_for_ROI_early(3,:),'LineWidth',2)
% hold on
% shadedErrorBar(window_start_list(8:end),mean(FEF_across_time(8:end,1,:),3),std(FEF_across_time(8:end,1,:),0,3),{'color',color_for_ROI_early(2,:),'LineWidth',2},2)
% hold on
% shadedErrorBar(window_start_list(8:end),mean(PFC_across_time(8:end,1,:),3),std(PFC_across_time(8:end,1,:),0,3),{'color',color_for_ROI_early(3,:),'LineWidth',2},2)
% hold on
% % for i=8:length(window_start_list)
% %     if p_PFC_FEF(i,1)<0.05/length(window_start_list(8:end)) %Bonferoni
% %         plot(window_start_list(i),0.8,'*','Color',(color_for_ROI_early(2,:)+color_for_ROI_early(3,:))/2,'LineWidth',2)
% %     end
% % end
% hold on
% yl=yline(1/3,'--',{'Chance level'},'Color',[0 0 0],'LineWidth',1);
%
% yl.LabelVerticalAlignment = 'middle';
% yl.LabelHorizontalAlignment = 'center';
% yl.FontSize=12;
% ylim([0.3 0.8])
% switch event_list{end}
%     case 'target'
%
% hold on
% xline(-300,'-','Color',[0.5 0.5 0.5],'LineWidth',1)
% hold on
% xl=xline(-400,'-',{'Fixation'},'Color',[0.5 0.5 0.5],'LineWidth',1);
% xl.LabelVerticalAlignment = 'bottom';
% xl.LabelHorizontalAlignment = 'left';
% xl.FontSize=12;
%
% hold on
% xline(300,'-','Color',[0.25 0.75 0.25],'LineWidth',1)
% hold on
% xl=xline(200,'-',{'Response'},'Color',[0.25 0.75 0.25],'LineWidth',1);
% xl.LabelVerticalAlignment = 'bottom';
% xl.LabelHorizontalAlignment = 'left';
% xl.FontSize=12;
% end
% ylabel('Expected template decoding accuracy','FontSize',14)
% xlabel(sprintf('Time to %s',event_list{end}'))
% ax=gca;
% ax.XAxis.FontSize=12;
% ax.YAxis.FontSize=12;
% box off
% %legend({'FEF early','PFC early'},'FontSize',16,'Location','South')
% ylim([0.18 0.82])
% xlim([-600 350])
%
% subplot(3,3,9)
%
% plot(window_start_list(8:end),mean(FEF_across_time(8:end,3,:),3),'Color',color_for_ROI_late(2,:),'LineWidth',2)
% hold on
% plot(window_start_list(8:end),mean(PFC_across_time(8:end,3,:),3),'Color',color_for_ROI_late(3,:),'LineWidth',2)
% hold on
% shadedErrorBar(window_start_list(8:end),mean(FEF_across_time(8:end,3,:),3),std(FEF_across_time(8:end,3,:),0,3),{'color',color_for_ROI_late(2,:),'LineWidth',2},2)
% hold on
% shadedErrorBar(window_start_list(8:end),mean(PFC_across_time(8:end,3,:),3),std(PFC_across_time(8:end,3,:),0,3),{'color',color_for_ROI_late(3,:),'LineWidth',2},2)
% hold on
% % for i=8:length(window_start_list)
% %     if p_PFC_FEF(i,3)<0.05/length(window_start_list(8:end)) %Bonferoni
% %         plot(window_start_list(i),0.8,'*','Color',(color_for_ROI_late(2,:)+color_for_ROI_late(3,:))/2,'LineWidth',2)
% %     end
% % end
% hold on
% yl=yline(1/3,'--',{'Chance level'},'Color',[0 0 0],'LineWidth',1);
% yl.LabelVerticalAlignment = 'middle';
% yl.LabelHorizontalAlignment = 'center';
% yl.FontSize=12;
% switch event_list{end}
%     case 'target'
%
% hold on
% xline(-300,'-','Color',[0.5 0.5 0.5],'LineWidth',1)
% hold on
% xl=xline(-400,'-',{'Fixation'},'Color',[0.5 0.5 0.5],'LineWidth',1);
% xl.LabelVerticalAlignment = 'bottom';
% xl.LabelHorizontalAlignment = 'left';
% xl.FontSize=12;
%
% hold on
% xline(300,'-','Color',[0.25 0.75 0.25],'LineWidth',1)
% hold on
% xl=xline(200,'-',{'Response'},'Color',[0.25 0.75 0.25],'LineWidth',1);
% xl.LabelVerticalAlignment = 'bottom';
% xl.LabelHorizontalAlignment = 'left';
% xl.FontSize=12;
% end
% ylabel('Expected template decoding accuracy','FontSize',14)
% xlabel(sprintf('Time to %s',event_list{end}'))
% ax=gca;
% ax.XAxis.FontSize=12;
% ax.YAxis.FontSize=12;
% box off
% %legend({'FEF late','PFC late'},'FontSize',16,'Location','South')
% ylim([0.18 0.82])
% xlim([-600 350])
%
%
% subplot(3,3,8)
% plot(window_start_list(8:end),mean(FEF_across_time(8:end,2,:),3),'Color',color_for_ROI(2,:),'LineWidth',2)
% hold on
% plot(window_start_list(8:end),mean(PFC_across_time(8:end,2,:),3),'Color',color_for_ROI(3,:),'LineWidth',2)
% hold on
% shadedErrorBar(window_start_list(8:end),mean(FEF_across_time(8:end,2,:),3),std(FEF_across_time(8:end,2,:),0,3),{'color',color_for_ROI(2,:),'LineWidth',2},2)
% hold on
% shadedErrorBar(window_start_list(8:end),mean(PFC_across_time(8:end,2,:),3),std(PFC_across_time(8:end,2,:),0,3),{'color',color_for_ROI(3,:),'LineWidth',2},2)
% hold on
% % for i=8:length(window_start_list)
% %     if p_PFC_FEF(i,2)<0.05/length(window_start_list(8:end)) %Bonferoni
% %         plot(window_start_list(i),0.8,'*','Color',(color_for_ROI(2,:)+color_for_ROI(3,:))/2,'LineWidth',2)
% %     end
% % end
% hold on
% yl=yline(1/3,'--',{'Chance level'},'Color',[0 0 0],'LineWidth',1);
% yl.LabelVerticalAlignment = 'middle';
% yl.LabelHorizontalAlignment = 'center';
% yl.FontSize=12;
% switch event_list{end}
%     case 'target'
%
% hold on
% xline(-300,'-','Color',[0.5 0.5 0.5],'LineWidth',1)
% hold on
% xl=xline(-400,'-',{'Fixation'},'Color',[0.5 0.5 0.5],'LineWidth',1);
% xl.LabelVerticalAlignment = 'bottom';
% xl.LabelHorizontalAlignment = 'left';
% xl.FontSize=12;
%
% hold on
% xline(300,'-','Color',[0.25 0.75 0.25],'LineWidth',1)
% hold on
% xl=xline(200,'-',{'Response'},'Color',[0.25 0.75 0.25],'LineWidth',1);
% xl.LabelVerticalAlignment = 'bottom';
% xl.LabelHorizontalAlignment = 'left';
% xl.FontSize=12;
% end
% ylabel('Expected template decoding accuracy','FontSize',14)
% xlabel(sprintf('Time to %s',event_list{end}'))
% ax=gca;
% ax.XAxis.FontSize=12;
% ax.YAxis.FontSize=12;
% box off
% %legend({'FEF mid','PFC mid'},'FontSize',16,'Location','South')
% ylim([0.18 0.82])
% xlim([-600 350])
%
%
%


%%
%%
% figure
% subplot(3,3,1)
% plot(window_start_list(1:7),mean(PFC_across_time(1:7,1,:),3),'Color',color_for_ROI_early(3,:),'LineWidth',2)
% hold on
% plot(window_start_list(1:7),mean(LIP_across_time(1:7,1,:),3),'Color',color_for_ROI_early(1,:),'LineWidth',2)
% hold on
% shadedErrorBar(window_start_list(1:7),mean(PFC_across_time(1:7,1,:),3),std(PFC_across_time(1:7,1,:),0,3),{'color',color_for_ROI_early(3,:),'LineWidth',2},2)
% hold on
% shadedErrorBar(window_start_list(1:7),mean(LIP_across_time(1:7,1,:),3),std(LIP_across_time(1:7,1,:),0,3),{'color',color_for_ROI_early(1,:),'LineWidth',2},2)
% hold on
% for i=1:7
%     if p_LIP_PFC(i,1)<0.05/length(window_start_list(1:7)) %Bonferoni
%         plot(window_start_list(i),0.8,'*','Color',(color_for_ROI_early(3,:)+color_for_ROI_early(1,:))/2,'LineWidth',2)
%     end
% end
% hold on
% yl=yline(1/3,'--',{'Chance level'},'Color',[0 0 0],'LineWidth',1);
%
% yl.LabelVerticalAlignment = 'middle';
% yl.LabelHorizontalAlignment = 'center';
% yl.FontSize=12;
% ylim([0.3 0.8])
% switch event_list{1}
%     case 'target'
%
% hold on
% xline(-300,'-','Color',[0.5 0.5 0.5],'LineWidth',1)
% hold on
% xl=xline(-400,'-',{'Fixation'},'Color',[0.5 0.5 0.5],'LineWidth',1);
% xl.LabelVerticalAlignment = 'bottom';
% xl.LabelHorizontalAlignment = 'left';
% xl.FontSize=12;
%
% hold on
% xline(300,'-','Color',[0.25 0.75 0.25],'LineWidth',1)
% hold on
% xl=xline(200,'-',{'Response'},'Color',[0.25 0.75 0.25],'LineWidth',1);
% xl.LabelVerticalAlignment = 'bottom';
% xl.LabelHorizontalAlignment = 'left';
% xl.FontSize=12;
% end
% ylabel('Expected template decoding accuracy','FontSize',14)
% xlabel(sprintf('Time to end of reward'))
% ax=gca;
% ax.XAxis.FontSize=12;
% ax.YAxis.FontSize=12;
% box off
% %legend({'PFC early','LIP early'},'FontSize',16,'Location','South')
% ylim([0.18 0.82])
%
% subplot(3,3,3)
% plot(window_start_list(1:7),mean(PFC_across_time(1:7,3,:),3),'Color',color_for_ROI_late(3,:),'LineWidth',2)
% hold on
% plot(window_start_list(1:7),mean(LIP_across_time(1:7,3,:),3),'Color',color_for_ROI_late(1,:),'LineWidth',2)
% hold on
% shadedErrorBar(window_start_list(1:7),mean(PFC_across_time(1:7,3,:),3),std(PFC_across_time(1:7,3,:),0,3),{'color',color_for_ROI_late(3,:),'LineWidth',2},2)
% hold on
% shadedErrorBar(window_start_list(1:7),mean(LIP_across_time(1:7,3,:),3),std(LIP_across_time(1:7,3,:),0,3),{'color',color_for_ROI_late(1,:),'LineWidth',2},2)
% hold on
% for i=1:7
%     if p_LIP_PFC(i,3)<0.05/length(window_start_list(1:7)) %Bonferoni
%         plot(window_start_list(i),0.8,'*','Color',(color_for_ROI_late(3,:)+color_for_ROI_late(1,:))/2,'LineWidth',2)
%     end
% end
% hold on
% yl=yline(1/3,'--',{'Chance level'},'Color',[0 0 0],'LineWidth',1);
% yl.LabelVerticalAlignment = 'middle';
% yl.LabelHorizontalAlignment = 'center';
% yl.FontSize=12;
% switch event_list{1}
%     case 'target'
%
% hold on
% xline(-300,'-','Color',[0.5 0.5 0.5],'LineWidth',1)
% hold on
% xl=xline(-400,'-',{'Fixation'},'Color',[0.5 0.5 0.5],'LineWidth',1);
% xl.LabelVerticalAlignment = 'bottom';
% xl.LabelHorizontalAlignment = 'left';
% xl.FontSize=12;
%
% hold on
% xline(300,'-','Color',[0.25 0.75 0.25],'LineWidth',1)
% hold on
% xl=xline(200,'-',{'Response'},'Color',[0.25 0.75 0.25],'LineWidth',1);
% xl.LabelVerticalAlignment = 'bottom';
% xl.LabelHorizontalAlignment = 'left';
% xl.FontSize=12;
% end
% ylabel('Expected template decoding accuracy','FontSize',14)
% xlabel(sprintf('Time to end of reward'))
% ax=gca;
% ax.XAxis.FontSize=12;
% ax.YAxis.FontSize=12;
% box off
% %legend({'PFC late','LIP late'},'FontSize',16,'Location','South')
% ylim([0.18 0.82])
%
%
% subplot(3,3,2)
% plot(window_start_list(1:7),mean(PFC_across_time(1:7,2,:),3),'Color',color_for_ROI(3,:),'LineWidth',2)
% hold on
% plot(window_start_list(1:7),mean(LIP_across_time(1:7,2,:),3),'Color',color_for_ROI(1,:),'LineWidth',2)
% hold on
% shadedErrorBar(window_start_list(1:7),mean(PFC_across_time(1:7,2,:),3),std(PFC_across_time(1:7,2,:),0,3),{'color',color_for_ROI(3,:),'LineWidth',2},2)
% hold on
% shadedErrorBar(window_start_list(1:7),mean(LIP_across_time(1:7,2,:),3),std(LIP_across_time(1:7,2,:),0,3),{'color',color_for_ROI(1,:),'LineWidth',2},2)
% hold on
% for i=1:7
%     if p_LIP_PFC(i,2)<0.05/length(window_start_list(1:7)) %Bonferoni
%         plot(window_start_list(i),0.8,'*','Color',(color_for_ROI(3,:)+color_for_ROI(1,:))/2,'LineWidth',2)
%     end
% end
% hold on
% yl=yline(1/3,'--',{'Chance level'},'Color',[0 0 0],'LineWidth',1);
% yl.LabelVerticalAlignment = 'middle';
% yl.LabelHorizontalAlignment = 'center';
% yl.FontSize=12;
% switch event_list{1}
%     case 'target'
%
% hold on
% xline(-300,'-','Color',[0.5 0.5 0.5],'LineWidth',1)
% hold on
% xl=xline(-400,'-',{'Fixation'},'Color',[0.5 0.5 0.5],'LineWidth',1);
% xl.LabelVerticalAlignment = 'bottom';
% xl.LabelHorizontalAlignment = 'left';
% xl.FontSize=12;
%
% hold on
% xline(300,'-','Color',[0.25 0.75 0.25],'LineWidth',1)
% hold on
% xl=xline(200,'-',{'Response'},'Color',[0.25 0.75 0.25],'LineWidth',1);
% xl.LabelVerticalAlignment = 'bottom';
% xl.LabelHorizontalAlignment = 'left';
% xl.FontSize=12;
% end
% ylabel('Expected template decoding accuracy','FontSize',14)
% xlabel(sprintf('Time to end of reward'))
% ax=gca;
% ax.XAxis.FontSize=12;
% ax.YAxis.FontSize=12;
% box off
% %legend({'PFC mid','LIP mid'},'FontSize',16,'Location','South')
% ylim([0.18 0.82])
%
%
% subplot(3,3,4)
% plot(window_start_list(1:7),mean(FEF_across_time(1:7,1,:),3),'Color',color_for_ROI_early(2,:),'LineWidth',2)
% hold on
% plot(window_start_list(1:7),mean(LIP_across_time(1:7,1,:),3),'Color',color_for_ROI_early(1,:),'LineWidth',2)
% hold on
% shadedErrorBar(window_start_list(1:7),mean(FEF_across_time(1:7,1,:),3),std(FEF_across_time(1:7,1,:),0,3),{'color',color_for_ROI_early(2,:),'LineWidth',2},2)
% hold on
% shadedErrorBar(window_start_list(1:7),mean(LIP_across_time(1:7,1,:),3),std(LIP_across_time(1:7,1,:),0,3),{'color',color_for_ROI_early(1,:),'LineWidth',2},2)
% hold on
% for i=1:7
%     if p_LIP_FEF(i,1)<0.05/length(window_start_list(1:7)) %Bonferoni
%         plot(window_start_list(i),0.8,'*','Color',(color_for_ROI_early(2,:)+color_for_ROI_early(1,:))/2,'LineWidth',2)
%     end
% end
% hold on
% yl=yline(1/3,'--',{'Chance level'},'Color',[0 0 0],'LineWidth',1);
%
% yl.LabelVerticalAlignment = 'middle';
% yl.LabelHorizontalAlignment = 'center';
% yl.FontSize=12;
% ylim([0.3 0.8])
% switch event_list{1}
%     case 'target'
%
% hold on
% xline(-300,'-','Color',[0.5 0.5 0.5],'LineWidth',1)
% hold on
% xl=xline(-400,'-',{'Fixation'},'Color',[0.5 0.5 0.5],'LineWidth',1);
% xl.LabelVerticalAlignment = 'bottom';
% xl.LabelHorizontalAlignment = 'left';
% xl.FontSize=12;
%
% hold on
% xline(300,'-','Color',[0.25 0.75 0.25],'LineWidth',1)
% hold on
% xl=xline(200,'-',{'Response'},'Color',[0.25 0.75 0.25],'LineWidth',1);
% xl.LabelVerticalAlignment = 'bottom';
% xl.LabelHorizontalAlignment = 'left';
% xl.FontSize=12;
% end
% ylabel('Expected template decoding accuracy','FontSize',14)
% xlabel(sprintf('Time to end of reward'))
% ax=gca;
% ax.XAxis.FontSize=12;
% ax.YAxis.FontSize=12;
% box off
% %legend({'FEF early','LIP early'},'FontSize',16,'Location','South')
% ylim([0.18 0.82])
%
% subplot(3,3,6)
% plot(window_start_list(1:7),mean(FEF_across_time(1:7,3,:),3),'Color',color_for_ROI_late(2,:),'LineWidth',2)
% hold on
% plot(window_start_list(1:7),mean(LIP_across_time(1:7,3,:),3),'Color',color_for_ROI_late(1,:),'LineWidth',2)
% hold on
% shadedErrorBar(window_start_list(1:7),mean(FEF_across_time(1:7,3,:),3),std(FEF_across_time(1:7,3,:),0,3),{'color',color_for_ROI_late(2,:),'LineWidth',2},2)
% hold on
% shadedErrorBar(window_start_list(1:7),mean(LIP_across_time(1:7,3,:),3),std(LIP_across_time(1:7,3,:),0,3),{'color',color_for_ROI_late(1,:),'LineWidth',2},2)
% hold on
% for i=1:7
%     if p_LIP_FEF(i,3)<0.05/length(window_start_list(1:7)) %Bonferoni
%         plot(window_start_list(i),0.8,'*','Color',(color_for_ROI_late(2,:)+color_for_ROI_late(1,:))/2,'LineWidth',2)
%     end
% end
% hold on
% yl=yline(1/3,'--',{'Chance level'},'Color',[0 0 0],'LineWidth',1);
% yl.LabelVerticalAlignment = 'middle';
% yl.LabelHorizontalAlignment = 'center';
% yl.FontSize=12;
% switch event_list{1}
%     case 'target'
%
% hold on
% xline(-300,'-','Color',[0.5 0.5 0.5],'LineWidth',1)
% hold on
% xl=xline(-400,'-',{'Fixation'},'Color',[0.5 0.5 0.5],'LineWidth',1);
% xl.LabelVerticalAlignment = 'bottom';
% xl.LabelHorizontalAlignment = 'left';
% xl.FontSize=12;
%
% hold on
% xline(300,'-','Color',[0.25 0.75 0.25],'LineWidth',1)
% hold on
% xl=xline(200,'-',{'Response'},'Color',[0.25 0.75 0.25],'LineWidth',1);
% xl.LabelVerticalAlignment = 'bottom';
% xl.LabelHorizontalAlignment = 'left';
% xl.FontSize=12;
% end
% ylabel('Expected template decoding accuracy','FontSize',14)
% xlabel(sprintf('Time to end of reward'))
% ax=gca;
% ax.XAxis.FontSize=12;
% ax.YAxis.FontSize=12;
% box off
% %legend({'FEF late','LIP late'},'FontSize',16,'Location','South')
% ylim([0.18 0.82])
%
%
% subplot(3,3,5)
% plot(window_start_list(1:7),mean(FEF_across_time(1:7,2,:),3),'Color',color_for_ROI(2,:),'LineWidth',2)
% hold on
% plot(window_start_list(1:7),mean(LIP_across_time(1:7,2,:),3),'Color',color_for_ROI(1,:),'LineWidth',2)
% hold on
% shadedErrorBar(window_start_list(1:7),mean(FEF_across_time(1:7,2,:),3),std(FEF_across_time(1:7,2,:),0,3),{'color',color_for_ROI(2,:),'LineWidth',2},2)
% hold on
% shadedErrorBar(window_start_list(1:7),mean(LIP_across_time(1:7,2,:),3),std(LIP_across_time(1:7,2,:),0,3),{'color',color_for_ROI(1,:),'LineWidth',2},2)
% hold on
% for i=1:7
%     if p_LIP_FEF(i,2)<0.05/length(window_start_list(1:7)) %Bonferoni
%         plot(window_start_list(i),0.8,'*','Color',(color_for_ROI(2,:)+color_for_ROI(1,:))/2,'LineWidth',2)
%     end
% end
% hold on
% yl=yline(1/3,'--',{'Chance level'},'Color',[0 0 0],'LineWidth',1);
% yl.LabelVerticalAlignment = 'middle';
% yl.LabelHorizontalAlignment = 'center';
% yl.FontSize=12;
% switch event_list{1}
%     case 'target'
%
% hold on
% xline(-300,'-','Color',[0.5 0.5 0.5],'LineWidth',1)
% hold on
% xl=xline(-400,'-',{'Fixation'},'Color',[0.5 0.5 0.5],'LineWidth',1);
% xl.LabelVerticalAlignment = 'bottom';
% xl.LabelHorizontalAlignment = 'left';
% xl.FontSize=12;
%
% hold on
% xline(300,'-','Color',[0.25 0.75 0.25],'LineWidth',1)
% hold on
% xl=xline(200,'-',{'Response'},'Color',[0.25 0.75 0.25],'LineWidth',1);
% xl.LabelVerticalAlignment = 'bottom';
% xl.LabelHorizontalAlignment = 'left';
% xl.FontSize=12;
% end
% ylabel('Expected template decoding accuracy','FontSize',14)
% xlabel(sprintf('Time to end of reward'))
% ax=gca;
% ax.XAxis.FontSize=12;
% ax.YAxis.FontSize=12;
% box off
% %legend({'FEF mid','LIP mid'},'FontSize',16,'Location','South')
% ylim([0.18 0.82])
%
%
% subplot(3,3,7)
% plot(window_start_list(1:7),mean(FEF_across_time(1:7,1,:),3),'Color',color_for_ROI_early(2,:),'LineWidth',2)
% hold on
% plot(window_start_list(1:7),mean(PFC_across_time(1:7,1,:),3),'Color',color_for_ROI_early(3,:),'LineWidth',2)
% hold on
% shadedErrorBar(window_start_list(1:7),mean(FEF_across_time(1:7,1,:),3),std(FEF_across_time(1:7,1,:),0,3),{'color',color_for_ROI_early(2,:),'LineWidth',2},2)
% hold on
% shadedErrorBar(window_start_list(1:7),mean(PFC_across_time(1:7,1,:),3),std(PFC_across_time(1:7,1,:),0,3),{'color',color_for_ROI_early(3,:),'LineWidth',2},2)
% hold on
% for i=1:7
%     if p_PFC_FEF(i,1)<0.05/length(window_start_list(1:7)) %Bonferoni
%         plot(window_start_list(i),0.8,'*','Color',(color_for_ROI_early(2,:)+color_for_ROI_early(3,:))/2,'LineWidth',2)
%     end
% end
% hold on
% yl=yline(1/3,'--',{'Chance level'},'Color',[0 0 0],'LineWidth',1);
%
% yl.LabelVerticalAlignment = 'middle';
% yl.LabelHorizontalAlignment = 'center';
% yl.FontSize=12;
% ylim([0.3 0.8])
% switch event_list{1}
%     case 'target'
%
% hold on
% xline(-300,'-','Color',[0.5 0.5 0.5],'LineWidth',1)
% hold on
% xl=xline(-400,'-',{'Fixation'},'Color',[0.5 0.5 0.5],'LineWidth',1);
% xl.LabelVerticalAlignment = 'bottom';
% xl.LabelHorizontalAlignment = 'left';
% xl.FontSize=12;
%
% hold on
% xline(300,'-','Color',[0.25 0.75 0.25],'LineWidth',1)
% hold on
% xl=xline(200,'-',{'Response'},'Color',[0.25 0.75 0.25],'LineWidth',1);
% xl.LabelVerticalAlignment = 'bottom';
% xl.LabelHorizontalAlignment = 'left';
% xl.FontSize=12;
% end
% ylabel('Expected template decoding accuracy','FontSize',14)
% xlabel(sprintf('Time to end of reward'))
% ax=gca;
% ax.XAxis.FontSize=12;
% ax.YAxis.FontSize=12;
% box off
% %legend({'FEF early','PFC early'},'FontSize',16,'Location','South')
% ylim([0.18 0.82])
%
% subplot(3,3,9)
%
% plot(window_start_list(1:7),mean(FEF_across_time(1:7,3,:),3),'Color',color_for_ROI_late(2,:),'LineWidth',2)
% hold on
% plot(window_start_list(1:7),mean(PFC_across_time(1:7,3,:),3),'Color',color_for_ROI_late(3,:),'LineWidth',2)
% hold on
% shadedErrorBar(window_start_list(1:7),mean(FEF_across_time(1:7,3,:),3),std(FEF_across_time(1:7,3,:),0,3),{'color',color_for_ROI_late(2,:),'LineWidth',2},2)
% hold on
% shadedErrorBar(window_start_list(1:7),mean(PFC_across_time(1:7,3,:),3),std(PFC_across_time(1:7,3,:),0,3),{'color',color_for_ROI_late(3,:),'LineWidth',2},2)
% hold on
% for i=1:7
%     if p_PFC_FEF(i,3)<0.05/length(window_start_list(1:7)) %Bonferoni
%         plot(window_start_list(i),0.8,'*','Color',(color_for_ROI_late(2,:)+color_for_ROI_late(3,:))/2,'LineWidth',2)
%     end
% end
% hold on
% yl=yline(1/3,'--',{'Chance level'},'Color',[0 0 0],'LineWidth',1);
% yl.LabelVerticalAlignment = 'middle';
% yl.LabelHorizontalAlignment = 'center';
% yl.FontSize=12;
% switch event_list{1}
%     case 'target'
%
% hold on
% xline(-300,'-','Color',[0.5 0.5 0.5],'LineWidth',1)
% hold on
% xl=xline(-400,'-',{'Fixation'},'Color',[0.5 0.5 0.5],'LineWidth',1);
% xl.LabelVerticalAlignment = 'bottom';
% xl.LabelHorizontalAlignment = 'left';
% xl.FontSize=12;
%
% hold on
% xline(300,'-','Color',[0.25 0.75 0.25],'LineWidth',1)
% hold on
% xl=xline(200,'-',{'Response'},'Color',[0.25 0.75 0.25],'LineWidth',1);
% xl.LabelVerticalAlignment = 'bottom';
% xl.LabelHorizontalAlignment = 'left';
% xl.FontSize=12;
% end
% ylabel('Expected template decoding accuracy','FontSize',14)
% xlabel(sprintf('Time to end of reward'))
% ax=gca;
% ax.XAxis.FontSize=12;
% ax.YAxis.FontSize=12;
% box off
% %legend({'FEF late','PFC late'},'FontSize',16,'Location','South')
% ylim([0.18 0.82])
%
%
% subplot(3,3,8)
% plot(window_start_list(1:7),mean(FEF_across_time(1:7,2,:),3),'Color',color_for_ROI(2,:),'LineWidth',2)
% hold on
% plot(window_start_list(1:7),mean(PFC_across_time(1:7,2,:),3),'Color',color_for_ROI(3,:),'LineWidth',2)
% hold on
% shadedErrorBar(window_start_list(1:7),mean(FEF_across_time(1:7,2,:),3),std(FEF_across_time(1:7,2,:),0,3),{'color',color_for_ROI(2,:),'LineWidth',2},2)
% hold on
% shadedErrorBar(window_start_list(1:7),mean(PFC_across_time(1:7,2,:),3),std(PFC_across_time(1:7,2,:),0,3),{'color',color_for_ROI(3,:),'LineWidth',2},2)
% hold on
% for i=1:7
%     if p_PFC_FEF(i,2)<0.05/length(window_start_list(1:7)) %Bonferoni
%         plot(window_start_list(i),0.8,'*','Color',(color_for_ROI(2,:)+color_for_ROI(3,:))/2,'LineWidth',2)
%     end
% end
% hold on
% yl=yline(1/3,'--',{'Chance level'},'Color',[0 0 0],'LineWidth',1);
% yl.LabelVerticalAlignment = 'middle';
% yl.LabelHorizontalAlignment = 'center';
% yl.FontSize=12;
% switch event_list{1}
%     case 'target'
%
% hold on
% xline(-300,'-','Color',[0.5 0.5 0.5],'LineWidth',1)
% hold on
% xl=xline(-400,'-',{'Fixation'},'Color',[0.5 0.5 0.5],'LineWidth',1);
% xl.LabelVerticalAlignment = 'bottom';
% xl.LabelHorizontalAlignment = 'left';
% xl.FontSize=12;
%
% hold on
% xline(300,'-','Color',[0.25 0.75 0.25],'LineWidth',1)
% hold on
% xl=xline(200,'-',{'Response'},'Color',[0.25 0.75 0.25],'LineWidth',1);
% xl.LabelVerticalAlignment = 'bottom';
% xl.LabelHorizontalAlignment = 'left';
% xl.FontSize=12;
% end
% ylabel('Expected template decoding accuracy','FontSize',14)
% xlabel(sprintf('Time to end of reward'))
% ax=gca;
% ax.XAxis.FontSize=12;
% ax.YAxis.FontSize=12;
% box off
% %legend({'FEF mid','PFC mid'},'FontSize',16,'Location','South')
% ylim([0.18 0.82])
%





% %%
% figure
% plot(window_start_list,mean(LIP_across_time(:,1,:),3),'Color',color_for_ROI_early(1,:),'LineWidth',2)
% hold on
% plot(window_start_list,mean(LIP_across_time(:,2,:),3),'Color',color_for_ROI(1,:),'LineWidth',2)
% hold on
% plot(window_start_list,mean(LIP_across_time(:,3,:),3),'Color',color_for_ROI_late(1,:),'LineWidth',2)
% hold on
% shadedErrorBar(window_start_list,mean(LIP_across_time(:,1,:),3),std(LIP_across_time(:,1,:),0,3),{'color',color_for_ROI_early(1,:),'LineWidth',2},2)
% hold on
% shadedErrorBar(window_start_list,mean(LIP_across_time(:,2,:),3),std(LIP_across_time(:,2,:),0,3),{'color',color_for_ROI(1,:),'LineWidth',2},2)
% hold on
% shadedErrorBar(window_start_list,mean(LIP_across_time(:,3,:),3),std(LIP_across_time(:,3,:),0,3),{'color',color_for_ROI_late(1,:),'LineWidth',2},2)
% hold on
% yl=yline(1/3,'--',{'Chance level'},'Color',[0 0 0],'LineWidth',1);
% hold on
%
% yl.LabelVerticalAlignment = 'middle';
% yl.LabelHorizontalAlignment = 'center';
% yl.FontSize=12;
% switch event
%     case 'target'
%
% hold on
% xline(-300,'-','Color',[0.5 0.5 0.5],'LineWidth',1)
% hold on
% xl=xline(-400,'-',{'Fixation'},'Color',[0.5 0.5 0.5],'LineWidth',1);
% xl.LabelVerticalAlignment = 'bottom';
% xl.LabelHorizontalAlignment = 'left';
% xl.FontSize=12;
%
% hold on
% xline(300,'-','Color',[0.25 0.75 0.25],'LineWidth',1)
% hold on
% xl=xline(200,'-',{'Response'},'Color',[0.25 0.75 0.25],'LineWidth',1);
% xl.LabelVerticalAlignment = 'bottom';
% xl.LabelHorizontalAlignment = 'left';
% xl.FontSize=12;
% end
% ylabel('Expected template decoding accuracy','FontSize',14)
% xlabel(sprintf('Time to %s',event_list{end}'))
% ax=gca;
% ax.XAxis.FontSize=12;
% ax.YAxis.FontSize=12;
% box off
% %legend({'Early in block','Mid block','Late in block'},'FontSize',16,'Location','South')
% ylim([0.18 0.82])
%
% figure
% plot(window_start_list,mean(FEF_across_time(:,1,:),3),'Color',color_for_ROI_early(2,:),'LineWidth',2)
% hold on
% plot(window_start_list,mean(FEF_across_time(:,2,:),3),'Color',color_for_ROI(2,:),'LineWidth',2)
% hold on
% plot(window_start_list,mean(FEF_across_time(:,3,:),3),'Color',color_for_ROI_late(2,:),'LineWidth',2)
% hold on
% shadedErrorBar(window_start_list,mean(FEF_across_time(:,1,:),3),std(FEF_across_time(:,1,:),0,3),{'color',color_for_ROI_early(2,:),'LineWidth',2},2)
% hold on
% shadedErrorBar(window_start_list,mean(FEF_across_time(:,2,:),3),std(FEF_across_time(:,2,:),0,3),{'color',color_for_ROI(2,:),'LineWidth',2},2)
% hold on
% shadedErrorBar(window_start_list,mean(FEF_across_time(:,3,:),3),std(FEF_across_time(:,3,:),0,3),{'color',color_for_ROI_late(2,:),'LineWidth',2},2)
% hold on
% yl=yline(1/3,'--',{'Chance level'},'Color',[0 0 0],'LineWidth',1);
% hold on
%
% yl.LabelVerticalAlignment = 'middle';
% yl.LabelHorizontalAlignment = 'center';
% yl.FontSize=12;
% switch event
%     case 'target'
%
% hold on
% xline(-300,'-','Color',[0.5 0.5 0.5],'LineWidth',1)
% hold on
% xl=xline(-400,'-',{'Fixation'},'Color',[0.5 0.5 0.5],'LineWidth',1);
% xl.LabelVerticalAlignment = 'bottom';
% xl.LabelHorizontalAlignment = 'left';
% xl.FontSize=12;
%
% hold on
% xline(300,'-','Color',[0.25 0.75 0.25],'LineWidth',1)
% hold on
% xl=xline(200,'-',{'Response'},'Color',[0.25 0.75 0.25],'LineWidth',1);
% xl.LabelVerticalAlignment = 'bottom';
% xl.LabelHorizontalAlignment = 'left';
% xl.FontSize=12;
% end
% ylabel('Expected template decoding accuracy','FontSize',14)
% xlabel(sprintf('Time to %s',event_list{end}'))
% ax=gca;
% ax.XAxis.FontSize=12;
% ax.YAxis.FontSize=12;
% box off
% %legend({'Early in block','Mid block','Late in block'},'FontSize',16,'Location','South')
% ylim([0.18 0.82])
%
%
%
% figure
% plot(window_start_list,mean(PFC_across_time(:,1,:),3),'Color',color_for_ROI_early(3,:),'LineWidth',2)
% hold on
% plot(window_start_list,mean(PFC_across_time(:,2,:),3),'Color',color_for_ROI(3,:),'LineWidth',2)
% hold on
% plot(window_start_list,mean(PFC_across_time(:,3,:),3),'Color',color_for_ROI_late(3,:),'LineWidth',2)
% hold on
% shadedErrorBar(window_start_list,mean(PFC_across_time(:,1,:),3),std(PFC_across_time(:,1,:),0,3),{'color',color_for_ROI_early(3,:),'LineWidth',2},2)
% hold on
% shadedErrorBar(window_start_list,mean(PFC_across_time(:,2,:),3),std(PFC_across_time(:,2,:),0,3),{'color',color_for_ROI(3,:),'LineWidth',2},2)
% hold on
% shadedErrorBar(window_start_list,mean(PFC_across_time(:,3,:),3),std(PFC_across_time(:,3,:),0,3),{'color',color_for_ROI_late(3,:),'LineWidth',2},2)
% hold on
% yl=yline(1/3,'--',{'Chance level'},'Color',[0 0 0],'LineWidth',1);
% hold on
%
% yl.LabelVerticalAlignment = 'middle';
% yl.LabelHorizontalAlignment = 'center';
% yl.FontSize=12;
% switch event
%     case 'target'
%
% hold on
% xline(-300,'-','Color',[0.5 0.5 0.5],'LineWidth',1)
% hold on
% xl=xline(-400,'-',{'Fixation'},'Color',[0.5 0.5 0.5],'LineWidth',1);
% xl.LabelVerticalAlignment = 'bottom';
% xl.LabelHorizontalAlignment = 'left';
% xl.FontSize=12;
%
% hold on
% xline(300,'-','Color',[0.25 0.75 0.25],'LineWidth',1)
% hold on
% xl=xline(200,'-',{'Response'},'Color',[0.25 0.75 0.25],'LineWidth',1);
% xl.LabelVerticalAlignment = 'bottom';
% xl.LabelHorizontalAlignment = 'left';
% xl.FontSize=12;
% end
% ylabel('Expected template decoding accuracy','FontSize',14)
% xlabel(sprintf('Time to %s',event_list{end}'))
% ax=gca;
% ax.XAxis.FontSize=12;
% ax.YAxis.FontSize=12;
% box off
% %legend({'Early in block','Mid block','Late in block'},'FontSize',16,'Location','South')
% ylim([0.18 0.82])




%% now load the combined

fsroot='/Volumes/buschman';
event='target';

task='Learning_Attentional_Templates';
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

task='Learning_Attentional_Templates';
subtask='exploreexploit/Restricted_FEF/Reset_RW_model';

subsubtask_classifier=sprintf('Classifier_%s_%d_to_%s_%d',event_list{1},window_start_list(1),event_list{end},window_start_list(end));

dirstem = fullfile('Projects',task,'Analysis','Electrophy_analysis',subtask,subsubtask_classifier);
data_path_clasifier = fullfile(fsroot,dirstem);

color_for_ROI=[0 0.447 0.741; 0.85 0.325 0.098; 0.929 0.694 0.125]; %LIP FEF PFC

load('colors')

color_belief(1,:)=colors(1,:);
color_belief(2,:)=colors(14,:);
color_belief(3,:)=colors(27,:);

this_combination=2;


for n_tt=1:100
    
    ROI='LIP';
    if this_combination==1
        results_name=sprintf('Pseudo_pop_peak_belief_prog_combined_time_results_more_trials_NN_%s_600ms_%d',ROI,n_tt);
    elseif this_combination==2
        results_name=sprintf('Pseudo_pop_peak_belief_prog_combined_time_results_more_trials_%s_900ms_NN_%d',ROI,n_tt);
    end
    
    %     if exist(sprintf('%s.mat',fullfile(data_path_clasifier,results_name)))==2
    
    load(fullfile(data_path_clasifier,results_name),'Classification_net_correct','Classification_net_correct_across_time');
    
    LIP.Classification_net_correct_comb(:,:,n_tt)=Classification_net_correct;
    LIP.Classification_net_correct_across_time(:,:,:,n_tt)=Classification_net_correct_across_time;
    %         LIP.W_combined(:,:,n_tt)=W;
    
    clear Classification_net_correct* 
    %
    %     else
    %         LIP.Classification_net_correct_combined(:,:,n_tt)=NaN;
    %         LIP.W_combined(:,:,n_tt)=NaN;
    %         sprintf('%s %d missing',ROI,n_tt)
    %
    %     end
    
    %FEF
    ROI='FEF';
    
    if this_combination==1
        results_name=sprintf('Pseudo_pop_peak_belief_prog_combined_time_results_more_trials_NN_%s_600ms_%d',ROI,n_tt);
    elseif this_combination==2
        results_name=sprintf('Pseudo_pop_peak_belief_prog_combined_time_results_more_trials_%s_900ms_NN_%d',ROI,n_tt);
    end
    
    %     if exist(sprintf('%s.mat',fullfile(data_path_clasifier,results_name)))==2
    
    load(fullfile(data_path_clasifier,results_name),'Classification_net_correct','Classification_net_correct_across_time');
    
    FEF.Classification_net_correct_comb(:,:,n_tt)=Classification_net_correct;
    FEF.Classification_net_correct_across_time(:,:,:,n_tt)=Classification_net_correct_across_time;
    %         FEF.W_combined(:,:,n_tt)=W;
    
    clear Classification_net_correct*
    
    %     else
    %         FEF.Classification_net_correct_combined(:,:,n_tt)=NaN;
    %         FEF.W_combined(:,:,:,n_tt)=NaN;
    %
    %         sprintf('%s %d missing',ROI,n_tt)
    %
    %     end
    
    %PFC
    ROI='PFC';
    if this_combination==1
        results_name=sprintf('Pseudo_pop_peak_belief_prog_combined_time_results_more_trials_NN_%s_600ms_%d',ROI,n_tt);
    elseif this_combination==2
        results_name=sprintf('Pseudo_pop_peak_belief_prog_combined_time_results_more_trials_%s_900ms_NN_%d',ROI,n_tt);
    end
    
    %     if exist(sprintf('%s.mat',fullfile(data_path_clasifier,results_name)))==2
    
    
    load(fullfile(data_path_clasifier,results_name),'Classification_net_correct','Classification_net_correct_across_time')
    
    PFC.Classification_net_correct_comb(:,:,n_tt)=Classification_net_correct;
    PFC.Classification_net_correct_across_time(:,:,:,n_tt)=Classification_net_correct_across_time;
    %         PFC.W_combined(:,:,n_tt)=W;
    
    clear Classification_net_correct*
    
    %     else
    %         PFC.Classification_net_correct_combined(:,:,n_tt)=NaN;
    %         PFC.W_combined(:,:,:,n_tt)=NaN;
    %         sprintf('%s %d missing',ROI,n_tt)
    %
    %     end
    
    
end

LIP_mean_classification_combined(:,:)=nanmean(LIP.Classification_net_correct_comb,2);
FEF_mean_classification_combined(:,:)=nanmean(FEF.Classification_net_correct_comb,2);
PFC_mean_classification_combined(:,:)=nanmean(PFC.Classification_net_correct_comb,2);

for tp=1:length(window_start_list)
    for n_tt=1:100
        LIP_across_time_dyn(tp,n_tt)=mean(mean(LIP.Classification_net_correct_across_time(:,:,tp,n_tt),1),2);
        FEF_across_time_dyn(tp,n_tt)=mean(mean(FEF.Classification_net_correct_across_time(:,:,tp,n_tt),1),2);
        PFC_across_time_dyn(tp,n_tt)=mean(mean(PFC.Classification_net_correct_across_time(:,:,tp,n_tt),1),2);
    end
end

LIP_across_time_diag(:,:)=mean(LIP_across_time,2);
FEF_across_time_diag(:,:)=mean(FEF_across_time,2);
PFC_across_time_diag(:,:)=mean(PFC_across_time,2);

%%

for tp=1:length(window_start_list)
    P_LIP_across_dyn(tp)=z_test_function_bootstrap(LIP_across_time_diag(tp,:)-LIP_across_time_dyn(tp,:),0);
    P_FEF_across_dyn(tp)=z_test_function_bootstrap(FEF_across_time_diag(tp,:)-FEF_across_time_dyn(tp,:),0);
    P_PFC_across_dyn(tp)=z_test_function_bootstrap(PFC_across_time_diag(tp,:)-PFC_across_time_dyn(tp,:),0);
    
    P_LIP_dyn(tp)=z_test_function_bootstrap(LIP_across_time_dyn(tp,:),1/3);
    P_FEF_dyn(tp)=z_test_function_bootstrap(FEF_across_time_dyn(tp,:),1/3);
    P_PFC_dyn(tp)=z_test_function_bootstrap(PFC_across_time_dyn(tp,:),1/3);
   
    P_LIP_diag(tp)=z_test_function_bootstrap(LIP_across_time_diag(tp,:),1/3);
    P_FEF_diag(tp)=z_test_function_bootstrap(FEF_across_time_diag(tp,:),1/3);
    P_PFC_diag(tp)=z_test_function_bootstrap(PFC_across_time_diag(tp,:),1/3);
end

%%
figure
subplot(3,1,1)
hold on
shadedErrorBar(window_start_list(8:20),mean(LIP_across_time_diag(8:20,:),2)',[prctile(LIP_across_time_diag(8:20,:),95,2),prctile(LIP_across_time_diag(8:20,:),5,2)]' ,{'color',color_for_ROI_early(1,:),'LineWidth',2},2)
shadedErrorBar(window_start_list(8:20),mean(LIP_across_time_dyn(8:20,:),2)',[prctile(LIP_across_time_dyn(8:20,:),95,2),prctile(LIP_across_time_dyn(8:20,:),5,2)]' ,{'--','color',color_for_ROI_late(1,:),'LineWidth',2},2)
plot_significance_level(window_start_list(8:20),P_LIP_diag(8:20),0.8,color_for_ROI_early(1,:),[0.01,0.05/13, 0.01/13])
plot_significance_level(window_start_list(8:20),P_LIP_dyn(8:20),0.85,color_for_ROI_late(1,:),[0.01,0.05/13, 0.01/13])

xline(-300,'-','Color',[0.5 0.5 0.5],'LineWidth',1)
xl=xline(-400,'-',{'Fixation'},'Color',[0.5 0.5 0.5],'LineWidth',1);
xl.LabelVerticalAlignment = 'bottom';
xl.LabelHorizontalAlignment = 'left';
xl.FontSize=12;

% xline(300,'-','Color',[0.25 0.75 0.25],'LineWidth',1)
% xl=xline(200,'-',{'Response'},'Color',[0.25 0.75 0.25],'LineWidth',1);
% xl.LabelVerticalAlignment = 'bottom';
% xl.LabelHorizontalAlignment = 'left';
xl.FontSize=12;
yline(1/3,'--')
ylabel({'Decoding accuracy'})
box off
xlabel({'Time to targets onset'})

subplot(3,1,2)
hold on
shadedErrorBar(window_start_list(8:20),mean(FEF_across_time_diag(8:20,:),2)',[prctile(FEF_across_time_diag(8:20,:),95,2),prctile(FEF_across_time_diag(8:20,:),5,2)]' ,{'color',color_for_ROI_early(2,:),'LineWidth',2},2)
shadedErrorBar(window_start_list(8:20),mean(FEF_across_time_dyn(8:20,:),2)',[prctile(FEF_across_time_dyn(8:20,:),95,2),prctile(FEF_across_time_dyn(8:20,:),5,2)]' ,{'--','color',color_for_ROI_late(2,:),'LineWidth',2},2)
plot_significance_level(window_start_list(8:20),P_FEF_diag(8:20),0.8,color_for_ROI_early(2,:),[0.01,0.05/13, 0.01/13])
plot_significance_level(window_start_list(8:20),P_FEF_dyn(8:20),0.85,color_for_ROI_late(2,:),[0.01,0.05/13, 0.01/13])


xline(-300,'-','Color',[0.5 0.5 0.5],'LineWidth',1)
xl=xline(-400,'-',{'Fixation'},'Color',[0.5 0.5 0.5],'LineWidth',1);
xl.LabelVerticalAlignment = 'bottom';
xl.LabelHorizontalAlignment = 'left';
xl.FontSize=12;

% xline(300,'-','Color',[0.25 0.75 0.25],'LineWidth',1)
% xl=xline(200,'-',{'Response'},'Color',[0.25 0.75 0.25],'LineWidth',1);
% xl.LabelVerticalAlignment = 'bottom';
% xl.LabelHorizontalAlignment = 'left';
xl.FontSize=12;
yline(1/3,'--')
ylabel({'Decoding accuracy'})
box off
xlabel({'Time to targets onset'})

subplot(3,1,3)
hold on
shadedErrorBar(window_start_list(8:20),mean(PFC_across_time_diag(8:20,:),2)',[prctile(PFC_across_time_diag(8:20,:),95,2),prctile(PFC_across_time_diag(8:20,:),5,2)]' ,{'color',color_for_ROI_early(3,:),'LineWidth',2},2)
shadedErrorBar(window_start_list(8:20),mean(PFC_across_time_dyn(8:20,:),2)',[prctile(PFC_across_time_dyn(8:20,:),95,2),prctile(PFC_across_time_dyn(8:20,:),5,2)]' ,{'--','color',color_for_ROI_late(3,:),'LineWidth',2},2)

plot_significance_level(window_start_list(8:20),P_PFC_diag(8:20),0.8,color_for_ROI_early(3,:),[0.01,0.05/13, 0.01/13])
plot_significance_level(window_start_list(8:20),P_PFC_dyn(8:20),0.85,color_for_ROI_late(3,:),[0.01,0.05/13, 0.01/13])

xline(-300,'-','Color',[0.5 0.5 0.5],'LineWidth',1)
xl=xline(-400,'-',{'Fixation'},'Color',[0.5 0.5 0.5],'LineWidth',1);
xl.LabelVerticalAlignment = 'bottom';
xl.LabelHorizontalAlignment = 'left';
xl.FontSize=12;

% xline(300,'-','Color',[0.25 0.75 0.25],'LineWidth',1)
% xl=xline(200,'-',{'Response'},'Color',[0.25 0.75 0.25],'LineWidth',1);
% xl.LabelVerticalAlignment = 'bottom';
% xl.LabelHorizontalAlignment = 'left';
xl.FontSize=12;
yline(1/3,'--')
ylabel({'Decoding accuracy'})
box off
xlabel({'Time to targets onset'})




%%
figure
subplot(3,1,1)
hold on
shadedErrorBar(window_start_list(8:end),mean(LIP_across_time_diag(8:end,:),2)',[prctile(LIP_across_time_diag(8:end,:),95,2),prctile(LIP_across_time_diag(8:end,:),5,2)]' ,{'color',color_for_ROI_early(1,:),'LineWidth',2},2)
shadedErrorBar(window_start_list(8:end),mean(LIP_across_time_dyn(8:end,:),2)',[prctile(LIP_across_time_dyn(8:end,:),95,2),prctile(LIP_across_time_dyn(8:end,:),5,2)]' ,{'--','color',color_for_ROI_late(1,:),'LineWidth',2},2)
plot_significance_level(window_start_list(8:end),P_LIP_diag(8:end),0.8,color_for_ROI_early(1,:),[0.01,0.05/length(window_start_list), 0.01/length(window_start_list)])
plot_significance_level(window_start_list(8:end),P_LIP_dyn(8:end),0.85,color_for_ROI_late(1,:),[0.01,0.05/length(window_start_list), 0.01/length(window_start_list)])

xline(-300,'-','Color',[0.5 0.5 0.5],'LineWidth',1)
xl=xline(-400,'-',{'Fixation'},'Color',[0.5 0.5 0.5],'LineWidth',1);
xl.LabelVerticalAlignment = 'bottom';
xl.LabelHorizontalAlignment = 'left';
xl.FontSize=12;

xline(300,'-','Color',[0.25 0.75 0.25],'LineWidth',1)
xl=xline(200,'-',{'Response'},'Color',[0.25 0.75 0.25],'LineWidth',1);
xl.LabelVerticalAlignment = 'bottom';
xl.LabelHorizontalAlignment = 'left';
xl.FontSize=12;
yline(1/3,'--')
ylabel({'Decoding accuracy'})
box off
xlabel({'Time to targets onset'})

subplot(3,1,2)
hold on
shadedErrorBar(window_start_list(8:end),mean(FEF_across_time_diag(8:end,:),2)',[prctile(FEF_across_time_diag(8:end,:),95,2),prctile(FEF_across_time_diag(8:end,:),5,2)]' ,{'color',color_for_ROI_early(2,:),'LineWidth',2},2)
shadedErrorBar(window_start_list(8:end),mean(FEF_across_time_dyn(8:end,:),2)',[prctile(FEF_across_time_dyn(8:end,:),95,2),prctile(FEF_across_time_dyn(8:end,:),5,2)]' ,{'--','color',color_for_ROI_late(2,:),'LineWidth',2},2)
plot_significance_level(window_start_list(8:end),P_FEF_diag(8:end),0.8,color_for_ROI_early(2,:),[0.01,0.05/length(window_start_list), 0.01/length(window_start_list)])
plot_significance_level(window_start_list(8:end),P_FEF_dyn(8:end),0.85,color_for_ROI_late(2,:),[0.01,0.05/length(window_start_list), 0.01/length(window_start_list)])


xline(-300,'-','Color',[0.5 0.5 0.5],'LineWidth',1)
xl=xline(-400,'-',{'Fixation'},'Color',[0.5 0.5 0.5],'LineWidth',1);
xl.LabelVerticalAlignment = 'bottom';
xl.LabelHorizontalAlignment = 'left';
xl.FontSize=12;

xline(300,'-','Color',[0.25 0.75 0.25],'LineWidth',1)
xl=xline(200,'-',{'Response'},'Color',[0.25 0.75 0.25],'LineWidth',1);
xl.LabelVerticalAlignment = 'bottom';
xl.LabelHorizontalAlignment = 'left';
xl.FontSize=12;
yline(1/3,'--')
ylabel({'Decoding accuracy'})
box off
xlabel({'Time to targets onset'})

subplot(3,1,3)
hold on
shadedErrorBar(window_start_list(8:end),mean(PFC_across_time_diag(8:end,:),2)',[prctile(PFC_across_time_diag(8:end,:),95,2),prctile(PFC_across_time_diag(8:end,:),5,2)]' ,{'color',color_for_ROI_early(3,:),'LineWidth',2},2)
shadedErrorBar(window_start_list(8:end),mean(PFC_across_time_dyn(8:end,:),2)',[prctile(PFC_across_time_dyn(8:end,:),95,2),prctile(PFC_across_time_dyn(8:end,:),5,2)]' ,{'--','color',color_for_ROI_late(3,:),'LineWidth',2},2)

plot_significance_level(window_start_list(8:end),P_PFC_diag(8:end),0.8,color_for_ROI_early(3,:),[0.01,0.05/length(window_start_list), 0.01/length(window_start_list)])
plot_significance_level(window_start_list(8:end),P_PFC_dyn(8:end),0.85,color_for_ROI_late(3,:),[0.01,0.05/length(window_start_list), 0.01/length(window_start_list)])

xline(-300,'-','Color',[0.5 0.5 0.5],'LineWidth',1)
xl=xline(-400,'-',{'Fixation'},'Color',[0.5 0.5 0.5],'LineWidth',1);
xl.LabelVerticalAlignment = 'bottom';
xl.LabelHorizontalAlignment = 'left';
xl.FontSize=12;

xline(300,'-','Color',[0.25 0.75 0.25],'LineWidth',1)
xl=xline(200,'-',{'Response'},'Color',[0.25 0.75 0.25],'LineWidth',1);
xl.LabelVerticalAlignment = 'bottom';
xl.LabelHorizontalAlignment = 'left';
xl.FontSize=12;
yline(1/3,'--')
ylabel({'Decoding accuracy'})
box off
xlabel({'Time to targets onset'})


% figure
% subplot(3,1,1)
% hold on
% shadedErrorBar(window_start_list(1:7),mean(LIP_across_time_diag(1:7,:),2)',[prctile(LIP_across_time_diag(1:7,:),95,2),prctile(LIP_across_time_diag(1:7,:),5,2)]' ,{'color',color_for_ROI_early(1,:),'LineWidth',2},2)
% shadedErrorBar(window_start_list(1:7),mean(LIP_across_time_dyn(1:7,:),2)',[prctile(LIP_across_time_dyn(1:7,:),95,2),prctile(LIP_across_time_dyn(1:7,:),5,2)]' ,{'--','color',color_for_ROI_late(1,:),'LineWidth',2},2)
% plot_significance_level(window_start_list(1:7),P_LIP_diag(1:7),0.8,color_for_ROI_early(1,:),[0.01,0.05/length(window_start_list), 0.01/length(window_start_list)])
% plot_significance_level(window_start_list(1:7),P_LIP_dyn(1:7),0.85,color_for_ROI_late(1,:),[0.01,0.05/length(window_start_list), 0.01/length(window_start_list)])
% 
% yline(1/3,'--')
% ylabel({'Decoding accuracy'})
% box off
% xlabel({'Time to reward end'})
% 
% subplot(3,1,2)
% hold on
% shadedErrorBar(window_start_list(1:7),mean(FEF_across_time_diag(1:7,:),2)',[prctile(FEF_across_time_diag(1:7,:),95,2),prctile(FEF_across_time_diag(1:7,:),5,2)]' ,{'color',color_for_ROI_early(2,:),'LineWidth',2},2)
% shadedErrorBar(window_start_list(1:7),mean(FEF_across_time_dyn(1:7,:),2)',[prctile(FEF_across_time_dyn(1:7,:),95,2),prctile(FEF_across_time_dyn(1:7,:),5,2)]' ,{'--','color',color_for_ROI_late(2,:),'LineWidth',2},2)
% plot_significance_level(window_start_list(1:7),P_FEF_diag(1:7),0.8,color_for_ROI_early(2,:),[0.01,0.05/length(window_start_list), 0.01/length(window_start_list)])
% plot_significance_level(window_start_list(1:7),P_FEF_dyn(1:7),0.85,color_for_ROI_late(2,:),[0.01,0.05/length(window_start_list), 0.01/length(window_start_list)])
% 
% 
% yline(1/3,'--')
% ylabel({'Decoding accuracy'})
% box off
% xlabel({'Time to reward end'})
% 
% subplot(3,1,3)
% hold on
% shadedErrorBar(window_start_list(1:7),mean(PFC_across_time_diag(1:7,:),2)',[prctile(PFC_across_time_diag(1:7,:),95,2),prctile(PFC_across_time_diag(1:7,:),5,2)]' ,{'color',color_for_ROI_early(3,:),'LineWidth',2},2)
% shadedErrorBar(window_start_list(1:7),mean(PFC_across_time_dyn(1:7,:),2)',[prctile(PFC_across_time_dyn(1:7,:),95,2),prctile(PFC_across_time_dyn(1:7,:),5,2)]' ,{'--','color',color_for_ROI_late(3,:),'LineWidth',2},2)
% 
% plot_significance_level(window_start_list(1:7),P_PFC_diag(1:7),0.8,color_for_ROI_early(3,:),[0.01,0.05/length(window_start_list), 0.01/length(window_start_list)])
% plot_significance_level(window_start_list(1:7),P_PFC_dyn(1:7),0.85,color_for_ROI_late(3,:),[0.01,0.05/length(window_start_list), 0.01/length(window_start_list)])
% 
% yline(1/3,'--')
% ylabel({'Decoding accuracy'})
% box off
% xlabel({'Time to reward end'})



%%
% for k=1:3
%     p_LIP_comb(k) = z_test_function_bootstrap(LIP_mean_classification_combined(k,:),1/3);
%     p_FEF_comb(k) = z_test_function_bootstrap(FEF_mean_classification_combined(k,:),1/3);
%     p_PFC_comb(k) = z_test_function_bootstrap(PFC_mean_classification_combined(k,:),1/3);
%     
% end
% %%
% for k=1:2
%     for j=k+1:3
%         p_LIP_comb2(k,j) = z_test_function_bootstrap_2_tail(LIP_mean_classification_combined(k,:)-LIP_mean_classification_combined(j,:),0);
%         p_FEF_comb2(k,j) = z_test_function_bootstrap_2_tail(FEF_mean_classification_combined(k,:)-FEF_mean_classification_combined(j,:),0);
%         p_PFC_comb2(k,j) = z_test_function_bootstrap_2_tail(PFC_mean_classification_combined(k,:)-PFC_mean_classification_combined(j,:),0);
%     end
% end
% 
% %%
% figure
% subplot(1,3,1)
% vs = violinplot(LIP_mean_classification_combined',{'1/3','2/3','3/3'},'ViolinColor',color_for_ROI(1,:),'ShowData',false)
% hold on
% yline(1/3,'--',{'Chance level'},'Color',[0.5 0.5 0.5],'LabelVerticalAlignment','bottom','LabelHorizontalAlignment','center')
% box off
% title('LIP')
% ylabel('Expected template classification accuracy')
% xlabel('Progression in block')
% ylim([0.25 1])
% subplot(1,3,2)
% 
% vs = violinplot(FEF_mean_classification_combined',{'1/3','2/3','3/3'},'ViolinColor',color_for_ROI(2,:),'ShowData',false)
% hold on
% yline(1/3,'--',{'Chance level'},'Color',[0.5 0.5 0.5],'LabelVerticalAlignment','bottom','LabelHorizontalAlignment','center')
% box off
% title('FEF')
% ylabel('Expected template classification accuracy')
% xlabel('Progression in block')
% ylim([0.25 1])
% 
% subplot(1,3,3)
% vs = violinplot(PFC_mean_classification_combined',{'1/3','2/3','3/3'},'ViolinColor',color_for_ROI(3,:),'ShowData',false)
% hold on
% yline(1/3,'--',{'Chance level'},'Color',[0.5 0.5 0.5],'LabelVerticalAlignment','bottom','LabelHorizontalAlignment','center')
% box off
% title('PFC')
% ylabel('Expected template classification accuracy')
% xlabel('Progression in block')
% ylim([0.25 1])
% 
% 
% %% now compare the two
% 
% for i=1:20
%     for j=1:3
%         for k=1:100
%             Compare_LIP(i,j,k)=LIP_across_time(i+7,j,k)-LIP_mean_classification_combined(j,k);
%             Compare_FEF(i,j,k)=FEF_across_time(i+7,j,k)-FEF_mean_classification_combined(j,k);
%             Compare_PFC(i,j,k)=PFC_across_time(i+7,j,k)-PFC_mean_classification_combined(j,k);
%             
%         end
%         p_compare_LIP(i,j)=z_test_function_bootstrap(Compare_LIP(i,j,:),0);
%         p_compare_LIP(i,j)=z_test_function_bootstrap(Compare_LIP(i,j,:),0);
%         p_compare_LIP(i,j)=z_test_function_bootstrap(Compare_LIP(i,j,:),0);
%     end
% end
% 
% 
% %%
% figure
% subplot(3,2,1)
% hold on
% shadedErrorBar(window_start_list(8:end),mean(Compare_LIP(:,1,:),3)',[prctile(Compare_LIP(:,1,:),95,3),prctile(Compare_LIP(:,1,:),5,3)]' ,{'color',color_for_ROI_early(1,:),'LineWidth',2},2)
% shadedErrorBar(window_start_list(8:end),mean(Compare_LIP(:,2,:),3)',[prctile(Compare_LIP(:,2,:),95,3),prctile(Compare_LIP(:,2,:),5,3)]',{'color',color_for_ROI(1,:),'LineWidth',2},2)
% shadedErrorBar(window_start_list(8:end),mean(Compare_LIP(:,3,:),3)',[prctile(Compare_LIP(:,3,:),95,3),prctile(Compare_LIP(:,3,:),5,3)]',{'color',color_for_ROI_late(1,:),'LineWidth',2},2)
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
% yline(0,'--')
% ylabel({'Decoding accuracy across time - combined'})
% box off
% xlabel({'Time to target'})
% ylim([-0.5 0.2])
% 
% subplot(3,2,3)
% hold on
% shadedErrorBar(window_start_list(8:end),mean(Compare_FEF(:,1,:),3)',[prctile(Compare_FEF(:,1,:),95,3),prctile(Compare_FEF(:,1,:),5,3)]',{'color',color_for_ROI_early(2,:),'LineWidth',2},2)
% shadedErrorBar(window_start_list(8:end),mean(Compare_FEF(:,2,:),3)',[prctile(Compare_FEF(:,2,:),95,3),prctile(Compare_FEF(:,2,:),5,3)]',{'color',color_for_ROI(2,:),'LineWidth',2},2)
% shadedErrorBar(window_start_list(8:end),mean(Compare_FEF(:,3,:),3)',[prctile(Compare_FEF(:,3,:),95,3),prctile(Compare_FEF(:,3,:),5,3)]',{'color',color_for_ROI_late(2,:),'LineWidth',2},2)
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
% 
% yline(0,'--')
% ylabel({'Decoding accuracy across time - combined'})
% box off
% xlabel({'Time to target'})
% ylim([-0.5 0.2])
% 
% subplot(3,2,5)
% hold on
% shadedErrorBar(window_start_list(8:end),mean(Compare_PFC(:,1,:),3)',[prctile(Compare_PFC(:,1,:),95,3),prctile(Compare_PFC(:,1,:),5,3)]',{'color',color_for_ROI_early(3,:),'LineWidth',2},2)
% shadedErrorBar(window_start_list(8:end),mean(Compare_PFC(:,2,:),3)',[prctile(Compare_PFC(:,2,:),95,3),prctile(Compare_PFC(:,2,:),5,3)]',{'color',color_for_ROI(3,:),'LineWidth',2},2)
% shadedErrorBar(window_start_list(8:end),mean(Compare_PFC(:,3,:),3)',[prctile(Compare_PFC(:,3,:),95,3),prctile(Compare_PFC(:,3,:),5,3)]',{'color',color_for_ROI_late(3,:),'LineWidth',2},2)
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
% yline(0,'--')
% ylabel({'Decoding accuracy across time - combined'})
% box off
% xlabel({'Time to target'})
% ylim([-0.5 0.2])
% 
% 
% 
% 
%%
function p = z_test_function_bootstrap(dist,null)

m = mean(dist);
s = std(dist);
z = (m-null)/s;
p=1-normcdf(z);

end

function p = z_test_function_bootstrap_2_tail(dist,null)

m = mean(dist);
s = std(dist);
z = (m-null)/s;
p=1-normcdf(abs(z))/2;

end

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





