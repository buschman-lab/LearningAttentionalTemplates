clear all;

fsroot='/Volumes/buschman';

event='target';
for i=1:27
    window_start_list(i)=-500+(i-1)*50;
end

window_size=200;
subsubtask_classifier=sprintf('Classifier_%s_%d_to_%s_%d',event,window_start_list(1),event,window_start_list(end));

N_boot=100;
task='Learning_Attentional_Templates';
subtask='exploreexploit/Restricted_FEF/Reset_RW_model';

dirstem = fullfile('Projects',task,'Analysis','Electrophy_analysis',subtask,subsubtask_classifier);
data_path_clasifier = fullfile(fsroot,dirstem);

save_name='Pseudo_pop_stim_color_no_prog_results_for_plot';

color_for_ROI=[0 0.447 0.741; 0.85 0.325 0.098; 0.929 0.694 0.125]; %LIP FEF PFC

load('colors')


color_belief(1,:)=colors(1,:);
color_belief(2,:)=colors(14,:);
color_belief(3,:)=colors(27,:);

N_loc=4;

N_prog=1;

LOAD=1;

if LOAD==1

for loc=1:N_loc
    
    %LIP
    ROI='LIP';
    
    for n_tt=1:N_boot
        
        for this_time=1:length(window_start_list)
            
            results_name=sprintf('Pseudo_pop_stim_color_no_prog_results_%s_%d_%d_%d',ROI,this_time,n_tt,loc);
            
            if exist(sprintf('%s.mat',fullfile(data_path_clasifier,results_name)))==2
                load(fullfile(data_path_clasifier,results_name),'Classification_correct','Classification_net_correct','W')
                
                
                LIP.Classification_correct(this_time,:,:,n_tt,loc)=Classification_correct(:,:);
                LIP.Classification_net_correct(this_time,:,:,n_tt,loc)=Classification_net_correct(:,:);
                
                LIP.W(this_time,:,n_tt,loc)=W(:);
                
                clear Classification* W*
                
            else
                sprintf('%s is missing',results_name)
                
            end
            
        end
    end
    
    %FEF
    ROI='FEF';
    
    for n_tt=1:N_boot
        
        for this_time=1:length(window_start_list)
            
            %         results_name=sprintf('Pseudo_pop_belief_results_%s_%d_%d_%d',ROI,this_time,n_tt,loc);
            results_name=sprintf('Pseudo_pop_stim_color_no_prog_results_%s_%d_%d_%d',ROI,this_time,n_tt,loc);
            
            if exist(sprintf('%s.mat',fullfile(data_path_clasifier,results_name)))==2
                
                load(fullfile(data_path_clasifier,results_name),'Classification_correct','Classification_net_correct','W')
                
                
                FEF.Classification_correct(this_time,:,:,n_tt,loc)=Classification_correct(:,:);
                FEF.Classification_net_correct(this_time,:,:,n_tt,loc)=Classification_net_correct(:,:);
                
                FEF.W(this_time,:,n_tt,loc)=W(:);
                
                clear Classification* W*
            else
                sprintf('%s is missing',results_name)
                
            end
            
        end
    end
    
    %PFC
    ROI='PFC';
    
    for n_tt=1:N_boot
        
        for this_time=1:length(window_start_list)
            
            %         results_name=sprintf('Pseudo_pop_belief_results_%s_%d_%d_%d',ROI,this_time,n_tt,loc);
            results_name=sprintf('Pseudo_pop_stim_color_no_prog_results_%s_%d_%d_%d',ROI,this_time,n_tt,loc);
            if exist(sprintf('%s.mat',fullfile(data_path_clasifier,results_name)))==2
                load(fullfile(data_path_clasifier,results_name),'Classification_correct','Classification_net_correct','W')
                
                
                PFC.Classification_correct(this_time,:,:,n_tt,loc)=Classification_correct(:,:);
                PFC.Classification_net_correct(this_time,:,:,n_tt,loc)=Classification_net_correct(:,:);
                
                PFC.W(this_time,:,n_tt,loc)=W(:);
                
                clear Classification* W*
            else
                sprintf('%s is missing',results_name)
                
            end
            
        end
    end
end


%%

for loc=1:N_loc
    for i=1:length(window_start_list)
        LIP_across_time(i,:,loc)=mean(LIP.Classification_correct(i,i,:,:,loc),3);
        FEF_across_time(i,:,loc)=mean(FEF.Classification_correct(i,i,:,:,loc),3);
        PFC_across_time(i,:,loc)=mean(PFC.Classification_correct(i,i,:,:,loc),3);
    end
end

for loc=1:N_loc
    for i=1:length(window_start_list)
        LIP_across_time_net(i,:,loc)=mean(LIP.Classification_net_correct(i,i,:,:,loc),3);
        FEF_across_time_net(i,:,loc)=mean(FEF.Classification_net_correct(i,i,:,:,loc),3);
        PFC_across_time_net(i,:,loc)=mean(PFC.Classification_net_correct(i,i,:,:,loc),3);
                
        p_LIP_across_time(i,loc)=z_test_function_bootstrap(mean(LIP.Classification_net_correct(i,i,:,:,loc),3),1/3);
        p_FEF_across_time(i,loc)=z_test_function_bootstrap(mean(FEF.Classification_net_correct(i,i,:,:,loc),3),1/3);
        p_PFC_across_time(i,loc)=z_test_function_bootstrap(mean(PFC.Classification_net_correct(i,i,:,:,loc),3),1/3);
    end
end

save(fullfile(data_path_clasifier,save_name))

else
    load(fullfile(data_path_clasifier,save_name))
end

%
for i=1:3
    color_for_ROI(i,:)=rgb2hsv(color_for_ROI(i,:));
    color_for_ROI(i,:)=hsv2rgb(color_for_ROI(i,1),color_for_ROI(i,2),1);
    color_for_ROI_late(i,:)=rgb2hsv(color_for_ROI(i,:));
    color_for_ROI_late(i,:)=hsv2rgb(color_for_ROI_late(i,1),color_for_ROI_late(i,2),0.55);
    
end

%%
figure

for loc=1:N_loc
    subplot(3,4,loc)
    plot(window_start_list(:),mean(LIP_across_time(:,:,loc),2),'Color',color_for_ROI(1,:),'LineWidth',2)
    hold on
    shadedErrorBar(window_start_list(:),mean(LIP_across_time(:,:,loc),2),std(LIP_across_time(:,:,loc),0,2),{'color',color_for_ROI(1,:),'LineWidth',2},2)
    hold on
    yl=yline(0.33,'--',{'Chance level'},'Color',[0 0 0],'LineWidth',1);
    hold on
    yl.LabelVerticalAlignment = 'bottom';
    yl.LabelHorizontalAlignment = 'center';
    yl.FontSize=12;
    
    hold on
    xline(-300,'-','Color',[0.5 0.5 0.5],'LineWidth',1)
    hold on
    xl=xline(-400,'-',{'Fixation'},'Color',[0.5 0.5 0.5],'LineWidth',1);
    xl.LabelVerticalAlignment = 'bottom';
    xl.LabelHorizontalAlignment = 'left';
    xl.FontSize=12;
    
    hold on
    xline(300,'-','Color',[0.25 0.75 0.25],'LineWidth',1)
    hold on
    xl=xline(200,'-',{'Response'},'Color',[0.25 0.75 0.25],'LineWidth',1);
    xl.LabelVerticalAlignment = 'bottom';
    xl.LabelHorizontalAlignment = 'left';
    xl.FontSize=12;
    ylabel('Accuracy','FontSize',14)
    xlabel(sprintf('Time to %s',event'))
    ax=gca;
    ax.XAxis.FontSize=12;
    ax.YAxis.FontSize=12;
    box off
    % legend({'Early in block','Late in block'},'FontSize',16,'Location','South')
    ylim([0.1 0.5])
    
    xlim([-600 350])
    
    
    subplot(3,4,loc+4)
    plot(window_start_list(:),mean(FEF_across_time(:,:,loc),2),'Color',color_for_ROI(2,:),'LineWidth',2)
    hold on
    shadedErrorBar(window_start_list(:),mean(FEF_across_time(:,:,loc),2),std(FEF_across_time(:,:,loc),0,2),{'color',color_for_ROI(2,:),'LineWidth',2},2)
    hold on
    yl=yline(0.33,'--',{'Chance level'},'Color',[0 0 0],'LineWidth',1);
    hold on
    
    yl.LabelVerticalAlignment = 'bottom';
    yl.LabelHorizontalAlignment = 'center';
    yl.FontSize=12;
    
    hold on
    xline(-300,'-','Color',[0.5 0.5 0.5],'LineWidth',1)
    hold on
    xl=xline(-400,'-',{'Fixation'},'Color',[0.5 0.5 0.5],'LineWidth',1);
    xl.LabelVerticalAlignment = 'bottom';
    xl.LabelHorizontalAlignment = 'left';
    xl.FontSize=12;
    
    hold on
    xline(300,'-','Color',[0.25 0.75 0.25],'LineWidth',1)
    hold on
    xl=xline(200,'-',{'Response'},'Color',[0.25 0.75 0.25],'LineWidth',1);
    xl.LabelVerticalAlignment = 'bottom';
    xl.LabelHorizontalAlignment = 'left';
    xl.FontSize=12;
    ylabel('Accuracy','FontSize',14)
    xlabel(sprintf('Time to %s',event'))
    ax=gca;
    ax.XAxis.FontSize=12;
    ax.YAxis.FontSize=12;
    box off
    % legend({'Early in block','Late in block'},'FontSize',16,'Location','South')
    ylim([0.1 0.5])
    
    xlim([-600 350])
    
    subplot(3,4,loc+8)
    plot(window_start_list(:),mean(PFC_across_time(:,:,loc),2),'Color',color_for_ROI(3,:),'LineWidth',2)
    hold on
    shadedErrorBar(window_start_list(:),mean(PFC_across_time(:,:,loc),2),std(PFC_across_time(:,:,loc),0,2),{'color',color_for_ROI(3,:),'LineWidth',2},2)
    hold on
    yl=yline(0.33,'--',{'Chance level'},'Color',[0 0 0],'LineWidth',1);
    hold on
    
    yl.LabelVerticalAlignment = 'bottom';
    yl.LabelHorizontalAlignment = 'center';
    yl.FontSize=12;
    
    hold on
    xline(-300,'-','Color',[0.5 0.5 0.5],'LineWidth',1)
    hold on
    xl=xline(-400,'-',{'Fixation'},'Color',[0.5 0.5 0.5],'LineWidth',1);
    xl.LabelVerticalAlignment = 'bottom';
    xl.LabelHorizontalAlignment = 'left';
    xl.FontSize=12;
    
    hold on
    xline(300,'-','Color',[0.25 0.75 0.25],'LineWidth',1)
    hold on
    xl=xline(200,'-',{'Response'},'Color',[0.25 0.75 0.25],'LineWidth',1);
    xl.LabelVerticalAlignment = 'bottom';
    xl.LabelHorizontalAlignment = 'left';
    xl.FontSize=12;
    ylabel('Accuracy','FontSize',14)
    xlabel(sprintf('Time to %s',event'))
    ax=gca;
    ax.XAxis.FontSize=12;
    ax.YAxis.FontSize=12;
    box off
    % legend({'Early in block','Late in block'},'FontSize',16,'Location','South')
    ylim([0.1 0.5])
    
    xlim([-600 350])
end


%%
figure

for loc=1:N_loc
    subplot(3,4,loc)
    hold on
    shadedErrorBar(window_start_list(:),mean(LIP_across_time_net(:,:,loc),2),[prctile(LIP_across_time_net(:,:,loc),95,2),prctile(LIP_across_time_net(:,:,loc),5,2)]',{'color',color_for_ROI(1,:),'LineWidth',2},2)
    yl=yline(0.33,'--',{'Chance level'},'Color',[0 0 0],'LineWidth',1); 
    yl.LabelVerticalAlignment = 'bottom';
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
    ylabel('Stim color decoding accuracy','FontSize',14)
    xlabel(sprintf('Time to %s',event'))
    ax=gca;
    ax.XAxis.FontSize=12;
    ax.YAxis.FontSize=12;
    box off
    % legend({'Early in block','Late in block'},'FontSize',16,'Location','South')
    ylim([0.1 0.5])
    
xlim([-600 window_start_list(end)+150])
xticks([-600:200:window_start_list(end)+150])
    
    
    subplot(3,4,loc+4)
    hold on
    shadedErrorBar(window_start_list(:),mean(FEF_across_time_net(:,:,loc),2),[prctile(FEF_across_time_net(:,:,loc),95,2),prctile(FEF_across_time_net(:,:,loc),5,2)]',{'color',color_for_ROI(2,:),'LineWidth',2},2)
    yl=yline(0.33,'--',{'Chance level'},'Color',[0 0 0],'LineWidth',1);
    yl.LabelVerticalAlignment = 'bottom';
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
    ylabel('Stim color decoding accuracy','FontSize',14)
    xlabel(sprintf('Time to %s',event'))
    ax=gca;
    ax.XAxis.FontSize=12;
    ax.YAxis.FontSize=12;
    box off
    % legend({'Early in block','Late in block'},'FontSize',16,'Location','South')
    ylim([0.1 0.5])
    
xlim([-600 window_start_list(end)+150])
xticks([-600:200:window_start_list(end)+150])
    
    subplot(3,4,loc+8)
        hold on
    shadedErrorBar(window_start_list(:),mean(PFC_across_time_net(:,:,loc),2),[prctile(PFC_across_time_net(:,:,loc),95,2),prctile(PFC_across_time_net(:,:,loc),5,2)]',{'color',color_for_ROI(3,:),'LineWidth',2},2)
    yl=yline(0.33,'--',{'Chance level'},'Color',[0 0 0],'LineWidth',1);    
    yl.LabelVerticalAlignment = 'bottom';
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
    ylabel('Stim color decoding accuracy','FontSize',14)
    xlabel(sprintf('Time to %s',event'))
    ax=gca;
    ax.XAxis.FontSize=12;
    ax.YAxis.FontSize=12;
    box off
    % legend({'Early in block','Late in block'},'FontSize',16,'Location','South')
    ylim([0.1 0.5])
    
xlim([-600 window_start_list(end)+150])
xticks([-600:200:window_start_list(end)+150])
end




%%

%
%
%
% for loc=1:N_loc
% figure
% for k=1:1
% subplot(3,3,1+3*(k-1))
% imagesc(mean(mean(LIP.Classification_correct(:,:,:,:,loc),3),4),[0.3 1])
% colorbar
% colormap hot
% ax=gca;
%
% ax.XAxis.FontSize=12;
% ax.YAxis.FontSize=12;
%
% title('LIP')
%
% subplot(3,3,2+3*(k-1))
% imagesc(mean(mean(FEF.Classification_correct(:,:,:,:,loc),3),4),[0.3 1])
% colorbar
% colormap hot
%
% title('FEF')
% ax=gca;
% ax.XAxis.FontSize=12;
% ax.YAxis.FontSize=12;
%
% subplot(3,3,3+3*(k-1))
% imagesc(mean(mean(PFC.Classification_correct(:,:,:,:,loc),3),4),[0.3 1])
% colorbar
% colormap hot
%
% title('PFC')
% ax=gca;
% ax.XAxis.FontSize=12;
% ax.YAxis.FontSize=12;
% end
% end


%%


function p = z_test_function_bootstrap(dist,null)

m = mean(dist);
s = std(dist);
z = (m-null)/s;
p=1-normcdf(z);

end