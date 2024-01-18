clear all;

fsroot='/Volumes/buschman';

for i=1:23
    initial_window=-500;
    event_list{i}='response';
    window_start_list(i)=initial_window+(i-1)*50;
end

event='response';

task='Learning_Attentional_Templates';
subtask='exploreexploit/Restricted_FEF/Reset_RW_model';

window_size=200;

subsubtask_classifier=sprintf('Classifier_%s_%d_to_%s_%d',event_list{1},window_start_list(1),event_list{end},window_start_list(end));
% subsubtask_classifier=sprintf('Classifier_%s_%d_to_%s_%d',event,window_start_list(1),event,window_start_list(end));

dirstem = fullfile('Projects',task,'Analysis','Electrophy_analysis',subtask,subsubtask_classifier);
data_path_clasifier = fullfile(fsroot,dirstem);

data_name='Pseudo_pop_choice_in_context_classifier';


dirstem = fullfile('Projects',task,'Analysis','Electrophy_analysis',subtask,subsubtask_classifier);
data_path_clasifier = fullfile(fsroot,dirstem);

save_name='Pseudo_pop_choice_in_context_results_for_plot';

color_for_ROI=[0 0.447 0.741; 0.85 0.325 0.098; 0.929 0.694 0.125]; %LIP FEF PFC

load('colors')

N_boot=100;

color_belief(1,:)=colors(1,:);
color_belief(2,:)=colors(14,:);
color_belief(3,:)=colors(27,:);

N_context=4;

N_prog=1;

LOAD=1;

if LOAD==1

    
    %LIP
    ROI='LIP';
    
    for n_tt=1:N_boot
        
        for this_time=1:length(window_start_list)
            
            results_name=sprintf('Pseudo_pop_choice_in_context_results_%s_%d_%d',ROI,this_time,n_tt);

            if exist(sprintf('%s.mat',fullfile(data_path_clasifier,results_name)))==2
                load(fullfile(data_path_clasifier,results_name),'Classification_correct')
                
                
                LIP.Classification_correct(this_time,:,:,:,n_tt)=Classification_correct;
                                
                clear Classification*
                
            else
                sprintf('%s is missing',results_name)
                
            end
            
        end
    end
    
    %FEF
    ROI='FEF';
    
    for n_tt=1:N_boot
        
        for this_time=1:length(window_start_list)
            
            results_name=sprintf('Pseudo_pop_choice_in_context_results_%s_%d_%d',ROI,this_time,n_tt);
            
            if exist(sprintf('%s.mat',fullfile(data_path_clasifier,results_name)))==2
                
                load(fullfile(data_path_clasifier,results_name),'Classification_correct')
                
                
                FEF.Classification_correct(this_time,:,:,:,n_tt)=Classification_correct;
                
                clear Classification*
            else
                sprintf('%s is missing',results_name)
                
            end
            
        end
    end
    
    %PFC
    ROI='PFC';
    
    for n_tt=1:N_boot
        
        for this_time=1:length(window_start_list)
            
            results_name=sprintf('Pseudo_pop_choice_in_context_results_%s_%d_%d',ROI,this_time,n_tt);
            
            if exist(sprintf('%s.mat',fullfile(data_path_clasifier,results_name)))==2
                load(fullfile(data_path_clasifier,results_name),'Classification_correct')
                
                PFC.Classification_correct(this_time,:,:,:,n_tt)=Classification_correct;
                
                clear Classification*
            else
                sprintf('%s is missing',results_name)
                
            end
            
        end
    end


%%

for loc=1:N_context
    for i=1:length(window_start_list)
        LIP_across_time(i,:,loc)=mean(LIP.Classification_correct(i,i,:,loc,:),3);
        FEF_across_time(i,:,loc)=mean(FEF.Classification_correct(i,i,:,loc,:),3);
        PFC_across_time(i,:,loc)=mean(PFC.Classification_correct(i,i,:,loc,:),3);
    end
end

for loc=1:N_context
    for i=1:length(window_start_list)
        p_LIP_across_time(i,loc)=z_test_function_bootstrap(mean(LIP.Classification_correct(i,i,:,loc,:),3),1/3);
        p_FEF_across_time(i,loc)=z_test_function_bootstrap(mean(FEF.Classification_correct(i,i,:,loc,:),3),1/3);
        p_PFC_across_time(i,loc)=z_test_function_bootstrap(mean(PFC.Classification_correct(i,i,:,loc,:),3),1/3);
    end
end

save(fullfile(data_path_clasifier,save_name))

else
    load(fullfile(data_path_clasifier,save_name))
end

%%
for i=1:3
    color_for_ROI(i,:)=rgb2hsv(color_for_ROI(i,:));
    color_for_ROI(i,:)=hsv2rgb(color_for_ROI(i,1),color_for_ROI(i,2),1);
    color_for_ROI_late(i,:)=rgb2hsv(color_for_ROI(i,:));
    color_for_ROI_late(i,:)=hsv2rgb(color_for_ROI_late(i,1),color_for_ROI_late(i,2),0.55);
    
end

%%
figure

for loc=1:N_context
    subplot(3,4,loc)
    hold on
    shadedErrorBar(window_start_list(:),mean(LIP_across_time(:,:,loc),2),std(LIP_across_time(:,:,loc),0,2),{'color',color_for_ROI(1,:),'LineWidth',2},2)
    plot_significance_level(window_start_list,p_LIP_across_time(:,loc),1,color_for_ROI(1,:),[0.01,0.05/(length(window_start_list))/4, 0.01/(length(window_start_list))/4])
    yl=yline(1/3,'--',{'Chance level'},'Color',[0 0 0],'LineWidth',1);
    yl.LabelVerticalAlignment = 'bottom';
    yl.LabelHorizontalAlignment = 'center';
    yl.FontSize=12;
    
    xl.FontSize=12;
    ylabel('Accuracy','FontSize',14)
    xlabel(sprintf('Time to %s',event'))
    ax=gca;
    ax.XAxis.FontSize=12;
    ax.YAxis.FontSize=12;
    box off
    % legend({'Early in block','Late in block'},'FontSize',16,'Location','South')
%     ylim([0.1 0.5])
    
    xlim([-600 700])
    
    
    subplot(3,4,loc+4)
    hold on
    shadedErrorBar(window_start_list(:),mean(FEF_across_time(:,:,loc),2),std(FEF_across_time(:,:,loc),0,2),{'color',color_for_ROI(2,:),'LineWidth',2},2)
    plot_significance_level(window_start_list,p_FEF_across_time(:,loc),1,color_for_ROI(2,:),[0.01,0.05/(length(window_start_list))/4, 0.01/(length(window_start_list))/4])
        yl=yline(1/3,'--',{'Chance level'},'Color',[0 0 0],'LineWidth',1);
    yl.LabelVerticalAlignment = 'bottom';
    yl.LabelHorizontalAlignment = 'center';
    yl.FontSize=12;

    ylabel('Accuracy','FontSize',14)
    xlabel(sprintf('Time to %s',event'))
    ax=gca;
    ax.XAxis.FontSize=12;
    ax.YAxis.FontSize=12;
    box off
    % legend({'Early in block','Late in block'},'FontSize',16,'Location','South')
%     ylim([0.1 0.5])
    
    xlim([-600 700])
    
    subplot(3,4,loc+8)
    shadedErrorBar(window_start_list(:),mean(PFC_across_time(:,:,loc),2),std(PFC_across_time(:,:,loc),0,2),{'color',color_for_ROI(3,:),'LineWidth',2},2)
    plot_significance_level(window_start_list,p_PFC_across_time(:,loc),1,color_for_ROI(3,:),[0.01,0.05/(length(window_start_list))/4, 0.01/(length(window_start_list))/4])
        yl=yline(1/3,'--',{'Chance level'},'Color',[0 0 0],'LineWidth',1);
    yl.LabelVerticalAlignment = 'bottom';
    yl.LabelHorizontalAlignment = 'center';
    yl.FontSize=12;

    ylabel('Accuracy','FontSize',14)
    xlabel(sprintf('Time to %s',event'))
    ax=gca;
    ax.XAxis.FontSize=12;
    ax.YAxis.FontSize=12;
    box off
    % legend({'Early in block','Late in block'},'FontSize',16,'Location','South')
%     ylim([0.1 0.5])
    
    xlim([-600 700])
end




%%

%
%
%
% for loc=1:N_loc
% figure
% for k=1:1
% subplot(3,3,1+3*(k-1))
% imagesc(mean(mean(LIP.Classification_correct(:,:,:,:),3),4),[0.3 1])
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
% imagesc(mean(mean(FEF.Classification_correct(:,:,:,:),3),4),[0.3 1])
% colorbar
% colormap hot
%
% title('FEF')
% ax=gca;
% ax.XAxis.FontSize=12;
% ax.YAxis.FontSize=12;
%
% subplot(3,3,3+3*(k-1))
% imagesc(mean(mean(PFC.Classification_correct(:,:,:,:),3),4),[0.3 1])
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
    if p(end-1)>=this_thr && p(end)<this_thr
        plot(x(end),a,'.','Color',c,'LineWidth',b)
    end
end


end