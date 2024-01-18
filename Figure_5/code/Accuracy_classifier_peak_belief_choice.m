clear all;

fsroot='/Volumes/buschman';

% event_list={'reward_end','reward_end','fixation','target','target','response'};
% event='target';

% switch event
%     case 'target'
%         for i=1:20
%             window_start_list(i)=-600+(i-1)*50;
%         end
%     case 'reward_end'
%         for i=1:7
%             window_start_list(i)=0+(i-1)*50;
%         end
% end
%
% window_size=300;
%
% subsubtask_classifier=sprintf('Classifier_%s_%d_to_%s_%d',event,window_start_list(1),event,window_start_list(end));
%

event='response';

for i=1:23
    initial_window=-500;
    event_list{i}='response';
    window_start_list(i)=initial_window+(i-1)*50;
end

subsubtask_classifier=sprintf('Classifier_%s_%d_to_%s_%d',event,window_start_list(1),event,window_start_list(end));


N_boot=100;
task='Learning_Attentional_Templates';
subtask='exploreexploit/Restricted_FEF/Reset_RW_model';

dirstem = fullfile('Projects',task,'Analysis','Electrophy_analysis',subtask,subsubtask_classifier);
data_path_clasifier = fullfile(fsroot,dirstem);

color_for_ROI=[0 0.447 0.741; 0.85 0.325 0.098; 0.929 0.694 0.125]; %LIP FEF PFC

load('colors')

color_belief(1,:)=colors(1,:);
color_belief(2,:)=colors(14,:);
color_belief(3,:)=colors(27,:);

N_loc=4;

N_prog=1;


save_name='Choice_peak_template_classifier_results_plot';

LOAD=0;

if LOAD==1

%%
for loc=1:N_loc
    
    %LIP
    ROI='LIP';
    
    for n_tt=1:N_boot
        
        for this_time=1:length(window_start_list)
            
            %         results_name=sprintf('Pseudo_pop_belief_choice_results_%s_%d_%d_%d',ROI,this_time,n_tt,loc);
            %         results_name=sprintf('Pseudo_pop_belief_choice_no_prog_results_%s_%d_%d_%d',ROI,this_time,n_tt,loc);
            results_name=sprintf('Pseudo_pop_choice_peak_template_classifier_results_%s_%d_%d_%d',ROI,this_time,n_tt,loc);
            if exist(sprintf('%s.mat',fullfile(data_path_clasifier,results_name)))==2
                
                load(fullfile(data_path_clasifier,results_name),'Classification_belief_correct','Classification_choice_each_belief_correct','Classification_choice_correct','W_choice','W','W_choice_each_belief')
                
                for k=1:N_prog
                    LIP.Classification_correct_belief(this_time,:,k,n_tt,loc)=Classification_belief_correct(:,k);
                    
                    LIP.Classification_correct_choice(this_time,:,k,n_tt,loc)=Classification_choice_correct(:,k);
                    
                    LIP.Classification_correct_choice_each_belief(this_time,:,:,:,k,n_tt,loc)=Classification_choice_each_belief_correct(:,:,:,k);
                end
                LIP.W_template(this_time,:,:,n_tt,loc)=W(:,:);
                LIP.W_choice(this_time,:,n_tt,loc)=W_choice(:);
                LIP.W_choice_each_belief(this_time,:,:,n_tt,loc)=W_choice_each_belief(:,:);
                
                clear Classification_belief_correct Classification_choice_correct Classification_choice_each_belief_correct W*
            else
                sprintf('%s is missing',results_name)
                
            end
            
        end
    end
    
    %FEF
    ROI='FEF';
    
    for n_tt=1:N_boot
        
        for this_time=1:length(window_start_list)
            
            %         results_name=sprintf('Pseudo_pop_belief_choice_results_%s_%d_%d_%d',ROI,this_time,n_tt,loc);
            %         results_name=sprintf('Pseudo_pop_belief_choice_no_prog_results_%s_%d_%d_%d',ROI,this_time,n_tt,loc);
            results_name=sprintf('Pseudo_pop_choice_peak_template_classifier_results_%s_%d_%d_%d',ROI,this_time,n_tt,loc);
            if exist(sprintf('%s.mat',fullfile(data_path_clasifier,results_name)))==2
                load(fullfile(data_path_clasifier,results_name),'Classification_belief_correct','Classification_choice_each_belief_correct','Classification_choice_correct','W_choice','W','W_choice_each_belief')
                
                for k=1:N_prog
                    FEF.Classification_correct_belief(this_time,:,k,n_tt,loc)=Classification_belief_correct(:,k);
                    
                    FEF.Classification_correct_choice(this_time,:,k,n_tt,loc)=Classification_choice_correct(:,k);
                    
                    FEF.Classification_correct_choice_each_belief(this_time,:,:,:,k,n_tt,loc)=Classification_choice_each_belief_correct(:,:,:,k);
                end
                FEF.W_template(this_time,:,:,n_tt,loc)=W(:,:);
                FEF.W_choice(this_time,:,n_tt,loc)=W_choice(:);
                FEF.W_choice_each_belief(this_time,:,:,n_tt,loc)=W_choice_each_belief(:,:);
                
                clear Classification_belief_correct Classification_choice_correct Classification_choice_each_belief_correct W*
            else
                sprintf('%s is missing',results_name)
                
            end
        end
    end
    
    %PFC
    ROI='PFC';
    
    for n_tt=1:N_boot
        
        for this_time=1:length(window_start_list)
            
            %         results_name=sprintf('Pseudo_pop_belief_choice_results_%s_%d_%d_%d',ROI,this_time,n_tt,loc);
            %         results_name=sprintf('Pseudo_pop_belief_choice_no_prog_results_%s_%d_%d_%d',ROI,this_time,n_tt,loc);
            results_name=sprintf('Pseudo_pop_choice_peak_template_classifier_results_%s_%d_%d_%d',ROI,this_time,n_tt,loc);
            if exist(sprintf('%s.mat',fullfile(data_path_clasifier,results_name)))==2
                load(fullfile(data_path_clasifier,results_name),'Classification_belief_correct','Classification_choice_each_belief_correct','Classification_choice_correct','W_choice','W','W_choice_each_belief')
                
                for k=1:N_prog
                    PFC.Classification_correct_belief(this_time,:,k,n_tt,loc)=Classification_belief_correct(:,k);
                    
                    PFC.Classification_correct_choice(this_time,:,k,n_tt,loc)=Classification_choice_correct(:,k);
                    
                    PFC.Classification_correct_choice_each_belief(this_time,:,:,:,k,n_tt,loc)=Classification_choice_each_belief_correct(:,:,:,k);
                end
                PFC.W_template(this_time,:,:,n_tt,loc)=W(:,:);
                PFC.W_choice(this_time,:,n_tt,loc)=W_choice(:);
                PFC.W_choice_each_belief(this_time,:,:,n_tt,loc)=W_choice_each_belief(:,:);
                
                clear Classification_belief_correct Classification_choice_correct Classification_choice_each_belief_correct W*
            else
                sprintf('%s is missing',results_name)
                
            end
        end
    end
end

%%
for loc=1:N_loc
    for i=1:length(window_start_list)
        for k=1:N_prog
            p_belief_LIP(i,k,loc)=z_test_function_bootstrap(mean(LIP.Classification_correct_belief(i,:,k,:,loc),2),1/3);
            p_belief_FEF(i,k,loc)=z_test_function_bootstrap(mean(FEF.Classification_correct_belief(i,:,k,:,loc),2),1/3);
            p_belief_PFC(i,k,loc)=z_test_function_bootstrap(mean(PFC.Classification_correct_belief(i,:,k,:,loc),2),1/3);
            
            p_choice_LIP(i,k,loc)=z_test_function_bootstrap(mean(LIP.Classification_correct_choice(i,:,k,:,loc),2),1/2);
            p_choice_FEF(i,k,loc)=z_test_function_bootstrap(mean(FEF.Classification_correct_choice(i,:,k,:,loc),2),1/2);
            p_choice_PFC(i,k,loc)=z_test_function_bootstrap(mean(PFC.Classification_correct_choice(i,:,k,:,loc),2),1/2);
            for l=1:3
                for m=1:3
                    p_choice_each_belief_LIP(i,l,m,k,loc)=z_test_function_bootstrap(mean(LIP.Classification_correct_choice_each_belief(i,:,l,m,k,:,loc),2),1/2);
                    p_choice_each_belief_FEF(i,l,m,k,loc)=z_test_function_bootstrap(mean(FEF.Classification_correct_choice_each_belief(i,:,l,m,k,:,loc),2),1/2);
                    p_choice_each_belief_PFC(i,l,m,k,loc)=z_test_function_bootstrap(mean(PFC.Classification_correct_choice_each_belief(i,:,l,m,k,:,loc),2),1/2);
                end
            end
        end
    end
end
%%

LIP.Same_belief_choice=horzcat(reshape(LIP.Classification_correct_choice_each_belief(:,:,1,1,1,:,:),size(LIP.Classification_correct_choice_each_belief,1),size(LIP.Classification_correct_choice_each_belief,2),size(LIP.Classification_correct_choice_each_belief,6),size(LIP.Classification_correct_choice_each_belief,7)),...
    reshape(LIP.Classification_correct_choice_each_belief(:,:,2,2,1,:,:),size(LIP.Classification_correct_choice_each_belief,1),size(LIP.Classification_correct_choice_each_belief,2),size(LIP.Classification_correct_choice_each_belief,6),size(LIP.Classification_correct_choice_each_belief,7)),...
    reshape(LIP.Classification_correct_choice_each_belief(:,:,3,3,1,:,:),size(LIP.Classification_correct_choice_each_belief,1),size(LIP.Classification_correct_choice_each_belief,2),size(LIP.Classification_correct_choice_each_belief,6),size(LIP.Classification_correct_choice_each_belief,7)));

LIP.Other_belief_choice=horzcat(reshape(LIP.Classification_correct_choice_each_belief(:,:,1,2,1,:,:),size(LIP.Classification_correct_choice_each_belief,1),size(LIP.Classification_correct_choice_each_belief,2),size(LIP.Classification_correct_choice_each_belief,6),size(LIP.Classification_correct_choice_each_belief,7)),...
    reshape(LIP.Classification_correct_choice_each_belief(:,:,1,3,1,:,:),size(LIP.Classification_correct_choice_each_belief,1),size(LIP.Classification_correct_choice_each_belief,2),size(LIP.Classification_correct_choice_each_belief,6),size(LIP.Classification_correct_choice_each_belief,7)),...
    reshape(LIP.Classification_correct_choice_each_belief(:,:,2,1,1,:,:),size(LIP.Classification_correct_choice_each_belief,1),size(LIP.Classification_correct_choice_each_belief,2),size(LIP.Classification_correct_choice_each_belief,6),size(LIP.Classification_correct_choice_each_belief,7)),...
    reshape(LIP.Classification_correct_choice_each_belief(:,:,2,3,1,:,:),size(LIP.Classification_correct_choice_each_belief,1),size(LIP.Classification_correct_choice_each_belief,2),size(LIP.Classification_correct_choice_each_belief,6),size(LIP.Classification_correct_choice_each_belief,7)),...
    reshape(LIP.Classification_correct_choice_each_belief(:,:,3,1,1,:,:),size(LIP.Classification_correct_choice_each_belief,1),size(LIP.Classification_correct_choice_each_belief,2),size(LIP.Classification_correct_choice_each_belief,6),size(LIP.Classification_correct_choice_each_belief,7)),...
    reshape(LIP.Classification_correct_choice_each_belief(:,:,3,2,1,:,:),size(LIP.Classification_correct_choice_each_belief,1),size(LIP.Classification_correct_choice_each_belief,2),size(LIP.Classification_correct_choice_each_belief,6),size(LIP.Classification_correct_choice_each_belief,7)));

FEF.Same_belief_choice=horzcat(reshape(FEF.Classification_correct_choice_each_belief(:,:,1,1,1,:,:),size(FEF.Classification_correct_choice_each_belief,1),size(FEF.Classification_correct_choice_each_belief,2),size(FEF.Classification_correct_choice_each_belief,6),size(FEF.Classification_correct_choice_each_belief,7)),...
    reshape(FEF.Classification_correct_choice_each_belief(:,:,2,2,1,:,:),size(FEF.Classification_correct_choice_each_belief,1),size(FEF.Classification_correct_choice_each_belief,2),size(FEF.Classification_correct_choice_each_belief,6),size(FEF.Classification_correct_choice_each_belief,7)),...
    reshape(FEF.Classification_correct_choice_each_belief(:,:,3,3,1,:,:),size(FEF.Classification_correct_choice_each_belief,1),size(FEF.Classification_correct_choice_each_belief,2),size(FEF.Classification_correct_choice_each_belief,6),size(FEF.Classification_correct_choice_each_belief,7)));

FEF.Other_belief_choice=horzcat(reshape(FEF.Classification_correct_choice_each_belief(:,:,1,2,1,:,:),size(FEF.Classification_correct_choice_each_belief,1),size(FEF.Classification_correct_choice_each_belief,2),size(FEF.Classification_correct_choice_each_belief,6),size(FEF.Classification_correct_choice_each_belief,7)),...
    reshape(FEF.Classification_correct_choice_each_belief(:,:,1,3,1,:,:),size(FEF.Classification_correct_choice_each_belief,1),size(FEF.Classification_correct_choice_each_belief,2),size(FEF.Classification_correct_choice_each_belief,6),size(FEF.Classification_correct_choice_each_belief,7)),...
    reshape(FEF.Classification_correct_choice_each_belief(:,:,2,1,1,:,:),size(FEF.Classification_correct_choice_each_belief,1),size(FEF.Classification_correct_choice_each_belief,2),size(FEF.Classification_correct_choice_each_belief,6),size(FEF.Classification_correct_choice_each_belief,7)),...
    reshape(FEF.Classification_correct_choice_each_belief(:,:,2,3,1,:,:),size(FEF.Classification_correct_choice_each_belief,1),size(FEF.Classification_correct_choice_each_belief,2),size(FEF.Classification_correct_choice_each_belief,6),size(FEF.Classification_correct_choice_each_belief,7)),...
    reshape(FEF.Classification_correct_choice_each_belief(:,:,3,1,1,:,:),size(FEF.Classification_correct_choice_each_belief,1),size(FEF.Classification_correct_choice_each_belief,2),size(FEF.Classification_correct_choice_each_belief,6),size(FEF.Classification_correct_choice_each_belief,7)),...
    reshape(FEF.Classification_correct_choice_each_belief(:,:,3,2,1,:,:),size(FEF.Classification_correct_choice_each_belief,1),size(FEF.Classification_correct_choice_each_belief,2),size(FEF.Classification_correct_choice_each_belief,6),size(FEF.Classification_correct_choice_each_belief,7)));

PFC.Same_belief_choice=horzcat(reshape(PFC.Classification_correct_choice_each_belief(:,:,1,1,1,:,:),size(PFC.Classification_correct_choice_each_belief,1),size(PFC.Classification_correct_choice_each_belief,2),size(PFC.Classification_correct_choice_each_belief,6),size(PFC.Classification_correct_choice_each_belief,7)),...
    reshape(PFC.Classification_correct_choice_each_belief(:,:,2,2,1,:,:),size(PFC.Classification_correct_choice_each_belief,1),size(PFC.Classification_correct_choice_each_belief,2),size(PFC.Classification_correct_choice_each_belief,6),size(PFC.Classification_correct_choice_each_belief,7)),...
    reshape(PFC.Classification_correct_choice_each_belief(:,:,3,3,1,:,:),size(PFC.Classification_correct_choice_each_belief,1),size(PFC.Classification_correct_choice_each_belief,2),size(PFC.Classification_correct_choice_each_belief,6),size(PFC.Classification_correct_choice_each_belief,7)));

PFC.Other_belief_choice=horzcat(reshape(PFC.Classification_correct_choice_each_belief(:,:,1,2,1,:,:),size(PFC.Classification_correct_choice_each_belief,1),size(PFC.Classification_correct_choice_each_belief,2),size(PFC.Classification_correct_choice_each_belief,6),size(PFC.Classification_correct_choice_each_belief,7)),...
    reshape(PFC.Classification_correct_choice_each_belief(:,:,1,3,1,:,:),size(PFC.Classification_correct_choice_each_belief,1),size(PFC.Classification_correct_choice_each_belief,2),size(PFC.Classification_correct_choice_each_belief,6),size(PFC.Classification_correct_choice_each_belief,7)),...
    reshape(PFC.Classification_correct_choice_each_belief(:,:,2,1,1,:,:),size(PFC.Classification_correct_choice_each_belief,1),size(PFC.Classification_correct_choice_each_belief,2),size(PFC.Classification_correct_choice_each_belief,6),size(PFC.Classification_correct_choice_each_belief,7)),...
    reshape(PFC.Classification_correct_choice_each_belief(:,:,2,3,1,:,:),size(PFC.Classification_correct_choice_each_belief,1),size(PFC.Classification_correct_choice_each_belief,2),size(PFC.Classification_correct_choice_each_belief,6),size(PFC.Classification_correct_choice_each_belief,7)),...
    reshape(PFC.Classification_correct_choice_each_belief(:,:,3,1,1,:,:),size(PFC.Classification_correct_choice_each_belief,1),size(PFC.Classification_correct_choice_each_belief,2),size(PFC.Classification_correct_choice_each_belief,6),size(PFC.Classification_correct_choice_each_belief,7)),...
    reshape(PFC.Classification_correct_choice_each_belief(:,:,3,2,1,:,:),size(PFC.Classification_correct_choice_each_belief,1),size(PFC.Classification_correct_choice_each_belief,2),size(PFC.Classification_correct_choice_each_belief,6),size(PFC.Classification_correct_choice_each_belief,7)));

%%
for loc=1:4
    for t=1:size(LIP.Same_belief_choice,1)
        p_LIP_same(t,loc)=z_test_function_bootstrap(reshape(mean(LIP.Same_belief_choice(t,:,:,loc),2),size(LIP.Same_belief_choice,3),1),1/2);
        p_LIP_other(t,loc)=z_test_function_bootstrap(reshape(mean(LIP.Other_belief_choice(t,:,:,loc),2),size(LIP.Same_belief_choice,3),1),1/2);
        this_dist(:)=mean(LIP.Same_belief_choice(t,:,:,loc),2)-mean(LIP.Other_belief_choice(t,:,:,loc),2);
        p_LIP_same_other(t,loc) = z_test_function_bootstrap(this_dist,0);
        
        
        p_FEF_same(t,loc)=z_test_function_bootstrap(reshape(mean(FEF.Same_belief_choice(t,:,:,loc),2),size(FEF.Same_belief_choice,3),1),1/2);
        p_FEF_other(t,loc)=z_test_function_bootstrap(reshape(mean(FEF.Other_belief_choice(t,:,:,loc),2),size(FEF.Same_belief_choice,3),1),1/2);
        this_dist(:)=mean(FEF.Same_belief_choice(t,:,:,loc),2)-mean(FEF.Other_belief_choice(t,:,:,loc),2);
        p_FEF_same_other(t,loc) = z_test_function_bootstrap(this_dist,0);
        
        p_PFC_same(t,loc)=z_test_function_bootstrap(reshape(mean(PFC.Same_belief_choice(t,:,:,loc),2),size(PFC.Same_belief_choice,3),1),1/2);
        p_PFC_other(t,loc)=z_test_function_bootstrap(reshape(mean(PFC.Other_belief_choice(t,:,:,loc),2),size(PFC.Same_belief_choice,3),1),1/2);
        this_dist(:)=mean(PFC.Same_belief_choice(t,:,:,loc),2)-mean(PFC.Other_belief_choice(t,:,:,loc),2);
        p_PFC_same_other(t,loc) = z_test_function_bootstrap(this_dist,0);
        
     end
end



save(fullfile(data_path_clasifier,save_name))

else
    
load(fullfile(data_path_clasifier,save_name))

end

%%

for loc=1:4
    for t=1:size(LIP.Same_belief_choice,1)
        this_dist(:)=mean(LIP.Same_belief_choice(t,:,:,loc),2)-mean(LIP.Other_belief_choice(t,:,:,loc),2);
        [bf_LIP_same_other(t,loc), ~] = bf.ttest(this_dist);
        
        
        this_dist(:)=mean(FEF.Same_belief_choice(t,:,:,loc),2)-mean(FEF.Other_belief_choice(t,:,:,loc),2);
        [bf_FEF_same_other(t,loc), ~] = bf.ttest(this_dist);
        
        this_dist(:)=mean(PFC.Same_belief_choice(t,:,:,loc),2)-mean(PFC.Other_belief_choice(t,:,:,loc),2);
        [bf_PFC_same_other(t,loc), ~] = bf.ttest(this_dist,0);
        
     end
end

%%
for i=1:3
    color_for_ROI_early(i,:)=rgb2hsv(color_for_ROI(i,:));
    color_for_ROI_early(i,:)=hsv2rgb(color_for_ROI_early(i,1),color_for_ROI_early(i,2),1);
    color_for_ROI_late(i,:)=rgb2hsv(color_for_ROI(i,:));
    color_for_ROI_late(i,:)=hsv2rgb(color_for_ROI_late(i,1),color_for_ROI_late(i,2),0.55);
    
end

%%
figure

for loc=1:N_loc
    subplot(3,4,loc)
    hold on
    
    plot(window_start_list,mean(mean(LIP.Classification_correct_choice(:,:,1,:,loc),2),4),'Color',color_for_ROI(1,:),'LineWidth',1)
    
    shadedErrorBar(window_start_list,mean(mean(LIP.Classification_correct_choice(:,:,1,:,loc),2),4),[prctile(mean(LIP.Classification_correct_choice(:,:,1,:,loc),2),95,4),prctile(mean(LIP.Classification_correct_choice(:,:,1,:,loc),2),5,4)]',{'Color',color_for_ROI(1,:),'LineWidth',2},2)
    
    plot_significance_level(window_start_list,p_choice_LIP(:,loc),1.1,color_for_ROI(1,:),[0.01,0.05/(length(window_start_list))/4, 0.01/(length(window_start_list))/4])
    
    yl=yline(0.5,'--',{'Chance level'},'Color',[0 0 0],'LineWidth',1);
    xline(-200,'-','Color',[0.5 0.5 0.5],'LineWidth',1)
    xl=xline(-300,'-',{'Targets on'},'Color',[0.5 0.5 0.5],'LineWidth',1);
    xl.LabelVerticalAlignment = 'middle';
    xl.LabelHorizontalAlignment = 'left';
   xl.FontSize=12;
    yl.LabelVerticalAlignment = 'middle';
    yl.LabelHorizontalAlignment = 'right';
    yl.FontSize=12;
    ylim([0.35 1.15])
%     text(-575,1,'Choice')
    box off
    xlabel('Time to response')
    ylabel('Choice decoding accuracy')
    xlim([-600 window_start_list(end)+150])
    xticks([-600:200:window_start_list(end)+150])
    
    
    
    subplot(3,4,loc+4)
    hold on
    plot(window_start_list,mean(mean(FEF.Classification_correct_choice(:,:,1,:,loc),2),4),'Color',color_for_ROI(2,:),'LineWidth',1)
    
    shadedErrorBar(window_start_list,mean(mean(FEF.Classification_correct_choice(:,:,1,:,loc),2),4),[prctile(mean(FEF.Classification_correct_choice(:,:,1,:,loc),2),95,4),prctile(mean(FEF.Classification_correct_choice(:,:,1,:,loc),2),5,4)]',{'Color',color_for_ROI(2,:),'LineWidth',2},2)
    
    plot_significance_level(window_start_list,p_choice_FEF(:,loc),1.1,color_for_ROI(2,:),[0.01,0.05/(length(window_start_list))/4, 0.01/(length(window_start_list))/4])
    
    yl=yline(0.5,'--',{'Chance level'},'Color',[0 0 0],'LineWidth',1);
    xline(-200,'-','Color',[0.5 0.5 0.5],'LineWidth',1)
    xl=xline(-300,'-',{'Targets on'},'Color',[0.5 0.5 0.5],'LineWidth',1);
    xl.LabelVerticalAlignment = 'middle';
    xl.LabelHorizontalAlignment = 'left';
    xl.FontSize=12;
    yl.LabelVerticalAlignment = 'middle';
    yl.LabelHorizontalAlignment = 'right';
    yl.FontSize=12;
    ylim([0.35 1.15])
%     text(-575,1,'Within - across')
%     text(-575,1.05,'Within')
%     text(-575,1.1,'Across')
    box off
    xlabel('Time to response')
    ylabel('Choice decoding accuracy')
xlim([-600 window_start_list(end)+150])
xticks([-600:200:window_start_list(end)+150])
    

subplot(3,4,loc+8)
    hold on
    plot(window_start_list,mean(mean(PFC.Classification_correct_choice(:,:,1,:,loc),2),4),'Color',color_for_ROI(3,:),'LineWidth',1)
    
    shadedErrorBar(window_start_list,mean(mean(PFC.Classification_correct_choice(:,:,1,:,loc),2),4),[prctile(mean(PFC.Classification_correct_choice(:,:,1,:,loc),2),95,4),prctile(mean(PFC.Classification_correct_choice(:,:,1,:,loc),2),5,4)]',{'Color',color_for_ROI(3,:),'LineWidth',2},2)
    
    plot_significance_level(window_start_list,p_choice_PFC(:,loc),1.1,color_for_ROI(3,:),[0.01,0.05/(length(window_start_list))/4, 0.01/(length(window_start_list))/4])
    yl=yline(0.5,'--',{'Chance level'},'Color',[0 0 0],'LineWidth',1);
    xline(-200,'-','Color',[0.5 0.5 0.5],'LineWidth',1)
    xl=xline(-300,'-',{'Targets on'},'Color',[0.5 0.5 0.5],'LineWidth',1);
   xl.LabelVerticalAlignment = 'middle';
    xl.LabelHorizontalAlignment = 'left';
    xl.FontSize=12;
    yl.LabelVerticalAlignment = 'middle';
    yl.LabelHorizontalAlignment = 'right';
    yl.FontSize=12;
    ylim([0.35 1.15])
%     text(-575,1,'Within - across')
%     text(-575,1.05,'Within')
%     text(-575,1.1,'Across')
    box off
    xlabel('Time to response')
    ylabel('Choice decoding accuracy')
xlim([-600 window_start_list(end)+150])
xticks([-600:200:window_start_list(end)+150])
    
end



%%




figure

for loc=1:N_loc
    subplot(3,4,loc)
    hold on
    
    shadedErrorBar(window_start_list,mean(mean(LIP.Same_belief_choice(:,:,:,loc),2),3),[prctile(mean(LIP.Same_belief_choice(:,:,:,loc),2),95,3),prctile(mean(LIP.Same_belief_choice(:,:,:,loc),2),5,3)]',{'Color',color_for_ROI(1,:),'LineWidth',2},2)
    shadedErrorBar(window_start_list,mean(mean(LIP.Other_belief_choice(:,:,:,loc),2),3),[prctile(mean(LIP.Other_belief_choice(:,:,:,loc),2),95,3),prctile(mean(LIP.Other_belief_choice(:,:,:,loc),2),5,3)]',{'--','Color',color_for_ROI(1,:),'LineWidth',2},2)
    
    plot_significance_level(window_start_list,p_LIP_same_other(:,loc),1,color_for_ROI(1,:),[0.01,0.05/(length(window_start_list))/4, 0.01/(length(window_start_list))/4])
    plot_significance_level(window_start_list,p_LIP_other(:,loc),1.05,color_for_ROI(1,:),[0.01,0.05/(length(window_start_list))/4, 0.01/(length(window_start_list))/4])
    plot_significance_level(window_start_list,p_LIP_same(:,loc),1.1,color_for_ROI(1,:),[0.01,0.05/(length(window_start_list))/4, 0.01/(length(window_start_list))/4])
    
    yl=yline(0.5,'--',{'Chance level'},'Color',[0 0 0],'LineWidth',1);
    xline(-200,'-','Color',[0.5 0.5 0.5],'LineWidth',1)
    xl=xline(-300,'-',{'Targets on'},'Color',[0.5 0.5 0.5],'LineWidth',1);
    xl.LabelVerticalAlignment = 'middle';
    xl.LabelHorizontalAlignment = 'left';
   xl.FontSize=12;
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
xlim([-600 window_start_list(end)+150])
xticks([-600:200:window_start_list(end)+150])
    subplot(3,4,loc+4)
    hold on
    
    shadedErrorBar(window_start_list,mean(mean(FEF.Same_belief_choice(:,:,:,loc),2),3),[prctile(mean(FEF.Same_belief_choice(:,:,:,loc),2),95,3),prctile(mean(FEF.Same_belief_choice(:,:,:,loc),2),5,3)]',{'Color',color_for_ROI(2,:),'LineWidth',2},2)
    shadedErrorBar(window_start_list,mean(mean(FEF.Other_belief_choice(:,:,:,loc),2),3),[prctile(mean(FEF.Other_belief_choice(:,:,:,loc),2),95,3),prctile(mean(FEF.Other_belief_choice(:,:,:,loc),2),5,3)]',{'--','Color',color_for_ROI(2,:),'LineWidth',2},2)
    
    plot_significance_level(window_start_list,p_FEF_same_other(:,loc),1,color_for_ROI(2,:),[0.01,0.05/(length(window_start_list))/4, 0.01/(length(window_start_list))/4])
    plot_significance_level(window_start_list,p_FEF_other(:,loc),1.05,color_for_ROI(2,:),[0.01,0.05/(length(window_start_list))/4, 0.01/(length(window_start_list))/4])
    plot_significance_level(window_start_list,p_FEF_same(:,loc),1.1,color_for_ROI(2,:),[0.01,0.05/(length(window_start_list))/4, 0.01/(length(window_start_list))/4])
    
    yl=yline(0.5,'--',{'Chance level'},'Color',[0 0 0],'LineWidth',1);
    xline(-200,'-','Color',[0.5 0.5 0.5],'LineWidth',1)
    xl=xline(-300,'-',{'Targets on'},'Color',[0.5 0.5 0.5],'LineWidth',1);
    xl.LabelVerticalAlignment = 'middle';
    xl.LabelHorizontalAlignment = 'left';
    xl.FontSize=12;
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
xlim([-600 window_start_list(end)+150])
xticks([-600:200:window_start_list(end)+150])
    subplot(3,4,loc+8)
    hold on
    
    shadedErrorBar(window_start_list,mean(mean(PFC.Same_belief_choice(:,:,:,loc),2),3),[prctile(mean(PFC.Same_belief_choice(:,:,:,loc),2),95,3),prctile(mean(PFC.Same_belief_choice(:,:,:,loc),2),5,3)]',{'Color',color_for_ROI(3,:),'LineWidth',2},2)
    shadedErrorBar(window_start_list,mean(mean(PFC.Other_belief_choice(:,:,:,loc),2),3),[prctile(mean(PFC.Other_belief_choice(:,:,:,loc),2),95,3),prctile(mean(PFC.Other_belief_choice(:,:,:,loc),2),5,3)]',{'--','Color',color_for_ROI(3,:),'LineWidth',2},2)
    
    plot_significance_level(window_start_list,p_PFC_same_other(:,loc),1,color_for_ROI(3,:),[0.01,0.05/(length(window_start_list))/4, 0.01/(length(window_start_list))/4])
    plot_significance_level(window_start_list,p_PFC_other(:,loc),1.05,color_for_ROI(3,:),[0.01,0.05/(length(window_start_list))/4, 0.01/(length(window_start_list))/4])
    plot_significance_level(window_start_list,p_PFC_same(:,loc),1.1,color_for_ROI(3,:),[0.01,0.05/(length(window_start_list))/4, 0.01/(length(window_start_list))/4])
    
    yl=yline(0.5,'--',{'Chance level'},'Color',[0 0 0],'LineWidth',1);
    xline(-200,'-','Color',[0.5 0.5 0.5],'LineWidth',1)
    xl=xline(-300,'-',{'Targets on'},'Color',[0.5 0.5 0.5],'LineWidth',1);
   xl.LabelVerticalAlignment = 'middle';
    xl.LabelHorizontalAlignment = 'left';
    xl.FontSize=12;
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
xlim([-600 window_start_list(end)+150])
xticks([-600:200:window_start_list(end)+150])
    
end




%%

for loc=1:N_loc
    
    figure
    for m=1:3
        subplot(3,3,m)
        for n=1:3
            
            plot(window_start_list,mean(mean(LIP.Classification_correct_choice_each_belief(:,:,m,n,1,:,loc),2),6),'Color',color_for_ROI_early(1,:),'LineWidth',1)
            % hold on
            % plot(window_start_list,mean(mean(LIP.Classification_correct_choice(:,:,3,:,loc),2),4),'Color',color_for_ROI_late(1,:),'LineWidth',2)
            hold on
            shadedErrorBar(window_start_list,mean(mean(LIP.Classification_correct_choice_each_belief(:,:,m,n,1,:,loc),2),6),std(mean(LIP.Classification_correct_choice_each_belief(:,:,m,n,1,:,loc),2),0,6),{'Color',color_for_ROI_early(1,:),'LineWidth',(m==n)*2+1},2)
            % hold on
            % shadedErrorBar(window_start_list,mean(mean(LIP.Classification_correct_choice(:,:,3,:,loc),2),4),std(mean(LIP.Classification_correct_choice(:,:,3,:,loc),2),0,4),{'Color',color_for_ROI_late(1,:),'LineWidth',2},2)
            hold on
            %             for i=1:length(window_start_list)
            %                 if p_choice_each_belief_LIP(i,m,n,1,loc)<0.05%%/length(window_start_list) %Bonferoni
            %                     plot(window_start_list(i),1,'*','Color',color_for_ROI_early(1,:),'LineWidth',2)
            %                 end
            %                 %     if p_choice_LIP(i,3,loc)<0.05%/length(window_start_list) %Bonferoni
            %                 %         plot(window_start_list(i),0.85,'*','Color',color_for_ROI_late(1,:),'LineWidth',2)
            %                 %     end
            %             end
        end
        yl=yline(0.5,'--',{'Chance level'},'Color',[0 0 0],'LineWidth',1);
        
        yl.LabelVerticalAlignment = 'middle';
        yl.LabelHorizontalAlignment = 'right';
        yl.FontSize=16;
        
        
        switch event
            case 'target'
                hold on
                xline(-300,'-','Color',[0.5 0.5 0.5],'LineWidth',1)
                hold on
                xl=xline(-400,'-',{'Fixation'},'Color',[0.5 0.5 0.5],'LineWidth',1);
                xl.LabelVerticalAlignment = 'bottom';
                xl.LabelHorizontalAlignment = 'left';
                xl.FontSize=16;
                
                hold on
                xline(300,'-','Color',[0.25 0.75 0.25],'LineWidth',1)
                hold on
                xl=xline(200,'-',{'Response'},'Color',[0.25 0.75 0.25],'LineWidth',1);
                xl.LabelVerticalAlignment = 'bottom';
                xl.LabelHorizontalAlignment = 'left';
                xl.FontSize=16;
        end
        ylabel('Decoder accuracy','FontSize',14)
        xlabel(sprintf('Time to %s',event'))
        ax=gca;
        ax.XAxis.FontSize=16;
        ax.YAxis.FontSize=16;
        box off
        ylim([0.18 1])
    end
    
    for m=1:3
        subplot(3,3,m+3)
        for n=1:3
            plot(window_start_list,mean(mean(FEF.Classification_correct_choice_each_belief(:,:,m,n,1,:,loc),2),6),'Color',color_for_ROI_early(2,:),'LineWidth',1)
            % hold on
            % plot(window_start_list,mean(mean(FEF.Classification_correct_choice(:,:,3,:,loc),2),4),'Color',color_for_ROI_late(1,:),'LineWidth',2)
            hold on
            shadedErrorBar(window_start_list,mean(mean(FEF.Classification_correct_choice_each_belief(:,:,m,n,1,:,loc),2),6),std(mean(FEF.Classification_correct_choice_each_belief(:,:,m,n,1,:,loc),2),0,6),{'Color',color_for_ROI_early(2,:),'LineWidth',(m==n)*2+1},2)
            % hold on
            % shadedErrorBar(window_start_list,mean(mean(FEF.Classification_correct_choice(:,:,3,:,loc),2),4),std(mean(FEF.Classification_correct_choice(:,:,3,:,loc),2),0,4),{'Color',color_for_ROI_late(1,:),'LineWidth',2},2)
            hold on
            hold on
        end
        %             for i=1:length(window_start_list)
        %                 if p_choice_each_belief_FEF(i,m,n,1,loc)<0.05%%/length(window_start_list) %Bonferoni
        %                     plot(window_start_list(i),1,'*','Color',color_for_ROI_early(2,:),'LineWidth',2)
        %                 end
        %                 %     if p_choice_FEF(i,3,loc)<0.05%/length(window_start_list) %Bonferoni
        %                 %         plot(window_start_list(i),0.85,'*','Color',color_for_ROI_late(1,:),'LineWidth',2)
        %                 %     end
        %             end
        yl=yline(0.5,'--',{'Chance level'},'Color',[0 0 0],'LineWidth',1);
        
        yl.LabelVerticalAlignment = 'middle';
        yl.LabelHorizontalAlignment = 'right';
        yl.FontSize=16;
        
        
        switch event
            case 'target'
                hold on
                xline(-300,'-','Color',[0.5 0.5 0.5],'LineWidth',1)
                hold on
                xl=xline(-400,'-',{'Fixation'},'Color',[0.5 0.5 0.5],'LineWidth',1);
                xl.LabelVerticalAlignment = 'bottom';
                xl.LabelHorizontalAlignment = 'left';
                xl.FontSize=16;
                
                hold on
                xline(300,'-','Color',[0.25 0.75 0.25],'LineWidth',1)
                hold on
                xl=xline(200,'-',{'Response'},'Color',[0.25 0.75 0.25],'LineWidth',1);
                xl.LabelVerticalAlignment = 'bottom';
                xl.LabelHorizontalAlignment = 'left';
                xl.FontSize=16;
        end
        ylabel('Decoder accuracy','FontSize',14)
        xlabel(sprintf('Time to %s',event'))
        ax=gca;
        ax.XAxis.FontSize=16;
        ax.YAxis.FontSize=16;
        box off
        ylim([0.18 1])
    end
    
    
    for m=1:3
        subplot(3,3,m+6)
        for n=1:3
            plot(window_start_list,mean(mean(PFC.Classification_correct_choice_each_belief(:,:,m,n,1,:,loc),2),6),'Color',color_for_ROI_early(3,:),'LineWidth',1)
            % hold on
            % plot(window_start_list,mean(mean(PFC.Classification_correct_choice(:,:,3,:,loc),2),4),'Color',color_for_ROI_late(1,:),'LineWidth',2)
            hold on
            shadedErrorBar(window_start_list,mean(mean(PFC.Classification_correct_choice_each_belief(:,:,m,n,1,:,loc),2),6),std(mean(PFC.Classification_correct_choice_each_belief(:,:,m,n,1,:,loc),2),0,6),{'Color',color_for_ROI_early(3,:),'LineWidth',(m==n)*2+1},2)
            % hold on
            % shadedErrorBar(window_start_list,mean(mean(PFC.Classification_correct_choice(:,:,3,:,loc),2),4),std(mean(PFC.Classification_correct_choice(:,:,3,:,loc),2),0,4),{'Color',color_for_ROI_late(1,:),'LineWidth',2},2)
            hold on
            hold on
            %             for i=1:length(window_start_list)
            %                 if p_choice_each_belief_PFC(i,m,n,1,loc)<0.05%%/length(window_start_list) %Bonferoni
            %                     plot(window_start_list(i),1,'*','Color',color_for_ROI_early(3,:),'LineWidth',2)
            %                 end
            %                 %     if p_choice_PFC(i,3,loc)<0.05%/length(window_start_list) %Bonferoni
            %                 %         plot(window_start_list(i),0.85,'*','Color',color_for_ROI_late(1,:),'LineWidth',2)
            %                 %     end
            %             end
        end
        yl=yline(0.5,'--',{'Chance level'},'Color',[0 0 0],'LineWidth',1);
        
        yl.LabelVerticalAlignment = 'middle';
        yl.LabelHorizontalAlignment = 'right';
        yl.FontSize=16;
        
        switch event
            case 'target'
                hold on
                xline(-300,'-','Color',[0.5 0.5 0.5],'LineWidth',1)
                hold on
                xl=xline(-400,'-',{'Fixation'},'Color',[0.5 0.5 0.5],'LineWidth',1);
                xl.LabelVerticalAlignment = 'bottom';
                xl.LabelHorizontalAlignment = 'left';
                xl.FontSize=16;
                
                hold on
                xline(300,'-','Color',[0.25 0.75 0.25],'LineWidth',1)
                hold on
                xl=xline(200,'-',{'Response'},'Color',[0.25 0.75 0.25],'LineWidth',1);
                xl.LabelVerticalAlignment = 'bottom';
                xl.LabelHorizontalAlignment = 'left';
                xl.FontSize=16;
        end
        ylabel('Decoder accuracy','FontSize',14)
        xlabel(sprintf('Time to %s',event'))
        ax=gca;
        ax.XAxis.FontSize=16;
        ax.YAxis.FontSize=16;
        box off
        ylim([0.18 1])
    end
    
    
    
end



%%
figure
for loc=1:N_loc
    subplot(3,4,loc)
    buffer(:,:)=mean(mean(LIP.Classification_correct_choice_each_belief(13,:,:,:,1,:,loc),2),6);
    imagesc(buffer,[0.5 0.9])
    title(sprintf('%s location %d','LIP',loc))
    colorbar
    xlabel('Template color train')
    ylabel('Template color test')
    
    subplot(3,4,loc+4)
    buffer(:,:)=mean(mean(FEF.Classification_correct_choice_each_belief(13,:,:,:,1,:,loc),2),6);
    imagesc(buffer,[0.5 0.8])
    title(sprintf('%s location %d','FEF',loc))
    colorbar
    xlabel('Template color train')
    ylabel('Template color test')
    
    subplot(3,4,loc+8)
    buffer(:,:)=mean(mean(PFC.Classification_correct_choice_each_belief(13,:,:,:,1,:,loc),2),6);
    imagesc(buffer,[0.5 0.8])
    title(sprintf('%s location %d','PFC',loc))
    colorbar
    xlabel('Template color train')
    ylabel('Template color test')
    
end

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
