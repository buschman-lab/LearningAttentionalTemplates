clear all

ROI='LIP';

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
save_path=fullfile(fsroot,dirstem2);

save_name=sprintf('glm_explained_variance_split_%s',ROI);

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

switch ROI
    case 'LIP'
        color_for_plot=color_for_ROI(1,:);
        n_color=1;
    case 'FEF'
        color_for_plot=color_for_ROI(2,:);
        n_color=2;
    case 'PFC'
        color_for_plot=color_for_ROI(3,:);
        n_color=3;
end



%%

for i=1:length(arrayID_list)
    arrayID=arrayID_list(i);
    
    if (arrayID==34 && strcmp(ROI,'LIP')) || (arrayID==30 && strcmp(ROI,'LIP')) || (arrayID==14 && strcmp(ROI,'FEF'))
        
    else
        for this_time=1:length(window_start_list)
            load_name=sprintf('glm_explained_variance_split_%s_%d_%d',ROI,arrayID,this_time);
            load(fullfile(save_path,load_name),'Mean_z_R_square_within_loc_across_split','Mean_z_R_square_explainable_across_loc','Mean_z_R_square_across_loc','Mean_z_R_square_within_loc')
            
            for n=1:size(Mean_z_R_square_within_loc_across_split,1)
                
                All_Mean_z_R_square_within_loc_across_split(count_ROI(this_time),this_time,:)=Mean_z_R_square_within_loc_across_split(n,:);
                All_Mean_z_R_square_explainable_across_loc(count_ROI(this_time),this_time,:,:) = Mean_z_R_square_explainable_across_loc(n,:,:);
                All_Mean_z_R_square_across_loc(count_ROI(this_time),this_time,:,:)=(Mean_z_R_square_across_loc(n,:,:,1,2)+Mean_z_R_square_across_loc(n,:,:,2,1))./2;
                All_Mean_z_R_square_within_loc_within_split(count_ROI(this_time),this_time,:) = (Mean_z_R_square_within_loc(n,:,1,1)+Mean_z_R_square_within_loc(n,:,2,2))/2;
                
                count_ROI(this_time)=count_ROI(this_time)+1;
            end
            
            clear Mean_z_R_square_within_loc_across_split Mean_z_R_square_explainable_across_loc Mean_z_R_square_across_loc Mean_z_R_square_within_loc
        end
    end
end


%%

contra_z_R_square_across_loc=vertcat(All_Mean_z_R_square_across_loc(:,:,1,1),All_Mean_z_R_square_across_loc(:,:,2,2));
contra_contra_z_R_square_across_loc=vertcat(All_Mean_z_R_square_across_loc(:,:,1,2),All_Mean_z_R_square_across_loc(:,:,2,1));
contra_ipsi_z_R_square_across_loc=vertcat(All_Mean_z_R_square_across_loc(:,:,1,4),All_Mean_z_R_square_across_loc(:,:,2,3));
contra_ipsi_diag_z_R_square_across_loc=vertcat(All_Mean_z_R_square_across_loc(:,:,1,3),All_Mean_z_R_square_across_loc(:,:,2,4));

%%

for t=1:length(window_start_list)
    [~, p_contra_z_R_square_across_loc(t)]=ttest(contra_z_R_square_across_loc(:,t),0,'Tail','right');
    [~, p_contra_contra_z_R_square_across_loc(t)]=ttest(contra_contra_z_R_square_across_loc(:,t),0,'Tail','right');
    [~, p_contra_ipsi_z_R_square_across_loc(t)]=ttest(contra_ipsi_z_R_square_across_loc(:,t),0,'Tail','right');
    [~, p_contra_ipsi_diag_z_R_square_across_loc(t)]=ttest(contra_ipsi_diag_z_R_square_across_loc(:,t),0,'Tail','right');
    
    [~, p_contra_to_contra_contra(t)]=ttest(contra_z_R_square_across_loc(:,t)-contra_contra_z_R_square_across_loc(:,t),0,'Tail','right');
    [~, p_contra_to_contra_ipsi(t)]=ttest(contra_z_R_square_across_loc(:,t)-contra_ipsi_z_R_square_across_loc(:,t),0,'Tail','right');
    [~, p_contra_to_contra_ipsi_diag(t)]=ttest(contra_z_R_square_across_loc(:,t)-contra_ipsi_diag_z_R_square_across_loc(:,t),0,'Tail','right');
end

%%

save(fullfile(save_path,save_name))
%%
figure
subplot(1,3,1)
hold on

shadedErrorBar(window_start_list(4:end),nanmean(contra_z_R_square_across_loc(:,4:end),1),nanstd(contra_z_R_square_across_loc(:,4:end),0,1)./sqrt(sum(~isnan(contra_z_R_square_across_loc(:,4:end)),1)-1),{'-','color',color_each_loc(n_color,:,4),'LineWidth',2},2)
shadedErrorBar(window_start_list(4:end),nanmean(contra_contra_z_R_square_across_loc(:,4:end),1),nanstd(contra_contra_z_R_square_across_loc(:,4:end),0,1)./sqrt(sum(~isnan(contra_contra_z_R_square_across_loc(:,4:end)),1)-1),{'--','color',color_each_loc(n_color,:,3),'LineWidth',2},2)

plot_significance_level(window_start_list,p_contra_z_R_square_across_loc,1.05,color_each_loc(n_color,:,4),[0.01,0.05/(length(window_start_list)-4), 0.01/(length(window_start_list)-4)])
plot_significance_level(window_start_list,p_contra_contra_z_R_square_across_loc,1,color_each_loc(n_color,:,3),[0.01,0.05/(length(window_start_list)-4), 0.01/(length(window_start_list)-4)])
plot_significance_level(window_start_list,p_contra_to_contra_contra,0.95,'k',[0.01,0.05/(length(window_start_list)-4), 0.01/(length(window_start_list)-4)])

subplot(1,3,2)
hold on
shadedErrorBar(window_start_list(4:end),nanmean(contra_z_R_square_across_loc(:,4:end),1),nanstd(contra_z_R_square_across_loc(:,4:end),0,1)./sqrt(sum(~isnan(contra_z_R_square_across_loc(:,4:end)),1)-1),{'-','color',color_each_loc(n_color,:,4),'LineWidth',2},2)
shadedErrorBar(window_start_list(4:end),nanmean(contra_ipsi_z_R_square_across_loc(:,4:end),1),nanstd(contra_ipsi_z_R_square_across_loc(:,4:end),0,1)./sqrt(sum(~isnan(contra_ipsi_z_R_square_across_loc(:,4:end)),1)-1),{'--','color',color_each_loc(n_color,:,2),'LineWidth',2},2)

plot_significance_level(window_start_list,p_contra_z_R_square_across_loc,1.05,color_each_loc(n_color,:,4),[0.01,0.05/(length(window_start_list)-4), 0.01/(length(window_start_list)-4)])
plot_significance_level(window_start_list,p_contra_ipsi_z_R_square_across_loc,1,color_each_loc(n_color,:,2),[0.01,0.05/(length(window_start_list)-4), 0.01/(length(window_start_list)-4)])
plot_significance_level(window_start_list,p_contra_to_contra_ipsi,0.95,'k',[0.01,0.05/(length(window_start_list)-4), 0.01/(length(window_start_list)-4)])

subplot(1,3,3)
hold on
shadedErrorBar(window_start_list(4:end),nanmean(contra_z_R_square_across_loc(:,4:end),1),nanstd(contra_z_R_square_across_loc(:,4:end),0,1)./sqrt(sum(~isnan(contra_z_R_square_across_loc(:,4:end)),1)-1),{'-','color',color_each_loc(n_color,:,4),'LineWidth',2},2)
shadedErrorBar(window_start_list(4:end),nanmean(contra_ipsi_diag_z_R_square_across_loc(:,4:end),1),nanstd(contra_ipsi_diag_z_R_square_across_loc(:,4:end),0,1)./sqrt(sum(~isnan(contra_ipsi_diag_z_R_square_across_loc(:,4:end)),1)-1),{'--','color',color_each_loc(n_color,:,1),'LineWidth',2},2)

plot_significance_level(window_start_list,p_contra_z_R_square_across_loc,1.05,color_each_loc(n_color,:,4),[0.01,0.05/(length(window_start_list)-4), 0.01/(length(window_start_list)-4)])
plot_significance_level(window_start_list,p_contra_ipsi_diag_z_R_square_across_loc,1,color_each_loc(n_color,:,1),[0.01,0.05/(length(window_start_list)-4), 0.01/(length(window_start_list)-4)])
plot_significance_level(window_start_list,p_contra_to_contra_ipsi_diag,0.95,'k',[0.01,0.05/(length(window_start_list)-4), 0.01/(length(window_start_list)-4)])







%% Functions

function plot_significance_level(x,p,a,c,thr)
hold on
for b=1:length(thr)
    this_thr=thr(b);
    for i=4:length(p)-1
        if p(i)<=this_thr && p(i+1)<=this_thr
            plot(x(i):x(i+1),a*ones(1,x(i+1)-x(i)+1),'-','Color',c,'LineWidth',b)
        end
        if i>1
            if p(i-1)>=this_thr && p(i)<=this_thr && p(i+1)>=this_thr
                plot(x(i),a,'.','Color',c,'LineWidth',b)
            end
        end
    end
    if p(end-1)>=this_thr && p(end)<=this_thr
        plot(x(end),a,'.','Color',c,'LineWidth',b)
    end
end

end


