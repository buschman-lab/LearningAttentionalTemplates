%load and plot model stim belief cos

fsroot='/Volumes/buschman';

arrayID_list=[3 8 12 14 18 20 22 26 28 30 32 34 36 38 40 42 44]; % 42

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
% subtask='exploreexploit/Reset_RW_model';
subtask='exploreexploit/Restricted_FEF/Reset_RW_model';

window_size=300;


color_for_ROI=[0 0.447 0.741; 0.85 0.325 0.098; 0.929 0.694 0.125]; %LIP FEF PFC

count_LIP=0;
count_FEF=0;
count_PFC=0;

Wv_thr=12;


%%

for arrayID_ind=1:length(arrayID_list)
    
    arrayID=arrayID_list(arrayID_ind);
    
    if arrayID~=34
        
        %LIP
        data_name=sprintf('ID_%d.mat',arrayID);
        subsubtask=sprintf('FR_%s_%d_window_from_%d_to_%d','target',window_size,-300,-300+window_size);
        dirstem = fullfile('Projects',task,'Analysis','Electrophy_analysis',subtask,subsubtask);
        data_path = fullfile(fsroot,dirstem);
        
        load(fullfile(data_path,data_name),'LIP_wv')
        
        for n=1:size(LIP_wv,1)
            Wv_LIP(count_LIP+n,:)=LIP_wv(n,:);
            Min_to_Max_LIP(count_LIP+n)=find(Wv_LIP(count_LIP+n,:)==max(Wv_LIP(count_LIP+n,:)))-find(Wv_LIP(count_LIP+n,:)==min(Wv_LIP(count_LIP+n,:)));
        end
        clear LIP_wv
        
        for this_time=1:length(window_start_list)
            
            event=event_list{this_time};
            
            subsubtask=sprintf('FR_%s_%d_window_from_%d_to_%d',event,window_size,window_start_list(this_time),window_start_list(this_time)+window_size);
            
            dirstem = fullfile('Projects',task,'Analysis','Electrophy_analysis',subtask,subsubtask);
            data_path = fullfile(fsroot,dirstem);
            
            %             data_name=sprintf('models_all_types_belief_vm_baseline_%s_%d','LIP',arrayID);
            data_name=sprintf('models_FR_figure2_value_function_%s_%d','LIP',arrayID);
            
            str=[fullfile(data_path,data_name) '.mat'];
            if exist(str)>0
                
                load(fullfile(data_path,data_name),'R2','param','Neuron_id')
                
                if exist('R2')
                    
                    for n=1:size(R2,1)
                        
                        R2_LIP(count_LIP+Neuron_id(n),:,:,this_time)=R2(n,:,:);
                        param_LIP(count_LIP+Neuron_id(n),:,:,:,:,this_time)=param(n,:,:,:,:);
                        LIP_sess_all(count_LIP+Neuron_id(n),this_time)=arrayID_ind;
                        LIP_neur_all(count_LIP+Neuron_id(n),this_time)=Neuron_id(n);
                    end
                    
                    clear R2 param Neuron_id
                    
                end
                
            end
            
        end
        
        count_LIP=size(Wv_LIP,1);
        
    end
    
    %FEF
    data_name=sprintf('ID_%d.mat',arrayID);
    subsubtask=sprintf('FR_%s_%d_window_from_%d_to_%d','target',window_size,-300,-300+window_size);
    dirstem = fullfile('Projects',task,'Analysis','Electrophy_analysis',subtask,subsubtask);
    data_path = fullfile(fsroot,dirstem);
    
    load(fullfile(data_path,data_name),'FEF_wv')
    
    for n=1:size(FEF_wv,1)
        Wv_FEF(count_FEF+n,:)=FEF_wv(n,:);
        Min_to_Max_FEF(count_FEF+n)=find(Wv_FEF(count_FEF+n,:)==max(Wv_FEF(count_FEF+n,:)))-find(Wv_FEF(count_FEF+n,:)==min(Wv_FEF(count_FEF+n,:)));
    end
    clear FEF_wv
    
    for this_time=1:length(window_start_list)
        
        event=event_list{this_time};
        
        
        subsubtask=sprintf('FR_%s_%d_window_from_%d_to_%d',event,window_size,window_start_list(this_time),window_start_list(this_time)+window_size);
        
        dirstem = fullfile('Projects',task,'Analysis','Electrophy_analysis',subtask,subsubtask);
        data_path = fullfile(fsroot,dirstem);
        
        %         data_name=sprintf('models_all_types_belief_vm_baseline_%s_%d','FEF',arrayID);
        data_name=sprintf('models_FR_figure2_value_function_%s_%d','FEF',arrayID);
        
        str=[fullfile(data_path,data_name) '.mat'];
        if exist(str)>0
            
            load(fullfile(data_path,data_name),'R2','param','Neuron_id')
            
            
            if exist('R2')
                
                for n=1:size(R2,1)
                    
                    R2_FEF(count_FEF+Neuron_id(n),:,:,this_time)=R2(n,:,:);
                    param_FEF(count_FEF+Neuron_id(n),:,:,:,:,this_time)=param(n,:,:,:,:);
                    FEF_sess_all(count_FEF+Neuron_id(n),this_time)=arrayID_ind;
                    FEF_neur_all(count_FEF+Neuron_id(n),this_time)=Neuron_id(n);
                    
                end
                
                clear R2 param Neuron_id
                
            end
            
        end
        
    end
    
    count_FEF=size(Wv_FEF,1);
    
    %PFC
    data_name=sprintf('ID_%d.mat',arrayID);
    subsubtask=sprintf('FR_%s_%d_window_from_%d_to_%d','target',window_size,-300,-300+window_size);
    dirstem = fullfile('Projects',task,'Analysis','Electrophy_analysis',subtask,subsubtask);
    data_path = fullfile(fsroot,dirstem);
    
    load(fullfile(data_path,data_name),'PFC_wv')
    
    for n=1:size(PFC_wv,1)
        Wv_PFC(count_PFC+n,:)=PFC_wv(n,:);
        Min_to_Max_PFC(count_PFC+n)=find(Wv_PFC(count_PFC+n,:)==max(Wv_PFC(count_PFC+n,:)))-find(Wv_PFC(count_PFC+n,:)==min(Wv_PFC(count_PFC+n,:)));
    end
    clear PFC_wv
    for this_time=1:length(window_start_list)
        
        event=event_list{this_time};
        
        subsubtask=sprintf('FR_%s_%d_window_from_%d_to_%d',event,window_size,window_start_list(this_time),window_start_list(this_time)+window_size);
        
        dirstem = fullfile('Projects',task,'Analysis','Electrophy_analysis',subtask,subsubtask);
        data_path = fullfile(fsroot,dirstem);
        
        %         data_name=sprintf('models_all_types_belief_vm_baseline_%s_%d','PFC',arrayID);
        data_name=sprintf('models_FR_figure2_value_function_%s_%d','PFC',arrayID);
        
        
        str=[fullfile(data_path,data_name) '.mat'];
        if exist(str)>0
            
            
            load(fullfile(data_path,data_name),'R2','param','Neuron_id')
            
            
            if exist('R2')
                
                for n=1:size(R2,1)
                    
                    R2_PFC(count_PFC+Neuron_id(n),:,:,this_time)=R2(n,:,:);
                    param_PFC(count_PFC+Neuron_id(n),:,:,:,:,this_time)=param(n,:,:,:,:);
                    PFC_sess_all(count_PFC+Neuron_id(n),this_time)=arrayID_ind;
                    PFC_neur_all(count_PFC+Neuron_id(n),this_time)=Neuron_id(n);
                end
                
                clear R2 param Neuron_id
                
            end
        end
    end
    
    count_PFC=size(Wv_PFC,1);
    
end

%% remove wrong Wv
Min_to_Max_LIP(Min_to_Max_LIP>=25 | Min_to_Max_LIP<=0)=NaN;
Min_to_Max_FEF(Min_to_Max_FEF>=25 | Min_to_Max_FEF<=0)=NaN;
Min_to_Max_PFC(Min_to_Max_PFC>=25 | Min_to_Max_PFC<=0)=NaN;

%% padding
R2_LIP(size(R2_LIP,1)+1:size(Min_to_Max_LIP,2),:,:,:)=NaN;
R2_FEF(size(R2_FEF,1)+1:size(Min_to_Max_FEF,2),:,:,:)=NaN;
R2_PFC(size(R2_PFC,1)+1:size(Min_to_Max_PFC,2),:,:,:)=NaN;

param_LIP(size(param_LIP,1)+1:size(Min_to_Max_LIP,2),:,:,:,:,:)=NaN;
param_FEF(size(param_FEF,1)+1:size(Min_to_Max_FEF,2),:,:,:,:,:)=NaN;
param_PFC(size(param_PFC,1)+1:size(Min_to_Max_PFC,2),:,:,:,:,:)=NaN;

%%
clear h* p_*

h_LIP=NaN(size(R2_LIP,2),length(window_start_list),size(Min_to_Max_LIP,2));
h_FEF=NaN(size(R2_FEF,2),length(window_start_list),size(Min_to_Max_FEF,2));
h_PFC=NaN(size(R2_PFC,2),length(window_start_list),size(Min_to_Max_PFC,2));

p_LIP=NaN(size(R2_LIP,2),length(window_start_list),size(Min_to_Max_LIP,2));
p_FEF=NaN(size(R2_FEF,2),length(window_start_list),size(Min_to_Max_FEF,2));
p_PFC=NaN(size(R2_PFC,2),length(window_start_list),size(Min_to_Max_PFC,2));


for this_time=1:length(window_start_list)
    for j=1:size(R2_LIP,2)
        for i=1:size(Min_to_Max_LIP,2)
            if ~isnan(mean(R2_LIP(i,j,:,this_time)))
                if sum(R2_LIP(i,j,:,this_time)==0,3)==size(R2_LIP,3)
                    h_LIP(j,this_time,i)=NaN;
                    p_LIP(j,this_time,i)=NaN;
                else
                    h_LIP(j,this_time,i)=(mean(R2_LIP(i,j,:,this_time),3)>0);
                    p_LIP(j,this_time,i)=(mean(R2_LIP(i,j,:,this_time),3)<=0);
                end
            end
        end
        for i=1:size(Min_to_Max_FEF,2)
            if ~isnan(mean(R2_FEF(i,j,:,this_time)))
                if sum(R2_FEF(i,j,:,this_time)==0,3)==size(R2_FEF,3)
                    h_FEF(j,this_time,i)=NaN;
                    p_FEF(j,this_time,i)=NaN;
                else
                    h_FEF(j,this_time,i)=(mean(R2_FEF(i,j,:,this_time),3)>0);
                    p_FEF(j,this_time,i)=(mean(R2_FEF(i,j,:,this_time),3)<=0);
                end
            end
        end
        for i=1:size(Min_to_Max_PFC,2)
            if ~isnan(mean(R2_PFC(i,j,:,this_time)))
                if sum(R2_PFC(i,j,:,this_time)==0,3)==size(R2_FEF,3)
                    h_PFC(j,this_time,i)=NaN;
                    p_PFC(j,this_time,i)=NaN;
                else
                    h_PFC(j,this_time,i)=(mean(R2_PFC(i,j,:,this_time),3)>0);
                    p_PFC(j,this_time,i)=(mean(R2_PFC(i,j,:,this_time),3)<=0);
                end
            end
        end
    end
end

%%
for this_time=1:length(window_start_list)
    total_cells_LIP(this_time)=sum(~isnan(h_LIP(1,this_time,:)),3);
    total_cells_FEF(this_time)=sum(~isnan(h_FEF(1,this_time,:)),3);
    total_cells_PFC(this_time)=sum(~isnan(h_PFC(1,this_time,:)),3);
end

%%
h_LIP_best_model=zeros(4,length(window_start_list),size(h_LIP,3));
h_FEF_best_model=zeros(4,length(window_start_list),size(h_FEF,3));
h_PFC_best_model=zeros(4,length(window_start_list),size(h_PFC,3));

for this_time=1:length(window_start_list)
    for i=1:size(Min_to_Max_LIP,2)
        if sum(h_LIP(:,this_time,i))>=1
            [~, ind]=max(mean(R2_LIP(i,:,:,this_time),3),[],2);
            h_LIP_best_model(ind,this_time,i)=1;
        end
    end
    for i=1:size(Min_to_Max_FEF,2)
        if sum(h_FEF(:,this_time,i))>=1
            [~, ind]=max(mean(R2_FEF(i,:,:,this_time),3),[],2);
            h_FEF_best_model(ind,this_time,i)=1;
        end
    end
    for i=1:size(Min_to_Max_PFC,2)
        if sum(h_PFC(:,this_time,i))>=1
            [~, ind]=max(mean(R2_PFC(i,:,:,this_time),3),[],2);
            h_PFC_best_model(ind,this_time,i)=1;
        end
    end
end


%% Plot

figure
for i=1:4
    subplot(4,2,1+(i-1)*2)
    plot(window_start_list(1:7),nansum(h_LIP_best_model(i,1:7,:),3)./total_cells_LIP(1:7),'Color',color_for_ROI(1,:),'LineWidth',2)
    hold on
    plot(window_start_list(1:7),nansum(h_FEF_best_model(i,1:7,:),3)./total_cells_FEF(1:7),'Color',color_for_ROI(2,:),'LineWidth',2)
    hold on
    plot(window_start_list(1:7),nansum(h_PFC_best_model(i,1:7,:),3)./total_cells_PFC(1:7),'Color',color_for_ROI(3,:),'LineWidth',2)
    
    box off
    xlabel('time to end of reward')
    ylabel('Proportion of significant neurons')
    
    xlim([window_start_list(1) window_start_list(7)])
    if i==1
        title('Expected template')
        ylim([0 0.35])
    elseif i==2
        title('Expected value')
        ylim([0 0.25])
    elseif i==3
        title('Mean value')
        ylim([0 0.15])
    elseif i==4
        title('Previous chosen color')
        ylim([0 0.1])
    end
    
    subplot(4,2,2+(i-1)*2)
    plot(window_start_list(8:end),nansum(h_LIP_best_model(i,8:end,:),3)./total_cells_LIP(8:end),'Color',color_for_ROI(1,:),'LineWidth',2)
    hold on
    plot(window_start_list(8:end),nansum(h_FEF_best_model(i,8:end,:),3)./total_cells_FEF(8:end),'Color',color_for_ROI(2,:),'LineWidth',2)
    hold on
    plot(window_start_list(8:end),nansum(h_PFC_best_model(i,8:end,:),3)./total_cells_PFC(8:end),'Color',color_for_ROI(3,:),'LineWidth',2)
    
    xline(-300,'-','Color',[0.5 0.5 0.5],'LineWidth',1)
    xl=xline(-400,'-',{'Fixation'},'Color',[0.5 0.5 0.5],'LineWidth',1);
    if i>2
        xl.LabelVerticalAlignment = 'top';
    else
        xl.LabelVerticalAlignment = 'bottom';
    end
    xl.LabelHorizontalAlignment = 'left';
    xl.FontSize=12;
    
    xline(300,'-','Color',[0.25 0.75 0.25],'LineWidth',1)
    xl=xline(200,'-',{'Response'},'Color',[0.25 0.75 0.25],'LineWidth',1);
    if i>2
        xl.LabelVerticalAlignment = 'top';
    else
        xl.LabelVerticalAlignment = 'bottom';
    end
    xl.LabelHorizontalAlignment = 'left';
    xl.FontSize=12;
    box off
    xlabel('time to targets on')
    ylabel('Proportion of significant neurons')
    ylim([0 0.35])
    xlim([window_start_list(8) window_start_list(end)])
    if i==1
        title('Expected template')
        ylim([0 0.35])
    elseif i==2
        title('Expected value')
        ylim([0 0.25])
    elseif i==3
        title('Mean value')
        ylim([0 0.15])
    elseif i==4
        title('chosen color')
        ylim([0 0.1])
    end
    
    
end



%%
% figure
% for i=1:4
%     subplot(4,2,1+(i-1)*2)
%     plot(window_start_list(1:7),nansum(h_LIP(i,1:7,:),3)./total_cells_LIP(1:7),'Color',color_for_ROI(1,:),'LineWidth',2)
%     hold on
%     plot(window_start_list(1:7),nansum(h_FEF(i,1:7,:),3)./total_cells_FEF(1:7),'Color',color_for_ROI(2,:),'LineWidth',2)
%     hold on
%     plot(window_start_list(1:7),nansum(h_PFC(i,1:7,:),3)./total_cells_PFC(1:7),'Color',color_for_ROI(3,:),'LineWidth',2)
%     
%     box off
%     xlabel('time to end of reward')
%     ylabel('Proportion of significant neurons')
%     
%     xlim([window_start_list(1) window_start_list(7)])
%     if i==1
%         title('Expected template')
%         ylim([0 0.35])
%     elseif i==2
%         title('Expected value')
%         ylim([0 0.25])
%     elseif i==3
%         title('Mean value')
%         ylim([0 0.15])
%     elseif i==4
%         title('Previous chosen color')
%         ylim([0 0.1])
%     end
%     
%     subplot(4,2,2+(i-1)*2)
%     plot(window_start_list(8:end),nansum(h_LIP(i,8:end,:),3)./total_cells_LIP(8:end),'Color',color_for_ROI(1,:),'LineWidth',2)
%     hold on
%     plot(window_start_list(8:end),nansum(h_FEF(i,8:end,:),3)./total_cells_FEF(8:end),'Color',color_for_ROI(2,:),'LineWidth',2)
%     hold on
%     plot(window_start_list(8:end),nansum(h_PFC(i,8:end,:),3)./total_cells_PFC(8:end),'Color',color_for_ROI(3,:),'LineWidth',2)
%     
%     xline(-300,'-','Color',[0.5 0.5 0.5],'LineWidth',1)
%     xl=xline(-400,'-',{'Fixation'},'Color',[0.5 0.5 0.5],'LineWidth',1);
%     if i>2
%         xl.LabelVerticalAlignment = 'top';
%     else
%         xl.LabelVerticalAlignment = 'bottom';
%     end
%     xl.LabelHorizontalAlignment = 'left';
%     xl.FontSize=12;
%     
%     xline(300,'-','Color',[0.25 0.75 0.25],'LineWidth',1)
%     xl=xline(200,'-',{'Response'},'Color',[0.25 0.75 0.25],'LineWidth',1);
%     if i>2
%         xl.LabelVerticalAlignment = 'top';
%     else
%         xl.LabelVerticalAlignment = 'bottom';
%     end
%     xl.LabelHorizontalAlignment = 'left';
%     xl.FontSize=12;
%     box off
%     xlabel('time to targets on')
%     ylabel('Proportion of significant neurons')
%     ylim([0 0.35])
%     xlim([window_start_list(8) window_start_list(end)])
%     if i==1
%         title('Expected template')
%         ylim([0 0.35])
%     elseif i==2
%         title('Expected value')
%         ylim([0 0.25])
%     elseif i==3
%         title('Mean value')
%         ylim([0 0.15])
%     elseif i==4
%         title('chosen color')
%         ylim([0 0.1])
%     end
%     
%     
% end
