%load and bar model stim belief cos

fsroot='/Volumes/buschman';

arrayID_list=[3 8 12 14 18 20 22 26 28 30 32 34 36 38 40 42 44]; % 42

event='target';

window_start_list=-600;
window_size=900;

task='Learning_Attentional_Templates';
% subtask='exploreexploit/Reset_RW_model';
subtask='exploreexploit/Restricted_FEF/Reset_RW_model';

color_for_ROI=[0 0.447 0.741; 0.85 0.325 0.098; 0.929 0.694 0.125]; %LIP FEF PFC

count_LIP=1;
count_FEF=1;
count_PFC=1;

this_time=1;

%% Load

for arrayID_ind=1:length(arrayID_list)
    
    arrayID=arrayID_list(arrayID_ind);
    
    if arrayID~=34
        
        
        subsubtask=sprintf('FR_%s_%d_window_from_%d_to_%d',event,window_size,window_start_list(this_time),window_start_list(this_time)+window_size);
        
        dirstem = fullfile('Projects',task,'Analysis','Electrophy_analysis',subtask,subsubtask);
        data_path = fullfile(fsroot,dirstem);
        
        data_name=sprintf('models_FR_RR_value_function_%s_%d','LIP',arrayID);
        null_dist_name=sprintf('null_distribution_models_FR_RR_value_function_%s_%d','LIP',arrayID);
        
        str=[fullfile(data_path,data_name) '.mat'];
        if exist(str)>0
            
            load(fullfile(data_path,data_name),'R2','param','Neuron_id')
%             load(fullfile(data_path,null_dist_name),'R2_null')
            
            if exist('R2')
                
                for n=1:size(R2,1)
                    
                    R2_LIP(count_LIP,:,:)=R2(n,:,:);
%                     R2_null_LIP(count_LIP,:,:,:)=R2_null(n,:,:,:);
                    param_LIP(count_LIP,:,:,:,:)=param(n,:,:,:,:);
                    LIP_sess_all(count_LIP)=arrayID_ind;
                    LIP_neur_all(count_LIP)=Neuron_id(n);
                    count_LIP=count_LIP+1;
                end
                
                clear R2 param Neuron_id R2_null
                
            end
            
        end
        
    end
    
    subsubtask=sprintf('FR_%s_%d_window_from_%d_to_%d',event,window_size,window_start_list(this_time),window_start_list(this_time)+window_size);
    
    dirstem = fullfile('Projects',task,'Analysis','Electrophy_analysis',subtask,subsubtask);
    data_path = fullfile(fsroot,dirstem);
    
    data_name=sprintf('models_FR_RR_value_function_%s_%d','FEF',arrayID);
    null_dist_name=sprintf('null_distribution_models_FR_RR_value_function_%s_%d','FEF',arrayID);
    
    str=[fullfile(data_path,data_name) '.mat'];
    if exist(str)>0
        
        load(fullfile(data_path,data_name),'R2','param','Neuron_id')
%             load(fullfile(data_path,null_dist_name),'R2_null')
            
        if exist('R2')
            
            for n=1:size(R2,1)
                
                R2_FEF(count_FEF,:,:)=R2(n,:,:);
%                 R2_null_FEF(count_FEF,:,:,:)=R2_null(n,:,:,:);
                param_FEF(count_FEF,:,:,:,:)=param(n,:,:,:,:);
                FEF_sess_all(count_FEF)=arrayID_ind;
                FEF_neur_all(count_FEF)=Neuron_id(n);
                count_FEF=count_FEF+1;
            end
            
            clear R2 param Neuron_id R2_null
            
        end
        
    end
    
    subsubtask=sprintf('FR_%s_%d_window_from_%d_to_%d',event,window_size,window_start_list(this_time),window_start_list(this_time)+window_size);
    
    dirstem = fullfile('Projects',task,'Analysis','Electrophy_analysis',subtask,subsubtask);
    data_path = fullfile(fsroot,dirstem);
    
    data_name=sprintf('models_FR_RR_value_function_%s_%d','PFC',arrayID);
    null_dist_name=sprintf('null_distribution_models_FR_RR_value_function_%s_%d','PFC',arrayID);
    
    str=[fullfile(data_path,data_name) '.mat'];
    if exist(str)>0
        
        load(fullfile(data_path,data_name),'R2','param','Neuron_id')
%             load(fullfile(data_path,null_dist_name),'R2_null')
            
        if exist('R2')
            
            for n=1:size(R2,1)
                
                R2_PFC(count_PFC,:,:)=R2(n,:,:);
%                 R2_null_PFC(count_PFC,:,:,:)=R2_null(n,:,:,:);
                param_PFC(count_PFC,:,:,:,:)=param(n,:,:,:,:);
                PFC_sess_all(count_PFC)=arrayID_ind;
                PFC_neur_all(count_PFC)=Neuron_id(n);
                count_PFC=count_PFC+1;
            end
            
            clear R2 param Neuron_id R2_null
            
        end
        
    end
    
    
end



%% True dist

clear h* p_*

h_LIP=NaN(size(R2_LIP,2),size(R2_LIP,1));
h_FEF=NaN(size(R2_FEF,2),size(R2_FEF,1));
h_PFC=NaN(size(R2_PFC,2),size(R2_PFC,1));

p_LIP=NaN(size(R2_LIP,2),size(R2_LIP,1));
p_FEF=NaN(size(R2_FEF,2),size(R2_FEF,1));
p_PFC=NaN(size(R2_PFC,2),size(R2_PFC,1));

p_LIP_null=NaN(size(R2_LIP,2),size(R2_LIP,1));
p_FEF_null=NaN(size(R2_FEF,2),size(R2_FEF,1));
p_PFC_null=NaN(size(R2_PFC,2),size(R2_PFC,1));

h_LIP_null=NaN(size(R2_LIP,2),size(R2_LIP,1));
h_FEF_null=NaN(size(R2_FEF,2),size(R2_FEF,1));
h_PFC_null=NaN(size(R2_PFC,2),size(R2_PFC,1));

for j=1:size(R2_LIP,2)
    for i=1:size(R2_LIP,1)
        if ~isnan(mean(R2_LIP(i,j,:)))
            if sum(R2_LIP(i,j,:)==0)==size(R2_LIP,3)
                h_LIP(j,i)=NaN;
                p_LIP(j,i)=NaN;
            else
                h_LIP(j,i)=(mean(R2_LIP(i,j,:),3)>0);
                p_LIP(j,i)=(mean(R2_LIP(i,j,:),3)<=0);
%                 p_LIP_null(j,i)=sum(mean(R2_LIP(i,j,:),3)<mean(R2_null_LIP(i,j,:,:),3),4)/size(R2_null_LIP,4);
%                 h_LIP_null(j,i)=(p_LIP_null(j,i)<0.05);
            end
        end
    end
    for i=1:size(R2_FEF,1)
        if ~isnan(mean(R2_FEF(i,j,:)))
            if sum(R2_FEF(i,j,:)==0)==size(R2_FEF,3)
                h_FEF(j,i)=NaN;
                p_FEF(j,i)=NaN;
            else
                h_FEF(j,i)=(mean(R2_FEF(i,j,:),3)>0);
                p_FEF(j,i)=(mean(R2_FEF(i,j,:),3)<=0);
%                 p_FEF_null(j,i)=sum(mean(R2_FEF(i,j,:),3)<mean(R2_null_FEF(i,j,:,:),3),4)/size(R2_null_FEF,4);
%                 h_FEF_null(j,i)=(p_FEF_null(j,i)<0.05);
            end
        end
    end
    for i=1:size(R2_PFC,1)
        if ~isnan(mean(R2_PFC(i,j,:)))
            if sum(R2_PFC(i,j,:)==0)==size(R2_PFC,3)
                h_PFC(j,i)=NaN;
                p_PFC(j,i)=NaN;
            else
                h_PFC(j,i)=(mean(R2_PFC(i,j,:),3)>0);
                p_PFC(j,i)=(mean(R2_PFC(i,j,:),3)<=0);
%                 p_PFC_null(j,i)=sum(mean(R2_PFC(i,j,:),3)<mean(R2_null_PFC(i,j,:,:),3),4)/size(R2_null_PFC,4);
%                 h_PFC_null(j,i)=(p_PFC_null(j,i)<0.05);
            end
        end
    end
end

total_cells_LIP=sum(~isnan(h_LIP(1,:)));
total_cells_FEF=sum(~isnan(h_FEF(1,:)));
total_cells_PFC=sum(~isnan(h_PFC(1,:)));

h_LIP_best_model=zeros(4,size(h_LIP,3));
h_FEF_best_model=zeros(4,size(h_FEF,3));
h_PFC_best_model=zeros(4,size(h_PFC,3));

for i=1:size(R2_LIP,1)
    if sum(h_LIP(:,i))>=1
        [~, ind]=max(mean(R2_LIP(i,:,:),3),[],2);
        h_LIP_best_model(ind,i)=1;
    end
end
for i=1:size(R2_FEF,1)
    if sum(h_FEF(:,i))>=1
        [~, ind]=max(mean(R2_FEF(i,:,:),3),[],2);
        h_FEF_best_model(ind,i)=1;
    end
end
for i=1:size(R2_PFC,1)
    if sum(h_PFC(:,i))>=1
        [~, ind]=max(mean(R2_PFC(i,:,:),3),[],2);
        h_PFC_best_model(ind,i)=1;
    end
end

% h_LIP_null_best_model=zeros(4,size(h_LIP_null,3));
% h_FEF_null_best_model=zeros(4,size(h_FEF_null,3));
% h_PFC_null_best_model=zeros(4,size(h_PFC_null,3));
% 
% for i=1:size(R2_LIP,1)
%     if sum(h_LIP_null(:,i))>=1
%         [~, ind]=max(mean(R2_LIP(i,:,:),3),[],2);
%         h_LIP_null_best_model(ind,i)=1;
%     end
% end
% for i=1:size(R2_FEF,1)
%     if sum(h_FEF_null(:,i))>=1
%         [~, ind]=max(mean(R2_FEF(i,:,:),3),[],2);
%         h_FEF_null_best_model(ind,i)=1;
%     end
% end
% for i=1:size(R2_PFC,1)
%     if sum(h_PFC_null(:,i))>=1
%         [~, ind]=max(mean(R2_PFC(i,:,:),3),[],2);
%         h_PFC_null_best_model(ind,i)=1;
%     end
% end


%% Null dist
% 
% h_LIP_null_dist=NaN(size(R2_LIP,2),size(R2_LIP,1),size(R2_null_PFC,4));
% h_FEF_null_dist=NaN(size(R2_FEF,2),size(R2_FEF,1),size(R2_null_PFC,4));
% h_PFC_null_dist=NaN(size(R2_PFC,2),size(R2_PFC,1),size(R2_null_PFC,4));
% 
% for np=1:size(R2_null_PFC,4)
%     for j=1:size(R2_LIP,2)
%         for i=1:size(R2_LIP,1)
%             if ~isnan(mean(R2_null_LIP(i,j,:,np)))
%                 if sum(R2_null_LIP(i,j,:,np)==0)==size(R2_LIP,3)
%                     h_LIP_null_dist(j,i,np)=NaN;
%                 else
%                     h_LIP_null_dist(j,i,np)=(mean(R2_null_LIP(i,j,:,np),3)>0);
%                 end
%             end
%         end
%         for i=1:size(R2_FEF,1)
%             if ~isnan(mean(R2_null_FEF(i,j,:,np)))
%                 if sum(R2_null_FEF(i,j,:,np)==0)==size(R2_FEF,3)
%                     h_FEF_null_dist(j,i,np)=NaN;
%                 else
%                     h_FEF_null_dist(j,i,np)=(mean(R2_null_FEF(i,j,:,np),3)>0);
%                 end
%             end
%         end
%         for i=1:size(R2_PFC,1)
%             if ~isnan(mean(R2_null_PFC(i,j,:,np)))
%                 if sum(R2_null_PFC(i,j,:,np)==0)==size(R2_PFC,3)
%                     h_PFC_null_dist(j,i,np)=NaN;
%                 else
%                     h_PFC_null_dist(j,i,np)=(mean(R2_null_PFC(i,j,:,np),3)>0);
%                 end
%             end
%         end
%     end
% end
% 
% h_LIP_best_model_null_dist=zeros(4,size(h_LIP,3),size(R2_null_PFC,4));
% h_FEF_best_model_null_dist=zeros(4,size(h_FEF,3),size(R2_null_PFC,4));
% h_PFC_best_model_null_dist=zeros(4,size(h_PFC,3),size(R2_null_PFC,4));
% 
% for np=1:size(R2_null_PFC,4)
%     for i=1:size(R2_LIP,1)
%         if sum(h_LIP_null_dist(:,i,np))>=1
%             [~, ind]=max(mean(R2_null_LIP(i,:,:,np),3),[],2);
%             h_LIP_best_model_null_dist(ind,i,np)=1;
%         end
%     end
%     for i=1:size(R2_FEF,1)
%         if sum(h_FEF_null_dist(:,i,np))>=1
%             [~, ind]=max(mean(R2_null_FEF(i,:,:,np),3),[],2);
%             h_FEF_best_model_null_dist(ind,i,np)=1;
%         end
%     end
%     for i=1:size(R2_PFC,1)
%         if sum(h_PFC_null_dist(:,i,np))>=1
%             [~, ind]=max(mean(R2_null_PFC(i,:,:,np),3),[],2);
%             h_PFC_best_model_null_dist(ind,i,np)=1;
%         end
%     end
% end
% 
% Proba_LIP_best_model_null_dist(:,:)=nansum(h_LIP_best_model_null_dist,2)./total_cells_LIP;
% Proba_FEF_best_model_null_dist(:,:)=nansum(h_FEF_best_model_null_dist,2)./total_cells_FEF;
% Proba_PFC_best_model_null_dist(:,:)=nansum(h_PFC_best_model_null_dist,2)./total_cells_PFC;
% 
% for i=1:4
%     Null_dist_LIP_best_model(i,1)=prctile(Proba_LIP_best_model_null_dist(i,:),5,2);
%     Null_dist_LIP_best_model(i,2)=prctile(Proba_LIP_best_model_null_dist(i,:),95,2);
%     
%     Null_dist_FEF_best_model(i,1)=prctile(Proba_FEF_best_model_null_dist(i,:),5,2);
%     Null_dist_FEF_best_model(i,2)=prctile(Proba_FEF_best_model_null_dist(i,:),95,2);
%     
%     Null_dist_PFC_best_model(i,1)=prctile(Proba_PFC_best_model_null_dist(i,:),5,2);
%     Null_dist_PFC_best_model(i,2)=prctile(Proba_PFC_best_model_null_dist(i,:),95,2);
% end
% 
% 
% 
% 
% %% Best model
% 
% % figure
% % subplot(3,1,1)
% % bar(nansum(h_LIP_best_model,2)./total_cells_LIP,'FaceColor',color_for_ROI(1,:),'LineWidth',1)
% % yline(0.05,'--')
% % xticks([1 2 3 4])
% % xticklabels({'ET','EV','MV','CC'})
% % box off
% % ylabel('Proportion of significant neurons')
% % ylim([0 0.5])
% % subplot(3,1,2)
% % bar(nansum(h_FEF_best_model,2)./total_cells_FEF,'FaceColor',color_for_ROI(2,:),'LineWidth',1)
% % yline(0.05,'--')
% % xticks([1 2 3 4])
% % xticklabels({'ET','EV','MV','CC'})
% % box off
% % ylabel('Proportion of significant neurons')
% % ylim([0 0.5])
% % subplot(3,1,3)
% % bar(nansum(h_PFC_best_model,2)./total_cells_PFC,'FaceColor',color_for_ROI(3,:),'LineWidth',1)
% % yline(0.05,'--')
% % xticks([1 2 3 4])
% % xticklabels({'ET','EV','MV','CC'})
% % box off
% % ylabel('Proportion of significant neurons')
% % ylim([0 0.5])
% 
% % %% Best model with stats
% % 
% % figure
% % subplot(3,1,1)
% % bar(nansum(h_LIP_null_best_model,2)./total_cells_LIP,'FaceColor',color_for_ROI(1,:),'LineWidth',1)
% % yline(0.05,'--')
% % xticks([1 2 3 4])
% % xticklabels({'ET','EV','MV','CC'})
% % box off
% % ylabel('Proportion of significant neurons')
% % ylim([0 0.5])
% % subplot(3,1,2)
% % bar(nansum(h_FEF_null_best_model,2)./total_cells_FEF,'FaceColor',color_for_ROI(2,:),'LineWidth',1)
% % yline(0.05,'--')
% % xticks([1 2 3 4])
% % xticklabels({'ET','EV','MV','CC'})
% % box off
% % ylabel('Proportion of significant neurons')
% % ylim([0 0.5])
% % subplot(3,1,3)
% % bar(nansum(h_PFC_null_best_model,2)./total_cells_PFC,'FaceColor',color_for_ROI(3,:),'LineWidth',1)
% % yline(0.05,'--')
% % xticks([1 2 3 4])
% % xticklabels({'ET','EV','MV','CC'})
% % box off
% % ylabel('Proportion of significant neurons')
% % ylim([0 0.5])

%% Best model with null dist

% p_LIP_test_against_null=(1+sum(nansum(h_LIP_best_model,2)./total_cells_LIP<Proba_LIP_best_model_null_dist,2))/(1+size(Proba_LIP_best_model_null_dist,2));
% p_FEF_test_against_null=(1+sum(nansum(h_FEF_best_model,2)./total_cells_FEF<Proba_FEF_best_model_null_dist,2))/(1+size(Proba_FEF_best_model_null_dist,2));
% p_PFC_test_against_null=(1+sum(nansum(h_PFC_best_model,2)./total_cells_PFC<Proba_PFC_best_model_null_dist,2))/(1+size(Proba_PFC_best_model_null_dist,2));

figure
subplot(1,3,1)
bar(nansum(h_LIP_best_model,2)./total_cells_LIP,'FaceColor',color_for_ROI(1,:),'LineWidth',1)
xticks([1 2 3 4 5])
xticklabels({'ET','EV','MV','PCC','WSLF'})
box off
ylabel('Proportion of significant neurons')
ylim([0 0.42])
subplot(1,3,2)
bar(nansum(h_FEF_best_model,2)./total_cells_FEF,'FaceColor',color_for_ROI(2,:),'LineWidth',1)
xticks([1 2 3 4 5])
xticklabels({'ET','EV','MV','PCC','WSLF'})
box off
ylabel('Proportion of significant neurons')
ylim([0 0.42])
subplot(1,3,3)
bar(nansum(h_PFC_best_model,2)./total_cells_PFC,'FaceColor',color_for_ROI(3,:),'LineWidth',1)
xticks([1 2 3 4 5])
xticklabels({'ET','EV','MV','PCC','WSLF'})
box off
ylabel('Proportion of significant neurons')
ylim([0 0.42])


% figure
% for i=1:4
%     subplot(1,4,i)
%     bar(1,nansum(h_LIP_best_model(i,:))./total_cells_LIP,'FaceColor',color_for_ROI(1,:),'LineWidth',2)
%     hold on
%     bar(2,nansum(h_FEF_best_model(i,:))./total_cells_FEF,'FaceColor',color_for_ROI(2,:),'LineWidth',2)
%     hold on
%     bar(3,nansum(h_PFC_best_model(i,:))./total_cells_PFC,'FaceColor',color_for_ROI(3,:),'LineWidth',2)
%     yline(0.05,'--')
%     xticks([1 2 3])
%     xticklabels({'LIP','FEF','PFC'})
%     box off
%     ylabel('Proportion of significant neurons')
%     ylim([0 0.5])
%     if i==1
%         title('Expected template')
%     elseif i==2
%         title('Expected value')
%     elseif i==3
%         title('Mean value')
%     elseif i==4
%         title('Chosen color')
%     end
% 
% end




