arrayID_list=[3 8 12 14 18 20 22 26 28 30 32 34 36 38 40 42 44];

fsroot='/Volumes/buschman';

event='target';

for i=1:20
    window_start_list(i)=-600+(i-1)*50;
end

target_on=13;
response_on=17;

task='Learning_Attentional_Templates';
subtask='exploreexploit/Restricted_FEF/Reset_RW_model';

window_size=300;

subsubtask_classifier=sprintf('PCA_%s_%d_to_%s_%d',event,window_start_list(1),event,window_start_list(end));

dirstem = fullfile('Projects',task,'Analysis','Electrophy_analysis',subtask,subsubtask_classifier);
save_path_clasifier = fullfile(fsroot,dirstem);
mkdir(save_path_clasifier);

save_name='PCA_expected_template';

color_for_ROI=[0 0.447 0.741; 0.85 0.325 0.098; 0.929 0.694 0.125]; %LIP FEF PFC


N_channels_stim=3;
N_bins=N_channels_stim;

count_LIP=1;
count_FEF=1;
count_PFC=1;

N_sessions=length(arrayID_list);

for n_sess=1:N_sessions
    
    arrayID=arrayID_list(n_sess);
    
    %LIP
    if arrayID~=34 && arrayID~=30
        
        %         subsubtask_data=sprintf('FR_%s_%d_window_from_%d_to_%d',event_list{1},window_size,window_start_list(1),window_start_list(1)+window_size);
        subsubtask_data=sprintf('FR_%s_%d_window_from_%d_to_%d',event,window_size,window_start_list(1),window_start_list(1)+window_size);
        data_name=sprintf('ID_%d.mat',arrayID);
        dirstem = fullfile('Projects',task,'Analysis','Electrophy_analysis',subtask,subsubtask_data);
        data_path_clasifier = fullfile(fsroot,dirstem);
        
        load(fullfile(data_path_clasifier,data_name),'Conditions','LIP')
        
        %get the progression in trial
        for i=1:length(Conditions.block_nb)
            Conditions.Progression_in_block(i)=Conditions.trial_in_block(i)/sum(Conditions.block_nb==Conditions.block_nb(i));
        end
        %get the peak belief
        N_bins=100;
        color_binned =  0:(2*pi/N_bins):(2*pi-(2*pi/N_bins));
        
        for i=1:size(Conditions.stim,2)
            [~, peak_belief_index(i)]=max(Conditions.belief(:,i));
        end
        Conditions.Peak_belief=color_binned(peak_belief_index);
        
        for i=1:size(Conditions.stim,2)
            if sum(Conditions.belief(:,i)==0)==size(Conditions.belief,1)
                Conditions.Peak_belief(i)=NaN;
            end
        end
        for t=1:length(window_start_list)
            
            %             event=event_list{t};
            window_start=window_start_list(t);
            
            subsubtask_data=sprintf('FR_%s_%d_window_from_%d_to_%d',event,window_size,window_start,window_start+window_size);
            data_name=sprintf('ID_%d.mat',arrayID);
            dirstem = fullfile('Projects',task,'Analysis','Electrophy_analysis',subtask,subsubtask_data);
            data_path_clasifier = fullfile(fsroot,dirstem);
            
            load(fullfile(data_path_clasifier,data_name),'LIP')
            
            N_neurons=size(LIP,1);
            
            for nn=1:N_neurons
                for i=1:length(Conditions.block_nb)
                    Raw_data_LIP(nn,i,t)=LIP(nn,i);
                end
            end
            clear LIP
        end
        
        N_neurons=size(Raw_data_LIP,1);
        
        %remove NaN if neuron was active at some point
        for nn=1:N_neurons
            for i=1:size(Raw_data_LIP,2)
                if sum(~isnan(Raw_data_LIP(nn,i,:)))>0
                    Raw_data_LIP(nn,i,isnan(Raw_data_LIP(nn,i,:)))=0;
                end
            end
        end
        %remove partially active neurons
        data_LIP=Raw_data_LIP(~isnan(sum(Raw_data_LIP(:,:,1),2)),:,:); %remove partially active neurons
        clear Raw_data_LIP
        for nn=1:size(data_LIP,1)
            
            %first get the conditions
            for i=2:length(Conditions.block_nb)
                for t=1:length(window_start_list)
                    Pseudo_pop_LIP(count_LIP,i-1).Classifier_FR(t)=data_LIP(nn,i,t);
                end
                Pseudo_pop_LIP(count_LIP,i-1).Sess=n_sess;
                Pseudo_pop_LIP(count_LIP,i-1).Block_nb=Conditions.block_nb(i); %we need to take the updated belief here
                Pseudo_pop_LIP(count_LIP,i-1).Precision=Conditions.belief_precision(i);
                Pseudo_pop_LIP(count_LIP,i-1).Peak_belief=mod(Conditions.Peak_belief(i),2*pi);
                Pseudo_pop_LIP(count_LIP,i-1).progression_in_block=Conditions.Progression_in_block(i);
                Pseudo_pop_LIP(count_LIP,i-1).stim(:)=Conditions.stim(:,i);
                for n_c=1:N_channels_stim
                    if mod(Pseudo_pop_LIP(count_LIP,i-1).Peak_belief,2*pi)>=(n_c-1)*2*pi/N_channels_stim && mod(Pseudo_pop_LIP(count_LIP,i-1).Peak_belief,2*pi)<n_c*2*pi/N_channels_stim
                        Pseudo_pop_LIP(count_LIP,i-1).Belief_color=n_c;
                    end
                end
                
                for j=1:4
                    if isnan(Pseudo_pop_LIP(count_LIP,i-1).stim(j)) && Conditions.choice(i)~=j
                        Pseudo_pop_LIP(count_LIP,i-1).Stim_color(j)=NaN;
                    else
                        for n_c=1:N_channels_stim
                            if mod(Pseudo_pop_LIP(count_LIP,i-1).stim(j),2*pi)>=(n_c-1)*2*pi/N_channels_stim && mod(Pseudo_pop_LIP(count_LIP,i-1).stim(j),2*pi)<n_c*2*pi/N_channels_stim
                                Pseudo_pop_LIP(count_LIP,i-1).Stim_color(j)=n_c;
                            end
                        end
                    end
                end
                
            end
            
            count_LIP=count_LIP+1;
        end
        
        clear Conditions
        
    end
    
    %FEF
    if arrayID~=14
        %         subsubtask_data=sprintf('FR_%s_%d_window_from_%d_to_%d',event_list{1},window_size,window_start_list(1),window_start_list(1)+window_size);
        subsubtask_data=sprintf('FR_%s_%d_window_from_%d_to_%d',event,window_size,window_start_list(1),window_start_list(1)+window_size);
        data_name=sprintf('ID_%d.mat',arrayID);
        dirstem = fullfile('Projects',task,'Analysis','Electrophy_analysis',subtask,subsubtask_data);
        data_path_clasifier = fullfile(fsroot,dirstem);
        
        load(fullfile(data_path_clasifier,data_name),'Conditions','FEF')
        
        %get the progression in trial
        for i=1:length(Conditions.block_nb)
            Conditions.Progression_in_block(i)=Conditions.trial_in_block(i)/sum(Conditions.block_nb==Conditions.block_nb(i));
        end
        %get the peak belief
        N_bins=100;
        color_binned =  0:(2*pi/N_bins):(2*pi-(2*pi/N_bins));
        
        for i=1:size(Conditions.stim,2)
            [~, peak_belief_index(i)]=max(Conditions.belief(:,i));
        end
        Conditions.Peak_belief=color_binned(peak_belief_index);
        
        for i=1:size(Conditions.stim,2)
            if sum(Conditions.belief(:,i)==0)==size(Conditions.belief,1)
                Conditions.Peak_belief(i)=NaN;
            end
        end
        for t=1:length(window_start_list)
            
            %         event=event_list{t};
            window_start=window_start_list(t);
            
            subsubtask_data=sprintf('FR_%s_%d_window_from_%d_to_%d',event,window_size,window_start,window_start+window_size);
            data_name=sprintf('ID_%d.mat',arrayID);
            dirstem = fullfile('Projects',task,'Analysis','Electrophy_analysis',subtask,subsubtask_data);
            data_path_clasifier = fullfile(fsroot,dirstem);
            
            load(fullfile(data_path_clasifier,data_name),'FEF')
            
            N_neurons=size(FEF,1);
            
            for nn=1:N_neurons
                for i=1:length(Conditions.block_nb)
                    Raw_data_FEF(nn,i,t)=FEF(nn,i);
                end
            end
            
            clear FEF
        end
        
        N_neurons=size(Raw_data_FEF,1);
        
        %remove NaN if neuron was active at some point
        for nn=1:N_neurons
            for i=1:size(Raw_data_FEF,2)
                if sum(~isnan(Raw_data_FEF(nn,i,:)))>0
                    Raw_data_FEF(nn,i,isnan(Raw_data_FEF(nn,i,:)))=0;
                end
            end
        end
        %remove partially active neurons
        data_FEF=Raw_data_FEF(~isnan(sum(Raw_data_FEF(:,:,1),2)),:,:); %remove partially active neurons
        clear Raw_data_FEF
        for nn=1:size(data_FEF,1)
            
            %first get the conditions
            for i=2:length(Conditions.block_nb)
                for t=1:length(window_start_list)
                    Pseudo_pop_FEF(count_FEF,i-1).Classifier_FR(t)=data_FEF(nn,i,t);
                end
                Pseudo_pop_FEF(count_FEF,i-1).Sess=n_sess;
                Pseudo_pop_FEF(count_FEF,i-1).Block_nb=Conditions.block_nb(i); %we need to take the updated belief here
                Pseudo_pop_FEF(count_FEF,i-1).Precision=Conditions.belief_precision(i);
                Pseudo_pop_FEF(count_FEF,i-1).Peak_belief=mod(Conditions.Peak_belief(i),2*pi);
                Pseudo_pop_FEF(count_FEF,i-1).progression_in_block=Conditions.Progression_in_block(i);
                Pseudo_pop_FEF(count_FEF,i-1).stim(:)=Conditions.stim(:,i);
                for n_c=1:N_channels_stim
                    if mod(Pseudo_pop_FEF(count_FEF,i-1).Peak_belief,2*pi)>=(n_c-1)*2*pi/N_channels_stim && mod(Pseudo_pop_FEF(count_FEF,i-1).Peak_belief,2*pi)<n_c*2*pi/N_channels_stim
                        Pseudo_pop_FEF(count_FEF,i-1).Belief_color=n_c;
                    end
                end
                
                for j=1:4
                    if isnan(Pseudo_pop_FEF(count_FEF,i-1).stim(j))  && Conditions.choice(i)~=j
                        Pseudo_pop_FEF(count_FEF,i-1).Stim_color(j)=NaN;
                    else
                        for n_c=1:N_channels_stim
                            if mod(Pseudo_pop_FEF(count_FEF,i-1).stim(j),2*pi)>=(n_c-1)*2*pi/N_channels_stim && mod(Pseudo_pop_FEF(count_FEF,i-1).stim(j),2*pi)<n_c*2*pi/N_channels_stim
                                Pseudo_pop_FEF(count_FEF,i-1).Stim_color(j)=n_c;
                            end
                        end
                    end
                end
                
            end
            
            count_FEF=count_FEF+1;
        end
        
        clear Conditions
    end
    
    %PFC
    %         subsubtask_data=sprintf('FR_%s_%d_window_from_%d_to_%d',event_list{1},window_size,window_start_list(1),window_start_list(1)+window_size);
    subsubtask_data=sprintf('FR_%s_%d_window_from_%d_to_%d',event,window_size,window_start_list(1),window_start_list(1)+window_size);
    data_name=sprintf('ID_%d.mat',arrayID);
    dirstem = fullfile('Projects',task,'Analysis','Electrophy_analysis',subtask,subsubtask_data);
    data_path_clasifier = fullfile(fsroot,dirstem);
    
    load(fullfile(data_path_clasifier,data_name),'Conditions','PFC')
    
    %get the progression in trial
    for i=1:length(Conditions.block_nb)
        Conditions.Progression_in_block(i)=Conditions.trial_in_block(i)/sum(Conditions.block_nb==Conditions.block_nb(i));
    end
    %get the peak belief
    N_bins=100;
    color_binned =  0:(2*pi/N_bins):(2*pi-(2*pi/N_bins));
    
    for i=1:size(Conditions.stim,2)
        [~, peak_belief_index(i)]=max(Conditions.belief(:,i));
    end
    Conditions.Peak_belief=color_binned(peak_belief_index);
    
    for i=1:size(Conditions.stim,2)
        if sum(Conditions.belief(:,i)==0)==size(Conditions.belief,1)
            Conditions.Peak_belief(i)=NaN;
        end
    end
    for t=1:length(window_start_list)
        
        %         event=event_list{t};
        window_start=window_start_list(t);
        
        subsubtask_data=sprintf('FR_%s_%d_window_from_%d_to_%d',event,window_size,window_start,window_start+window_size);
        data_name=sprintf('ID_%d.mat',arrayID);
        dirstem = fullfile('Projects',task,'Analysis','Electrophy_analysis',subtask,subsubtask_data);
        data_path_clasifier = fullfile(fsroot,dirstem);
        
        load(fullfile(data_path_clasifier,data_name),'PFC')
        
        N_neurons=size(PFC,1);
        
        start_ind=2;
        end_ind=length(Conditions.block_nb);
        
        for nn=1:N_neurons
            for i=1:length(Conditions.block_nb)
                Raw_data_PFC(nn,i,t)=PFC(nn,i);
            end
        end
        
        clear PFC
    end
    
    N_neurons=size(Raw_data_PFC,1);
    
    %remove NaN if neuron was active at some point
    for nn=1:N_neurons
        for i=1:size(Raw_data_PFC,2)
            if sum(~isnan(Raw_data_PFC(nn,i,:)))>0
                Raw_data_PFC(nn,i,isnan(Raw_data_PFC(nn,i,:)))=0;
            end
        end
    end
    %remove partially active neurons
    data_PFC=Raw_data_PFC(~isnan(sum(Raw_data_PFC(:,:,1),2)),:,:); %remove partially active neurons
    clear Raw_data_PFC
    for nn=1:size(data_PFC,1)
        
        %first get the conditions
        for i=2:length(Conditions.block_nb)
            for t=1:length(window_start_list)
                Pseudo_pop_PFC(count_PFC,i-1).Classifier_FR(t)=data_PFC(nn,i,t);
            end
            Pseudo_pop_PFC(count_PFC,i-1).Sess=n_sess;
            Pseudo_pop_PFC(count_PFC,i-1).Block_nb=Conditions.block_nb(i); %we need to take the updated belief here
            Pseudo_pop_PFC(count_PFC,i-1).Precision=Conditions.belief_precision(i);
            Pseudo_pop_PFC(count_PFC,i-1).Peak_belief=mod(Conditions.Peak_belief(i),2*pi);
            Pseudo_pop_PFC(count_PFC,i-1).progression_in_block=Conditions.Progression_in_block(i);
            Pseudo_pop_PFC(count_PFC,i-1).stim(:)=Conditions.stim(:,i);
            for n_c=1:N_channels_stim
                if mod(Pseudo_pop_PFC(count_PFC,i-1).Peak_belief,2*pi)>=(n_c-1)*2*pi/N_channels_stim && mod(Pseudo_pop_PFC(count_PFC,i-1).Peak_belief,2*pi)<n_c*2*pi/N_channels_stim
                    Pseudo_pop_PFC(count_PFC,i-1).Belief_color=n_c;
                end
            end
            
            for j=1:4
                if isnan(Pseudo_pop_PFC(count_PFC,i-1).stim(j)) && Conditions.choice(i)~=j
                    Pseudo_pop_PFC(count_PFC,i-1).Stim_color(j)=NaN;
                else
                    for n_c=1:N_channels_stim
                        if mod(Pseudo_pop_PFC(count_PFC,i-1).stim(j),2*pi)>=(n_c-1)*2*pi/N_channels_stim && mod(Pseudo_pop_PFC(count_PFC,i-1).stim(j),2*pi)<n_c*2*pi/N_channels_stim
                            Pseudo_pop_PFC(count_PFC,i-1).Stim_color(j)=n_c;
                        end
                    end
                end
            end
            
        end
        
        count_PFC=count_PFC+1;
    end
    
    clear Conditions
    
end

%% Now create the matrix N x 9 conditions x time x 3 prog x 4 loc

N_prog=1;


for nn=1:size(Pseudo_pop_LIP,1)
    for i=1:size(Pseudo_pop_LIP(nn,:),2)
        if ~isempty(Pseudo_pop_LIP(nn,i).progression_in_block) && ~isnan(Pseudo_pop_LIP(nn,i).progression_in_block)
            if Pseudo_pop_LIP(nn,i).progression_in_block<=1/N_prog
                selected_LIP(nn,i,Pseudo_pop_LIP(nn,i).Belief_color)=1;
            end
        end
    end
    Sess_LIP(nn)= Pseudo_pop_LIP(nn,1).Sess;
end

for nn=1:size(Pseudo_pop_FEF,1)
    for i=1:size(Pseudo_pop_FEF(nn,:),2)
        if ~isempty(Pseudo_pop_FEF(nn,i).progression_in_block) && ~isnan(Pseudo_pop_FEF(nn,i).progression_in_block)
            if Pseudo_pop_FEF(nn,i).progression_in_block<=1/N_prog
                selected_FEF(nn,i,Pseudo_pop_FEF(nn,i).Belief_color)=1;
            end
        end
    end
    Sess_FEF(nn)= Pseudo_pop_FEF(nn,1).Sess;
end

for nn=1:size(Pseudo_pop_PFC,1)
    for i=1:size(Pseudo_pop_PFC(nn,:),2)
        if ~isempty(Pseudo_pop_PFC(nn,i).progression_in_block) && ~isnan(Pseudo_pop_PFC(nn,i).progression_in_block)
            
            if Pseudo_pop_PFC(nn,i).progression_in_block<=1/N_prog
                selected_PFC(nn,i,Pseudo_pop_PFC(nn,i).Belief_color)=1;
            end
        end
    end
    Sess_PFC(nn)= Pseudo_pop_PFC(nn,1).Sess;
end
%% Now create the matrix N x 9 conditions x time x 3 prog x 4 loc


for n_c=1:N_channels_stim
    for t=1:length(window_start_list)
        %LIP
        for nn=1:size(Pseudo_pop_LIP,1)
            if sum(selected_LIP(nn,:,n_c))>1
                list=find(selected_LIP(nn,:,n_c)==1);
                for i=1:length(list)
                    buffer(i)=Pseudo_pop_LIP(nn,list(i)).Classifier_FR(t);
                end
                LIP_PCA_belief(nn,n_c,t)=mean(buffer);
                clear buffer list
            else
                LIP_PCA_belief(nn,n_c,t)=NaN;
            end
        end
        %FEF
        for nn=1:size(Pseudo_pop_FEF,1)
            if sum(selected_FEF(nn,:,n_c))>1
                list=find(selected_FEF(nn,:,n_c)==1);
                for i=1:length(list)
                    buffer(i)=Pseudo_pop_FEF(nn,list(i)).Classifier_FR(t);
                end
                FEF_PCA_belief(nn,n_c,t)=mean(buffer);
                clear buffer list
            else
                FEF_PCA_belief(nn,n_c,t)=NaN;
            end
        end
        %PFC
        for nn=1:size(Pseudo_pop_PFC,1)
            if sum(selected_PFC(nn,:,n_c))>1
                list=find(selected_PFC(nn,:,n_c)==1);
                for i=1:length(list)
                    buffer(i)=Pseudo_pop_PFC(nn,list(i)).Classifier_FR(t);
                end
                PFC_PCA_belief(nn,n_c,t)=mean(buffer);
                clear buffer list
            else
                PFC_PCA_belief(nn,n_c,t)=NaN;
            end
        end
    end
end

%% remove NaN

LIP_PCA_belief=LIP_PCA_belief(~isnan(sum(sum(LIP_PCA_belief,2),3)),:,:);
FEF_PCA_belief=FEF_PCA_belief(~isnan(sum(sum(FEF_PCA_belief,2),3)),:,:);
PFC_PCA_belief=PFC_PCA_belief(~isnan(sum(sum(PFC_PCA_belief,2),3)),:,:);

%% reshape for the PCA


LIP_PCA_belief_for_pca(:,:)=reshape(LIP_PCA_belief(:,:,:),size(LIP_PCA_belief,1),size(LIP_PCA_belief,2)*size(LIP_PCA_belief,3));
FEF_PCA_belief_for_pca(:,:)=reshape(FEF_PCA_belief(:,:,:),size(FEF_PCA_belief,1),size(FEF_PCA_belief,2)*size(FEF_PCA_belief,3));
PFC_PCA_belief_for_pca(:,:)=reshape(PFC_PCA_belief(:,:,:),size(PFC_PCA_belief,1),size(PFC_PCA_belief,2)*size(PFC_PCA_belief,3));


%% do the pca

[coeff_LIP(:,:), score_LIP(:,:), eigenvalues_LIP,~,explained_LIP(:),~] = pca(LIP_PCA_belief_for_pca(:,:)');
[coeff_FEF(:,:), score_FEF(:,:), eigenvalues_FEF,~,explained_FEF(:),~] = pca(FEF_PCA_belief_for_pca(:,:)');
[coeff_PFC(:,:), score_PFC(:,:), eigenvalues_PFC,~,explained_PFC(:),~] = pca(PFC_PCA_belief_for_pca(:,:)');

%% reshape or platting

for k=1:4
    for pc=1:size(score_LIP,2)
        PC_LIP(:,:,pc)=reshape(score_LIP(:,pc),size(LIP_PCA_belief,2),size(LIP_PCA_belief,3));
        PC_FEF(:,:,pc)=reshape(score_FEF(:,pc),size(FEF_PCA_belief,2),size(FEF_PCA_belief,3));
        PC_PFC(:,:,pc)=reshape(score_PFC(:,pc),size(PFC_PCA_belief,2),size(PFC_PCA_belief,3));
    end
end

%% now plot

load('colors.mat')
color_bin(1,:)=colors(1,:);
color_bin(2,:)=colors(14,:);
color_bin(3,:)=colors(27,:);

figure

% PC_1=2;
% PC_2=3;
% PC_3=4;

PC_1=2;
PC_2=3;
PC_3=4;


subplot(3,1,1)
hold on
for n_c=1:N_channels_stim
    
    x(:)=PC_LIP(n_c,:,PC_1);
    y(:)=PC_LIP(n_c,:,PC_2);
    z(:)=PC_LIP(n_c,:,PC_3);
    if n_c==1
        plot3(x,y,z,'-','MarkerFaceColor',color_bin(n_c,:),'MarkerEdgeColor',color_bin(n_c,:),'Color',color_bin(n_c,:),'LineWidth',2)
    elseif n_c==2
        plot3(x,y,z,'-','MarkerFaceColor',color_bin(n_c,:),'MarkerEdgeColor',color_bin(n_c,:),'Color',color_bin(n_c,:),'LineWidth',2)
    else
        plot3(x,y,z,'-','MarkerFaceColor',color_bin(n_c,:),'MarkerEdgeColor',color_bin(n_c,:),'Color',color_bin(n_c,:),'LineWidth',2)
    end
    plot3(PC_LIP(n_c,target_on,PC_1),PC_LIP(n_c,target_on,PC_2),PC_LIP(n_c,target_on,PC_3),'^','MarkerSize',7,'MarkerFaceColor',[0.75 0.75 0.75],'MarkerEdgeColor',[0.75 0.75 0.75])
    hold on
    plot3(PC_LIP(n_c,response_on,PC_1),PC_LIP(n_c,response_on,PC_2),PC_LIP(n_c,response_on,PC_3),'^','MarkerSize',7,'MarkerFaceColor',[0.25 0.75 0.25],'MarkerEdgeColor',[0.25 0.75 0.25])
    grid on
    clear x y z
    
end
title(sprintf('LIP '))
xlabel(sprintf('PC %d',PC_1))
ylabel(sprintf('PC %d',PC_2))
zlabel(sprintf('PC %d',PC_3))
% view([-1 -1 1])
view([-1 -1 1])

subplot(3,1,2)
hold on
for n_c=1:N_channels_stim
    x(:)=PC_FEF(n_c,:,PC_1);
    y(:)=PC_FEF(n_c,:,PC_2);
    z(:)=PC_FEF(n_c,:,PC_3);
    if n_c==1
        plot3(x,y,z,'-','MarkerFaceColor',color_bin(n_c,:),'MarkerEdgeColor',color_bin(n_c,:),'Color',color_bin(n_c,:),'LineWidth',2)
    elseif n_c==2
        plot3(x,y,z,'-','MarkerFaceColor',color_bin(n_c,:),'MarkerEdgeColor',color_bin(n_c,:),'Color',color_bin(n_c,:),'LineWidth',2)
    else
        plot3(x,y,z,'-','MarkerFaceColor',color_bin(n_c,:),'MarkerEdgeColor',color_bin(n_c,:),'Color',color_bin(n_c,:),'LineWidth',2)
    end
    plot3(PC_FEF(n_c,target_on,PC_1),PC_FEF(n_c,target_on,PC_2),PC_FEF(n_c,target_on,PC_3),'^','MarkerSize',7,'MarkerFaceColor',[0.75 0.75 0.75],'MarkerEdgeColor',[0.75 0.75 0.75])
    hold on
    plot3(PC_FEF(n_c,response_on,PC_1),PC_FEF(n_c,response_on,PC_2),PC_FEF(n_c,response_on,PC_3),'^','MarkerSize',7,'MarkerFaceColor',[0.25 0.75 0.25],'MarkerEdgeColor',[0.25 0.75 0.25])
    
    clear x y z
    grid on
end
title(sprintf('FEF '))
xlabel(sprintf('PC %d',PC_1))
ylabel(sprintf('PC %d',PC_2))
zlabel(sprintf('PC %d',PC_3))
view([-1 1 1])



subplot(3,1,3)
hold on

for n_c=1:N_channels_stim
    x(:)=PC_PFC(n_c,:,PC_1);
    y(:)=PC_PFC(n_c,:,PC_2);
    z(:)=PC_PFC(n_c,:,PC_3);
    if n_c==1
        plot3(x,y,z,'-','MarkerFaceColor',color_bin(n_c,:),'MarkerEdgeColor',color_bin(n_c,:),'Color',color_bin(n_c,:),'LineWidth',2)
    elseif n_c==2
        plot3(x,y,z,'-','MarkerFaceColor',color_bin(n_c,:),'MarkerEdgeColor',color_bin(n_c,:),'Color',color_bin(n_c,:),'LineWidth',2)
    else
        plot3(x,y,z,'-','MarkerFaceColor',color_bin(n_c,:),'MarkerEdgeColor',color_bin(n_c,:),'Color',color_bin(n_c,:),'LineWidth',2)
    end
    plot3(PC_PFC(n_c,target_on,PC_1),PC_PFC(n_c,target_on,PC_2),PC_PFC(n_c,target_on,PC_3),'^','MarkerSize',7,'MarkerFaceColor',[0.75 0.75 0.75],'MarkerEdgeColor',[0.75 0.75 0.75])
    hold on
    plot3(PC_PFC(n_c,response_on,PC_1),PC_PFC(n_c,response_on,PC_2),PC_PFC(n_c,response_on,PC_3),'^','MarkerSize',7,'MarkerFaceColor',[0.25 0.75 0.25],'MarkerEdgeColor',[0.25 0.75 0.25])
    
    clear x y z
    
end
title(sprintf('PFC '))
xlabel(sprintf('PC %d',PC_1))
ylabel(sprintf('PC %d',PC_2))
zlabel(sprintf('PC %d',PC_3))
% view([-1 1 1])
view([-1 1 1])

grid on

%%

figure
hold on
plot(explained_LIP(1:10),'Color',color_for_ROI(1,:),'LineWidth',2)
plot(explained_FEF(1:10),'Color',color_for_ROI(2,:),'LineWidth',2)
plot(explained_PFC(1:10),'Color',color_for_ROI(3,:),'LineWidth',2)
xlabel('PC')
ylabel('% explained variance')

figure
hold on
bar(1,sum(eigenvalues_LIP)^2/sum(eigenvalues_LIP.^2),'FaceColor',color_for_ROI(1,:));
bar(2,sum(eigenvalues_FEF)^2/sum(eigenvalues_FEF.^2),'FaceColor',color_for_ROI(2,:));
bar(3,sum(eigenvalues_PFC)^2/sum(eigenvalues_PFC.^2),'FaceColor',color_for_ROI(3,:));
xlabel('ROI')
ylabel('Effective dimensionality')


%% bootstrap to get variance

for nb=1:100
    ind_LIP(:,nb)=randsample(size(LIP_PCA_belief_for_pca,1),size(LIP_PCA_belief_for_pca,1),'true');
    LIP_PCA_belief_for_pca_boot=LIP_PCA_belief_for_pca(ind_LIP,:);
    
    ind_FEF(:,nb)=randsample(size(FEF_PCA_belief_for_pca,1),size(FEF_PCA_belief_for_pca,1),'true');
    FEF_PCA_belief_for_pca_boot=FEF_PCA_belief_for_pca(ind_FEF,:);
    
    ind_PFC(:,nb)=randsample(size(PFC_PCA_belief_for_pca,1),size(PFC_PCA_belief_for_pca,1),'true');
    PFC_PCA_belief_for_pca_boot=PFC_PCA_belief_for_pca(ind_PFC,:);
    
    
    [~,~, eigenvalues_LIP_boot(:,nb),~,explained_LIP_boot(:,nb),~] = pca(LIP_PCA_belief_for_pca_boot');
    [~,~, eigenvalues_FEF_boot(:,nb),~,explained_FEF_boot(:,nb),~] = pca(FEF_PCA_belief_for_pca_boot');
    [~,~, eigenvalues_PFC_boot(:,nb),~,explained_PFC_boot(:,nb),~] = pca(PFC_PCA_belief_for_pca_boot');
    
end

%%
figure
hold on
shadedErrorBar([],mean(explained_LIP_boot(1:10,:),2)',[prctile(explained_LIP_boot(1:10,:),95,2),prctile(explained_LIP_boot(1:10,:),5,2)]',{'color',color_for_ROI(1,:),'LineWidth',2},2)
shadedErrorBar([],mean(explained_FEF_boot(1:10,:),2)',[prctile(explained_FEF_boot(1:10,:),95,2),prctile(explained_FEF_boot(1:10,:),5,2)]',{'color',color_for_ROI(2,:),'LineWidth',2},2)
shadedErrorBar([],mean(explained_PFC_boot(1:10,:),2)',[prctile(explained_PFC_boot(1:10,:),95,2),prctile(explained_PFC_boot(1:10,:),5,2)]',{'color',color_for_ROI(3,:),'LineWidth',2},2)
xlabel('PC')
ylabel('% explained variance')

figure
hold on
bar(1,mean(sum(eigenvalues_LIP_boot,1).^2./sum(eigenvalues_LIP_boot.^2,1)),'FaceColor',color_for_ROI(1,:));
bar(2,mean(sum(eigenvalues_FEF_boot,1).^2./sum(eigenvalues_FEF_boot.^2,1)),'FaceColor',color_for_ROI(2,:));
bar(3,mean(sum(eigenvalues_PFC_boot,1).^2./sum(eigenvalues_PFC_boot.^2,1)),'FaceColor',color_for_ROI(3,:));

errorbar(1,mean(sum(eigenvalues_LIP_boot,1).^2./sum(eigenvalues_LIP_boot.^2,1)),std(sum(eigenvalues_LIP_boot,1).^2./sum(eigenvalues_LIP_boot.^2,1)),'k');
errorbar(2,mean(sum(eigenvalues_FEF_boot,1).^2./sum(eigenvalues_FEF_boot.^2,1)),std(sum(eigenvalues_FEF_boot,1).^2./sum(eigenvalues_FEF_boot.^2,1)),'k');
errorbar(3,mean(sum(eigenvalues_PFC_boot,1).^2./sum(eigenvalues_PFC_boot.^2,1)),std(sum(eigenvalues_PFC_boot,1).^2./sum(eigenvalues_PFC_boot.^2,1)),'k');

xlabel('ROI')
ylabel('Effective dimensionality')


% %%
% figure
% 
% % PC_1=2;
% % PC_2=3;
% % PC_3=4;
% 
% PC_1=1;
% PC_2=2;
% 
% subplot(3,2,1)
% hold on
% for n_c=1:N_channels_stim
%     
%     x(:)=PC_LIP(n_c,:,PC_1);
%     y(:)=PC_LIP(n_c,:,PC_2);
%     plot(x,y,'-','MarkerFaceColor',color_bin(n_c,:),'MarkerEdgeColor',color_bin(n_c,:),'Color',color_bin(n_c,:),'LineWidth',2)
%     plot(PC_LIP(n_c,target_on,PC_1),PC_LIP(n_c,target_on,PC_2),'^','MarkerSize',7,'MarkerFaceColor',[0.75 0.75 0.75],'MarkerEdgeColor',[0.75 0.75 0.75])
%     hold on
%     plot(PC_LIP(n_c,response_on,PC_1),PC_LIP(n_c,response_on,PC_2),'^','MarkerSize',7,'MarkerFaceColor',[0.25 0.75 0.25],'MarkerEdgeColor',[0.25 0.75 0.25])
% %     grid on
%     clear x y z
%     
% end
% title(sprintf('LIP '))
% xlabel(sprintf('PC %d',PC_1))
% ylabel(sprintf('PC %d',PC_2))
% 
% PC_1=1;
% PC_2=3;
% 
% subplot(3,2,3)
% hold on
% for n_c=1:N_channels_stim
%     x(:)=PC_FEF(n_c,:,PC_1);
%     y(:)=PC_FEF(n_c,:,PC_2);
%         plot(x,y,'-','MarkerFaceColor',color_bin(n_c,:),'MarkerEdgeColor',color_bin(n_c,:),'Color',color_bin(n_c,:),'LineWidth',2)
%     plot(PC_FEF(n_c,target_on,PC_1),PC_FEF(n_c,target_on,PC_2),'^','MarkerSize',7,'MarkerFaceColor',[0.75 0.75 0.75],'MarkerEdgeColor',[0.75 0.75 0.75])
%     hold on
%     plot(PC_FEF(n_c,response_on,PC_1),PC_FEF(n_c,response_on,PC_2),'^','MarkerSize',7,'MarkerFaceColor',[0.25 0.75 0.25],'MarkerEdgeColor',[0.25 0.75 0.25])
%     
%     clear x y z
% %     grid on
% end
% title(sprintf('FEF '))
% xlabel(sprintf('PC %d',PC_1))
% ylabel(sprintf('PC %d',PC_2))
% 
% PC_1=1;
% PC_2=2;
% 
% subplot(3,2,5)
% hold on
% 
% for n_c=1:N_channels_stim
%     x(:)=PC_PFC(n_c,:,PC_1);
%     y(:)=PC_PFC(n_c,:,PC_2);
%     plot(x,y,'-','MarkerFaceColor',color_bin(n_c,:),'MarkerEdgeColor',color_bin(n_c,:),'Color',color_bin(n_c,:),'LineWidth',2)
%     plot(PC_PFC(n_c,target_on,PC_1),PC_PFC(n_c,target_on,PC_2),'^','MarkerSize',7,'MarkerFaceColor',[0.75 0.75 0.75],'MarkerEdgeColor',[0.75 0.75 0.75])
%     hold on
%     plot(PC_PFC(n_c,response_on,PC_1),PC_PFC(n_c,response_on,PC_2),'^','MarkerSize',7,'MarkerFaceColor',[0.25 0.75 0.25],'MarkerEdgeColor',[0.25 0.75 0.25])
%     
%     clear x y z
%     
% end
% title(sprintf('PFC '))
% xlabel(sprintf('PC %d',PC_1))
% ylabel(sprintf('PC %d',PC_2))
% 
% % grid on
% 
% 
% PC_1=3;
% PC_2=4;
% 
% subplot(3,2,2)
% hold on
% for n_c=1:N_channels_stim
%     
%     x(:)=PC_LIP(n_c,:,PC_1);
%     y(:)=PC_LIP(n_c,:,PC_2);
%     plot(x,y,'-','MarkerFaceColor',color_bin(n_c,:),'MarkerEdgeColor',color_bin(n_c,:),'Color',color_bin(n_c,:),'LineWidth',2)
%     plot(PC_LIP(n_c,target_on,PC_1),PC_LIP(n_c,target_on,PC_2),'^','MarkerSize',7,'MarkerFaceColor',[0.75 0.75 0.75],'MarkerEdgeColor',[0.75 0.75 0.75])
%     hold on
%     plot(PC_LIP(n_c,response_on,PC_1),PC_LIP(n_c,response_on,PC_2),'^','MarkerSize',7,'MarkerFaceColor',[0.25 0.75 0.25],'MarkerEdgeColor',[0.25 0.75 0.25])
% %     grid on
%     clear x y z
%     
% end
% title(sprintf('LIP '))
% xlabel(sprintf('PC %d',PC_1))
% ylabel(sprintf('PC %d',PC_2))
% 
% PC_1=2;
% PC_2=4;
% 
% 
% subplot(3,2,4)
% hold on
% for n_c=1:N_channels_stim
%     x(:)=PC_FEF(n_c,:,PC_1);
%     y(:)=PC_FEF(n_c,:,PC_2);
%         plot(x,y,'-','MarkerFaceColor',color_bin(n_c,:),'MarkerEdgeColor',color_bin(n_c,:),'Color',color_bin(n_c,:),'LineWidth',2)
%     plot(PC_FEF(n_c,target_on,PC_1),PC_FEF(n_c,target_on,PC_2),'^','MarkerSize',7,'MarkerFaceColor',[0.75 0.75 0.75],'MarkerEdgeColor',[0.75 0.75 0.75])
%     hold on
%     plot(PC_FEF(n_c,response_on,PC_1),PC_FEF(n_c,response_on,PC_2),'^','MarkerSize',7,'MarkerFaceColor',[0.25 0.75 0.25],'MarkerEdgeColor',[0.25 0.75 0.25])
%     
%     clear x y z
% %     grid on
% end
% title(sprintf('FEF '))
% xlabel(sprintf('PC %d',PC_1))
% ylabel(sprintf('PC %d',PC_2))
% 
% PC_1=3;
% PC_2=4;
% 
% 
% subplot(3,2,6)
% hold on
% 
% for n_c=1:N_channels_stim
%     x(:)=PC_PFC(n_c,:,PC_1);
%     y(:)=PC_PFC(n_c,:,PC_2);
%     plot(x,y,'-','MarkerFaceColor',color_bin(n_c,:),'MarkerEdgeColor',color_bin(n_c,:),'Color',color_bin(n_c,:),'LineWidth',2)
%     plot(PC_PFC(n_c,target_on,PC_1),PC_PFC(n_c,target_on,PC_2),'^','MarkerSize',7,'MarkerFaceColor',[0.75 0.75 0.75],'MarkerEdgeColor',[0.75 0.75 0.75])
%     hold on
%     plot(PC_PFC(n_c,response_on,PC_1),PC_PFC(n_c,response_on,PC_2),'^','MarkerSize',7,'MarkerFaceColor',[0.25 0.75 0.25],'MarkerEdgeColor',[0.25 0.75 0.25])
%     
%     clear x y z
%     
% end
% title(sprintf('PFC '))
% xlabel(sprintf('PC %d',PC_1))
% ylabel(sprintf('PC %d',PC_2))
% 
% % grid on
% 
% 
% 
% 
% %% remove time
% 
% for t=1:length(window_start_list)
%     %LIP
%     for nn=1:size(LIP_PCA_belief,1)
%         LIP_PCA_no_time(nn,:,t)=LIP_PCA_belief(nn,:,t)-mean(LIP_PCA_belief(nn,:,t),'all');
%     end
%     %FEF
%     for nn=1:size(FEF_PCA_belief,1)
%         FEF_PCA_no_time(nn,:,t)=FEF_PCA_belief(nn,:,t)-mean(FEF_PCA_belief(nn,:,t),'all');
%     end
%     %PFC
%     for nn=1:size(PFC_PCA_belief,1)
%         PFC_PCA_no_time(nn,:,t)=PFC_PCA_belief(nn,:,t)-mean(PFC_PCA_belief(nn,:,t),'all');
%     end
% end
% 
% %% reshape for the PCA
% 
% 
% LIP_PCA_no_time_for_pca(:,:)=reshape(LIP_PCA_no_time(:,:,:),size(LIP_PCA_belief,1),size(LIP_PCA_belief,2)*size(LIP_PCA_belief,3));
% FEF_PCA_no_time_for_pca(:,:)=reshape(FEF_PCA_no_time(:,:,:),size(FEF_PCA_belief,1),size(FEF_PCA_belief,2)*size(FEF_PCA_belief,3));
% PFC_PCA_no_time_for_pca(:,:)=reshape(PFC_PCA_no_time(:,:,:),size(PFC_PCA_belief,1),size(PFC_PCA_belief,2)*size(PFC_PCA_belief,3));
% 
% 
% %% do the pca
% 
% [coeff_LIP_no_time(:,:), score_LIP_no_time(:,:),~,~,explained_LIP_no_time(:),~] = pca(LIP_PCA_no_time_for_pca(:,:)');
% [coeff_FEF_no_time(:,:), score_FEF_no_time(:,:),~,~,explained_FEF_no_time(:),~] = pca(FEF_PCA_no_time_for_pca(:,:)');
% [coeff_PFC_no_time(:,:), score_PFC_no_time(:,:),~,~,explained_PFC_no_time(:),~] = pca(PFC_PCA_no_time_for_pca(:,:)');
% 
% %% reshape or platting
% 
% for k=1:4
%     for pc=1:8
%         PC_LIP_no_time(:,:,pc)=reshape(score_LIP_no_time(:,pc),size(LIP_PCA_belief,2),size(LIP_PCA_belief,3));
%         PC_FEF_no_time(:,:,pc)=reshape(score_FEF_no_time(:,pc),size(FEF_PCA_belief,2),size(FEF_PCA_belief,3));
%         PC_PFC_no_time(:,:,pc)=reshape(score_PFC_no_time(:,pc),size(PFC_PCA_belief,2),size(PFC_PCA_belief,3));
%     end
% end
% 
% 
% 
% %% now plot
% 
% N_bins_color=N_bins;
% 
% load('colors.mat')
% color_bin(1,:)=colors(1,:);
% color_bin(2,:)=colors(14,:);
% color_bin(3,:)=colors(27,:);
% 
% % PC_1=1;
% % PC_2=2;
% % PC_3=5;
% 
% N_t=length(window_start_list);
% 
% j=1;
% figure
% 
% PC_1=1;
% PC_2=2;
% PC_3=3;
% 
% subplot(3,1,1)
% hold on
% for n_c=1:N_channels_stim
%     
%     x(:)=PC_LIP_no_time(n_c,:,PC_1);
%     y(:)=PC_LIP_no_time(n_c,:,PC_2);
%     z(:)=PC_LIP_no_time(n_c,:,PC_3);
%     if n_c==1
%         plot3(x,y,z,'-o','MarkerFaceColor',color_bin(n_c,:),'MarkerEdgeColor',color_bin(n_c,:),'Color',color_bin(n_c,:),'LineWidth',(n_c==n_c)+2)
%     elseif n_c==2
%         plot3(x,y,z,'-o','MarkerFaceColor',color_bin(n_c,:),'MarkerEdgeColor',color_bin(n_c,:),'Color',color_bin(n_c,:),'LineWidth',(n_c==n_c)+2)
%     else
%         plot3(x,y,z,'-o','MarkerFaceColor',color_bin(n_c,:),'MarkerEdgeColor',color_bin(n_c,:),'Color',color_bin(n_c,:),'LineWidth',(n_c==n_c)+2)
%     end
%     plot3(PC_LIP_no_time(n_c,target_on,PC_1),PC_LIP_no_time(n_c,target_on,PC_2),PC_LIP_no_time(n_c,target_on,PC_3),'^','MarkerSize',15,'MarkerFaceColor',[0.75 0.75 0.75],'MarkerEdgeColor',[0.75 0.75 0.75])
%     hold on
%     plot3(PC_LIP_no_time(n_c,response_on,PC_1),PC_LIP_no_time(n_c,response_on,PC_2),PC_LIP_no_time(n_c,response_on,PC_3),'^','MarkerSize',15,'MarkerFaceColor',[0.25 0.75 0.25],'MarkerEdgeColor',[0.25 0.75 0.25])
%     grid on
%     clear x y z
%     
% end
% title(sprintf('LIP '))
% xlabel(sprintf('PC %d',PC_1))
% ylabel(sprintf('PC %d',PC_2))
% zlabel(sprintf('PC %d',PC_3))
% view([-1 -1 1])
% 
% subplot(3,1,2)
% hold on
% for n_c=1:N_channels_stim
%     x(:)=PC_FEF_no_time(n_c,:,PC_1);
%     y(:)=PC_FEF_no_time(n_c,:,PC_2);
%     z(:)=PC_FEF_no_time(n_c,:,PC_3);
%     if n_c==1
%         plot3(x,y,z,'-o','MarkerFaceColor',color_bin(n_c,:),'MarkerEdgeColor',color_bin(n_c,:),'Color',color_bin(n_c,:),'LineWidth',(n_c==n_c)+2)
%     elseif n_c==2
%         plot3(x,y,z,'-o','MarkerFaceColor',color_bin(n_c,:),'MarkerEdgeColor',color_bin(n_c,:),'Color',color_bin(n_c,:),'LineWidth',(n_c==n_c)+2)
%     else
%         plot3(x,y,z,'-o','MarkerFaceColor',color_bin(n_c,:),'MarkerEdgeColor',color_bin(n_c,:),'Color',color_bin(n_c,:),'LineWidth',(n_c==n_c)+2)
%     end
%     plot3(PC_FEF_no_time(n_c,target_on,PC_1),PC_FEF_no_time(n_c,target_on,PC_2),PC_FEF_no_time(n_c,target_on,PC_3),'^','MarkerSize',15,'MarkerFaceColor',[0.75 0.75 0.75],'MarkerEdgeColor',[0.75 0.75 0.75])
%     hold on
%     plot3(PC_FEF_no_time(n_c,response_on,PC_1),PC_FEF_no_time(n_c,response_on,PC_2),PC_FEF_no_time(n_c,response_on,PC_3),'^','MarkerSize',15,'MarkerFaceColor',[0.25 0.75 0.25],'MarkerEdgeColor',[0.25 0.75 0.25])
%     
%     clear x y z
%     grid on
% end
% title(sprintf('FEF '))
% xlabel(sprintf('PC %d',PC_1))
% ylabel(sprintf('PC %d',PC_2))
% zlabel(sprintf('PC %d',PC_3))
% view([1 1 1])
% 
% subplot(3,1,3)
% hold on
% 
% for n_c=1:N_channels_stim
%     x(:)=PC_PFC_no_time(n_c,:,PC_1);
%     y(:)=PC_PFC_no_time(n_c,:,PC_2);
%     z(:)=PC_PFC_no_time(n_c,:,PC_3);
%     if n_c==1
%         plot3(x,y,z,'-o','MarkerFaceColor',color_bin(n_c,:),'MarkerEdgeColor',color_bin(n_c,:),'Color',color_bin(n_c,:),'LineWidth',(n_c==n_c)+2)
%     elseif n_c==2
%         plot3(x,y,z,'-o','MarkerFaceColor',color_bin(n_c,:),'MarkerEdgeColor',color_bin(n_c,:),'Color',color_bin(n_c,:),'LineWidth',(n_c==n_c)+2)
%     else
%         plot3(x,y,z,'-o','MarkerFaceColor',color_bin(n_c,:),'MarkerEdgeColor',color_bin(n_c,:),'Color',color_bin(n_c,:),'LineWidth',(n_c==n_c)+2)
%     end
%     plot3(PC_PFC_no_time(n_c,target_on,PC_1),PC_PFC_no_time(n_c,target_on,PC_2),PC_PFC_no_time(n_c,target_on,PC_3),'^','MarkerSize',15,'MarkerFaceColor',[0.75 0.75 0.75],'MarkerEdgeColor',[0.75 0.75 0.75])
%     hold on
%     plot3(PC_PFC_no_time(n_c,response_on,PC_1),PC_PFC_no_time(n_c,response_on,PC_2),PC_PFC_no_time(n_c,response_on,PC_3),'^','MarkerSize',15,'MarkerFaceColor',[0.25 0.75 0.25],'MarkerEdgeColor',[0.25 0.75 0.25])
%     
%     clear x y z
%     
% end
% title(sprintf('PFC '))
% xlabel(sprintf('PC %d',PC_1))
% ylabel(sprintf('PC %d',PC_2))
% zlabel(sprintf('PC %d',PC_3))
% view([-1 1 1])
% grid on
% 
% 
% 
% 
% 
% 
% 
% 
