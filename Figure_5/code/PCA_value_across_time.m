arrayID_list=[3 8 12 14 18 20 22 26 28 30 32 34 36 38 40 42 44];

fsroot='/Volumes/buschman';

event='target';
window_size=200;
for i=1:25
    initial_window=-400;
    event_list{i}='target';
    window_start_list(i)=initial_window+(i-1)*50;
end

task='Learning_Attentional_Templates';
subtask='exploreexploit/Restricted_FEF/Reset_RW_model';

subsubtask_classifier=sprintf('Classifier_%s_%d_to_%s_%d',event,window_start_list(1),event,window_start_list(end));

dirstem = fullfile('Projects',task,'Analysis','Electrophy_analysis',subtask,subsubtask_classifier);
save_path_clasifier = fullfile(fsroot,dirstem);
mkdir(save_path_clasifier);

save_name='PCA_value';

N_value_bins=2;
N_channels_stim=3;

N_bins=100;
color_binned =  0:(2*pi/N_bins):(2*pi-(2*pi/N_bins));


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
            Conditions.Chosen_value_bin(i)=Conditions.trial_in_block(i)/sum(Conditions.block_nb==Conditions.block_nb(i));
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
        
        Conditions.value=zeros(size(Conditions.stim));
        Conditions.chosen_value=zeros(size(Conditions.chosen_color));
        for i=1:size(Conditions.stim,2)
            for j=1:size(Conditions.stim,1)
                if isnan(Conditions.stim(j,i))
                    Conditions.value(j,i)=NaN;
                else
                    if Conditions.stim(j,i)>=color_binned(N_bins)
                        this_bin=N_bins;
                    else
                        for k=1:N_bins-1
                            if Conditions.stim(j,i)>=color_binned(k) && Conditions.stim(j,i)<color_binned(k+1)
                                this_bin=k;
                            end
                        end
                    end
                    Conditions.value(j,i)=Conditions.belief(this_bin,i);
                end
            end
            if Conditions.chosen_color(i)>=color_binned(N_bins)
                this_bin=N_bins;
            else
                for k=1:N_bins-1
                    if Conditions.chosen_color(i)>=color_binned(k) && Conditions.chosen_color(i)<color_binned(k+1)
                        this_bin=k;
                    end
                end
            end
            Conditions.chosen_value(i)=Conditions.belief(this_bin,i);
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
            for i=1:length(Conditions.block_nb)
                for t=1:length(window_start_list)
                    Pseudo_pop_LIP(count_LIP,i).Classifier_FR(t)=data_LIP(nn,i,t);
                end
                Pseudo_pop_LIP(count_LIP,i).Sess=n_sess;
                Pseudo_pop_LIP(count_LIP,i).Block_nb=Conditions.block_nb(i); %we need to take the updated belief here
                Pseudo_pop_LIP(count_LIP,i).choice=Conditions.choice(i);
                Pseudo_pop_LIP(count_LIP,i).stim_on(:)=~isnan(Conditions.stim(:,i));
                Pseudo_pop_LIP(count_LIP,i).chosen_value=Conditions.chosen_value(i);
                Pseudo_pop_LIP(count_LIP,i).Peak_belief=mod(Conditions.Peak_belief(i),2*pi);
                prc=prctile(Conditions.chosen_value,[0 50 100]);
                for np=1:N_value_bins
                    if Pseudo_pop_LIP(count_LIP,i).chosen_value>=prc(np) && Pseudo_pop_LIP(count_LIP,i).chosen_value<prc(np+1)
                        Pseudo_pop_LIP(count_LIP,i).Chosen_value_bin=np;
                    end
                end
                for n_c=1:N_channels_stim
                    if mod(Pseudo_pop_LIP(count_LIP,i).Peak_belief,2*pi)>=(n_c-1)*2*pi/N_channels_stim && mod(Pseudo_pop_LIP(count_LIP,i).Peak_belief,2*pi)<n_c*2*pi/N_channels_stim
                        Pseudo_pop_LIP(count_LIP,i).Belief_color=n_c;
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
        
        %get the progression in trial
        for i=1:length(Conditions.block_nb)
            Conditions.Chosen_value_bin(i)=Conditions.trial_in_block(i)/sum(Conditions.block_nb==Conditions.block_nb(i));
        end
        Conditions.value=zeros(size(Conditions.stim));
        Conditions.chosen_value=zeros(size(Conditions.chosen_color));
        for i=1:size(Conditions.stim,2)
            for j=1:size(Conditions.stim,1)
                if isnan(Conditions.stim(j,i))
                    Conditions.value(j,i)=NaN;
                else
                    if Conditions.stim(j,i)>=color_binned(N_bins)
                        this_bin=N_bins;
                    else
                        for k=1:N_bins-1
                            if Conditions.stim(j,i)>=color_binned(k) && Conditions.stim(j,i)<color_binned(k+1)
                                this_bin=k;
                            end
                        end
                    end
                    Conditions.value(j,i)=Conditions.belief(this_bin,i);
                end
            end
            if Conditions.chosen_color(i)>=color_binned(N_bins)
                this_bin=N_bins;
            else
                for k=1:N_bins-1
                    if Conditions.chosen_color(i)>=color_binned(k) && Conditions.chosen_color(i)<color_binned(k+1)
                        this_bin=k;
                    end
                end
            end
            Conditions.chosen_value(i)=Conditions.belief(this_bin,i);
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
            for i=1:length(Conditions.block_nb)
                for t=1:length(window_start_list)
                    Pseudo_pop_FEF(count_FEF,i).Classifier_FR(t)=data_FEF(nn,i,t);
                end
                Pseudo_pop_FEF(count_FEF,i).Sess=n_sess;
                Pseudo_pop_FEF(count_FEF,i).Block_nb=Conditions.block_nb(i); %we need to take the updated belief here
                Pseudo_pop_FEF(count_FEF,i).choice=Conditions.choice(i);
                Pseudo_pop_FEF(count_FEF,i).stim_on(:)=~isnan(Conditions.stim(:,i));
                Pseudo_pop_FEF(count_FEF,i).chosen_value=Conditions.chosen_value(i);
                Pseudo_pop_FEF(count_FEF,i).Peak_belief=mod(Conditions.Peak_belief(i),2*pi);
                prc=prctile(Conditions.chosen_value,[0 50 100]);
                for np=1:N_value_bins
                    if Pseudo_pop_FEF(count_FEF,i).chosen_value>=prc(np) && Pseudo_pop_FEF(count_FEF,i).chosen_value<prc(np+1)
                        Pseudo_pop_FEF(count_FEF,i).Chosen_value_bin=np;
                    end
                end
                for n_c=1:N_channels_stim
                    if mod(Pseudo_pop_FEF(count_FEF,i).Peak_belief,2*pi)>=(n_c-1)*2*pi/N_channels_stim && mod(Pseudo_pop_FEF(count_FEF,i).Peak_belief,2*pi)<n_c*2*pi/N_channels_stim
                        Pseudo_pop_FEF(count_FEF,i).Belief_color=n_c;
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
    
    %get the progression in trial
    for i=1:length(Conditions.block_nb)
        Conditions.Chosen_value_bin(i)=Conditions.trial_in_block(i)/sum(Conditions.block_nb==Conditions.block_nb(i));
    end
    Conditions.value=zeros(size(Conditions.stim));
    Conditions.chosen_value=zeros(size(Conditions.chosen_color));
    for i=1:size(Conditions.stim,2)
        for j=1:size(Conditions.stim,1)
            if isnan(Conditions.stim(j,i))
                Conditions.value(j,i)=NaN;
            else
                if Conditions.stim(j,i)>=color_binned(N_bins)
                    this_bin=N_bins;
                else
                    for k=1:N_bins-1
                        if Conditions.stim(j,i)>=color_binned(k) && Conditions.stim(j,i)<color_binned(k+1)
                            this_bin=k;
                        end
                    end
                end
                Conditions.value(j,i)=Conditions.belief(this_bin,i);
            end
        end
        if Conditions.chosen_color(i)>=color_binned(N_bins)
            this_bin=N_bins;
        else
            for k=1:N_bins-1
                if Conditions.chosen_color(i)>=color_binned(k) && Conditions.chosen_color(i)<color_binned(k+1)
                    this_bin=k;
                end
            end
        end
        Conditions.chosen_value(i)=Conditions.belief(this_bin,i);
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
        for i=1:length(Conditions.block_nb)
            for t=1:length(window_start_list)
                Pseudo_pop_PFC(count_PFC,i).Classifier_FR(t)=data_PFC(nn,i,t);
            end
            Pseudo_pop_PFC(count_PFC,i).Sess=n_sess;
            Pseudo_pop_PFC(count_PFC,i).Block_nb=Conditions.block_nb(i); %we need to take the updated belief here
            Pseudo_pop_PFC(count_PFC,i).choice=Conditions.choice(i);
            Pseudo_pop_PFC(count_PFC,i).stim_on(:)=~isnan(Conditions.stim(:,i));
            Pseudo_pop_PFC(count_PFC,i).chosen_value=Conditions.chosen_value(i);
            Pseudo_pop_PFC(count_PFC,i).Peak_belief=mod(Conditions.Peak_belief(i),2*pi);
                prc=prctile(Conditions.chosen_value,[0 50 100]);
            for np=1:N_value_bins
                if Pseudo_pop_PFC(count_PFC,i).chosen_value>=prc(np) && Pseudo_pop_PFC(count_PFC,i).chosen_value<prc(np+1)
                    Pseudo_pop_PFC(count_PFC,i).Chosen_value_bin=np;
                end
            end
            for n_c=1:N_channels_stim
                if mod(Pseudo_pop_PFC(count_PFC,i).Peak_belief,2*pi)>=(n_c-1)*2*pi/N_channels_stim && mod(Pseudo_pop_PFC(count_PFC,i).Peak_belief,2*pi)<n_c*2*pi/N_channels_stim
                    Pseudo_pop_PFC(count_PFC,i).Belief_color=n_c;
                end
            end
            
        end
        
        count_PFC=count_PFC+1;
    end
    
    clear Conditions
    
end

%% Now create the matrix N x 9 conditions x time x 3 prog x 4 loc


for nn=1:size(Pseudo_pop_LIP,1)
    for i=1:size(Pseudo_pop_LIP(nn,:),2)
        if ~isempty(Pseudo_pop_LIP(nn,i).Chosen_value_bin) && ~isnan(Pseudo_pop_LIP(nn,i).Chosen_value_bin)
            selected_LIP_value(nn,i,Pseudo_pop_LIP(nn,i).Belief_color,Pseudo_pop_LIP(nn,i).Chosen_value_bin)=1;
        end
    end
end
for nn=1:size(Pseudo_pop_FEF,1)
    for i=1:size(Pseudo_pop_FEF(nn,:),2)
        if ~isempty(Pseudo_pop_FEF(nn,i).Chosen_value_bin) && ~isnan(Pseudo_pop_FEF(nn,i).Chosen_value_bin)
            selected_FEF_value(nn,i,Pseudo_pop_FEF(nn,i).Belief_color,Pseudo_pop_FEF(nn,i).Chosen_value_bin)=1;
        end
    end
end
for nn=1:size(Pseudo_pop_PFC,1)
    for i=1:size(Pseudo_pop_PFC(nn,:),2)
        if ~isempty(Pseudo_pop_PFC(nn,i).Chosen_value_bin) && ~isnan(Pseudo_pop_PFC(nn,i).Chosen_value_bin)
            selected_PFC_value(nn,i,Pseudo_pop_PFC(nn,i).Belief_color,Pseudo_pop_PFC(nn,i).Chosen_value_bin)=1;
        end
    end
end

%% Now create the matrix N x (4x5) conditions x time

for k=1:N_channels_stim %for each templte
    for j=1:N_value_bins
        for t=1:length(window_start_list)
            %LIP
            for nn=1:size(Pseudo_pop_LIP,1)
                if sum(selected_LIP_value(nn,:,k,j))>1
                    list=find(selected_LIP_value(nn,:,k,j)==1);
                    for i=1:length(list)
                        buffer(i)=Pseudo_pop_LIP(nn,list(i)).Classifier_FR(t);
                    end
                    LIP_PCA_chosen_value(nn,t,k,j)=mean(buffer);
                    clear buffer list
                else
                    LIP_PCA_chosen_value(nn,t,k,j)=NaN;
                end
            end
            %FEF
            for nn=1:size(Pseudo_pop_FEF,1)
                if sum(selected_FEF_value(nn,:,k,j))>1
                    list=find(selected_FEF_value(nn,:,k,j)==1);
                    for i=1:length(list)
                        buffer(i)=Pseudo_pop_FEF(nn,list(i)).Classifier_FR(t);
                    end
                    FEF_PCA_chosen_value(nn,t,k,j)=mean(buffer);
                    clear buffer list
                else
                    FEF_PCA_chosen_value(nn,t,k,j)=NaN;
                end
            end
            %PFC
            for nn=1:size(Pseudo_pop_PFC,1)
                if sum(selected_PFC_value(nn,:,k,j))>1
                    list=find(selected_PFC_value(nn,:,k,j)==1);
                    for i=1:length(list)
                        buffer(i)=Pseudo_pop_PFC(nn,list(i)).Classifier_FR(t);
                    end
                    PFC_PCA_chosen_value(nn,t,k,j)=mean(buffer);
                    clear buffer list
                else
                    PFC_PCA_chosen_value(nn,t,k,j)=NaN;
                end
            end
        end
    end
end

%% remove NaN

LIP_PCA_chosen_value=LIP_PCA_chosen_value(~isnan(sum(sum(sum(LIP_PCA_chosen_value,2),3),4)),:,:,:,:,:);
FEF_PCA_chosen_value=FEF_PCA_chosen_value(~isnan(sum(sum(sum(FEF_PCA_chosen_value,2),3),4)),:,:,:,:,:);
PFC_PCA_chosen_value=PFC_PCA_chosen_value(~isnan(sum(sum(sum(PFC_PCA_chosen_value,2),3),4)),:,:,:,:,:);

% remove time

for t=1:length(window_start_list)
    %LIP
    for nn=1:size(LIP_PCA_chosen_value,1)
        LIP_PCA_chosen_value(nn,t,:,:)=LIP_PCA_chosen_value(nn,t,:,:)-mean(LIP_PCA_chosen_value(nn,t,:,:),'all');
    end
    %FEF
    for nn=1:size(FEF_PCA_chosen_value,1)
        FEF_PCA_chosen_value(nn,t,:,:)=FEF_PCA_chosen_value(nn,t,:,:)-mean(FEF_PCA_chosen_value(nn,t,:,:),'all');
    end
    %PFC
    for nn=1:size(PFC_PCA_chosen_value,1)
        PFC_PCA_chosen_value(nn,t,:,:)=PFC_PCA_chosen_value(nn,t,:,:)-mean(PFC_PCA_chosen_value(nn,t,:,:),'all');
    end
end

%% reshape for the PCA


LIP_PCA_chosen_value_for_pca=reshape(LIP_PCA_chosen_value,size(LIP_PCA_chosen_value,1),size(LIP_PCA_chosen_value,2)*size(LIP_PCA_chosen_value,3)*size(LIP_PCA_chosen_value,4));
FEF_PCA_chosen_value_for_pca=reshape(FEF_PCA_chosen_value,size(FEF_PCA_chosen_value,1),size(FEF_PCA_chosen_value,2)*size(FEF_PCA_chosen_value,3)*size(FEF_PCA_chosen_value,4));
PFC_PCA_chosen_value_for_pca=reshape(PFC_PCA_chosen_value,size(PFC_PCA_chosen_value,1),size(PFC_PCA_chosen_value,2)*size(PFC_PCA_chosen_value,3)*size(PFC_PCA_chosen_value,4));



%% do the pca

[coeff_LIP, score_LIP,~,~,explained_LIP,~] = pca(LIP_PCA_chosen_value_for_pca');
[coeff_FEF, score_FEF,~,~,explained_FEF,~] = pca(FEF_PCA_chosen_value_for_pca');
[coeff_PFC, score_PFC,~,~,explained_PFC,~] = pca(PFC_PCA_chosen_value_for_pca');

%% reshape fpr plotting

for pc=1:10
    PC_LIP(:,:,:,pc)=reshape(score_LIP(:,pc),size(LIP_PCA_chosen_value,2),size(LIP_PCA_chosen_value,3),size(LIP_PCA_chosen_value,4));
    PC_FEF(:,:,:,pc)=reshape(score_FEF(:,pc),size(FEF_PCA_chosen_value,2),size(FEF_PCA_chosen_value,3),size(FEF_PCA_chosen_value,4));
    PC_PFC(:,:,:,pc)=reshape(score_PFC(:,pc),size(PFC_PCA_chosen_value,2),size(PFC_PCA_chosen_value,3),size(PFC_PCA_chosen_value,4));
end

%% now plot

N_value_bins_color=N_value_bins;

load('colors.mat')
color_for_plot(1,:)=colors(1,:);
color_for_plot(2,:)=colors(round(40/3),:);
color_for_plot(3,:)=colors(round(2*40/3),:);

for i=1:3
color_for_plot2(i,:)=rgb2hsv(color_for_plot(i,:));
end
color_for_plot2(:,2)=color_for_plot2(:,2)-0.3;
for i=1:3
color_for_plot2(i,:)=hsv2rgb(color_for_plot2(i,:));
end

%%
% 
% for PC_1=1:9
%     for PC_2=PC_1+1:10
%             figure
%             for k=1:N_channels_stim
%                 for j=1:N_value_bins
%                     x(:)=PC_PFC(:,k,j,PC_1);
%                     y(:)=PC_PFC(:,k,j,PC_2);
%                     if j==1
%                     plot(x,y,'-o','MarkerFaceColor',color_for_plot(k,:),'MarkerEdgeColor',color_for_plot(k,:),'Color',color_for_plot(k,:),'LineWidth',2)
%                     elseif j==2
%                     plot(x,y,'--o','MarkerFaceColor',color_for_plot(k,:),'MarkerEdgeColor',color_for_plot(k,:),'Color',color_for_plot(k,:),'LineWidth',2)
%                     end
%                         
%                     hold on
%                     clear x y z
%                 end
%             end
%             title(sprintf('PFC'))
%             xlabel(sprintf('PC %d',PC_1))
%             ylabel(sprintf('PC %d',PC_2))
%         end
% end
% 



%%
            figure
subplot(3,1,1)
for PC_1=1
    for PC_2=2
        for PC_3=3
            for k=1:N_channels_stim
                for j=1:N_value_bins
                    x(:)=PC_LIP(:,k,j,PC_1);
                    y(:)=PC_LIP(:,k,j,PC_2);
                    z(:)=PC_LIP(:,k,j,PC_3);
                    if j==1
                    plot3(x,y,z,'-o','MarkerFaceColor',color_for_plot(k,:),'MarkerEdgeColor',color_for_plot(k,:),'Color',color_for_plot(k,:),'LineWidth',2)
                    elseif j==2
                    plot3(x,y,z,'--o','MarkerFaceColor',color_for_plot2(k,:),'MarkerEdgeColor',color_for_plot2(k,:),'Color',color_for_plot2(k,:),'LineWidth',2)
                    end
                        
                    hold on
                    clear x y z
                end
            end
            title(sprintf('LIP'))
            xlabel(sprintf('PC %d',PC_1))
            ylabel(sprintf('PC %d',PC_2))
            zlabel(sprintf('PC %d',PC_3))
        end
    end
end
grid on

subplot(3,1,2)

for PC_1=1
    for PC_2=2
        for PC_3=3

            for k=1:N_channels_stim
                for j=1:N_value_bins
                    x(:)=PC_FEF(:,k,j,PC_1);
                    y(:)=PC_FEF(:,k,j,PC_2);
                    z(:)=PC_FEF(:,k,j,PC_3);
                    if j==1
                    plot3(x,y,z,'-o','MarkerFaceColor',color_for_plot(k,:),'MarkerEdgeColor',color_for_plot(k,:),'Color',color_for_plot(k,:),'LineWidth',2)
                    elseif j==2
                    plot3(x,y,z,'--o','MarkerFaceColor',color_for_plot2(k,:),'MarkerEdgeColor',color_for_plot2(k,:),'Color',color_for_plot2(k,:),'LineWidth',2)
                    end
                    hold on
                    clear x y z
                end
            end
            title(sprintf('FEF'))
            xlabel(sprintf('PC %d',PC_1))
            ylabel(sprintf('PC %d',PC_2))
            zlabel(sprintf('PC %d',PC_3))
        end
    end
end
grid on

subplot(3,1,3)
for PC_1=1
    for PC_2=2
        for PC_3=3
            for k=1:N_channels_stim
                for j=1:N_value_bins
                    x(:)=PC_PFC(:,k,j,PC_1);
                    y(:)=PC_PFC(:,k,j,PC_2);
                    z(:)=PC_PFC(:,k,j,PC_3);
                    if j==1
                    plot3(x,y,z,'-o','MarkerFaceColor',color_for_plot(k,:),'MarkerEdgeColor',color_for_plot(k,:),'Color',color_for_plot(k,:),'LineWidth',2)
                    elseif j==2
                    plot3(x,y,z,'--o','MarkerFaceColor',color_for_plot2(k,:),'MarkerEdgeColor',color_for_plot2(k,:),'Color',color_for_plot2(k,:),'LineWidth',2)
                    end
                    hold on
                    clear x y z
                end
            end
            title(sprintf('PFC'))
            xlabel(sprintf('PC %d',PC_1))
            ylabel(sprintf('PC %d',PC_2))
            zlabel(sprintf('PC %d',PC_3))
        end
    end
end
grid on




