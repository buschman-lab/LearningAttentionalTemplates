%prepare_data_pseudopop

%here we want to put together separate events
%we want to combine the end of the reward with the next trial
%We need to be caredul about the belief, we
%want the updated belief after the the reward
%We can't take the first trial as it doesn't have a previous reward end

arrayID_list=[3 8 12 14 18 20 22 26 28 30 32 34 36 38 40 42 44];

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

window_size=300;

subsubtask_classifier=sprintf('Classifier_%s_%d_to_%s_%d',event_list{1},window_start_list(1),event_list{end},window_start_list(end));
% subsubtask_classifier=sprintf('Classifier_%s_%d_to_%s_%d',event,window_start_list(1),event,window_start_list(end));

dirstem = fullfile('Projects',task,'Analysis','Electrophy_analysis',subtask,subsubtask_classifier);
save_path_clasifier = fullfile(fsroot,dirstem);
mkdir(save_path_clasifier);

save_name='Pseudo_pop_peak_belief_prog_classifier';

N_channels_stim=3;

count_LIP=1;
count_FEF=1;
count_PFC=1;

N_sessions=length(arrayID_list);

for n_sess=1:N_sessions
    
    arrayID=arrayID_list(n_sess);
    
    
    %LIP
    if arrayID~=34 && arrayID~=30
        
        subsubtask_data=sprintf('FR_%s_%d_window_from_%d_to_%d',event_list{1},window_size,window_start_list(1),window_start_list(1)+window_size);
        %         subsubtask_data=sprintf('FR_%s_%d_window_from_%d_to_%d',event,window_size,window_start_list(1),window_start_list(1)+window_size);
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
            
            event=event_list{t};
            window_start=window_start_list(t);
            
            subsubtask_data=sprintf('FR_%s_%d_window_from_%d_to_%d',event,window_size,window_start,window_start+window_size);
            data_name=sprintf('ID_%d.mat',arrayID);
            dirstem = fullfile('Projects',task,'Analysis','Electrophy_analysis',subtask,subsubtask_data);
            data_path_clasifier = fullfile(fsroot,dirstem);
            
            load(fullfile(data_path_clasifier,data_name),'LIP')
            
            N_neurons=size(LIP,1);
            
            start_ind=2;
            end_ind=length(Conditions.block_nb);
            
            switch event
                case 'reward_end'
                    start_ind=1;
                    end_ind=length(Conditions.block_nb)-1;
            end
            for nn=1:N_neurons
                for i=start_ind:end_ind
                    Raw_data_LIP(nn,i-start_ind+1,t)=LIP(nn,i);
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
            for i=1:length(Conditions.block_nb)-1
                for t=1:length(window_start_list)
                    Pseudo_pop_LIP(count_LIP,i).Classifier_FR(t)=data_LIP(nn,i,t);
                end
                Pseudo_pop_LIP(count_LIP,i).Sess=n_sess;
                Pseudo_pop_LIP(count_LIP,i).Block_nb=Conditions.block_nb(i+1); %we need to take the updated belief here
                Pseudo_pop_LIP(count_LIP,i).Precision=Conditions.belief_precision(i+1);
                Pseudo_pop_LIP(count_LIP,i).Peak_belief=mod(Conditions.Peak_belief(i+1),2*pi);
                Pseudo_pop_LIP(count_LIP,i).progression_in_block=Conditions.Progression_in_block(i+1);
                
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
        subsubtask_data=sprintf('FR_%s_%d_window_from_%d_to_%d',event_list{1},window_size,window_start_list(1),window_start_list(1)+window_size);
        %     subsubtask_data=sprintf('FR_%s_%d_window_from_%d_to_%d',event,window_size,window_start_list(1),window_start_list(1)+window_size);
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
            
            event=event_list{t};
            window_start=window_start_list(t);
            
            subsubtask_data=sprintf('FR_%s_%d_window_from_%d_to_%d',event,window_size,window_start,window_start+window_size);
            data_name=sprintf('ID_%d.mat',arrayID);
            dirstem = fullfile('Projects',task,'Analysis','Electrophy_analysis',subtask,subsubtask_data);
            data_path_clasifier = fullfile(fsroot,dirstem);
            
            load(fullfile(data_path_clasifier,data_name),'FEF')
            
            N_neurons=size(FEF,1);
            
            start_ind=2;
            end_ind=length(Conditions.block_nb);
            
            switch event
                case 'reward_end'
                    start_ind=1;
                    end_ind=length(Conditions.block_nb)-1;
            end
            for nn=1:N_neurons
                for i=start_ind:end_ind
                    Raw_data_FEF(nn,i-start_ind+1,t)=FEF(nn,i);
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
            for i=1:length(Conditions.block_nb)-1
                for t=1:length(window_start_list)
                    Pseudo_pop_FEF(count_FEF,i).Classifier_FR(t)=data_FEF(nn,i,t);
                end
                Pseudo_pop_FEF(count_FEF,i).Sess=n_sess;
                Pseudo_pop_FEF(count_FEF,i).Block_nb=Conditions.block_nb(i+1); %we need to take the updated belief here
                Pseudo_pop_FEF(count_FEF,i).Precision=Conditions.belief_precision(i+1);
                Pseudo_pop_FEF(count_FEF,i).Peak_belief=mod(Conditions.Peak_belief(i+1),2*pi);
                Pseudo_pop_FEF(count_FEF,i).progression_in_block=Conditions.Progression_in_block(i+1);
                
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
    subsubtask_data=sprintf('FR_%s_%d_window_from_%d_to_%d',event_list{1},window_size,window_start_list(1),window_start_list(1)+window_size);
    %     subsubtask_data=sprintf('FR_%s_%d_window_from_%d_to_%d',event,window_size,window_start_list(1),window_start_list(1)+window_size);
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
        
        event=event_list{t};
        window_start=window_start_list(t);
        
        subsubtask_data=sprintf('FR_%s_%d_window_from_%d_to_%d',event,window_size,window_start,window_start+window_size);
        data_name=sprintf('ID_%d.mat',arrayID);
        dirstem = fullfile('Projects',task,'Analysis','Electrophy_analysis',subtask,subsubtask_data);
        data_path_clasifier = fullfile(fsroot,dirstem);
        
        load(fullfile(data_path_clasifier,data_name),'PFC')
        
        N_neurons=size(PFC,1);
        
        start_ind=2;
        end_ind=length(Conditions.block_nb);
        
        switch event
            case 'reward_end'
                start_ind=1;
                end_ind=length(Conditions.block_nb)-1;
        end
        for nn=1:N_neurons
            for i=start_ind:end_ind
                Raw_data_PFC(nn,i-start_ind+1,t)=PFC(nn,i);
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
        for i=1:length(Conditions.block_nb)-1
            for t=1:length(window_start_list)
                Pseudo_pop_PFC(count_PFC,i).Classifier_FR(t)=data_PFC(nn,i,t);
            end
            Pseudo_pop_PFC(count_PFC,i).Sess=n_sess;
            Pseudo_pop_PFC(count_PFC,i).Block_nb=Conditions.block_nb(i+1); %we need to take the updated belief here
            Pseudo_pop_PFC(count_PFC,i).Precision=Conditions.belief_precision(i+1);
            Pseudo_pop_PFC(count_PFC,i).Peak_belief=mod(Conditions.Peak_belief(i+1),2*pi);
            Pseudo_pop_PFC(count_PFC,i).progression_in_block=Conditions.Progression_in_block(i+1);
            
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


save(fullfile(save_path_clasifier,save_name),'Pseudo_pop_LIP','Pseudo_pop_FEF','Pseudo_pop_PFC');




















