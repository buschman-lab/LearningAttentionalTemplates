%prepare_bootstrap_peak_belief_classifier


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
% subtask='exploreexploit/Reset_RW_model';
subtask='exploreexploit/Restricted_FEF/Reset_RW_model';

window_size=300;

subsubtask_classifier=sprintf('Classifier_%s_%d_to_%s_%d',event_list{1},window_start_list(1),event_list{end},window_start_list(end));
% subsubtask_classifier=sprintf('Classifier_%s_%d_to_%s_%d',event,window_start_list(1),event,window_start_list(end));
dirstem = fullfile('Projects',task,'Analysis','Electrophy_analysis',subtask,subsubtask_classifier);
data_path_clasifier = fullfile(fsroot,dirstem);

data_name='Pseudo_pop_peak_belief_prog_classifier'; %added peak

load(fullfile(data_path_clasifier,data_name),'Pseudo_pop_LIP','Pseudo_pop_FEF','Pseudo_pop_PFC');

save_name='Pseudo_pop_peak_belief_prog_classifier_bootstrap_rotated_60'; %added peak


%% Progression

N_prog=3;
N_bins=3;
N_channels_stim=3;  
rotation_offset=pi/3;


for nn=1:size(Pseudo_pop_LIP,1)
    for i=1:size(Pseudo_pop_LIP(nn,:),2)
        if ~isempty(Pseudo_pop_LIP(nn,i).progression_in_block) && ~isnan(Pseudo_pop_LIP(nn,i).progression_in_block)
            if logical(Pseudo_pop_LIP(nn,i).Belief_color>0)
                %first we add 2*pi to angles below pi/3
                Pseudo_pop_LIP(nn,i).Peak_belief=mod(Pseudo_pop_LIP(nn,i).Peak_belief,2*pi);
                if Pseudo_pop_LIP(nn,i).Peak_belief<rotation_offset
                    Pseudo_pop_LIP(nn,i).Peak_belief=Pseudo_pop_LIP(nn,i).Peak_belief+2*pi; %get angles smaller than pi/3 to move to x + 2*pi
                end
                %change the binning to offset
                for n_c=1:N_channels_stim
                    if Pseudo_pop_LIP(nn,i).Peak_belief>=(n_c-1)*2*pi/N_channels_stim + rotation_offset && Pseudo_pop_LIP(nn,i).Peak_belief<n_c*2*pi/N_channels_stim + rotation_offset
                        Pseudo_pop_LIP(nn,i).Belief_color=n_c;
                    end
                end
                
                if Pseudo_pop_LIP(nn,i).progression_in_block<=1/N_prog
                    selected_LIP_progression(nn,i,1,Pseudo_pop_LIP(nn,i).Belief_color)=1;
                end
                for j=2:N_prog-1
                    if Pseudo_pop_LIP(nn,i).progression_in_block>(j-1)/N_prog && Pseudo_pop_LIP(nn,i).progression_in_block<=j/N_prog
                        selected_LIP_progression(nn,i,j,Pseudo_pop_LIP(nn,i).Belief_color)=1;
                    end
                end
                if Pseudo_pop_LIP(nn,i).progression_in_block>=(N_prog-1)/N_prog
                    selected_LIP_progression(nn,i,N_prog,Pseudo_pop_LIP(nn,i).Belief_color)=1;
                end
            end
        end
    end
    Sess_LIP(nn)= Pseudo_pop_LIP(nn,1).Sess;
end

for nn=1:size(Pseudo_pop_FEF,1)
    for i=1:size(Pseudo_pop_FEF(nn,:),2)
        if ~isempty(Pseudo_pop_FEF(nn,i).progression_in_block) && ~isnan(Pseudo_pop_FEF(nn,i).progression_in_block)
            if Pseudo_pop_FEF(nn,i).Belief_color>0
                %first we add 2*pi to angles below pi/3
                Pseudo_pop_FEF(nn,i).Peak_belief=mod(Pseudo_pop_FEF(nn,i).Peak_belief,2*pi);
                if Pseudo_pop_FEF(nn,i).Peak_belief<rotation_offset
                    Pseudo_pop_FEF(nn,i).Peak_belief=Pseudo_pop_FEF(nn,i).Peak_belief+2*pi; %get angles smaller than pi/3 to move to x + 2*pi
                end
                %change the binning to offset
                for n_c=1:N_channels_stim
                    if Pseudo_pop_FEF(nn,i).Peak_belief>=(n_c-1)*2*pi/N_channels_stim + rotation_offset && Pseudo_pop_FEF(nn,i).Peak_belief<n_c*2*pi/N_channels_stim + rotation_offset
                        Pseudo_pop_FEF(nn,i).Belief_color=n_c;
                    end
                end

                if Pseudo_pop_FEF(nn,i).progression_in_block<=1/N_prog
                    selected_FEF_progression(nn,i,1,Pseudo_pop_FEF(nn,i).Belief_color)=1;
                end
                for j=2:N_prog-1
                    if Pseudo_pop_FEF(nn,i).progression_in_block>(j-1)/N_prog && Pseudo_pop_FEF(nn,i).progression_in_block<=j/N_prog
                        selected_FEF_progression(nn,i,j,Pseudo_pop_FEF(nn,i).Belief_color)=1;
                    end
                end
                if Pseudo_pop_FEF(nn,i).progression_in_block>=(N_prog-1)/N_prog
                    selected_FEF_progression(nn,i,N_prog,Pseudo_pop_FEF(nn,i).Belief_color)=1;
                end
            end
        end
    end
    Sess_FEF(nn)= Pseudo_pop_FEF(nn,1).Sess;
end

for nn=1:size(Pseudo_pop_PFC,1)
    for i=1:size(Pseudo_pop_PFC(nn,:),2)
        if ~isempty(Pseudo_pop_PFC(nn,i).progression_in_block) && ~isnan(Pseudo_pop_PFC(nn,i).progression_in_block)
            if Pseudo_pop_PFC(nn,i).Belief_color>0
                %first we add 2*pi to angles below pi/3
                Pseudo_pop_PFC(nn,i).Peak_belief=mod(Pseudo_pop_PFC(nn,i).Peak_belief,2*pi);
                if Pseudo_pop_PFC(nn,i).Peak_belief<rotation_offset
                    Pseudo_pop_PFC(nn,i).Peak_belief=Pseudo_pop_PFC(nn,i).Peak_belief+2*pi; %get angles smaller than pi/3 to move to x + 2*pi
                end
                %change the binning to offset
                for n_c=1:N_channels_stim
                    if Pseudo_pop_PFC(nn,i).Peak_belief>=(n_c-1)*2*pi/N_channels_stim + rotation_offset && Pseudo_pop_PFC(nn,i).Peak_belief<n_c*2*pi/N_channels_stim + rotation_offset
                        Pseudo_pop_PFC(nn,i).Belief_color=n_c;
                    end
                end
                if Pseudo_pop_PFC(nn,i).progression_in_block<=1/N_prog
                    selected_PFC_progression(nn,i,1,Pseudo_pop_PFC(nn,i).Belief_color)=1;
                end
                for j=2:N_prog-1
                    if Pseudo_pop_PFC(nn,i).progression_in_block>(j-1)/N_prog && Pseudo_pop_PFC(nn,i).progression_in_block<=j/N_prog
                        selected_PFC_progression(nn,i,j,Pseudo_pop_PFC(nn,i).Belief_color)=1;
                    end
                end
                if Pseudo_pop_PFC(nn,i).progression_in_block>=(N_prog-1)/N_prog
                    selected_PFC_progression(nn,i,N_prog,Pseudo_pop_PFC(nn,i).Belief_color)=1;
                end
            end
        end
    end
    Sess_PFC(nn)= Pseudo_pop_PFC(nn,1).Sess;
end

%%
for nn=1:size(Pseudo_pop_LIP,1)
    if min(sum(selected_LIP_progression(nn,:,:,:),2),[],'all')>=60
        I_LIP(nn)=1;
        Min_LIP(nn)=min(sum(selected_LIP_progression(nn,:,:,:),2),[],'all');
    else
        I_LIP(nn)=0;
        Min_LIP(nn)=NaN;
    end
end
for nn=1:size(Pseudo_pop_FEF,1)
    if min(sum(selected_FEF_progression(nn,:,:,:),2),[],'all')>=60
        I_FEF(nn)=1;
        Min_FEF(nn)=min(sum(selected_FEF_progression(nn,:,:,:),2),[],'all');
    else
        I_FEF(nn)=0;
        Min_FEF(nn)=NaN;
    end
end
for nn=1:size(Pseudo_pop_PFC,1)
    if min(sum(selected_PFC_progression(nn,:,:,:),2),[],'all')>=60
        I_PFC(nn)=1;
        Min_PFC(nn)=min(sum(selected_PFC_progression(nn,:,:,:),2),[],'all');
    else
        I_PFC(nn)=0;
        Min_PFC(nn)=NaN;
    end
end


N_progression=min([Min_LIP,Min_FEF,Min_PFC],[],'omitnan');
N_neurons=min([sum(I_LIP),sum(I_FEF),sum(I_PFC)]);

%% now create the bootstrap
N_bootstrap=100;

for n=1:N_bootstrap
    
    for nn=1:N_neurons
        Pseudo_pop_bootstrap.LIP.Neurons(n,nn)=find(cumsum(I_LIP)==randsample(1:sum(I_LIP),1,'true'),1,'first');
        Pseudo_pop_bootstrap.FEF.Neurons(n,nn)=find(cumsum(I_FEF)==randsample(1:sum(I_FEF),1,'true'),1,'first');
        Pseudo_pop_bootstrap.PFC.Neurons(n,nn)=find(cumsum(I_PFC)==randsample(1:sum(I_PFC),1,'true'),1,'first');
    end
    
    %progression
    for k=1:N_prog
        for nn=1:N_neurons
            for n_c=1:N_bins
                Pseudo_pop_bootstrap.LIP.progression.Trials(n,nn,:,k,n_c)=randsample(find(selected_LIP_progression(Pseudo_pop_bootstrap.LIP.Neurons(n,nn),:,k,n_c)==1),N_progression,'false');
                Pseudo_pop_bootstrap.FEF.progression.Trials(n,nn,:,k,n_c)=randsample(find(selected_FEF_progression(Pseudo_pop_bootstrap.FEF.Neurons(n,nn),:,k,n_c)==1),N_progression,'false');
                Pseudo_pop_bootstrap.PFC.progression.Trials(n,nn,:,k,n_c)=randsample(find(selected_PFC_progression(Pseudo_pop_bootstrap.PFC.Neurons(n,nn),:,k,n_c)==1),N_progression,'false');
                
                h_progression=cvpartition(N_progression,'Holdout',0.2);
                Pseudo_pop_bootstrap.LIP.progression.Val_trials(n,nn,:,k,n_c)=Pseudo_pop_bootstrap.LIP.progression.Trials(n,nn,test(h_progression),k,n_c);
                Pseudo_pop_bootstrap.LIP.progression.TT_trials(n,nn,:,k,n_c)=Pseudo_pop_bootstrap.LIP.progression.Trials(n,nn,training(h_progression),k,n_c);
                
                Pseudo_pop_bootstrap.FEF.progression.Val_trials(n,nn,:,k,n_c)=Pseudo_pop_bootstrap.FEF.progression.Trials(n,nn,test(h_progression),k,n_c);
                Pseudo_pop_bootstrap.FEF.progression.TT_trials(n,nn,:,k,n_c)=Pseudo_pop_bootstrap.FEF.progression.Trials(n,nn,training(h_progression),k,n_c);
                
                Pseudo_pop_bootstrap.PFC.progression.Val_trials(n,nn,:,k,n_c)=Pseudo_pop_bootstrap.PFC.progression.Trials(n,nn,test(h_progression),k,n_c);
                Pseudo_pop_bootstrap.PFC.progression.TT_trials(n,nn,:,k,n_c)=Pseudo_pop_bootstrap.PFC.progression.Trials(n,nn,training(h_progression),k,n_c);
            end
        end
    end
    
end


%% Blocks the trials were taken from
n=1;
    for k=1:N_prog
        for nn=1:N_neurons
            for n_c=1:N_bins
                for i=1:size(Pseudo_pop_bootstrap.LIP.progression.TT_trials,3)
                    Block_TT(i,nn,n_c,k)=Pseudo_pop_LIP(Pseudo_pop_bootstrap.LIP.Neurons(n,nn),Pseudo_pop_bootstrap.LIP.progression.TT_trials(n,nn,i,k,n_c)).Block_nb;
                end
            end
        end
    end
    
    for nn=1:N_neurons
        for n_c=1:N_bins
            for k=1:N_prog
            Unique_block(nn,n_c,k)=length(unique(Block_TT(:,nn,n_c,k)));
            
            end
        end
    end
    
    figure
    for i=1:3
    subplot(1,3,i)
    histogram(Unique_block(:,i,1))
    end


%% save

Pseudo_pop_bootstrap.LIP.Pseudo_pop=Pseudo_pop_LIP;
Pseudo_pop_bootstrap.FEF.Pseudo_pop=Pseudo_pop_FEF;
Pseudo_pop_bootstrap.PFC.Pseudo_pop=Pseudo_pop_PFC;


save(fullfile(data_path_clasifier,save_name),'Pseudo_pop_bootstrap','-v7.3')






























