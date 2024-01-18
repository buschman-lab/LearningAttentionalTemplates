%prepare_bootstrap_data_chosen_color_classifier

fsroot='/Volumes/buschman';

for i=1:27
    initial_window=-500;
    event_list{i}='target';
    window_start_list(i)=initial_window+(i-1)*50;
end

task='Learning_Attentional_Templates';
subtask='exploreexploit/Restricted_FEF/Reset_RW_model';

window_size=200;

subsubtask_classifier=sprintf('Classifier_%s_%d_to_%s_%d',event_list{1},window_start_list(1),event_list{end},window_start_list(end));
% subsubtask_classifier=sprintf('Classifier_%s_%d_to_%s_%d',event,window_start_list(1),event,window_start_list(end));

dirstem = fullfile('Projects',task,'Analysis','Electrophy_analysis',subtask,subsubtask_classifier);
data_path_clasifier = fullfile(fsroot,dirstem);

data_name='Pseudo_pop_peak_belief_chosen_color_classifier_2bins';

%if you don't run the pseudo-pop just before then load it
% load(fullfile(data_path_clasifier,data_name),'Pseudo_pop_LIP','Pseudo_pop_FEF','Pseudo_pop_PFC');

save_name='Pseudo_pop_chosen_color_peak_template_classifier_bootstrap_2bins';


%% Peak belief

N_belief=2;

N_bins=2;

for nn=1:size(Pseudo_pop_LIP,1)
    for i=1:size(Pseudo_pop_LIP(nn,:),2)
        if ~isempty(Pseudo_pop_LIP(nn,i).progression_in_block) && ~isempty(Pseudo_pop_LIP(nn,i).Belief_color)
            if ~isnan(Pseudo_pop_LIP(nn,i).progression_in_block) && ~isnan(Pseudo_pop_LIP(nn,i).Belief_color) && ~isnan(Pseudo_pop_LIP(nn,i).Chosen_color_bin)
                selected_LIP_chosen_color(nn,i,Pseudo_pop_LIP(nn,i).Chosen_color_bin,Pseudo_pop_LIP(nn,i).Belief_color)=1;
            end
        end
    end
    Sess_LIP(nn)= Pseudo_pop_LIP(nn,1).Sess;
end
for nn=1:size(Pseudo_pop_FEF,1)
    for i=1:size(Pseudo_pop_FEF(nn,:),2)
        if ~isempty(Pseudo_pop_FEF(nn,i).progression_in_block) && ~isempty(Pseudo_pop_FEF(nn,i).Belief_color)
            if ~isnan(Pseudo_pop_FEF(nn,i).progression_in_block) && ~isnan(Pseudo_pop_FEF(nn,i).Belief_color) && ~isnan(Pseudo_pop_FEF(nn,i).Chosen_color_bin)
                selected_FEF_chosen_color(nn,i,Pseudo_pop_FEF(nn,i).Chosen_color_bin,Pseudo_pop_FEF(nn,i).Belief_color)=1;
            end
        end
    end
    Sess_FEF(nn)= Pseudo_pop_FEF(nn,1).Sess;
end
for nn=1:size(Pseudo_pop_PFC,1)
    for i=1:size(Pseudo_pop_PFC(nn,:),2)
        if ~isempty(Pseudo_pop_PFC(nn,i).progression_in_block) && ~isempty(Pseudo_pop_PFC(nn,i).Belief_color)
            if ~isnan(Pseudo_pop_PFC(nn,i).progression_in_block) && ~isnan(Pseudo_pop_PFC(nn,i).Belief_color) && ~isnan(Pseudo_pop_PFC(nn,i).Chosen_color_bin)
                selected_PFC_chosen_color(nn,i,Pseudo_pop_PFC(nn,i).Chosen_color_bin,Pseudo_pop_PFC(nn,i).Belief_color)=1;
            end
        end
    end
    Sess_PFC(nn)= Pseudo_pop_PFC(nn,1).Sess;
end


%%

for nn=1:size(Pseudo_pop_LIP,1)
    if min(sum(selected_LIP_chosen_color(nn,:,:,:),2),[],'all')>=60
        I_LIP(nn)=1;
        Min_LIP(nn)=min(sum(selected_LIP_chosen_color(nn,:,:,:),2),[],'all');
    else
        I_LIP(nn)=0;
        Min_LIP(nn)=NaN;
    end
end
for nn=1:size(Pseudo_pop_FEF,1)
    if min(sum(selected_FEF_chosen_color(nn,:,:,:),2),[],'all')>=60
        I_FEF(nn)=1;
        Min_FEF(nn)=min(sum(selected_FEF_chosen_color(nn,:,:,:),2),[],'all');
    else
        I_FEF(nn)=0;
        Min_FEF(nn)=NaN;
    end
end
for nn=1:size(Pseudo_pop_PFC,1)
    if min(sum(selected_PFC_chosen_color(nn,:,:,:),2),[],'all')>=60
        I_PFC(nn)=1;
        Min_PFC(nn)=min(sum(selected_PFC_chosen_color(nn,:,:,:),2),[],'all');
    else
        I_PFC(nn)=0;
        Min_PFC(nn)=NaN;
    end
end

N_chosen_color=min([Min_LIP,Min_FEF,Min_PFC],[],'omitnan');
N_neurons=min([sum(I_LIP),sum(I_FEF),sum(I_PFC)]);

%% now create the bootstrap
N_bootstrap=100;

for n=1:N_bootstrap
    
    for nn=1:N_neurons
        Pseudo_pop_bootstrap.LIP.Neurons(n,nn)=find(cumsum(I_LIP)==randsample(1:sum(I_LIP),1,'true'),1,'first');
        Pseudo_pop_bootstrap.FEF.Neurons(n,nn)=find(cumsum(I_FEF)==randsample(1:sum(I_FEF),1,'true'),1,'first');
        Pseudo_pop_bootstrap.PFC.Neurons(n,nn)=find(cumsum(I_PFC)==randsample(1:sum(I_PFC),1,'true'),1,'first');
    end
    
    for k=1:N_belief
        for nn=1:N_neurons
            for n_sc=1:N_bins
                Pseudo_pop_bootstrap.LIP.chosen_color.Trials(n,nn,:,n_sc,k)=randsample(find(selected_LIP_chosen_color(Pseudo_pop_bootstrap.LIP.Neurons(n,nn),:,n_sc,k)==1),N_chosen_color,'false');
                Pseudo_pop_bootstrap.FEF.chosen_color.Trials(n,nn,:,n_sc,k)=randsample(find(selected_FEF_chosen_color(Pseudo_pop_bootstrap.FEF.Neurons(n,nn),:,n_sc,k)==1),N_chosen_color,'false');
                Pseudo_pop_bootstrap.PFC.chosen_color.Trials(n,nn,:,n_sc,k)=randsample(find(selected_PFC_chosen_color(Pseudo_pop_bootstrap.PFC.Neurons(n,nn),:,n_sc,k)==1),N_chosen_color,'false');
                
                h_chosen_color=cvpartition(N_chosen_color,'Holdout',0.2);
                Pseudo_pop_bootstrap.LIP.chosen_color.TT_trials(n,nn,:,n_sc,k)=Pseudo_pop_bootstrap.LIP.chosen_color.Trials(n,nn,training(h_chosen_color),n_sc,k);
                Pseudo_pop_bootstrap.LIP.chosen_color.Val_trials(n,nn,:,n_sc,k)=Pseudo_pop_bootstrap.LIP.chosen_color.Trials(n,nn,test(h_chosen_color),n_sc,k);
                
                Pseudo_pop_bootstrap.FEF.chosen_color.TT_trials(n,nn,:,n_sc,k)=Pseudo_pop_bootstrap.FEF.chosen_color.Trials(n,nn,training(h_chosen_color),n_sc,k);
                Pseudo_pop_bootstrap.FEF.chosen_color.Val_trials(n,nn,:,n_sc,k)=Pseudo_pop_bootstrap.FEF.chosen_color.Trials(n,nn,test(h_chosen_color),n_sc,k);
                
                Pseudo_pop_bootstrap.PFC.chosen_color.TT_trials(n,nn,:,n_sc,k)=Pseudo_pop_bootstrap.PFC.chosen_color.Trials(n,nn,training(h_chosen_color),n_sc,k);
                Pseudo_pop_bootstrap.PFC.chosen_color.Val_trials(n,nn,:,n_sc,k)=Pseudo_pop_bootstrap.PFC.chosen_color.Trials(n,nn,test(h_chosen_color),n_sc,k);
            end
        end
    end
end

%% save

Pseudo_pop_bootstrap.LIP.Pseudo_pop=Pseudo_pop_LIP;
Pseudo_pop_bootstrap.FEF.Pseudo_pop=Pseudo_pop_FEF;
Pseudo_pop_bootstrap.PFC.Pseudo_pop=Pseudo_pop_PFC;

 
save(fullfile(data_path_clasifier,save_name),'Pseudo_pop_bootstrap','-v7.3')






























