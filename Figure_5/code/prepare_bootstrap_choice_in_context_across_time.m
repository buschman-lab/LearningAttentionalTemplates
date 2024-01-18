%prepare_bootstrap_data_choice_classifier

%%choice context to choices offered
% 1 = 1 2 3
% 2 = 1 2 4
% 3 = 1 3 4
% 4 = 2 3 4

fsroot='/Volumes/buschman';

for i=1:23
    initial_window=-500;
    event_list{i}='response';
    window_start_list(i)=initial_window+(i-1)*50;
end

task='Learning_Attentional_Templates';
subtask='exploreexploit/Restricted_FEF/Reset_RW_model';

window_size=200;

subsubtask_classifier=sprintf('Classifier_%s_%d_to_%s_%d',event_list{1},window_start_list(1),event_list{end},window_start_list(end));
% subsubtask_classifier=sprintf('Classifier_%s_%d_to_%s_%d',event,window_start_list(1),event,window_start_list(end));

dirstem = fullfile('Projects',task,'Analysis','Electrophy_analysis',subtask,subsubtask_classifier);
data_path_clasifier = fullfile(fsroot,dirstem);

data_name='Pseudo_pop_choice_in_context_classifier';

% load(fullfile(data_path_clasifier,data_name),'Pseudo_pop_LIP','Pseudo_pop_FEF','Pseudo_pop_PFC');

save_name='Pseudo_pop_choice_in_context_classifier_bootstrap';


%% Peak belief

N_context=4;
N_choices=3;

for nn=1:size(Pseudo_pop_LIP,1)
    for i=1:size(Pseudo_pop_LIP(nn,:),2)
        if ~isempty(Pseudo_pop_LIP(nn,i).progression_in_block) && ~isempty(Pseudo_pop_LIP(nn,i).Choice_context)
            if ~isnan(Pseudo_pop_LIP(nn,i).progression_in_block) && ~isnan(Pseudo_pop_LIP(nn,i).Choice_context)
                selected_LIP_choice(nn,i,Pseudo_pop_LIP(nn,i).Choice_in_context,Pseudo_pop_LIP(nn,i).Choice_context)=1;
            end
        end
    end
    Sess_LIP(nn)= Pseudo_pop_LIP(nn,1).Sess;
end

for nn=1:size(Pseudo_pop_FEF,1)
    for i=1:size(Pseudo_pop_FEF(nn,:),2)
        if ~isempty(Pseudo_pop_FEF(nn,i).progression_in_block) && ~isempty(Pseudo_pop_FEF(nn,i).Choice_context) 
            if ~isnan(Pseudo_pop_FEF(nn,i).progression_in_block)  && ~isnan(Pseudo_pop_FEF(nn,i).Choice_context)
            selected_FEF_choice(nn,i,Pseudo_pop_FEF(nn,i).Choice_in_context,Pseudo_pop_FEF(nn,i).Choice_context)=1;
            end
        end
    end
    Sess_FEF(nn)= Pseudo_pop_FEF(nn,1).Sess;
end

for nn=1:size(Pseudo_pop_PFC,1)
    for i=1:size(Pseudo_pop_PFC(nn,:),2)
        if ~isempty(Pseudo_pop_PFC(nn,i).progression_in_block) && ~isempty(Pseudo_pop_PFC(nn,i).Choice_context) 
            if ~isnan(Pseudo_pop_PFC(nn,i).progression_in_block)  && ~isnan(Pseudo_pop_PFC(nn,i).Choice_context)
            selected_PFC_choice(nn,i,Pseudo_pop_PFC(nn,i).Choice_in_context,Pseudo_pop_PFC(nn,i).Choice_context)=1;
            end
        end
    end
    Sess_PFC(nn)= Pseudo_pop_PFC(nn,1).Sess;
end

%%

for j=1:N_context
    for nn=1:size(Pseudo_pop_LIP,1)
        if min(sum(selected_LIP_choice(nn,:,:,j),2),[],'all')>=60
            I_LIP(nn,j)=1;
            Min_LIP(nn,j)=min(sum(selected_LIP_choice(nn,:,:,j),2),[],'all');
        else
            I_LIP(nn,j)=0;
            Min_LIP(nn,j)=NaN;
        end
    end
    for nn=1:size(Pseudo_pop_FEF,1)
        if min(sum(selected_FEF_choice(nn,:,:,j),2),[],'all')>=60
            I_FEF(nn,j)=1;
            Min_FEF(nn,j)=min(sum(selected_FEF_choice(nn,:,:,j),2),[],'all');
        else
            I_FEF(nn,j)=0;
            Min_FEF(nn,j)=NaN;
        end
    end
    for nn=1:size(Pseudo_pop_PFC,1)
        if min(sum(selected_PFC_choice(nn,:,:,j),2),[],'all')>=60
            I_PFC(nn,j)=1;
            Min_PFC(nn,j)=min(sum(selected_PFC_choice(nn,:,:,j),2),[],'all');
        else
            I_PFC(nn,j)=0;
            Min_PFC(nn,j)=NaN;
        end
    end
    
end

for nn=1:size(Pseudo_pop_LIP,1)
    if isnan(sum(Min_LIP(nn,:)))
        Min_LIP(nn,:)=NaN;
        I_LIP(nn,:)=0;
    end
end
for nn=1:size(Pseudo_pop_FEF,1)
    if isnan(sum(Min_FEF(nn,:)))
        Min_FEF(nn,:)=NaN;
        I_FEF(nn,:)=0;
    end
end
for nn=1:size(Pseudo_pop_PFC,1)
    if isnan(sum(Min_PFC(nn,:)))
        Min_PFC(nn,:)=NaN;
        I_PFC(nn,:)=0;
    end
end

for j=1:4
    N_choice_loc(j)=min([Min_LIP(:,j);Min_FEF(:,j);Min_PFC(:,j)],[],'omitnan');
    N_neurons(j)=min([sum(I_LIP(:,j)),sum(I_FEF(:,j)),sum(I_PFC(:,j))]);
end

N_choice=min(N_choice_loc(:),[],'omitnan');

%% now create the bootstrap
N_bootstrap=100;

for n=1:N_bootstrap
    
    for nn=1:N_neurons
        Pseudo_pop_bootstrap.LIP.Neurons(n,nn)=find(cumsum(I_LIP)==randsample(1:sum(I_LIP),1,'true'),1,'first');
        Pseudo_pop_bootstrap.FEF.Neurons(n,nn)=find(cumsum(I_FEF)==randsample(1:sum(I_FEF),1,'true'),1,'first');
        Pseudo_pop_bootstrap.PFC.Neurons(n,nn)=find(cumsum(I_PFC)==randsample(1:sum(I_PFC),1,'true'),1,'first');
    end
    
    %for each location
        for k=1:N_context
            for nn=1:N_neurons
                for n_c=1:N_choices
                    Pseudo_pop_bootstrap.LIP.Choice_context.Trials(n,nn,:,n_c,k)=randsample(find(selected_LIP_choice(Pseudo_pop_bootstrap.LIP.Neurons(n,nn),:,n_c,k)==1),N_choice,'false');
                    Pseudo_pop_bootstrap.FEF.Choice_context.Trials(n,nn,:,n_c,k)=randsample(find(selected_FEF_choice(Pseudo_pop_bootstrap.FEF.Neurons(n,nn),:,n_c,k)==1),N_choice,'false');
                    Pseudo_pop_bootstrap.PFC.Choice_context.Trials(n,nn,:,n_c,k)=randsample(find(selected_PFC_choice(Pseudo_pop_bootstrap.PFC.Neurons(n,nn),:,n_c,k)==1),N_choice,'false');
                    
                    h_choice=cvpartition(N_choice,'Holdout',0.2);
                    Pseudo_pop_bootstrap.LIP.Choice_context.TT_trials(n,nn,:,n_c,k)=Pseudo_pop_bootstrap.LIP.Choice_context.Trials(n,nn,training(h_choice),n_c,k);
                    Pseudo_pop_bootstrap.LIP.Choice_context.Val_trials(n,nn,:,n_c,k)=Pseudo_pop_bootstrap.LIP.Choice_context.Trials(n,nn,test(h_choice),n_c,k);
                    
                    Pseudo_pop_bootstrap.FEF.Choice_context.TT_trials(n,nn,:,n_c,k)=Pseudo_pop_bootstrap.FEF.Choice_context.Trials(n,nn,training(h_choice),n_c,k);
                    Pseudo_pop_bootstrap.FEF.Choice_context.Val_trials(n,nn,:,n_c,k)=Pseudo_pop_bootstrap.FEF.Choice_context.Trials(n,nn,test(h_choice),n_c,k);
                    
                    Pseudo_pop_bootstrap.PFC.Choice_context.TT_trials(n,nn,:,n_c,k)=Pseudo_pop_bootstrap.PFC.Choice_context.Trials(n,nn,training(h_choice),n_c,k);
                    Pseudo_pop_bootstrap.PFC.Choice_context.Val_trials(n,nn,:,n_c,k)=Pseudo_pop_bootstrap.PFC.Choice_context.Trials(n,nn,test(h_choice),n_c,k);
                end
            end
        end
end

%% save

Pseudo_pop_bootstrap.LIP.Pseudo_pop=Pseudo_pop_LIP;
Pseudo_pop_bootstrap.FEF.Pseudo_pop=Pseudo_pop_FEF;
Pseudo_pop_bootstrap.PFC.Pseudo_pop=Pseudo_pop_PFC;

 
save(fullfile(data_path_clasifier,save_name),'Pseudo_pop_bootstrap','-v7.3')






























