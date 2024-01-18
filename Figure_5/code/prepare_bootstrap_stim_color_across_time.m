%stim color half chosen and half unchosen

fsroot='/Volumes/buschman';

event='target';

for i=1:27
    window_start_list(i)=-500+(i-1)*50;
end

task='Learning_Attentional_Templates';
subtask='exploreexploit/Restricted_FEF/Reset_RW_model';

window_size=200;

subsubtask_classifier=sprintf('Classifier_%s_%d_to_%s_%d',event,window_start_list(1),event,window_start_list(end));
dirstem = fullfile('Projects',task,'Analysis','Electrophy_analysis',subtask,subsubtask_classifier);
data_path_clasifier = fullfile(fsroot,dirstem);

data_name='Pseudo_pop_peak_belief_stim_color_prog_classifier';

%load if you didn't run it just before
load(fullfile(data_path_clasifier,data_name),'Pseudo_pop_LIP','Pseudo_pop_FEF','Pseudo_pop_PFC');

save_name='Pseudo_pop_stim_color_no_prog_bootstrap';


%% Equal mix of chosen and unchosen

N_bins=3;

for nn=1:size(Pseudo_pop_LIP,1)
    for i=1:size(Pseudo_pop_LIP(nn,:),2)
        if ~isempty(Pseudo_pop_LIP(nn,i).progression_in_block) && ~isnan(Pseudo_pop_LIP(nn,i).progression_in_block) && ~isempty(Pseudo_pop_LIP(nn,i).Stim_color)
            for k=1:4 %loc
                if ~isnan(Pseudo_pop_LIP(nn,i).Stim_color(k))
                    if Pseudo_pop_LIP(nn,i).choice==k
                        selected_LIP_stim_color(nn,i,Pseudo_pop_LIP(nn,i).Stim_color(k),1,k)=1;
                    elseif Pseudo_pop_LIP(nn,i).choice~=k
                        selected_LIP_stim_color(nn,i,Pseudo_pop_LIP(nn,i).Stim_color(k),2,k)=1;
                    end
                end
            end
        end
    end
    Sess_LIP(nn)= Pseudo_pop_LIP(nn,1).Sess;
end

for nn=1:size(Pseudo_pop_FEF,1)
    for i=1:size(Pseudo_pop_FEF(nn,:),2)
        if ~isempty(Pseudo_pop_FEF(nn,i).progression_in_block) && ~isnan(Pseudo_pop_FEF(nn,i).progression_in_block) && ~isempty(Pseudo_pop_FEF(nn,i).Stim_color)
            for k=1:4 %loc
                if ~isnan(Pseudo_pop_FEF(nn,i).Stim_color(k))
                    if Pseudo_pop_FEF(nn,i).choice==k
                        selected_FEF_stim_color(nn,i,Pseudo_pop_FEF(nn,i).Stim_color(k),1,k)=1;
                    elseif Pseudo_pop_FEF(nn,i).choice~=k
                        selected_FEF_stim_color(nn,i,Pseudo_pop_FEF(nn,i).Stim_color(k),2,k)=1;
                    end
                end
            end
        end
    end
    Sess_FEF(nn)= Pseudo_pop_FEF(nn,1).Sess;
end

for nn=1:size(Pseudo_pop_PFC,1)
    for i=1:size(Pseudo_pop_PFC(nn,:),2)
        if ~isempty(Pseudo_pop_PFC(nn,i).progression_in_block) && ~isnan(Pseudo_pop_PFC(nn,i).progression_in_block) && ~isempty(Pseudo_pop_PFC(nn,i).Stim_color)
            for k=1:4 %loc
                if ~isnan(Pseudo_pop_PFC(nn,i).Stim_color(k))
                    if Pseudo_pop_PFC(nn,i).choice==k
                        selected_PFC_stim_color(nn,i,Pseudo_pop_PFC(nn,i).Stim_color(k),1,k)=1;
                    elseif Pseudo_pop_PFC(nn,i).choice~=k
                        selected_PFC_stim_color(nn,i,Pseudo_pop_PFC(nn,i).Stim_color(k),2,k)=1;
                    end
                end
            end
        end
    end
    Sess_PFC(nn)= Pseudo_pop_PFC(nn,1).Sess;
end

%%

for j=1:4
    for nn=1:size(Pseudo_pop_LIP,1)
        if min(sum(selected_LIP_stim_color(nn,:,:,:,j),2),[],'all')>=60
            I_LIP(nn,j)=1;
            Min_LIP(nn,j)=min(sum(selected_LIP_stim_color(nn,:,:,:,j),2),[],'all');
        else
            I_LIP(nn,j)=0;
            Min_LIP(nn,j)=NaN;
        end
    end
    for nn=1:size(Pseudo_pop_FEF,1)
        if min(sum(selected_FEF_stim_color(nn,:,:,:,j),2),[],'all')>=60
            I_FEF(nn,j)=1;
            Min_FEF(nn,j)=min(sum(selected_FEF_stim_color(nn,:,:,:,j),2),[],'all');
        else
            I_FEF(nn,j)=0;
            Min_FEF(nn,j)=NaN;
        end
    end
    for nn=1:size(Pseudo_pop_PFC,1)
        if min(sum(selected_PFC_stim_color(nn,:,:,:,j),2),[],'all')>=60
            I_PFC(nn,j)=1;
            Min_PFC(nn,j)=min(sum(selected_PFC_stim_color(nn,:,:,:,j),2),[],'all');
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
    N_stim_color_loc(j)=min([Min_LIP(:,j);Min_FEF(:,j);Min_PFC(:,j)],[],'omitnan');
    N_neurons(j)=min([sum(I_LIP(:,j)),sum(I_FEF(:,j)),sum(I_PFC(:,j))]);
end

N_stim_color=min(N_stim_color_loc(:),[],'omitnan');

%% now create the bootstrap
N_bootstrap=100;
%
for n=1:N_bootstrap
    
    for nn=1:N_neurons
        Pseudo_pop_bootstrap.LIP.Neurons(n,nn)=find(cumsum(I_LIP)==randsample(1:sum(I_LIP),1,'true'),1,'first');
        Pseudo_pop_bootstrap.FEF.Neurons(n,nn)=find(cumsum(I_FEF)==randsample(1:sum(I_FEF),1,'true'),1,'first');
        Pseudo_pop_bootstrap.PFC.Neurons(n,nn)=find(cumsum(I_PFC)==randsample(1:sum(I_PFC),1,'true'),1,'first');
    end
    
    %for each location
    for j=1:4  %loc
        for k=1:2 %chosen_unchosen
            for nn=1:N_neurons
                for n_sc=1:N_bins
                    Pseudo_pop_bootstrap.LIP.Stim_color.Trials(n,nn,:,n_sc,k,j)=randsample(find(selected_LIP_stim_color(Pseudo_pop_bootstrap.LIP.Neurons(n,nn),:,n_sc,k,j)==1),N_stim_color,'false');
                    Pseudo_pop_bootstrap.FEF.Stim_color.Trials(n,nn,:,n_sc,k,j)=randsample(find(selected_FEF_stim_color(Pseudo_pop_bootstrap.FEF.Neurons(n,nn),:,n_sc,k,j)==1),N_stim_color,'false');
                    Pseudo_pop_bootstrap.PFC.Stim_color.Trials(n,nn,:,n_sc,k,j)=randsample(find(selected_PFC_stim_color(Pseudo_pop_bootstrap.PFC.Neurons(n,nn),:,n_sc,k,j)==1),N_stim_color,'false');
                    
                    h_stim_color=cvpartition(N_stim_color,'Holdout',0.2);
                    Pseudo_pop_bootstrap.LIP.Stim_color.TT_trials(n,nn,:,n_sc,k,j)=Pseudo_pop_bootstrap.LIP.Stim_color.Trials(n,nn,training(h_stim_color),n_sc,k,j);
                    Pseudo_pop_bootstrap.LIP.Stim_color.Val_trials(n,nn,:,n_sc,k,j)=Pseudo_pop_bootstrap.LIP.Stim_color.Trials(n,nn,test(h_stim_color),n_sc,k,j);
                    
                    Pseudo_pop_bootstrap.FEF.Stim_color.TT_trials(n,nn,:,n_sc,k,j)=Pseudo_pop_bootstrap.FEF.Stim_color.Trials(n,nn,training(h_stim_color),n_sc,k,j);
                    Pseudo_pop_bootstrap.FEF.Stim_color.Val_trials(n,nn,:,n_sc,k,j)=Pseudo_pop_bootstrap.FEF.Stim_color.Trials(n,nn,test(h_stim_color),n_sc,k,j);
                    
                    Pseudo_pop_bootstrap.PFC.Stim_color.TT_trials(n,nn,:,n_sc,k,j)=Pseudo_pop_bootstrap.PFC.Stim_color.Trials(n,nn,training(h_stim_color),n_sc,k,j);
                    Pseudo_pop_bootstrap.PFC.Stim_color.Val_trials(n,nn,:,n_sc,k,j)=Pseudo_pop_bootstrap.PFC.Stim_color.Trials(n,nn,test(h_stim_color),n_sc,k,j);
                end
            end
        end
    end
end


%% save

Pseudo_pop_bootstrap.LIP.Pseudo_pop=Pseudo_pop_LIP;
Pseudo_pop_bootstrap.FEF.Pseudo_pop=Pseudo_pop_FEF;
Pseudo_pop_bootstrap.PFC.Pseudo_pop=Pseudo_pop_PFC;


save(fullfile(data_path_clasifier,save_name),'Pseudo_pop_bootstrap','-v7.3')






























