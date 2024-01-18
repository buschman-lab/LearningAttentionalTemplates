function get_projection_classifier_combined_time_prog_restricted_FEF(fsroot, ROI, n_tt, this_combination)

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

subsubtask_classifier=sprintf('Classifier_%s_%d_to_%s_%d',event_list{1},window_start_list(1),event_list{end},window_start_list(end));

dirstem = fullfile('Projects',task,'Analysis','Electrophy_analysis',subtask,subsubtask_classifier);
data_path_clasifier = fullfile(fsroot,dirstem);

% boot_name='Pseudo_pop_belief_prog_classifier_bootstrap_more_trials';
boot_name='Pseudo_pop_peak_belief_prog_classifier_bootstrap_more_trials'; %added peak
if this_combination==1
    class_name=sprintf('Pseudo_pop_peak_belief_prog_combined_time_results_more_trials_%s_600ms_%d',ROI,n_tt);
    save_name=sprintf('Pseudo_pop_peak_belief_prog_combined_time_results_more_trials_%s_600ms_%d_proj',ROI,n_tt);
elseif this_combination==2
    class_name=sprintf('Pseudo_pop_peak_belief_prog_combined_time_results_more_trials_%s_900ms_%d',ROI,n_tt);
    save_name=sprintf('Pseudo_pop_peak_belief_prog_combined_time_results_more_trials_%s_900ms_%d_proj',ROI,n_tt);
end

load(fullfile(data_path_clasifier,boot_name),'Pseudo_pop_bootstrap') %data
load(fullfile(data_path_clasifier,class_name),'SVMModel') %classifiers

switch ROI
    case 'LIP'
        this_ROI=Pseudo_pop_bootstrap.LIP;
        clear Pseudo_pop_bootstrap
    case 'FEF'
        this_ROI=Pseudo_pop_bootstrap.FEF;
        clear Pseudo_pop_bootstrap
    case 'PFC'
        this_ROI=Pseudo_pop_bootstrap.PFC;
        clear Pseudo_pop_bootstrap
end

N_neurons=size(this_ROI.progression.TT_trials,2);
N_prog=size(this_ROI.progression.TT_trials,4);
N_belief=size(this_ROI.progression.TT_trials,5);


%%

%now apply the classfier to the validation data
for l=1:N_belief
    for t=1:size(this_ROI.progression.Val_trials,3)
        for n=1:N_neurons
            for q=1:N_prog
                if this_combination==1
                    data_val_this_belief(n,t+(l-1)*size(this_ROI.progression.Val_trials,3)+(q-1)*size(this_ROI.progression.Val_trials,3)*N_belief)=(this_ROI.Pseudo_pop(this_ROI.Neurons(n_tt,n),this_ROI.progression.Val_trials(n_tt,n,t,q,l)).Classifier_FR(8)+this_ROI.Pseudo_pop(this_ROI.Neurons(n_tt,n),this_ROI.progression.Val_trials(n_tt,n,t,q,l)).Classifier_FR(14))/2;  %%error Val and TT are inverted
                elseif this_combination==2
                    data_val_this_belief(n,t+(l-1)*size(this_ROI.progression.Val_trials,3)+(q-1)*size(this_ROI.progression.Val_trials,3)*N_belief)=(this_ROI.Pseudo_pop(this_ROI.Neurons(n_tt,n),this_ROI.progression.Val_trials(n_tt,n,t,q,l)).Classifier_FR(8)+this_ROI.Pseudo_pop(this_ROI.Neurons(n_tt,n),this_ROI.progression.Val_trials(n_tt,n,t,q,l)).Classifier_FR(14)+this_ROI.Pseudo_pop(this_ROI.Neurons(n_tt,n),this_ROI.progression.Val_trials(n_tt,n,t,q,l)).Classifier_FR(20))/3;  %%error Val and TT are inverted
                end
            end
            groups_val_this_belief(t+(l-1)*size(this_ROI.progression.Val_trials,3)+(q-1)*size(this_ROI.progression.Val_trials,3)*N_belief)=l;
        end
    end
end
for k=1:N_belief
    [~, probas,~] = predict(SVMModel(k,1).mdl,data_val_this_belief');
    for l=1:N_belief
        Classifier_proba_val(:,l,k)=probas(groups_val_this_belief==l,2);
    end
    clear probas
end

%%

save(fullfile(data_path_clasifier,save_name),'Classifier_proba_val')

end