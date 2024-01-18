function Run_belief_classifier_combined_time_each_prog_NN_restricted_FEF(fsroot, ROI, n_tt, this_combination)

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

subsubtask_classifier=sprintf('Classifier_%s_%d_to_%s_%d',event_list{1},window_start_list(1),event_list{end},window_start_list(end));

par_param=0; %can't use par with cluster

dirstem = fullfile('Projects',task,'Analysis','Electrophy_analysis',subtask,subsubtask_classifier);
data_path_clasifier = fullfile(fsroot,dirstem);

boot_name='Pseudo_pop_peak_belief_prog_classifier_bootstrap_more_trials'; 
if this_combination==1
    class_name=sprintf('Pseudo_pop_peak_belief_each_prog_combined_time_results_more_trials_%s_600ms_%d',ROI,n_tt);
    save_name=sprintf('Pseudo_pop_peak_belief_each_prog_combined_time_results_more_trials_%s_600ms_NN_%d',ROI,n_tt);
elseif this_combination==2
    class_name=sprintf('Pseudo_pop_peak_belief_each_prog_combined_time_results_more_trials_%s_900ms_%d',ROI,n_tt);
    save_name=sprintf('Pseudo_pop_peak_belief_each_prog_combined_time_results_more_trials_%s_900ms_NN_%d',ROI,n_tt);
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

N_boot=1; 
N_neurons=size(this_ROI.progression.TT_trials,2);
N_data=size(this_ROI.progression.TT_trials,3); %
N_prog=size(this_ROI.progression.TT_trials,4);
N_belief=size(this_ROI.progression.TT_trials,5);

N_data_val=size(this_ROI.progression.Val_trials,3)*3;
Classification_net_correct=NaN(N_prog,N_data_val,N_boot);


%get the data
for p=1:N_prog
    for t=1:N_data
        for k=1:N_belief
            for n=1:N_neurons
                if this_combination==1
                    data_tt(n,t+(k-1)*N_data)=(this_ROI.Pseudo_pop(this_ROI.Neurons(n_tt,n),this_ROI.progression.TT_trials(n_tt,n,t,p,k)).Classifier_FR(8)+this_ROI.Pseudo_pop(this_ROI.Neurons(n_tt,n),this_ROI.progression.TT_trials(n_tt,n,t,p,k)).Classifier_FR(14))/2;
                elseif this_combination==2
                    data_tt(n,t+(k-1)*N_data)=(this_ROI.Pseudo_pop(this_ROI.Neurons(n_tt,n),this_ROI.progression.TT_trials(n_tt,n,t,p,k)).Classifier_FR(8)+this_ROI.Pseudo_pop(this_ROI.Neurons(n_tt,n),this_ROI.progression.TT_trials(n_tt,n,t,p,k)).Classifier_FR(14)+this_ROI.Pseudo_pop(this_ROI.Neurons(n_tt,n),this_ROI.progression.TT_trials(n_tt,n,t,p,k)).Classifier_FR(20))/3; %%error Val and TT are inverted
                end
                groups_tt(t+(k-1)*N_data)=k;
            end
        end
    end
    %apply the classifiers
    for k=1:N_belief
        [~, probas,~] = predict(SVMModel(k,p).mdl,data_tt');
        Classifier_proba_tt(:,k)=probas(:,2);
        clear probas
    end
    mdl_net(p).mdl=fitcnet(Classifier_proba_tt,groups_tt);
    
    %test ontriaing set
    Test_net(:,p) = predict(mdl_net(p).mdl,Classifier_proba_tt);
    accuracy_tt(p)=1-loss(mdl_net(p).mdl,Classifier_proba_tt,groups_tt);
    % figure
    % confusionchart(Test_net,groups_tt)
end

%%

for p=1:N_prog
    %now apply the classfier to the validation data
    for l=1:N_belief
        for t=1:size(this_ROI.progression.Val_trials,3)
            for n=1:N_neurons
                if this_combination==1
                    data_val_this_belief(n,t+(l-1)*size(this_ROI.progression.Val_trials,3))=(this_ROI.Pseudo_pop(this_ROI.Neurons(n_tt,n),this_ROI.progression.Val_trials(n_tt,n,t,p,l)).Classifier_FR(8)+this_ROI.Pseudo_pop(this_ROI.Neurons(n_tt,n),this_ROI.progression.Val_trials(n_tt,n,t,p,l)).Classifier_FR(14))/2;  %%error Val and TT are inverted
                elseif this_combination==2
                    data_val_this_belief(n,t+(l-1)*size(this_ROI.progression.Val_trials,3))=(this_ROI.Pseudo_pop(this_ROI.Neurons(n_tt,n),this_ROI.progression.Val_trials(n_tt,n,t,p,l)).Classifier_FR(8)+this_ROI.Pseudo_pop(this_ROI.Neurons(n_tt,n),this_ROI.progression.Val_trials(n_tt,n,t,p,l)).Classifier_FR(14)+this_ROI.Pseudo_pop(this_ROI.Neurons(n_tt,n),this_ROI.progression.Val_trials(n_tt,n,t,p,l)).Classifier_FR(20))/3;  %%error Val and TT are inverted
                end
            end
            groups_val_this_belief(t+(l-1)*size(this_ROI.progression.Val_trials,3))=l;
        end
    end
    for k=1:N_belief
        [~, probas,~] = predict(SVMModel(k,p).mdl,data_val_this_belief');
        Classifier_proba_val(:,k)=probas(:,2);
        clear probas
    end
    
    Classification_net_label(p,:)=predict(mdl_net(p).mdl,Classifier_proba_val);
    for i=1:length(groups_val_this_belief)
        Classification_net_correct(p,i)=(Classification_net_label(p,i)==groups_val_this_belief(i));
    end
    accuracy_val(p)=1-loss(mdl_net(p).mdl,Classifier_proba_val,groups_val_this_belief);
    
    %     figure
    %     confusionchart(Classification_net_label(q,:),groups_val_this_belief)
    clear Classifier_proba_val *_this_belief
end

%%

save(fullfile(data_path_clasifier,save_name),'Classification_net_label','Classification_net_correct','mdl_net','accuracy*','Test_net','Classifier_proba_tt')

end


