function Run_belief_classifier_combined_time_prog_NN_rest_FEF_rotated(fsroot, ROI, n_tt, this_combination)
% fsroot='/Volumes/buschman';
% ROI='LIP';
% n_tt= 1;
% this_combination=2;
% switch event
%     case 'target'
%         initial_window=-600;
%         N_time_points=20;
%     case 'reward_end'
%         initial_window=0;
%         N_time_points=7;
% end
%
% for i=1:N_time_points
%     window_start_list(i)=initial_window+(i-1)*50;
% end
%
% task='Learning_Attentional_Templates';
% subtask='exploreexploit/Reset_RW_model';
%
% % subsubtask_classifier=sprintf('Classifier_%s_%d_to_%s_%d',event_list{1},window_start_list(1),event_list{end},window_start_list(end));
% subsubtask_classifier=sprintf('Classifier_%s_%d_to_%s_%d',event,window_start_list(1),event,window_start_list(end));

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

par_param=0; %can't use par with cluster

dirstem = fullfile('Projects',task,'Analysis','Electrophy_analysis',subtask,subsubtask_classifier);
data_path_clasifier = fullfile(fsroot,dirstem);

% boot_name='Pseudo_pop_belief_prog_classifier_bootstrap_rotated';
boot_name='Pseudo_pop_peak_belief_prog_classifier_bootstrap_rotated_60'; %added peak
if this_combination==1
    class_name=sprintf('Pseudo_pop_peak_belief_prog_combined_time_results_rotated_%s_600ms_%d',ROI,n_tt);
    save_name=sprintf('Pseudo_pop_peak_belief_prog_combined_time_results_rotated_%s_600ms_NN_%d',ROI,n_tt);
elseif this_combination==2
    class_name=sprintf('Pseudo_pop_peak_belief_prog_combined_time_results_rotated_%s_900ms_%d',ROI,n_tt);
    save_name=sprintf('Pseudo_pop_peak_belief_prog_combined_time_results_rotated_%s_900ms_NN_%d',ROI,n_tt);
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

N_boot=1; %size(this_ROI.progression.TT_trials,1);
N_neurons=size(this_ROI.progression.TT_trials,2);
N_data=size(this_ROI.progression.TT_trials,3); %
N_prog=size(this_ROI.progression.TT_trials,4);
N_belief=size(this_ROI.progression.TT_trials,5);

N_data_val=size(this_ROI.progression.Val_trials,3)*3;
Classification_label=NaN(N_prog,N_data_val,N_boot);
Classification_this_label=NaN(N_prog,N_data_val,N_boot);
Classification_net_correct=NaN(N_prog,N_data_val,N_boot);


%get the data
for t=1:N_data
    for p=1:N_prog
        for k=1:N_belief
            for n=1:N_neurons
                if this_combination==1
                    data_tt(n,t+(p-1)*N_data+(k-1)*N_prog*N_data)=(this_ROI.Pseudo_pop(this_ROI.Neurons(n_tt,n),this_ROI.progression.TT_trials(n_tt,n,t,p,k)).Classifier_FR(8)+this_ROI.Pseudo_pop(this_ROI.Neurons(n_tt,n),this_ROI.progression.TT_trials(n_tt,n,t,p,k)).Classifier_FR(14))/2; %%error Val and TT are inverted
                elseif this_combination==2
                    data_tt(n,t+(p-1)*N_data+(k-1)*N_prog*N_data)=(this_ROI.Pseudo_pop(this_ROI.Neurons(n_tt,n),this_ROI.progression.TT_trials(n_tt,n,t,p,k)).Classifier_FR(8)+this_ROI.Pseudo_pop(this_ROI.Neurons(n_tt,n),this_ROI.progression.TT_trials(n_tt,n,t,p,k)).Classifier_FR(14)+this_ROI.Pseudo_pop(this_ROI.Neurons(n_tt,n),this_ROI.progression.TT_trials(n_tt,n,t,p,k)).Classifier_FR(20))/3; %%error Val and TT are inverted
                end
                groups_tt(t+(p-1)*N_data+(k-1)*N_prog*N_data)=k;
            end
        end
    end
end
%apply the classifiers
for k=1:N_belief
    [~, probas,~] = predict(SVMModel(k,1).mdl,data_tt');
    Classifier_proba_tt(:,k)=probas(:,2);
    clear probas
end
mdl_net=fitcnet(Classifier_proba_tt,groups_tt);

%test
Test_net = predict(mdl_net,Classifier_proba_tt);
accuracy_tt=1-loss(mdl_net,Classifier_proba_tt,groups_tt);
% figure
% confusionchart(Test_net,groups_tt)

%%

for q=1:N_prog
    %now apply the classfier to the validation data
    for l=1:N_belief
        for t=1:size(this_ROI.progression.Val_trials,3)
            for n=1:N_neurons
                if this_combination==1
                    data_val_this_belief(n,t+(l-1)*size(this_ROI.progression.Val_trials,3))=(this_ROI.Pseudo_pop(this_ROI.Neurons(n_tt,n),this_ROI.progression.Val_trials(n_tt,n,t,q,l)).Classifier_FR(8)+this_ROI.Pseudo_pop(this_ROI.Neurons(n_tt,n),this_ROI.progression.Val_trials(n_tt,n,t,q,l)).Classifier_FR(14))/2;  %%error Val and TT are inverted
                elseif this_combination==2
                    data_val_this_belief(n,t+(l-1)*size(this_ROI.progression.Val_trials,3))=(this_ROI.Pseudo_pop(this_ROI.Neurons(n_tt,n),this_ROI.progression.Val_trials(n_tt,n,t,q,l)).Classifier_FR(8)+this_ROI.Pseudo_pop(this_ROI.Neurons(n_tt,n),this_ROI.progression.Val_trials(n_tt,n,t,q,l)).Classifier_FR(14)+this_ROI.Pseudo_pop(this_ROI.Neurons(n_tt,n),this_ROI.progression.Val_trials(n_tt,n,t,q,l)).Classifier_FR(20))/3;  %%error Val and TT are inverted
                end
            end
            groups_val_this_belief(t+(l-1)*size(this_ROI.progression.Val_trials,3))=l;
        end
    end
    for k=1:N_belief
        [~, probas,~] = predict(SVMModel(k,1).mdl,data_val_this_belief');
        Classifier_proba_val(:,k)=probas(:,2);
        clear probas
    end
    
    Classification_net_label(q,:)=predict(mdl_net,Classifier_proba_val);
    for i=1:length(groups_val_this_belief)
        Classification_net_correct(q,i)=(Classification_net_label(q,i)==groups_val_this_belief(i));
    end
    accuracy_val(q)=1-loss(mdl_net,Classifier_proba_val,groups_val_this_belief);
    
%     figure
%     confusionchart(Classification_net_label(q,:),groups_val_this_belief)
    clear Classifier_proba_val *_this_belief
end


%% apply to all time points

for tp=1:length(window_start_list)
    for q=1:N_prog
        %now apply the classfier to the validation data
        for l=1:N_belief
            for t=1:size(this_ROI.progression.Val_trials,3)
                for n=1:N_neurons
                    if this_combination==1
                        data_val_this_belief(n,t+(l-1)*size(this_ROI.progression.Val_trials,3))=this_ROI.Pseudo_pop(this_ROI.Neurons(n_tt,n),this_ROI.progression.Val_trials(n_tt,n,t,q,l)).Classifier_FR(tp);  %%error Val and TT are inverted
                    elseif this_combination==2
                        data_val_this_belief(n,t+(l-1)*size(this_ROI.progression.Val_trials,3))=this_ROI.Pseudo_pop(this_ROI.Neurons(n_tt,n),this_ROI.progression.Val_trials(n_tt,n,t,q,l)).Classifier_FR(tp);  %%error Val and TT are inverted
                    end
                end
                groups_val_this_belief(t+(l-1)*size(this_ROI.progression.Val_trials,3))=l;
            end
        end
        for k=1:N_belief
            [~, probas,~] = predict(SVMModel(k,1).mdl,data_val_this_belief');
            Classifier_proba_val(:,k)=probas(:,2);
            clear probas
        end
        
        Classification_net_label_across_time(q,:,tp)=predict(mdl_net,Classifier_proba_val);
        for i=1:length(groups_val_this_belief)
            Classification_net_correct_across_time(q,i,tp)=(Classification_net_label_across_time(q,i,tp)==groups_val_this_belief(i));
        end
        accuracy_val_across_time(q,tp)=1-loss(mdl_net,Classifier_proba_val,groups_val_this_belief);
        
        %     figure
        %     confusionchart(Classification_net_label(q,:),groups_val_this_belief)
        clear Classifier_proba_val *_this_belief
    end
end


%%

save(fullfile(data_path_clasifier,save_name),'Classification_net_label*','Classification_net_correct*','mdl_net','accuracy*','Test_net','Classifier_proba_tt')

end


