function Run_belief_classifier_combined_time_prog_restricted_FEF_4bins(fsroot, ROI, n_tt, this_combination)


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

% boot_name='Pseudo_pop_belief_prog_classifier_bootstrap_more_trials';
% boot_name='Pseudo_pop_peak_belief_prog_classifier_bootstrap_more_trials'; %added peak
boot_name='Pseudo_pop_peak_belief_prog_classifier_bootstrap_4bins';
if this_combination==1
    save_name=sprintf('Pseudo_pop_peak_belief_prog_combined_time_results_4bins_%s_600ms_%d',ROI,n_tt);
elseif this_combination==2
    save_name=sprintf('Pseudo_pop_peak_belief_prog_combined_time_results_4bins_%s_900ms_%d',ROI,n_tt);
end
    
load(fullfile(data_path_clasifier,boot_name),'Pseudo_pop_bootstrap')

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

N_times=length(this_ROI.Pseudo_pop(1,1).Classifier_FR);
N_boot=1; %size(this_ROI.progression.TT_trials,1);
N_neurons=size(this_ROI.progression.TT_trials,2);
N_data=size(this_ROI.progression.TT_trials,3); %
%we need a number of trials that can be divided by 3 (to subsample equally
%for the 3 other belief bins)
if mod(N_data,3)==1
    N_data=N_data-1;
elseif mod(N_data,3)==2
    N_data=N_data-2;
end
N_prog=size(this_ROI.progression.TT_trials,4);
N_belief=size(this_ROI.progression.TT_trials,5);

N_data_val=size(this_ROI.progression.Val_trials,3)*3;
Classification_label=NaN(N_prog,N_data_val,N_boot);
Classification_this_label=NaN(N_prog,N_data_val,N_boot);
Classification_correct=NaN(N_prog,N_data_val,N_boot);

%%
for k=1:N_belief
    for n=1:N_neurons
        for t=1:N_data
            for p=1:N_prog
                if this_combination==1
                    data_this_belief(n,t+(p-1)*N_data)=(this_ROI.Pseudo_pop(this_ROI.Neurons(n_tt,n),this_ROI.progression.TT_trials(n_tt,n,t,p,k)).Classifier_FR(8)+this_ROI.Pseudo_pop(this_ROI.Neurons(n_tt,n),this_ROI.progression.TT_trials(n_tt,n,t,p,k)).Classifier_FR(14))/2; %%error Val and TT are inverted
                elseif this_combination==2
                    data_this_belief(n,t+(p-1)*N_data)=(this_ROI.Pseudo_pop(this_ROI.Neurons(n_tt,n),this_ROI.progression.TT_trials(n_tt,n,t,p,k)).Classifier_FR(8)+this_ROI.Pseudo_pop(this_ROI.Neurons(n_tt,n),this_ROI.progression.TT_trials(n_tt,n,t,p,k)).Classifier_FR(14)+this_ROI.Pseudo_pop(this_ROI.Neurons(n_tt,n),this_ROI.progression.TT_trials(n_tt,n,t,p,k)).Classifier_FR(20))/3; %%error Val and TT are inverted
                end
                
            end
        end
    end
    %   subsample to have the same number of chosen/unchosen
    other_belief_ind=randsample(N_data,N_data/(N_belief-1),'false');
    increment=0;
    for l=1:N_belief
        if l~=k
            for n=1:N_neurons
                for t=1:N_data/(N_belief-1)
                    for p=1:N_prog
                        if this_combination==1
                            data_other_belief(n,t+increment*round(N_data/(N_belief-1))+(p-1)*N_data)=(this_ROI.Pseudo_pop(this_ROI.Neurons(n_tt,n),this_ROI.progression.TT_trials(n_tt,n,other_belief_ind(t),p,l)).Classifier_FR(8)+this_ROI.Pseudo_pop(this_ROI.Neurons(n_tt,n),this_ROI.progression.TT_trials(n_tt,n,other_belief_ind(t),p,l)).Classifier_FR(14))/2;
                        elseif this_combination==2
                            data_other_belief(n,t+increment*round(N_data/(N_belief-1))+(p-1)*N_data)=(this_ROI.Pseudo_pop(this_ROI.Neurons(n_tt,n),this_ROI.progression.TT_trials(n_tt,n,other_belief_ind(t),p,l)).Classifier_FR(8)+this_ROI.Pseudo_pop(this_ROI.Neurons(n_tt,n),this_ROI.progression.TT_trials(n_tt,n,other_belief_ind(t),p,l)).Classifier_FR(14)+this_ROI.Pseudo_pop(this_ROI.Neurons(n_tt,n),this_ROI.progression.TT_trials(n_tt,n,other_belief_ind(t),p,l)).Classifier_FR(20))/3;
                        end
                    end
                end
                increment=1;
            end
        end
    end
    
    data_classifier=[data_this_belief,data_other_belief]';
    
    N_tt=size(data_classifier,1)/2;
    groups_classifier=ones(2*N_tt,1);
    groups_classifier(N_tt+1:end)=0;
    
    c = cvpartition(size(groups_classifier,1),'KFold',10);
    
    %optimize classfier
    
    opts = struct('Optimizer','bayesopt','ShowPlots',false,'CVPartition',c,...
        'AcquisitionFunctionName','expected-improvement-plus','Verbose',0,'UseParallel',par_param);
    
    svmmod = fitcsvm(data_classifier,groups_classifier,'KernelFunction','linear',...
        'OptimizeHyperparameters','auto','HyperparameterOptimizationOptions',opts);
    
    if ~isempty(svmmod.Beta)
        
        mdl=fitPosterior(svmmod,data_classifier,groups_classifier);
        %
        W(:,k,1)=mdl.Beta;
        Intercept(k,1)=mdl.Bias;
        
        SVMModel(k,1).mdl=mdl;
        
    else
        
        W(:,k,1)=zeros(N_neurons,1);
        Intercept(k,1)=svmmod.Bias;
        
        SVMModel(k,1).mdl=svmmod;
        
    end
    
    clear c cv *_this_belief *_classifier labels probas scores mdl *_other_belief*
end

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
        Classifier_proba(:,k)=probas(:,2);
        clear probas
    end
    
    for i=1:length(Classifier_proba)
        [Classification_this_label(q,i), Classification_label(q,i)]=max(Classifier_proba(i,:));
        Classification_correct(q,i)=(Classification_label(q,i)==groups_val_this_belief(i));
    end
    clear Classifier_proba *_this_belief
end

%%

save(fullfile(data_path_clasifier,save_name),'Classification_label','Classification_correct','W','Intercept','SVMModel')


end


