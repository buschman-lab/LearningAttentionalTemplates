function Run_choice_classifier_across_time_each_context_restricted_FEF(fsroot, ROI, n_tt, this_time)

%%choice context to choices offered
% 1 = 1 2 3
% 2 = 1 2 4
% 3 = 1 3 4
% 4 = 2 3 4

% %Test
% fsroot='/Volumes/buschman';
% ROI='LIP';
% n_tt=1;
% this_time=23;


for i=1:23
    initial_window=-500;
    event_list{i}='response';
    window_start_list(i)=initial_window+(i-1)*50;
end

task='Learning_Attentional_Templates';
subtask='exploreexploit/Restricted_FEF/Reset_RW_model';

window_size=200;
subsubtask_classifier=sprintf('Classifier_%s_%d_to_%s_%d',event_list{1},window_start_list(1),event_list{end},window_start_list(end));

par_param=0; %can't use par with cluster

dirstem = fullfile('Projects',task,'Analysis','Electrophy_analysis',subtask,subsubtask_classifier);
data_path_clasifier = fullfile(fsroot,dirstem);

boot_name='Pseudo_pop_choice_in_context_classifier_bootstrap';
save_name=sprintf('Pseudo_pop_choice_in_context_results_%s_%d_%d',ROI,this_time,n_tt);

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

%%

N_times=length(this_ROI.Pseudo_pop(1,1).Classifier_FR);
N_neurons=size(this_ROI.Choice_context.TT_trials,2);
N_data=size(this_ROI.Choice_context.TT_trials,3); %
%we need an even number to have the same number of other choices (for 3
%categories)
if mod(N_data,2)~=0
    N_data=N_data-1;
end
N_choices=size(this_ROI.Choice_context.TT_trials,4);
N_context=size(this_ROI.Choice_context.TT_trials,5);

N_data_val=size(this_ROI.Choice_context.Val_trials,3)*N_choices;
Classification_label=NaN(N_times,N_data_val,N_context);
Classification_this_label=NaN(N_times,N_data_val,N_context);
Classification_correct=NaN(N_times,N_data_val,N_context);

for p=1:N_context
    for k=1:N_choices
        for n=1:N_neurons
            for t=1:N_data
                data_this_choice(n,t)=this_ROI.Pseudo_pop(this_ROI.Neurons(n_tt,n),this_ROI.Choice_context.TT_trials(n_tt,n,t,k,p)).Classifier_FR(this_time); %%error Val and TT are inverted
            end
        end
        %   subsample to have the same number of chosen/unchosen
        other_choice_ind=randsample(N_data,N_data/(N_choices-1),'false');
        increment=0;
        for l=1:N_choices
            if l~=k
                for n=1:N_neurons
                    for t=1:N_data/(N_choices-1)
                            data_other_choice(n,t+increment*N_data/(N_choices-1))=this_ROI.Pseudo_pop(this_ROI.Neurons(n_tt,n),this_ROI.Choice_context.TT_trials(n_tt,n,other_choice_ind(t),l,p)).Classifier_FR(this_time);
                    end
                end
                increment=1;
            end
        end
        
        data_classifier=[data_this_choice,data_other_choice]';
        
        N_tt=size(data_classifier,1)/2;
        groups_classifier=ones(2*N_tt,1);
        groups_classifier(N_tt+1:end)=0;
        
        c = cvpartition(size(groups_classifier,1),'KFold',10);
        
        %% optimize classfier
        
        opts = struct('Optimizer','bayesopt','ShowPlots',false,'CVPartition',c,...
            'AcquisitionFunctionName','expected-improvement-plus','Verbose',0,'UseParallel',par_param);
        
        svmmod = fitcsvm(data_classifier,groups_classifier,'KernelFunction','linear',...
            'OptimizeHyperparameters','auto','HyperparameterOptimizationOptions',opts);
        
        if ~isempty(svmmod.Beta)
            
            mdl=fitPosterior(svmmod,data_classifier,groups_classifier);
            %
            W(:,k,p)=mdl.Beta;
            Intercept(k,p)=mdl.Bias;
            
            SVMModel(k,p).mdl=mdl;
            
        else
            
            W(:,k,p)=zeros(N_neurons,1);
            Intercept(k,p)=svmmod.Bias;
            
            SVMModel(k,p).mdl=svmmod;
            
        end
        
        clear c cv *_this_choice *_classifier labels probas scores mdl *_other_choice*
    end
end

for tp=1:N_times
    for p=1:N_context
        %now apply the classfier to the validation data
        for l=1:N_choices
            for t=1:size(this_ROI.Choice_context.Val_trials,3)
                for n=1:N_neurons
                    data_val_this_choice(n,t+(l-1)*size(this_ROI.Choice_context.Val_trials,3))=this_ROI.Pseudo_pop(this_ROI.Neurons(n_tt,n),this_ROI.Choice_context.Val_trials(n_tt,n,t,l,p)).Classifier_FR(tp);  %%error Val and TT are inverted
                end
                groups_val_this_choice(t+(l-1)*size(this_ROI.Choice_context.Val_trials,3))=l;
            end
        end
        for k=1:N_choices
            [~, probas,~] = predict(SVMModel(k,p).mdl,data_val_this_choice');
            Classifier_proba(:,k,p)=probas(:,2);
            clear probas
        end
        for i=1:length(Classifier_proba)
            [Classification_this_label(tp,i,p), Classification_label(tp,i,p)]=max(Classifier_proba(i,:,p));
            Classification_correct(tp,i,p)=(Classification_label(tp,i,p)==groups_val_this_choice(i));
        end
        clear Classifier_proba *_this_choice
    end
end

save(fullfile(data_path_clasifier,save_name),'Classification_label','Classification_correct','W','Intercept','SVMModel')

end


