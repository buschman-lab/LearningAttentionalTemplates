function Run_classifier_peak_belief_choice_restricted_FEF(fsroot, ROI, n_tt, loc, this_time)

par_param=0; %can't use par with cluster

task='Learning_Attentional_Templates';
subtask='exploreexploit/Restricted_FEF/Reset_RW_model';

event='response';

for i=1:23
    initial_window=-500;
    event_list{i}='response';
    window_start_list(i)=initial_window+(i-1)*50;
end

subsubtask_classifier=sprintf('Classifier_%s_%d_to_%s_%d',event,window_start_list(1),event,window_start_list(end));

dirstem = fullfile('Projects',task,'Analysis','Electrophy_analysis',subtask,subsubtask_classifier);
data_path_clasifier = fullfile(fsroot,dirstem);


boot_name='Pseudo_pop_choice_peak_template_classifier_bootstrap';
save_name=sprintf('Pseudo_pop_choice_peak_template_classifier_results_%s_%d_%d_%d',ROI,this_time,n_tt,loc);


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

N_neurons=size(this_ROI.Choice.TT_trials,2);
N_data=size(this_ROI.Choice.TT_trials,3); %
%we need an even number to have the same number of other beliefs (for 3
%categories)
if mod(N_data,2)~=0
    N_data=N_data-1;
end
N_choice=size(this_ROI.Choice.TT_trials,4);
N_belief=size(this_ROI.Choice.TT_trials,5);
% N_prog=size(this_ROI.Choice.TT_trials,6);

N_data_val=size(this_ROI.Choice.Val_trials,3)*N_choice;

Classification_belief_label=NaN(N_data_val*N_belief);
Classification_choice_label=NaN(N_data_val*N_belief);

Classification_belief_correct=NaN(N_data_val*N_belief,1);
Classification_choice_correct=NaN(N_data_val*N_belief,1);

Classification_choice_each_belief_label=NaN(N_data_val*N_belief,N_belief);
Classification_choice_each_belief_correct=NaN(N_data_val,N_belief,N_belief);

%% Belief classifier
for n_c=1:N_belief
    for n=1:N_neurons
        for t=1:N_data
            for n_sc=1:N_choice
                data_this_belief(n,t+(n_sc-1)*N_data)=this_ROI.Pseudo_pop(this_ROI.Neurons(n_tt,n),this_ROI.Choice.TT_trials(n_tt,n,t,n_sc,n_c,loc)).Classifier_FR(this_time);
            end
        end
    end
    %   subsample to have the same number of belief / non belief
    other_belief_ind=randsample(N_data,N_data/(N_belief-1),'false');
    increment=0;
    for l=1:N_belief
        if l~=n_c
            for n=1:N_neurons
                for t=1:N_data/(N_belief-1)
                    for n_sc=1:N_choice
                        data_other_belief(n,t+increment*N_data/(N_belief-1)+(n_sc-1)*N_data)=this_ROI.Pseudo_pop(this_ROI.Neurons(n_tt,n),this_ROI.Choice.TT_trials(n_tt,n,other_belief_ind(t),n_sc,l,loc)).Classifier_FR(this_time);
                    end
                end
            end
            increment=increment+1;
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
        W(:,n_c)=mdl.Beta;
        Intercept(n_c)=mdl.Bias;
        
        SVMModel(n_c).mdl=mdl;
        
    else
        
        W(:,n_c)=zeros(N_neurons,1);
        Intercept(n_c)=svmmod.Bias;
        
        SVMModel(n_c).mdl=svmmod;
        
    end
    
    clear c cv svmmod *_this_belief *_classifier labels probas scores mdl *_other_belief*
end


%% Choice classifier

for n=1:N_neurons
    for t=1:N_data
        for n_c=1:N_belief
            data_this_choice(n,t+(n_c-1)*N_data)=this_ROI.Pseudo_pop(this_ROI.Neurons(n_tt,n),this_ROI.Choice.TT_trials(n_tt,n,t,1,n_c,loc)).Classifier_FR(this_time); %%error Val and TT are inverted
        end
    end
end
for n=1:N_neurons
    for t=1:N_data/(N_choice-1)
        for n_c=1:N_belief
            data_other_choice(n,t+(n_c-1)*N_data)=this_ROI.Pseudo_pop(this_ROI.Neurons(n_tt,n),this_ROI.Choice.TT_trials(n_tt,n,t,2,n_c,loc)).Classifier_FR(this_time);
        end
    end
end

data_classifier_choice=[data_this_choice,data_other_choice]';

N_tt=size(data_classifier_choice,1)/2;
groups_classifier_choice=ones(2*N_tt,1);
groups_classifier_choice(N_tt+1:end)=0;

c = cvpartition(size(groups_classifier_choice,1),'KFold',10);

%optimize classfier

opts = struct('Optimizer','bayesopt','ShowPlots',false,'CVPartition',c,...
    'AcquisitionFunctionName','expected-improvement-plus','Verbose',0,'UseParallel',par_param);

svmmod = fitcsvm(data_classifier_choice,groups_classifier_choice,'KernelFunction','linear',...
    'OptimizeHyperparameters','auto','HyperparameterOptimizationOptions',opts);

if ~isempty(svmmod.Beta)
    
    mdl=fitPosterior(svmmod,data_classifier_choice,groups_classifier_choice);
    %
    W_choice(:)=mdl.Beta;
    Intercept_choice(:)=mdl.Bias;
    
    SVMModel_choice(:).mdl=mdl;
    
else
    
    W_choice(:)=zeros(N_neurons,1);
    Intercept_choice(:)=svmmod.Bias;
    
    SVMModel_choice(:).mdl=svmmod;
    
end

clear c cv svmmod *_this_choice *_classifier mdl *_other_choice*



%% Now decode the belief then the choice

%now apply the classfier to the validation data
%for p=1:N_prog
for l=1:N_belief
    for q=1:N_choice
        
        for t=1:size(this_ROI.Choice.Val_trials,3)
            for n=1:N_neurons
                data_val_this_belief(n,t+(q-1)*size(this_ROI.Choice.Val_trials,3)+(l-1)*N_data_val)=this_ROI.Pseudo_pop(this_ROI.Neurons(n_tt,n),this_ROI.Choice.Val_trials(n_tt,n,t,q,l,loc)).Classifier_FR(this_time);
            end
            groups_val_this_belief(t+(q-1)*size(this_ROI.Choice.Val_trials,3)+(l-1)*N_data_val)=l;
            groups_val_this_choice(t+(q-1)*size(this_ROI.Choice.Val_trials,3)+(l-1)*N_data_val)=q;
        end
    end
end

for n_c=1:N_belief %change for look of N_belief is not the same as N_choice
    [~, probas,~] = predict(SVMModel(n_c).mdl,data_val_this_belief');
    Classifier_belief_proba(:,n_c)=probas(:,2);
    clear probas
end
[class, ~,~] = predict(SVMModel_choice.mdl,data_val_this_belief');
class(class==0)=2;
Classification_choice_label=class;

clear class

%by belief classifier
for i=1:size(Classifier_belief_proba,1)
    [~, Classification_belief_label(i)]=max(Classifier_belief_proba(i,:));
    Classification_belief_correct(i)=(Classification_belief_label(i)==groups_val_this_belief(i));
end

%by choice classifier
for i=1:size(Classification_choice_label,1)
    Classification_choice_correct(i)=(Classification_choice_label(i)==groups_val_this_choice(i));
end
clear *_this_belief *_this_choice Classifier_*_proba
%end


%% Now do a choice decoder for each template
for n_c=1:N_belief
    for n=1:N_neurons
        for t=1:N_data
            data_this_choice(n,t)=this_ROI.Pseudo_pop(this_ROI.Neurons(n_tt,n),this_ROI.Choice.TT_trials(n_tt,n,t,1,n_c,loc)).Classifier_FR(this_time);
        end
    end
    for n=1:N_neurons
        for t=1:N_data/(N_choice-1)
            data_other_choice(n,t)=this_ROI.Pseudo_pop(this_ROI.Neurons(n_tt,n),this_ROI.Choice.TT_trials(n_tt,n,t,2,n_c,loc)).Classifier_FR(this_time);
        end
    end
    
    data_classifier_choice=[data_this_choice,data_other_choice]';
    
    N_tt=size(data_classifier_choice,1)/2;
    groups_classifier_choice=ones(2*N_tt,1);
    groups_classifier_choice(N_tt+1:end)=0;
    
    c = cvpartition(size(groups_classifier_choice,1),'KFold',10);
    
    %optimize classfier
    
    opts = struct('Optimizer','bayesopt','ShowPlots',false,'CVPartition',c,...
        'AcquisitionFunctionName','expected-improvement-plus','Verbose',0,'UseParallel',par_param);
    
    svmmod = fitcsvm(data_classifier_choice,groups_classifier_choice,'KernelFunction','linear',...
        'OptimizeHyperparameters','auto','HyperparameterOptimizationOptions',opts);
    
    if ~isempty(svmmod.Beta)
        
        mdl=fitPosterior(svmmod,data_classifier_choice,groups_classifier_choice);
        %
        W_choice_each_belief(:,n_c)=mdl.Beta;
        Intercept_choice_each_belief(n_c)=mdl.Bias;
        
        SVMModel_choice_each_belief(n_c).mdl=mdl;
        
    else
        
        W_choice_each_belief(:,n_c)=zeros(N_neurons,1);
        Intercept_choice_each_belief(n_c)=svmmod.Bias;
        
        SVMModel_choice_each_belief(n_c).mdl=svmmod;
        
    end
    
    clear c cv svmmod *_this_choice *_classifier mdl *_other_choice*
end

%% Now decode the choice for each template

%now apply the classfier to the validation data
% for p=1:N_prog

for n_c=1:N_belief
    for q=1:N_choice
        for t=1:size(this_ROI.Choice.Val_trials,3)
            for n=1:N_neurons
                data_val_this_belief(n,t+(q-1)*size(this_ROI.Choice.Val_trials,3)+(n_c-1)*N_data_val)=this_ROI.Pseudo_pop(this_ROI.Neurons(n_tt,n),this_ROI.Choice.Val_trials(n_tt,n,t,q,n_c,loc)).Classifier_FR(this_time);
            end
            groups_val_this_belief(t+(q-1)*size(this_ROI.Choice.Val_trials,3)+(n_c-1)*N_data_val)=n_c;
            groups_val_this_choice(t+(q-1)*size(this_ROI.Choice.Val_trials,3)+(n_c-1)*N_data_val)=q;
        end
    end
end

for n_c=1:N_belief
    [class, ~,~] = predict(SVMModel_choice_each_belief(n_c).mdl,data_val_this_belief');
    class(class==0)=2;
    Classification_choice_each_belief_label(:,n_c)=class;
    clear class
end

%by choice classifier
for n_c=1:N_belief
    for l=1:N_belief
        Classification_choice_each_belief_correct(:,l,n_c)=(Classification_choice_each_belief_label(groups_val_this_belief==l,n_c)'==groups_val_this_choice(groups_val_this_belief==l));
    end
end
clear *_this_belief *_this_choice Classifier_*_proba

% end


%%

save(fullfile(data_path_clasifier,save_name),'Classification_choice_label','Classification_belief_label','Classification_belief_correct','Classification_choice_correct','W','W_choice','Intercept','Intercept_choice','SVMModel','SVMModel_choice','Classification_choice_each_belief_label','Classification_choice_each_belief_correct','W_choice_each_belief','Intercept_choice_each_belief','SVMModel_choice_each_belief')


end


