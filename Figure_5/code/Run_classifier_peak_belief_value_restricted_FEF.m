function Run_classifier_peak_belief_value_restricted_FEF(fsroot, ROI, n_tt, this_time)

par_param=0; %can't use par with cluster

event='target';

for i=1:27
    initial_window=-500;
    event_list{i}='target';
    window_start_list(i)=initial_window+(i-1)*50;
end

task='Learning_Attentional_Templates';
subtask='exploreexploit/Restricted_FEF/Reset_RW_model';

window_size=200;

subsubtask_classifier=sprintf('Classifier_%s_%d_to_%s_%d',event,window_start_list(1),event,window_start_list(end));

dirstem = fullfile('Projects',task,'Analysis','Electrophy_analysis',subtask,subsubtask_classifier);
data_path_clasifier = fullfile(fsroot,dirstem);

boot_name='Pseudo_pop_value_peak_template_classifier_bootstrap';
save_name=sprintf('Pseudo_pop_value_peak_template_classifier_results_%s_%d_%d_%d',ROI,this_time,n_tt);


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

N_neurons=size(this_ROI.value.TT_trials,2);
N_data=size(this_ROI.value.TT_trials,3); %
%we need an even number to have the same number of other beliefs (for 3
%categories)
if mod(N_data,2)~=0
    N_data=N_data-1;
end
N_value=size(this_ROI.value.TT_trials,4);
N_belief=size(this_ROI.value.TT_trials,5);

N_data_val=size(this_ROI.value.Val_trials,3)*N_value;

Classification_belief_label=NaN(N_data_val*N_belief,1);

Classification_belief_correct=NaN(N_data_val*N_belief,1);
Classification_value_correct=NaN(N_data_val*N_belief,1);

Classification_value_each_belief_label=NaN(N_data_val*N_belief,N_belief);
Classification_value_each_belief_correct=NaN(N_data_val,N_belief,N_belief);

%% Belief classifier
for n_c=1:N_belief
    for n=1:N_neurons
        for t=1:N_data
            for n_sc=1:N_value
                data_this_belief(n,t+(n_sc-1)*N_data)=this_ROI.Pseudo_pop(this_ROI.Neurons(n_tt,n),this_ROI.value.TT_trials(n_tt,n,t,n_sc,n_c)).Classifier_FR(this_time);
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
                    for n_sc=1:N_value
                        data_other_belief(n,t+increment*N_data/(N_belief-1)+(n_sc-1)*N_data)=this_ROI.Pseudo_pop(this_ROI.Neurons(n_tt,n),this_ROI.value.TT_trials(n_tt,n,other_belief_ind(t),n_sc,l)).Classifier_FR(this_time);
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


%% value classifier

for n=1:N_neurons
    for t=1:N_data
        for n_c=1:N_belief
            data_this_value(n,t+(n_c-1)*N_data)=this_ROI.Pseudo_pop(this_ROI.Neurons(n_tt,n),this_ROI.value.TT_trials(n_tt,n,t,1,n_c)).Classifier_FR(this_time); %%error Val and TT are inverted
        end
    end
end
for n=1:N_neurons
    for t=1:N_data/(N_value-1)
        for n_c=1:N_belief
            data_other_value(n,t+(n_c-1)*N_data)=this_ROI.Pseudo_pop(this_ROI.Neurons(n_tt,n),this_ROI.value.TT_trials(n_tt,n,t,2,n_c)).Classifier_FR(this_time);
        end
    end
end

data_classifier_value=[data_this_value,data_other_value]';

N_tt=size(data_classifier_value,1)/2;
groups_classifier_value=ones(2*N_tt,1);
groups_classifier_value(N_tt+1:end)=0;

c = cvpartition(size(groups_classifier_value,1),'KFold',10);

%optimize classfier

opts = struct('Optimizer','bayesopt','ShowPlots',false,'CVPartition',c,...
    'AcquisitionFunctionName','expected-improvement-plus','Verbose',0,'UseParallel',par_param);

svmmod = fitcsvm(data_classifier_value,groups_classifier_value,'KernelFunction','linear',...
    'OptimizeHyperparameters','auto','HyperparameterOptimizationOptions',opts);

if ~isempty(svmmod.Beta)
    
    mdl=fitPosterior(svmmod,data_classifier_value,groups_classifier_value);
    %
    W_value(:)=mdl.Beta;
    Intercept_value(:)=mdl.Bias;
    
    SVMModel_value(:).mdl=mdl;
    
else
    
    W_value(:)=zeros(N_neurons,1);
    Intercept_value(:)=svmmod.Bias;
    
    SVMModel_value(:).mdl=svmmod;
    
end

clear c cv svmmod *_this_value *_classifier mdl *_other_value*



%% Now decode the belief then the value

%now apply the classfier to the validation data
%for p=1:N_prog
for l=1:N_belief
    for q=1:N_value
        
        for t=1:size(this_ROI.value.Val_trials,3)
            for n=1:N_neurons
                data_val_this_belief(n,t+(q-1)*size(this_ROI.value.Val_trials,3)+(l-1)*N_data_val)=this_ROI.Pseudo_pop(this_ROI.Neurons(n_tt,n),this_ROI.value.Val_trials(n_tt,n,t,q,l)).Classifier_FR(this_time);
            end
            groups_val_this_belief(t+(q-1)*size(this_ROI.value.Val_trials,3)+(l-1)*N_data_val)=l;
            groups_val_this_value(t+(q-1)*size(this_ROI.value.Val_trials,3)+(l-1)*N_data_val)=q;
        end
    end
end

for n_c=1:N_belief %change for look of N_belief is not the same as N_value
    [~, probas,~] = predict(SVMModel(n_c).mdl,data_val_this_belief');
    Classifier_belief_proba(:,n_c)=probas(:,2);
    clear probas
end
[class, ~,~] = predict(SVMModel_value.mdl,data_val_this_belief');
class(class==0)=2;
Classification_value_label=class;

clear class

%by belief classifier
for i=1:size(Classifier_belief_proba,1)
    [~, Classification_belief_label(i)]=max(Classifier_belief_proba(i,:));
    Classification_belief_correct(i)=(Classification_belief_label(i)==groups_val_this_belief(i));
end

%by value classifier
for i=1:size(Classification_value_label,1)
    Classification_value_correct(i)=(Classification_value_label(i)==groups_val_this_value(i));
end
clear *_this_belief *_this_value Classifier_*_proba
%end


%% Now do a value decoder for each template
for n_c=1:N_belief
    for n=1:N_neurons
        for t=1:N_data
            data_this_value(n,t)=this_ROI.Pseudo_pop(this_ROI.Neurons(n_tt,n),this_ROI.value.TT_trials(n_tt,n,t,1,n_c)).Classifier_FR(this_time);
        end
    end
    for n=1:N_neurons
        for t=1:N_data/(N_value-1)
            data_other_value(n,t)=this_ROI.Pseudo_pop(this_ROI.Neurons(n_tt,n),this_ROI.value.TT_trials(n_tt,n,t,2,n_c)).Classifier_FR(this_time);
        end
    end
    
    data_classifier_value=[data_this_value,data_other_value]';
    
    N_tt=size(data_classifier_value,1)/2;
    groups_classifier_value=ones(2*N_tt,1);
    groups_classifier_value(N_tt+1:end)=0;
    
    c = cvpartition(size(groups_classifier_value,1),'KFold',10);
    
    %optimize classfier
    
    opts = struct('Optimizer','bayesopt','ShowPlots',false,'CVPartition',c,...
        'AcquisitionFunctionName','expected-improvement-plus','Verbose',0,'UseParallel',par_param);
    
    svmmod = fitcsvm(data_classifier_value,groups_classifier_value,'KernelFunction','linear',...
        'OptimizeHyperparameters','auto','HyperparameterOptimizationOptions',opts);
    
    if ~isempty(svmmod.Beta)
        
        mdl=fitPosterior(svmmod,data_classifier_value,groups_classifier_value);
        %
        W_value_each_belief(:,n_c)=mdl.Beta;
        Intercept_value_each_belief(n_c)=mdl.Bias;
        
        SVMModel_value_each_belief(n_c).mdl=mdl;
        
    else
        
        W_value_each_belief(:,n_c)=zeros(N_neurons,1);
        Intercept_value_each_belief(n_c)=svmmod.Bias;
        
        SVMModel_value_each_belief(n_c).mdl=svmmod;
        
    end
    
    clear c cv svmmod *_this_value *_classifier mdl *_other_value*
end

%% Now decode the value for each template

%now apply the classfier to the validation data

for n_c=1:N_belief
    for q=1:N_value
        for t=1:size(this_ROI.value.Val_trials,3)
            for n=1:N_neurons
                data_val_this_belief(n,t+(q-1)*size(this_ROI.value.Val_trials,3)+(n_c-1)*N_data_val)=this_ROI.Pseudo_pop(this_ROI.Neurons(n_tt,n),this_ROI.value.Val_trials(n_tt,n,t,q,n_c)).Classifier_FR(this_time);
            end
            groups_val_this_belief(t+(q-1)*size(this_ROI.value.Val_trials,3)+(n_c-1)*N_data_val)=n_c;
            groups_val_this_value(t+(q-1)*size(this_ROI.value.Val_trials,3)+(n_c-1)*N_data_val)=q;
        end
    end
end

for n_c=1:N_belief
    [class, ~,~] = predict(SVMModel_value_each_belief(n_c).mdl,data_val_this_belief');
    class(class==0)=2;
    Classification_value_each_belief_label(:,n_c)=class;
    clear class
end

%by value classifier
for n_c=1:N_belief
    for l=1:N_belief
        Classification_value_each_belief_correct(:,l,n_c)=(Classification_value_each_belief_label(groups_val_this_belief==l,n_c)'==groups_val_this_value(groups_val_this_belief==l));
    end
end
clear *_this_belief *_this_value Classifier_*_proba



%%

save(fullfile(data_path_clasifier,save_name),'Classification_value_label','Classification_belief_label','Classification_belief_correct','Classification_value_correct','W','W_value','Intercept','Intercept_value','SVMModel','SVMModel_value','Classification_value_each_belief_label','Classification_value_each_belief_correct','W_value_each_belief','Intercept_value_each_belief','SVMModel_value_each_belief')


end


