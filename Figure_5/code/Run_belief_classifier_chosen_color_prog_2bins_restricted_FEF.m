function Run_belief_classifier_chosen_color_prog_2bins_restricted_FEF(fsroot, ROI, n_tt, this_time)

par_param=0; %can't use par with cluster

% %test
% fsroot='/Volumes/buschman';
% ROI='LIP';
% n_tt=1;
% this_time=25;
% par_param=1;
% 

event='target';

for i=1:27
    initial_window=-500;
    event_list{i}='target';
    window_start_list(i)=initial_window+(i-1)*50;
end

task='Learning_Attentional_Templates';
subtask='exploreexploit/Restricted_FEF/Reset_RW_model';

window_size=200;subsubtask_classifier=sprintf('Classifier_%s_%d_to_%s_%d',event,window_start_list(1),event,window_start_list(end));

dirstem = fullfile('Projects',task,'Analysis','Electrophy_analysis',subtask,subsubtask_classifier);
data_path_clasifier = fullfile(fsroot,dirstem);


% boot_name='Pseudo_pop_chosen_color_no_prog_classifier_bootstrap';
% save_name=sprintf('Pseudo_pop_belief_chosen_color_no_prog_results_%s_%d_%d_%d',ROI,this_time,n_tt);

boot_name='Pseudo_pop_chosen_color_peak_template_classifier_bootstrap_2bins';
save_name=sprintf('Pseudo_pop_chosen_color_peak_template_2bins_classifier_results_%s_%d_%d_%d',ROI,this_time,n_tt);


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

N_neurons=size(this_ROI.chosen_color.TT_trials,2);
N_data=size(this_ROI.chosen_color.TT_trials,3); %
N_chosen_color=size(this_ROI.chosen_color.TT_trials,4);
N_belief=size(this_ROI.chosen_color.TT_trials,5);
% N_prog=size(this_ROI.chosen_color.TT_trials,6);

N_data_val=size(this_ROI.chosen_color.Val_trials,3)*N_chosen_color;

Classification_belief_correct=NaN(N_data_val*N_belief,1);
Classification_chosen_color_correct=NaN(N_data_val*N_belief,1);

Classification_chosen_color_each_belief_label=NaN(N_data_val*N_belief,N_belief);
Classification_chosen_color_each_belief_correct=NaN(N_data_val,N_belief,N_belief);

%% Belief classifier

for n=1:N_neurons
    for t=1:N_data
        for n_c=1:N_chosen_color
            data_this_belief(n,t+(n_c-1)*N_data)=this_ROI.Pseudo_pop(this_ROI.Neurons(n_tt,n),this_ROI.chosen_color.TT_trials(n_tt,n,t,n_c,1)).Classifier_FR(this_time); %%error Val and TT are inverted
        end
    end
end
for n=1:N_neurons
    for t=1:N_data/(N_chosen_color-1)
        for n_c=1:N_chosen_color
            data_other_belief(n,t+(n_c-1)*N_data)=this_ROI.Pseudo_pop(this_ROI.Neurons(n_tt,n),this_ROI.chosen_color.TT_trials(n_tt,n,t,n_c,2)).Classifier_FR(this_time);
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
    W_belief=mdl.Beta;
    Intercept_belief=mdl.Bias;
    
    SVMModel_belief.mdl=mdl;
    
else
    
    W_belief=zeros(N_neurons,1);
    Intercept_belief=svmmod.Bias;
    
    SVMModel_belief.mdl=svmmod;
    
end

clear c cv svmmod *_this_belief *_classifier labels probas scores mdl *_other_belief*


%% chosen_color classifier

for n=1:N_neurons
    for t=1:N_data
        for n_c=1:N_belief
            data_this_chosen_color(n,t+(n_c-1)*N_data)=this_ROI.Pseudo_pop(this_ROI.Neurons(n_tt,n),this_ROI.chosen_color.TT_trials(n_tt,n,t,1,n_c)).Classifier_FR(this_time); %%error Val and TT are inverted
        end
    end
end
for n=1:N_neurons
    for t=1:N_data/(N_chosen_color-1)
        for n_c=1:N_belief
            data_other_chosen_color(n,t+(n_c-1)*N_data)=this_ROI.Pseudo_pop(this_ROI.Neurons(n_tt,n),this_ROI.chosen_color.TT_trials(n_tt,n,t,2,n_c)).Classifier_FR(this_time);
        end
    end
end

data_classifier_chosen_color=[data_this_chosen_color,data_other_chosen_color]';

N_tt=size(data_classifier_chosen_color,1)/2;
groups_classifier_chosen_color=ones(2*N_tt,1);
groups_classifier_chosen_color(N_tt+1:end)=0;

c = cvpartition(size(groups_classifier_chosen_color,1),'KFold',10);

%optimize classfier

opts = struct('Optimizer','bayesopt','ShowPlots',false,'CVPartition',c,...
    'AcquisitionFunctionName','expected-improvement-plus','Verbose',0,'UseParallel',par_param);

svmmod = fitcsvm(data_classifier_chosen_color,groups_classifier_chosen_color,'KernelFunction','linear',...
    'OptimizeHyperparameters','auto','HyperparameterOptimizationOptions',opts);

if ~isempty(svmmod.Beta)
    
    mdl=fitPosterior(svmmod,data_classifier_chosen_color,groups_classifier_chosen_color);
    %
    W_chosen_color(:)=mdl.Beta;
    Intercept_chosen_color(:)=mdl.Bias;
    
    SVMModel_chosen_color(:).mdl=mdl;
    
else
    
    W_chosen_color(:)=zeros(N_neurons,1);
    Intercept_chosen_color(:)=svmmod.Bias;
    
    SVMModel_chosen_color(:).mdl=svmmod;
    
end

clear c cv svmmod *_this_chosen_color *_classifier mdl *_other_chosen_color*



%% Now decode the belief then the chosen_color

%now apply the classfier to the validation data
%for p=1:N_prog
for l=1:N_belief
    for q=1:N_chosen_color
        
        for t=1:size(this_ROI.chosen_color.Val_trials,3)
            for n=1:N_neurons
                data_val_this_belief(n,t+(q-1)*size(this_ROI.chosen_color.Val_trials,3)+(l-1)*N_data_val)=this_ROI.Pseudo_pop(this_ROI.Neurons(n_tt,n),this_ROI.chosen_color.Val_trials(n_tt,n,t,q,l)).Classifier_FR(this_time);
            end
            groups_val_this_belief(t+(q-1)*size(this_ROI.chosen_color.Val_trials,3)+(l-1)*N_data_val)=l;
            groups_val_this_chosen_color(t+(q-1)*size(this_ROI.chosen_color.Val_trials,3)+(l-1)*N_data_val)=q;
        end
    end
end

[class, ~,~] = predict(SVMModel_belief.mdl,data_val_this_belief');
class(class==0)=2;
Classification_belief_label=class;

[class, ~,~] = predict(SVMModel_chosen_color.mdl,data_val_this_belief');
class(class==0)=2;
Classification_chosen_color_label=class;

clear class

%by belief classifier
for i=1:size(Classification_belief_label,1)
    Classification_belief_correct(i)=(Classification_belief_label(i)==groups_val_this_belief(i));
end

%by chosen_color classifier
for i=1:size(Classification_chosen_color_label,1)
    Classification_chosen_color_correct(i)=(Classification_chosen_color_label(i)==groups_val_this_chosen_color(i));
end
clear *_this_belief *_this_chosen_color Classifier_*_proba
%end


%% Now do a chosen_color decoder for each template
for n_c=1:N_belief
    for n=1:N_neurons
        for t=1:N_data
            data_this_chosen_color(n,t)=this_ROI.Pseudo_pop(this_ROI.Neurons(n_tt,n),this_ROI.chosen_color.TT_trials(n_tt,n,t,1,n_c)).Classifier_FR(this_time);
        end
    end
    for n=1:N_neurons
        for t=1:N_data/(N_chosen_color-1)
            data_other_chosen_color(n,t)=this_ROI.Pseudo_pop(this_ROI.Neurons(n_tt,n),this_ROI.chosen_color.TT_trials(n_tt,n,t,2,n_c)).Classifier_FR(this_time);
        end
    end
    
    data_classifier_chosen_color=[data_this_chosen_color,data_other_chosen_color]';
    
    N_tt=size(data_classifier_chosen_color,1)/2;
    groups_classifier_chosen_color=ones(2*N_tt,1);
    groups_classifier_chosen_color(N_tt+1:end)=0;
    
    c = cvpartition(size(groups_classifier_chosen_color,1),'KFold',10);
    
    %optimize classfier
    
    opts = struct('Optimizer','bayesopt','ShowPlots',false,'CVPartition',c,...
        'AcquisitionFunctionName','expected-improvement-plus','Verbose',0,'UseParallel',par_param);
    
    svmmod = fitcsvm(data_classifier_chosen_color,groups_classifier_chosen_color,'KernelFunction','linear',...
        'OptimizeHyperparameters','auto','HyperparameterOptimizationOptions',opts);
    
    if ~isempty(svmmod.Beta)
        
        mdl=fitPosterior(svmmod,data_classifier_chosen_color,groups_classifier_chosen_color);
        %
        W_chosen_color_each_belief(:,n_c)=mdl.Beta;
        Intercept_chosen_color_each_belief(n_c)=mdl.Bias;
        
        SVMModel_chosen_color_each_belief(n_c).mdl=mdl;
        
    else
        
        W_chosen_color_each_belief(:,n_c)=zeros(N_neurons,1);
        Intercept_chosen_color_each_belief(n_c)=svmmod.Bias;
        
        SVMModel_chosen_color_each_belief(n_c).mdl=svmmod;
        
    end
    
    clear c cv svmmod *_this_chosen_color *_classifier mdl *_other_chosen_color*
end

%% Now decode the chosen_color for each template

%now apply the classfier to the validation data
% for p=1:N_prog

for n_c=1:N_belief
    for q=1:N_chosen_color
        for t=1:size(this_ROI.chosen_color.Val_trials,3)
            for n=1:N_neurons
                data_val_this_belief(n,t+(q-1)*size(this_ROI.chosen_color.Val_trials,3)+(n_c-1)*N_data_val)=this_ROI.Pseudo_pop(this_ROI.Neurons(n_tt,n),this_ROI.chosen_color.Val_trials(n_tt,n,t,q,n_c)).Classifier_FR(this_time);
            end
            groups_val_this_belief(t+(q-1)*size(this_ROI.chosen_color.Val_trials,3)+(n_c-1)*N_data_val)=n_c;
            groups_val_this_chosen_color(t+(q-1)*size(this_ROI.chosen_color.Val_trials,3)+(n_c-1)*N_data_val)=q;
        end
    end
end

for n_c=1:N_belief
    [class, ~,~] = predict(SVMModel_chosen_color_each_belief(n_c).mdl,data_val_this_belief');
    class(class==0)=2;
    Classification_chosen_color_each_belief_label(:,n_c)=class;
    clear class
end

%by chosen_color classifier
for n_c=1:N_belief
    for l=1:N_belief
        Classification_chosen_color_each_belief_correct(:,l,n_c)=(Classification_chosen_color_each_belief_label(groups_val_this_belief==l,n_c)'==groups_val_this_chosen_color(groups_val_this_belief==l));
    end
end
clear *_this_belief *_this_chosen_color Classifier_*_proba

% end


%%

save(fullfile(data_path_clasifier,save_name),'Classification_chosen_color_label','Classification_belief_label','Classification_belief_correct','Classification_chosen_color_correct','W_belief','W_chosen_color','Intercept_belief','Intercept_chosen_color','SVMModel_belief','SVMModel_chosen_color','Classification_chosen_color_each_belief_label','Classification_chosen_color_each_belief_correct','W_chosen_color_each_belief','Intercept_chosen_color_each_belief','SVMModel_chosen_color_each_belief')


end


