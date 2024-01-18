function Run_stim_color_classifier_across_time_with_NN_restricted_FEF(fsroot, ROI, n_tt, this_time, loc)

par_param=0; %can't use par with cluster

task='Learning_Attentional_Templates';
subtask='exploreexploit/Restricted_FEF/Reset_RW_model';

% window_size=300;
% step_size=1;
% subsubtask_classifier=sprintf('stim_color_encoding_%s_time_analysis_learning_%d_step_%d',event,window_size,step_size);

event='target';

for i=1:27
    window_start_list(i)=-500+(i-1)*50;
end

window_size=200;
subsubtask_classifier=sprintf('Classifier_%s_%d_to_%s_%d',event,window_start_list(1),event,window_start_list(end));

dirstem = fullfile('Projects',task,'Analysis','Electrophy_analysis',subtask,subsubtask_classifier);
data_path_clasifier = fullfile(fsroot,dirstem);


boot_name='Pseudo_pop_stim_color_no_prog_bootstrap';
save_name=sprintf('Pseudo_pop_stim_color_no_prog_results_%s_%d_%d_%d',ROI,this_time,n_tt,loc);

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
N_boot=1;
N_neurons=size(this_ROI.Stim_color.TT_trials,2);
N_data=size(this_ROI.Stim_color.TT_trials,3); %
%we need an even number to have the same number of other stim_colors (for 3
%categories)
if mod(N_data,2)~=0
    N_data=N_data-1;
end
N_data_tt=size(this_ROI.Stim_color.TT_trials,3); %

N_stim=size(this_ROI.Stim_color.TT_trials,4);
N_chosen=size(this_ROI.Stim_color.TT_trials,5);

N_data_val=size(this_ROI.Stim_color.Val_trials,3)*N_stim*N_chosen;

Classification_label=NaN(N_times,N_data_val,N_boot);
Classification_this_label=NaN(N_times,N_data_val,N_boot);
Classification_correct=NaN(N_times,N_data_val,N_boot);
Classification_net_correct=NaN(N_times,N_data_val,N_boot);
Classification_net_label=NaN(N_times,N_data_val,N_boot);
Classifier_proba=NaN(N_data_val,N_stim,N_times);

%%
for k=1:N_stim
    for n=1:N_neurons
        for t=1:N_data
            for j=1:N_chosen
                data_this_stim_color(n,t+(j-1)*N_data)=this_ROI.Pseudo_pop(this_ROI.Neurons(n_tt,n),this_ROI.Stim_color.TT_trials(n_tt,n,t,k,j,loc)).Classifier_FR(this_time); 
            end
        end
    end
    %   subsample to have the same number of chosen/unchosen
    other_stim_color_ind=randsample(N_data,N_data/(N_stim-1),'false');
    increment=0;
    for l=1:N_stim
        if l~=k
            for n=1:N_neurons
                for t=1:N_data/(N_stim-1)
                    for j=1:N_chosen
                        data_other_stim_color(n,t+increment*N_data/(N_stim-1)+(j-1)*N_data)=this_ROI.Pseudo_pop(this_ROI.Neurons(n_tt,n),this_ROI.Stim_color.TT_trials(n_tt,n,other_stim_color_ind(t),l,loc)).Classifier_FR(this_time);
                    end
                end
            end
            increment=1;
        end
    end
    
    data_classifier=[data_this_stim_color,data_other_stim_color]';
    
    N_tt=size(data_classifier,1)/2;
    groups_classifier=ones(2*N_tt,1);
    groups_classifier(N_tt+1:end)=0;
    
    c = cvpartition(size(groups_classifier,1),'KFold',10);
    
    % optimize classfier
    
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
    
    clear c cv *_this_stim_color *_classifier labels probas scores mdl *_other_stim_color*
end

%% Now get a NN for combination of classifiers
for t=1:N_data_tt
    for k=1:N_stim
        for j=1:N_chosen
            for n=1:N_neurons
                data_tt(n,t+(k-1)*N_data_tt+(j-1)*N_data_tt*N_stim)=this_ROI.Pseudo_pop(this_ROI.Neurons(n_tt,n),this_ROI.Stim_color.TT_trials(n_tt,n,t,k,j,loc)).Classifier_FR(this_time);
            end
            groups_tt(t+(k-1)*N_data_tt+(j-1)*N_data_tt*N_stim)=k;
        end
    end
end
%apply the classifiers
for k=1:N_stim
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


%% Apply the classifier and NN to validation data

for tp=1:N_times
    for l=1:N_stim
        for t=1:size(this_ROI.Stim_color.Val_trials,3)
            for j=1:N_chosen
                for n=1:N_neurons
                    data_val_this_stim_color(n,t+(l-1)*size(this_ROI.Stim_color.Val_trials,3)+(j-1)*N_stim*size(this_ROI.Stim_color.Val_trials,3))=this_ROI.Pseudo_pop(this_ROI.Neurons(n_tt,n),this_ROI.Stim_color.Val_trials(n_tt,n,t,l,j,loc)).Classifier_FR(tp);  %%error Val and TT are inverted
                end
                groups_val_this_stim_color(t+(l-1)*size(this_ROI.Stim_color.Val_trials,3)+(j-1)*N_stim*size(this_ROI.Stim_color.Val_trials,3))=l;
            end
        end
    end
    for k=1:N_stim
        [~, probas,~] = predict(SVMModel(k,1).mdl,data_val_this_stim_color');
        Classifier_proba(:,k,tp)=probas(:,2);
        clear probas labels
    end
    
    %no NN
    for i=1:length(Classifier_proba)
        [Classification_this_label(tp,i), Classification_label(tp,i)]=max(Classifier_proba(i,:,tp));
        Classification_correct(tp,i)=(Classification_label(tp,i)==groups_val_this_stim_color(i));
    end
    %with NN
    Classification_net_label(tp,:)=predict(mdl_net,Classifier_proba(:,:,tp));
    for i=1:length(groups_val_this_stim_color)
        Classification_net_correct(tp,i)=(Classification_net_label(tp,i)==groups_val_this_stim_color(i));
    end
    accuracy_val(tp)=1-loss(mdl_net,Classifier_proba(:,:,tp),groups_val_this_stim_color);
    
%     figure
%     confusionchart(Classification_net_label(tp,:),groups_val_this_stim_color)
    
    clear *_this_stim_color
end

%%


save(fullfile(data_path_clasifier,save_name),'Classification*','W','Intercept','SVMModel','Test_net','accuracy*','Classifier*')


end


