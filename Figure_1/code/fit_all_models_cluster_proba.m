function fit_all_models_cluster_proba(fsroot,save_name,model_index,proba_idx)

save_name=sprintf('%s_%d.mat',save_name,proba_idx);

%load the test parameters
save_path=fullfile(fsroot,'Users/Caroline/Learning_attentional_templates/Behavioral_Analysis/Figure1/proba');
load(fullfile(save_path,save_name),'Test_parameters','data');

total_reset=length(Test_parameters.BF_kappa_list)*length(Test_parameters.learning_rate_list)*length(Test_parameters.Reset_list);
total_no_reset=length(Test_parameters.BF_kappa_list)*length(Test_parameters.learning_rate_list);
total_WSLS=length(Test_parameters.Reward_threshold_list);

%fit each model with each model and save in a different file
if model_index<=(total_reset+total_no_reset+total_WSLS) %first model to test is the reset model for all generative models
    if model_index<=total_reset
        fprintf('\nStarted fit reset with reset index %d\n',model_index)
        fit_each_model_proba(data, 'Reset', 'Reset', model_index, save_name, save_path,proba_idx)
    elseif model_index>total_reset && model_index<=total_reset+total_no_reset
        fprintf('\nStarted fit no reset with reset index %d\n',model_index-total_reset)
        fit_each_model_proba(data, 'No_reset', 'Reset', model_index-total_reset, save_name, save_path,proba_idx)
    else
        fprintf('\nStarted fit WSLS with reset index %d\n',model_index-total_reset-total_no_reset)
        fit_each_model_proba(data, 'WSLS', 'Reset', model_index-total_reset-total_no_reset, save_name, save_path,proba_idx)
    end
elseif  model_index>(total_reset+total_no_reset+total_WSLS) && model_index<=2*(total_reset+total_no_reset+total_WSLS)
    if model_index-(total_reset+total_no_reset+total_WSLS)<=total_reset
        fprintf('\nStarted fit reset with no reset index %d\n',model_index-(total_reset+total_no_reset+total_WSLS))
        fit_each_model_proba(data, 'Reset', 'No_reset', model_index-(total_reset+total_no_reset+total_WSLS), save_name, save_path,proba_idx)
    elseif model_index-(total_reset+total_no_reset+total_WSLS)>total_reset && model_index-(total_reset+total_no_reset+total_WSLS)<=total_reset+total_no_reset
        fprintf('\nStarted fit no reset with no reset index %d\n',model_index-(total_reset+total_no_reset+total_WSLS)-total_reset)
        fit_each_model_proba(data, 'No_reset', 'No_reset', model_index-(total_reset+total_no_reset+total_WSLS)-total_reset, save_name, save_path,proba_idx)
    else
        fprintf('\nStarted fit WSLS with no reset index %d\n',model_index-(total_reset+total_no_reset+total_WSLS)-total_reset-total_no_reset)
        fit_each_model_proba(data, 'WSLS', 'No_reset', model_index-(total_reset+total_no_reset+total_WSLS)-total_reset-total_no_reset, save_name, save_path,proba_idx)
    end
else
    if model_index-2*(total_reset+total_no_reset+total_WSLS)<=total_reset
        fprintf('\nStarted fit reset with WSLS index %d\n',model_index-2*(total_reset+total_no_reset+total_WSLS))
        fit_each_model_proba(data, 'Reset', 'WSLS', model_index-2*(total_reset+total_no_reset+total_WSLS), save_name, save_path,proba_idx)
    elseif model_index-2*(total_reset+total_no_reset+total_WSLS)>total_reset && model_index-2*(total_reset+total_no_reset+total_WSLS)<=total_reset+total_no_reset
        fprintf('\nStarted fit no reset with WSLS index %d\n',model_index-2*(total_reset+total_no_reset+total_WSLS)-total_reset)
        fit_each_model_proba(data, 'No_reset', 'WSLS', model_index-2*(total_reset+total_no_reset+total_WSLS)-total_reset, save_name, save_path,proba_idx)
    else
        fprintf('\nStarted fit WSLS with WSLS index %d\n',model_index-2*(total_reset+total_no_reset+total_WSLS)-total_reset-total_no_reset)
        fit_each_model_proba(data, 'WSLS', 'WSLS', model_index-2*(total_reset+total_no_reset+total_WSLS)-total_reset-total_no_reset, save_name, save_path,proba_idx)
    end
end
