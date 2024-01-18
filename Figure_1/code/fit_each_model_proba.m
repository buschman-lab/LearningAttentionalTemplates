function fit_each_model_proba(data, generative_model, tested_model, model_index, load_name, save_path, proba_idx)


%% Fit the model using the generated data 
save_name = sprintf('%s_%s_%d_%d_proba',generative_model, tested_model, proba_idx, model_index);

switch generative_model
    case 'Reset'
        load(fullfile(save_path,load_name),'Reset_models')
        u=Reset_models(model_index).u;
        y=Reset_models(model_index).y;
    case 'No_reset'
        load(fullfile(save_path,load_name),'No_reset_models')
        u=No_reset_models(model_index).u;
        y=No_reset_models(model_index).y;
    case 'WSLS'
        load(fullfile(save_path,load_name),'WSLS_models')
        u=WSLS_models(model_index).u;
        y=WSLS_models(model_index).y;
end
switch tested_model
    case 'Reset'
            [posterior, out] = QlearningFuncApproxReset_RW_VBMC_recovery (data, u, y);
            save(fullfile(save_path,save_name),'posterior','out')
            fprintf('\nDone fitting %s with %s index %d\n',generative_model, tested_model, model_index);
    case 'No_reset'
            [posterior, out] = QlearningFuncApproxNoReset_RW_VBMC_recovery (data, u, y);
            save(fullfile(save_path,save_name),'posterior','out')
            fprintf('\nDone fitting %s with %s index %d\n',generative_model, tested_model, model_index);
    case 'WSLS'
            [posterior, out] = WSLS_VBMC_recovery (data, u, y);
            save(fullfile(save_path,save_name),'posterior','out')
            fprintf('\nDone fitting %s with %s index %d\n',generative_model, tested_model, model_index);
end


