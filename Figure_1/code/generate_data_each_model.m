function generate_data_each_model(Test_parameters, data, feedback_each_option,save_name,model_type)


%% Generate the data for each model
for i=1:2
    for j=1:2
        for k=1:2
            %replace REset with better estimates for the reset
            %threshold
            Reset_models(k+2*(j-1)+4*(i-1)).Test_parameters=Test_parameters;
            Reset_models(k+2*(j-1)+4*(i-1)).Test_parameters.BF_kappa=Test_parameters.BF_kappa_list(i);
            Reset_models(k+2*(j-1)+4*(i-1)).Test_parameters.learning_rate=Test_parameters.learning_rate_list(j);
            Reset_models(k+2*(j-1)+4*(i-1)).Test_parameters.Reset=Test_parameters.Reset_list(k);
            switch model_type
                case 'argmax'
                    [Reset_models(k+2*(j-1)+4*(i-1)).Model_predictions, Reset_models(k+2*(j-1)+4*(i-1)).u, Reset_models(k+2*(j-1)+4*(i-1)).y] = make_QlearningFuncApproxReset_model_RW_VBMC(Reset_models(k+2*(j-1)+4*(i-1)).Test_parameters,data,feedback_each_option);
                case 'proba'
                    [Reset_models(k+2*(j-1)+4*(i-1)).Model_predictions, Reset_models(k+2*(j-1)+4*(i-1)).u, Reset_models(k+2*(j-1)+4*(i-1)).y] = make_Reset_proba(Reset_models(k+2*(j-1)+4*(i-1)).Test_parameters,data,feedback_each_option);
            end
        end
    end
end
for i=1:2
    for j=1:2
        %don't update those
        No_reset_models(j+2*(i-1)).Test_parameters=Test_parameters;
        No_reset_models(j+2*(i-1)).Test_parameters.BF_kappa=Test_parameters.BF_kappa_list(i);
        No_reset_models(j+2*(i-1)).Test_parameters.learning_rate=Test_parameters.learning_rate_list(j);
        switch model_type
            case 'argmax'
                [No_reset_models(j+2*(i-1)).Model_predictions, No_reset_models(j+2*(i-1)).u, No_reset_models(j+2*(i-1)).y] = make_QlearningFuncApproxNoReset_model_RW_VBMC(No_reset_models(j+2*(i-1)).Test_parameters,data,feedback_each_option);
            case 'proba'
                [No_reset_models(j+2*(i-1)).Model_predictions, No_reset_models(j+2*(i-1)).u, No_reset_models(j+2*(i-1)).y] = make_NoReset_proba(No_reset_models(j+2*(i-1)).Test_parameters,data,feedback_each_option);
        end
    end
end
for i=1:2
    %don't update those
    WSLS_models(i).Test_parameters=Test_parameters;
    WSLS_models(i).Test_parameters.Reward_threshold=Test_parameters.Reward_threshold_list(i);
    switch model_type
        case 'argmax'
            [WSLS_models(i).Model_predictions, WSLS_models(i).u, WSLS_models(i).y] = make_WSLS_VBMC(WSLS_models(i).Test_parameters,data,feedback_each_option);
        case 'proba'
            [WSLS_models(i).Model_predictions, WSLS_models(i).u, WSLS_models(i).y] = make_WSLS_proba(WSLS_models(i).Test_parameters,data,feedback_each_option);
    end
end

switch model_type
    case 'argmax'
        save(fullfile(pwd,save_name),"WSLS_models","No_reset_models","Reset_models",'-append')
    case 'proba'
        save(fullfile(pwd,'proba',save_name),"WSLS_models","No_reset_models","Reset_models",'-append')
end