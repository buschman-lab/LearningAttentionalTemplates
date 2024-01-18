%First generate the dat with a model
%Load the data for both monkeys
clear all

monkey='Both';

model_type='proba';
%model_type='argmax';

fsroot='/Volumes/buschman';
task='Learning_Attentional_Templates';

proba_iter_list=[1 2 3 4 5 6 7 8 9 10];
load_name='Fitting_pipeline_validation_results';

BIC_matrix=NaN(14,3,length(proba_iter_list));
BIC_all=NaN(12,length(proba_iter_list)); %preferes no reset (= lower values)

addpath(genpath('CircStat2012a'))

data=load_monkey_data_continuous(fsroot,monkey,6,100);

cd('/Volumes/buschman/Users/Caroline/Learning_attentional_templates/Behavioral_Analysis/Figure1')


LOAD=0;

%% Before doing part 2, you must have ran fit_all_models_cluster

%% Part 2 of test parameter recovery

if LOAD==1

for proba_iter=1:length(proba_iter_list)
    
    %% Load everything and put them in a structure so it's easier to manipulate
    save_name=fullfile('proba',sprintf('Fitting_pipeline_validation_proba_%d.mat',proba_iter_list(proba_iter)));
    load(save_name,"WSLS_models","No_reset_models","Reset_models") %,"data"

    All_models(proba_iter).Reset_models=Reset_models;
    All_models(proba_iter).No_reset_models=No_reset_models;
    All_models(proba_iter).WSLS_models=WSLS_models;

    for i=1:8
        %reset
        this_load_name = sprintf('%s_%s_%d_%d','Reset', 'Reset', proba_iter_list(proba_iter), i);
        switch model_type
            case 'proba'
                this_load_name=fullfile('proba',[this_load_name '_proba']);
        end

        All_models(proba_iter).Reset_models(i).Reset=load(this_load_name,'posterior','out');
        %no reset
        this_load_name = sprintf('%s_%s_%d_%d','Reset', 'No_reset', proba_iter_list(proba_iter), i);
        switch model_type
            case 'proba'
                this_load_name=fullfile('proba',[this_load_name '_proba']);
        end
        All_models(proba_iter).Reset_models(i).No_reset=load(this_load_name,'posterior','out');
        %WSLS
        this_load_name = sprintf('%s_%s_%d_%d','Reset', 'WSLS', proba_iter_list(proba_iter), i);
        switch model_type
            case 'proba'
                this_load_name=fullfile('proba',[this_load_name '_proba']);
        end
        All_models(proba_iter).Reset_models(i).WSLS=load(this_load_name,'posterior','out');
    end
    for i=1:4
        %reset
        this_load_name = sprintf('%s_%s_%d_%d','No_reset', 'Reset', proba_iter_list(proba_iter), i);
        switch model_type
            case 'proba'
                this_load_name=fullfile('proba',[this_load_name '_proba']);
        end
        All_models(proba_iter).No_reset_models(i).Reset=load(this_load_name,'posterior','out');
        %no reset
        this_load_name = sprintf('%s_%s_%d_%d','No_reset', 'No_reset', proba_iter_list(proba_iter), i);
        switch model_type
            case 'proba'
                this_load_name=fullfile('proba',[this_load_name '_proba']);
        end
        All_models(proba_iter).No_reset_models(i).No_reset=load(this_load_name,'posterior','out');
        %WSLS
        this_load_name = sprintf('%s_%s_%d_%d','No_reset', 'WSLS', proba_iter_list(proba_iter), i);
        switch model_type
            case 'proba'
                this_load_name=fullfile('proba',[this_load_name '_proba']);
        end
        All_models(proba_iter).No_reset_models(i).WSLS=load(this_load_name,'posterior','out');
    end
    for i=1:2
        %reset
        this_load_name = sprintf('%s_%s_%d_%d','WSLS', 'Reset', proba_iter_list(proba_iter), i);
        switch model_type
            case 'proba'
                this_load_name=fullfile('proba',[this_load_name '_proba']);
        end

        All_models(proba_iter).WSLS_models(i).Reset=load(this_load_name,'posterior','out');
        %no reset
        this_load_name = sprintf('%s_%s_%d_%d','WSLS', 'No_reset', proba_iter_list(proba_iter), i);
        switch model_type
            case 'proba'
                this_load_name=fullfile('proba',[this_load_name '_proba']);
        end
        All_models(proba_iter).WSLS_models(i).No_reset=load(this_load_name,'posterior','out');
        %WSLS
        this_load_name = sprintf('%s_%s_%d_%d','WSLS', 'WSLS', proba_iter_list(proba_iter), i);
        switch model_type
            case 'proba'
                this_load_name=fullfile('proba',[this_load_name '_proba']);
        end
        All_models(proba_iter).WSLS_models(i).WSLS=load(this_load_name,'posterior','out');
    end

    %% Rename and predict with recovered parameters
    for i=1:2
        for j=1:2
            for k=1:2
                All_models(proba_iter).Reset_models(k+2*(j-1)+4*(i-1)).Recovered_parameters.BF_kappa=All_models(proba_iter).Reset_models(k+2*(j-1)+4*(i-1)).Reset.posterior.mean(1);
                All_models(proba_iter).Reset_models(k+2*(j-1)+4*(i-1)).Recovered_parameters.learning_rate=All_models(proba_iter).Reset_models(k+2*(j-1)+4*(i-1)).Reset.posterior.mean(2);
                All_models(proba_iter).Reset_models(k+2*(j-1)+4*(i-1)).Recovered_parameters.Reset=All_models(proba_iter).Reset_models(k+2*(j-1)+4*(i-1)).Reset.posterior.mean(9);
                All_models(proba_iter).Reset_models(k+2*(j-1)+4*(i-1)).Recovered_parameters.Volatility=All_models(proba_iter).Reset_models(k+2*(j-1)+4*(i-1)).Reset.posterior.mean(10);

                All_models(proba_iter).Reset_models(k+2*(j-1)+4*(i-1)).Recovered_parameters.B_location(2:4)=All_models(proba_iter).Reset_models(k+2*(j-1)+4*(i-1)).Reset.posterior.mean(3:5);
                All_models(proba_iter).Reset_models(k+2*(j-1)+4*(i-1)).Recovered_parameters.B_popout(1:2)=All_models(proba_iter).Reset_models(k+2*(j-1)+4*(i-1)).Reset.posterior.mean(6:7);
                All_models(proba_iter).Reset_models(k+2*(j-1)+4*(i-1)).Recovered_parameters.B_previous_color=All_models(proba_iter).Reset_models(k+2*(j-1)+4*(i-1)).Reset.posterior.mean(8);

                All_models(proba_iter).Reset_models(k+2*(j-1)+4*(i-1)).Recovered_parameters.color_pref=All_models(proba_iter).Reset_models(k+2*(j-1)+4*(i-1)).Reset.posterior.mean(11);
                All_models(proba_iter).Reset_models(k+2*(j-1)+4*(i-1)).Recovered_parameters.B_color_pref=All_models(proba_iter).Reset_models(k+2*(j-1)+4*(i-1)).Reset.posterior.mean(12);

                All_models(proba_iter).Reset_models(k+2*(j-1)+4*(i-1)).Model_recovery = get_QlearningFuncApproxReset_model_RW_VBMC(All_models(proba_iter).Reset_models(k+2*(j-1)+4*(i-1)).Recovered_parameters,All_models(proba_iter).Reset_models(k+2*(j-1)+4*(i-1)).u,data,data.X_data(:,5));
            end
        end
    end

    %% Latent space recovery for reset model
    for i=1:2
        for j=1:2
            for k=1:2
                for n=1:size(All_models(proba_iter).Reset_models(k+2*(j-1)+4*(i-1)).Model_predictions.model_outputs.Value_for_choice,2)
                    if sum(All_models(proba_iter).Reset_models(k+2*(j-1)+4*(i-1)).Model_predictions.model_outputs.Value_for_choice(:,n)==0)<size(All_models(proba_iter).Reset_models(k+2*(j-1)+4*(i-1)).Model_predictions.model_outputs.Value_for_choice,1)
                        [~, All_models(proba_iter).Reset_models(k+2*(j-1)+4*(i-1)).AT_generative(n)]=max(All_models(proba_iter).Reset_models(k+2*(j-1)+4*(i-1)).Model_predictions.model_outputs.Value_for_choice(:,n));
                    else
                        All_models(proba_iter).Reset_models(k+2*(j-1)+4*(i-1)).AT_generative(n)=NaN;
                    end
                end
                for n=1:size(All_models(proba_iter).Reset_models(k+2*(j-1)+4*(i-1)).Model_recovery.model_outputs.Value_for_choice,2)
                    if sum(All_models(proba_iter).Reset_models(k+2*(j-1)+4*(i-1)).Model_recovery.model_outputs.Value_for_choice(:,n)==0)<size(All_models(proba_iter).Reset_models(k+2*(j-1)+4*(i-1)).Model_recovery.model_outputs.Value_for_choice,1)
                        [~, All_models(proba_iter).Reset_models(k+2*(j-1)+4*(i-1)).AT_recovered(n)]=max(All_models(proba_iter).Reset_models(k+2*(j-1)+4*(i-1)).Model_recovery.model_outputs.Value_for_choice(:,n));
                    else
                        All_models(proba_iter).Reset_models(k+2*(j-1)+4*(i-1)).AT_recovered(n)=NaN;
                    end
                end
                %in rad
                All_models(proba_iter).Reset_models(k+2*(j-1)+4*(i-1)).AT_recovered=All_models(proba_iter).Reset_models(k+2*(j-1)+4*(i-1)).AT_recovered/100*2*pi;
                All_models(proba_iter).Reset_models(k+2*(j-1)+4*(i-1)).AT_generative=All_models(proba_iter).Reset_models(k+2*(j-1)+4*(i-1)).AT_generative/100*2*pi;
            end
        end
    end


    %% Rename and predict with recovered parameters for the no reset model
% for i=1:8
%     Reset_models(i).Recovered_parameters_no_reset.BF_kappa=Reset_models(i).No_reset.posterior.mean(1);
%     Reset_models(i).Recovered_parameters_no_reset.learning_rate=Reset_models(i).No_reset.posterior.mean(2);
% 
%     Reset_models(i).Recovered_parameters_no_reset.B_location(2:4)=Reset_models(i).No_reset.posterior.mean(3:5);
%     Reset_models(i).Recovered_parameters_no_reset.B_popout(1:2)=Reset_models(i).No_reset.posterior.mean(6:7);
%     Reset_models(i).Recovered_parameters_no_reset.B_previous_color=Reset_models(i).No_reset.posterior.mean(8);
% 
%     Reset_models(i).Recovered_parameters_no_reset.color_pref=Reset_models(i).No_reset.posterior.mean(9);
%     Reset_models(i).Recovered_parameters_no_reset.B_color_pref=Reset_models(i).No_reset.posterior.mean(10);
% 
%     Reset_models(i).Model_recovery_no_reset = get_QlearningFuncApproxNoReset_model_RW_VBMC(Reset_models(i).Recovered_parameters_no_reset,Reset_models(i).u,data);
%     for n=1:size(Reset_models(i).Model_recovery_no_reset.model_outputs.Value_for_choice,2)
%         if sum(Reset_models(i).Model_recovery_no_reset.model_outputs.Value_for_choice(:,n)==0)<size(Reset_models(i).Model_recovery_no_reset.model_outputs.Value_for_choice,1)
%             [~, Reset_models(i).AT_recovered_no_reset(n)]=max(Reset_models(i).Model_recovery_no_reset.model_outputs.Value_for_choice(:,n));
%         else
%             Reset_models(i).AT_recovered_no_reset(n)=NaN;
%         end
%     end
%     %in rad
%     Reset_models(i).AT_recovered_no_reset=Reset_models(i).AT_recovered_no_reset/100*2*pi;
% end
% 


    %% Rename parameters, predict and latent space recovery for No reset model

    for i=1:2
        for j=1:2

            All_models(proba_iter).No_reset_models(j+2*(i-1)).Recovered_parameters.BF_kappa=All_models(proba_iter).No_reset_models(j+2*(i-1)).No_reset.posterior.mean(1);
            All_models(proba_iter).No_reset_models(j+2*(i-1)).Recovered_parameters.learning_rate=All_models(proba_iter).No_reset_models(j+2*(i-1)).No_reset.posterior.mean(2);

            All_models(proba_iter).No_reset_models(j+2*(i-1)).Recovered_parameters.B_location(2:4)=All_models(proba_iter).No_reset_models(j+2*(i-1)).No_reset.posterior.mean(3:5);
            All_models(proba_iter).No_reset_models(j+2*(i-1)).Recovered_parameters.B_popout(1:2)=All_models(proba_iter).No_reset_models(j+2*(i-1)).No_reset.posterior.mean(6:7);
            All_models(proba_iter).No_reset_models(j+2*(i-1)).Recovered_parameters.B_previous_color=All_models(proba_iter).No_reset_models(j+2*(i-1)).No_reset.posterior.mean(8);

            All_models(proba_iter).No_reset_models(j+2*(i-1)).Recovered_parameters.color_pref=All_models(proba_iter).No_reset_models(j+2*(i-1)).No_reset.posterior.mean(9);
            All_models(proba_iter).No_reset_models(j+2*(i-1)).Recovered_parameters.B_color_pref=All_models(proba_iter).No_reset_models(j+2*(i-1)).No_reset.posterior.mean(10);

            All_models(proba_iter).No_reset_models(j+2*(i-1)).Model_recovery = get_QlearningFuncApproxNoReset_model_RW_VBMC(All_models(proba_iter).No_reset_models(j+2*(i-1)).Recovered_parameters,All_models(proba_iter).No_reset_models(j+2*(i-1)).u,data);

            for n=1:size(All_models(proba_iter).No_reset_models(j+2*(i-1)).Model_predictions.model_outputs.Value_for_choice,2)
                if sum(All_models(proba_iter).No_reset_models(j+2*(i-1)).Model_predictions.model_outputs.Value_for_choice(:,n)==0)<size(All_models(proba_iter).No_reset_models(j+2*(i-1)).Model_predictions.model_outputs.Value_for_choice,1)
                    [~, All_models(proba_iter).No_reset_models(j+2*(i-1)).AT_generative(n)]=max(All_models(proba_iter).No_reset_models(j+2*(i-1)).Model_predictions.model_outputs.Value_for_choice(:,n));
                else
                    All_models(proba_iter).No_reset_models(j+2*(i-1)).AT_generative(n)=NaN;
                end
            end
            for n=1:size(All_models(proba_iter).No_reset_models(j+2*(i-1)).Model_recovery.model_outputs.Value_for_choice,2)
                if sum(All_models(proba_iter).No_reset_models(j+2*(i-1)).Model_recovery.model_outputs.Value_for_choice(:,n)==0)<size(All_models(proba_iter).No_reset_models(j+2*(i-1)).Model_recovery.model_outputs.Value_for_choice,1)
                    [~, All_models(proba_iter).No_reset_models(j+2*(i-1)).AT_recovered(n)]=max(All_models(proba_iter).No_reset_models(j+2*(i-1)).Model_recovery.model_outputs.Value_for_choice(:,n));
                else
                    All_models(proba_iter).No_reset_models(j+2*(i-1)).AT_recovered(n)=NaN;
                end
            end
            %in rad
            All_models(proba_iter).No_reset_models(j+2*(i-1)).AT_recovered=All_models(proba_iter).No_reset_models(j+2*(i-1)).AT_recovered/100*2*pi;
            All_models(proba_iter).No_reset_models(j+2*(i-1)).AT_generative=All_models(proba_iter).No_reset_models(j+2*(i-1)).AT_generative/100*2*pi;

        end
    end

    %% Mean accuracy of the genrative models
    for i=1:2
        for j=1:2
            for k=1:2
                All_models(proba_iter).Reset_models(k+2*(j-1)+4*(i-1)).mean_accuracy=mean(All_models(proba_iter).Reset_models(k+2*(j-1)+4*(i-1)).y(1,:));
            end
            All_models(proba_iter).No_reset_models(j+2*(i-1)).mean_accuracy=mean(All_models(proba_iter).No_reset_models(j+2*(i-1)).y(1,:));
        end
        All_models(proba_iter).WSLS_models(i).mean_accuracy=mean(All_models(proba_iter).WSLS_models(i).y(1,:));
    end


    %% BIC for model selection
    for i=1:2
        %generative model = WSLS
        BIC_matrix(i,1,proba_iter)=All_models(proba_iter).WSLS_models(i).WSLS.out.fit.BIC-All_models(proba_iter).WSLS_models(i).WSLS.out.fit.BIC;
        BIC_matrix(i,2,proba_iter)=All_models(proba_iter).WSLS_models(i).No_reset.out.fit.BIC-All_models(proba_iter).WSLS_models(i).WSLS.out.fit.BIC;
        BIC_matrix(i,3,proba_iter)=All_models(proba_iter).WSLS_models(i).Reset.out.fit.BIC-All_models(proba_iter).WSLS_models(i).WSLS.out.fit.BIC;
    end
    for i=1:2
        for j=1:2
            %generative model = No Reset
            BIC_matrix(2+j+2*(i-1),1,proba_iter)=All_models(proba_iter).No_reset_models(j+2*(i-1)).WSLS.out.fit.BIC-All_models(proba_iter).No_reset_models(j+2*(i-1)).No_reset.out.fit.BIC;
            BIC_matrix(2+j+2*(i-1),2,proba_iter)=All_models(proba_iter).No_reset_models(j+2*(i-1)).No_reset.out.fit.BIC-All_models(proba_iter).No_reset_models(j+2*(i-1)).No_reset.out.fit.BIC;
            BIC_matrix(2+j+2*(i-1),3,proba_iter)=All_models(proba_iter).No_reset_models(j+2*(i-1)).Reset.out.fit.BIC-All_models(proba_iter).No_reset_models(j+2*(i-1)).No_reset.out.fit.BIC;
        end
    end
    for i=1:8
        %generative model = Reset
        BIC_matrix(6+i,1,proba_iter)=All_models(proba_iter).Reset_models(i).WSLS.out.fit.BIC-All_models(proba_iter).Reset_models(i).Reset.out.fit.BIC;
        BIC_matrix(6+i,2,proba_iter)=All_models(proba_iter).Reset_models(i).No_reset.out.fit.BIC-All_models(proba_iter).Reset_models(i).Reset.out.fit.BIC;
        BIC_matrix(6+i,3,proba_iter)=All_models(proba_iter).Reset_models(i).Reset.out.fit.BIC-All_models(proba_iter).Reset_models(i).Reset.out.fit.BIC;
    end


    %% BIC reset and no reset only
    for i=1:4
        %generative model = No Reset
        BIC_all(i,proba_iter)=All_models(proba_iter).No_reset_models(i).No_reset.out.fit.BIC-All_models(proba_iter).No_reset_models(i).Reset.out.fit.BIC;
    end
    for i=1:8
        %generative model = Reset
        BIC_all(4+i,proba_iter)=All_models(proba_iter).Reset_models(i).No_reset.out.fit.BIC-All_models(proba_iter).Reset_models(i).Reset.out.fit.BIC;
    end

    %% Accuracy
    for i=1:4
        %generative model = No Reset
        Accuracy_all(i,proba_iter)=All_models(proba_iter).No_reset_models(i).mean_accuracy;
    end
    for i=1:8
        %generative model = Reset
        Accuracy_all(4+i,proba_iter)=All_models(proba_iter).Reset_models(i).mean_accuracy;
    end


        %% Accuracy recovered
    for i=1:4
        %generative model = No Reset
        Accuracy_recovered_all(i,proba_iter)=mean(All_models(proba_iter).No_reset_models(i).No_reset.out.suffStat.gx(1,:)) ;
    end
    for i=1:8
        %generative model = Reset
        Accuracy_recovered_all(4+i,proba_iter)=mean(All_models(proba_iter).Reset_models(i).Reset.out.suffStat.gx(1,:));
    end

    %% Reset parameter

    for i=1:4
        %generative model = No Reset
        Reset_all(i,proba_iter)=All_models(proba_iter).No_reset_models(i).Reset.posterior.mean(9);
    end
    for i=1:8
        %generative model = Reset
        Reset_all(4+i,proba_iter)=All_models(proba_iter).Reset_models(i).Reset.posterior.mean(9);
    end

    %% template recovery
    for i=1:4
        [r(i,proba_iter), p(i,proba_iter)]=circ_corrcc(All_models(proba_iter).No_reset_models(i).AT_recovered(~isnan(All_models(proba_iter).No_reset_models(i).AT_recovered+All_models(proba_iter).No_reset_models(i).AT_generative)),All_models(proba_iter).No_reset_models(i).AT_generative(~isnan(All_models(proba_iter).No_reset_models(i).AT_recovered+All_models(proba_iter).No_reset_models(i).AT_generative)));
    end
    for i=1:8
        [r(i+4,proba_iter), p(i+4,proba_iter)]=circ_corrcc(All_models(proba_iter).Reset_models(i).AT_recovered(~isnan(All_models(proba_iter).Reset_models(i).AT_recovered+All_models(proba_iter).Reset_models(i).AT_generative)),All_models(proba_iter).Reset_models(i).AT_generative(~isnan(All_models(proba_iter).Reset_models(i).AT_recovered+All_models(proba_iter).Reset_models(i).AT_generative)));
    end

    %% number of reset generative
    for i=1:4
        %generative model = No Reset
        Total_reset_all(i,proba_iter)=0;
    end
    for i=1:8
        %generative model = Reset
        Total_reset_all(4+i,proba_iter)=sum(All_models(proba_iter).Reset_models(i).Model_predictions.model_outputs.Switch);
    end


    %% reset frecovery
    for i=1:4
        %generative model = No Reset
        Total_reset_recovery_all(i,proba_iter)=sum(All_models(proba_iter).Reset_models(i).Model_recovery.model_outputs.Switch);
    end
    for i=1:8
        %generative model = Reset
        Total_reset_recovery_all(4+i,proba_iter)=sum(All_models(proba_iter).Reset_models(i).Model_recovery.model_outputs.Switch);
    end
    %% Preferred color
    for i=1:8
        %generative model = Reset
        Pref_color_all(i,proba_iter)=mod(All_models(proba_iter).Reset_models(i).Test_parameters.color_pref-All_models(proba_iter).Reset_models(i).Recovered_parameters.color_pref+pi,2*pi)-pi;
    end

    %% save

%    save(fullfile(pwd,save_name),"-append")

clear Reset_models No_reset_models WSLS_models

end

else

    load(load_name);

    for proba_iter=1:length(proba_iter_list)
    %% template recovery
    for i=1:4
        [r(i,proba_iter), p(i,proba_iter)]=circ_corrcc(All_models(proba_iter).No_reset_models(i).AT_recovered(~isnan(All_models(proba_iter).No_reset_models(i).AT_recovered+All_models(proba_iter).No_reset_models(i).AT_generative)),All_models(proba_iter).No_reset_models(i).AT_generative(~isnan(All_models(proba_iter).No_reset_models(i).AT_recovered+All_models(proba_iter).No_reset_models(i).AT_generative)));
    end
    for i=1:8
        [r(i+4,proba_iter), p(i+4,proba_iter)]=circ_corrcc(All_models(proba_iter).Reset_models(i).AT_recovered(~isnan(All_models(proba_iter).Reset_models(i).AT_recovered+All_models(proba_iter).Reset_models(i).AT_generative)),All_models(proba_iter).Reset_models(i).AT_generative(~isnan(All_models(proba_iter).Reset_models(i).AT_recovered+All_models(proba_iter).Reset_models(i).AT_generative)));
    end
    end
end


%% reset template recovery
% figure
% imagesc((BIC_matrix./max(BIC_matrix,[],1))')
% colorbar


%%
new_order=[1 3 5 7 2 4 6 8]+4; %sort by reset threshold value

figure
subplot(3,2,3)
hold on
for proba_iter=1:length(proba_iter_list)
    plot(1:4,BIC_all(1:4,proba_iter),'o','MarkerFaceColor',[0.75 0.75 0.75],'MarkerEdgeColor','none','MarkerSize',5) %r=0.3
    plot(5:8,BIC_all(new_order(1:4),proba_iter),'o','MarkerFaceColor',[0 0 0],'MarkerEdgeColor','none','MarkerSize',5) %r=0.3
    plot(9:12,BIC_all(new_order(5:8),proba_iter),'o','MarkerFaceColor',[0 0.5 0],'MarkerEdgeColor','none','MarkerSize',5) %r=0.3

end
yline(0,'--','Color',[0.5 0.5 0.5])
xlabel('model')
ylabel('Prefers reset model (when delta BIC>0)')
box off
xlim([0 13])

subplot(3,2,2)
hold on
for proba_iter=1:length(proba_iter_list)
    plot(1:4,log(Total_reset_all(1:4,proba_iter)/length(data.chosen_color)),'o','MarkerFaceColor',[0.75 0.75 0.75],'MarkerEdgeColor','none','MarkerSize',5) %r=0.3
    plot(5:8,log(Total_reset_all(new_order(1:4),proba_iter)/length(data.chosen_color)),'o','MarkerFaceColor',[0 0 0],'MarkerEdgeColor','none','MarkerSize',5) %r=0.3
    plot(9:12,log(Total_reset_all(new_order(5:8),proba_iter)/length(data.chosen_color)),'o','MarkerFaceColor',[0 0.5 0],'MarkerEdgeColor','none','MarkerSize',5) %r=0.3

    plot(1.4:4.4,log(Total_reset_recovery_all(1:4,proba_iter)/length(data.chosen_color)),'o','MarkerEdgeColor',[0.75 0.75 0.75],'MarkerFaceColor','none','MarkerSize',5) %r=0.3
    plot(5.4:8.4,log(Total_reset_recovery_all(new_order(1:4),proba_iter)/length(data.chosen_color)),'o','MarkerEdgeColor',[0 0 0],'MarkerFaceColor','none','MarkerSize',5) %r=0.3
    plot(9.4:12.4,log(Total_reset_recovery_all(new_order(5:8),proba_iter)/length(data.chosen_color)),'o','MarkerEdgeColor',[0 0.5 0],'MarkerFaceColor','none','MarkerSize',5) %r=0.3

end
yline(log(113/length(data.chosen_color)),'--','Color',[0.5 0.5 0.5])
xlabel('model')
ylabel('Log proportion trials with reset')
box off
xlim([0 13])

subplot(3,2,1)
hold on
for proba_iter=1:length(proba_iter_list)
    plot(1:4,Accuracy_all(1:4,proba_iter),'o','MarkerFaceColor',[0.75 0.75 0.75],'MarkerEdgeColor','none','MarkerSize',5) %r=0.3
    plot(5:8,Accuracy_all(new_order(1:4),proba_iter),'o','MarkerFaceColor',[0 0 0],'MarkerEdgeColor','none','MarkerSize',5) %r=0.3
    plot(9:12,Accuracy_all(new_order(5:8),proba_iter),'o','MarkerFaceColor',[0 0.5 0],'MarkerEdgeColor','none','MarkerSize',5) %r=0.3

    plot(1.4:4.4,Accuracy_recovered_all(1:4,proba_iter),'o','MarkerEdgeColor',[0.75 0.75 0.75],'MarkerFaceColor','none','MarkerSize',5) %r=0.3
    plot(5.4:8.4,Accuracy_recovered_all(new_order(1:4),proba_iter),'o','MarkerEdgeColor',[0 0 0],'MarkerFaceColor','none','MarkerSize',5) %r=0.3
    plot(9.4:12.4,Accuracy_recovered_all(new_order(5:8),proba_iter),'o','MarkerEdgeColor',[0 0.5 0],'MarkerFaceColor','none','MarkerSize',5) %r=0.3
end
yline(1/3,'--','Color',[0.5 0.5 0.5])
xlabel('model')
ylabel('Accuracy')
box off
xlim([0 13])

subplot(3,2,4)
hold on
for proba_iter=1:length(proba_iter_list)
    plot(1:4,Reset_all(1:4,proba_iter),'o','MarkerFaceColor',[0.75 0.75 0.75],'MarkerEdgeColor','none','MarkerSize',5) %r=0.3
    plot(5:8,Reset_all(new_order(1:4),proba_iter),'o','MarkerFaceColor',[0 0 0],'MarkerEdgeColor','none','MarkerSize',5) %r=0.3
    plot(9:12,Reset_all(new_order(5:8),proba_iter),'o','MarkerFaceColor',[0 0.5 0],'MarkerEdgeColor','none','MarkerSize',5) %r=0.3

end
yline(0.3,'--','Color',[0 0 0])
yline(0.5,'--','Color',[0 0.5 0])
xlabel('model')
ylabel('Fitted reset threshold')
box off
xlim([0 13])


subplot(3,2,5)
hold on
for proba_iter=1:length(proba_iter_list)
    plot(1:4,r(1:4,proba_iter),'o','MarkerFaceColor',[0.75 0.75 0.75],'MarkerEdgeColor','none','MarkerSize',5) %r=0.3
    plot(5:8,r(new_order(1:4),proba_iter),'o','MarkerFaceColor',[0 0 0],'MarkerEdgeColor','none','MarkerSize',5) %r=0.3
    plot(9:12,r(new_order(5:8),proba_iter),'o','MarkerFaceColor',[0 0.5 0],'MarkerEdgeColor','none','MarkerSize',5) %r=0.3

end
xlabel('model')
ylabel('Template recovery')
box off
xlim([0 13])

% 
%% set color
new_order=[1 3 5 7 2 4 6 8]; %sort by reset threshold value
color_for_plot=repmat([0 0 0],4,1);
color_for_plot=[color_for_plot; repmat([0 0.5 0],4,1)];

figure

subplot(2,2,1)
for proba_iter=1:length(proba_iter_list)
for i=1:length(new_order)
    hold on
    plot(All_models(proba_iter).Reset_models(new_order(i)).Test_parameters.learning_rate/exp(All_models(proba_iter).Reset_models(new_order(i)).Test_parameters.BF_kappa),All_models(proba_iter).Reset_models(new_order(i)).Recovered_parameters.learning_rate/exp(All_models(proba_iter).Reset_models(new_order(i)).Recovered_parameters.BF_kappa),'o','MarkerFaceColor',color_for_plot(i,:),'MarkerEdgeColor','none','MarkerSize',5)
end
end
plot(-0.1:0.01:0.5,-0.1:0.01:0.5,'k--')
xlim([-0.1 0.5])
box off
xlabel('Generative learning rate')
ylabel('Recovered learning rate')

% figure
subplot(2,2,2)
for proba_iter=1:length(proba_iter_list)
for i=1:length(new_order)
    for j=2:4
    hold on
    plot(All_models(proba_iter).Reset_models(new_order(i)).Test_parameters.B_location(j),All_models(proba_iter).Reset_models(new_order(i)).Recovered_parameters.B_location(j),'o','MarkerFaceColor',color_for_plot(i,:),'MarkerEdgeColor','none','MarkerSize',5)
    end
end
end
plot(-0.1:0.01:0.15,-0.1:0.01:0.15,'k--')
xlim([-0.1 0.15])
box off
xlabel('Generative location bias')
ylabel('Recovered locations bias')

subplot(2,2,3)
for proba_iter=1:length(proba_iter_list)
for i=1:length(new_order)
    for j=1:2
    hold on
    plot(All_models(proba_iter).Reset_models(new_order(i)).Test_parameters.B_popout(j),All_models(proba_iter).Reset_models(new_order(i)).Recovered_parameters.B_popout(j),'o','MarkerFaceColor',color_for_plot(i,:),'MarkerEdgeColor','none','MarkerSize',5)
    end
end
end
plot(-1:0.01:0.4,-1:0.01:0.4,'k--')
xlim([-1 0.4])
box off
xlabel('Generative stim size bias')
ylabel('Recovered stim size bias')

subplot(2,2,4)
for proba_iter=1:length(proba_iter_list)
for i=1:length(new_order)
    hold on
    plot(All_models(proba_iter).Reset_models(new_order(i)).Test_parameters.B_previous_color,All_models(proba_iter).Reset_models(new_order(i)).Recovered_parameters.B_previous_color,'o','MarkerFaceColor',color_for_plot(i,:),'MarkerEdgeColor','none','MarkerSize',5)
        plot(All_models(proba_iter).Reset_models(new_order(i)).Test_parameters.B_color_pref,All_models(proba_iter).Reset_models(new_order(i)).Recovered_parameters.B_color_pref,'o','MarkerFaceColor',color_for_plot(i,:),'MarkerEdgeColor','none','MarkerSize',5)
    end
end
plot(-0.1:0.01:0.7,-0.1:0.01:0.7,'k--')
xlim([-0.1 0.7])
box off
xlabel('Generative color bias')
ylabel('Recovered color bias')


