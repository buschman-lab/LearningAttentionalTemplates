
%First generate the dat with a model
%Load the data for both monkeys
monkey='Both';

fsroot='/Volumes/buschman';
task='Learning_Attentional_Templates';


%if you want to load a monkey's fit (to check the param values, we use the bias of monkey B (Beaker))
% save_path=fullfile(fsroot,'Projects', task, 'Analysed_data',sprintf('%s','Beaker'));
% load_name=sprintf('Reset_%s_surprise_RW_%d_channels_VBMC','Beaker',6); %load the best model
% load(fullfile(save_path,load_name),'Reset')

data=load_monkey_data_continuous(fsroot,monkey,6,100);

cd('/Volumes/buschman/Users/Caroline/Learning_attentional_templates/Behavioral_Analysis/Figure1')


%% Make a fake data set
%we need to genrate what reward would have been obtained for each
%option (since the genrative model can choose any option)
true_template=data.X_data(:,1);
color_each_option=data.options_color;
location_each_option=data.locations;
for i=1:size(color_each_option,2)
    for l=1:3
        feedback_each_option(l,i)=round(vonmisespdf(color_each_option(l,i),true_template(i),2.5)*data.Reward_max(i))/data.Reward_max(i);
    end
end

%Make test parameters
%we don't change those parameters
Test_parameters.Volatility=0.5;

Test_parameters.B_location(2:4)=[0.1 0.05 -0.02];
Test_parameters.B_popout(1:2)=[-0.8 0.1];
Test_parameters.B_previous_color=0.5;

Test_parameters.color_pref=-pi/2;
Test_parameters.B_color_pref=0.05;

% we change those parameters
Test_parameters.BF_kappa_list=[1 1.7];
Test_parameters.learning_rate_list=[0.5 1];
Test_parameters.Reset_list=[0.3 0.5];
Test_parameters.Reward_threshold_list=[0 0.2];

for proba_idx=1:10
    %generate 10 models for each type

    save_name=sprintf('Fitting_pipeline_validation_proba_%d',proba_idx);

    %save the Test parameters and other variables
    save(fullfile(pwd,'proba',save_name))

    %generate data for each model and save in workspace
    %make the generative models and save them in the save file (that will
    %be appended)
    generate_data_each_model(Test_parameters, data, feedback_each_option,save_name,'proba')

end
