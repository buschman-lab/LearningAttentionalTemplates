%save the model predictions for each session

clear all

fsroot='/Volumes/buschman';
task='Learning_attentional_templates';
save_path = '/Volumes/buschman/Projects/Learning_Attentional_Templates/Analysed_data';

%% choose the model
monkey='Beaker';
N_channels=6;
model_type='RW';
switch_type='Reset_no_bias_prev_cc';
N_bins=100;

load_path=fullfile(fsroot,'Projects', task, 'Analysed_data', sprintf('%s',monkey));

load_name=sprintf('%s_%s_surprise_%s_%d_channels_VBMC',switch_type,monkey,model_type,N_channels);

switch switch_type
    case 'Reset'
        load(fullfile(load_path,load_name),'Reset');
        Parameters=Reset.Mean;
    case 'No_reset'
        load(fullfile(load_path,load_name),'NoReset');
        Parameters=NoReset.Mean;
    case 'Reset_no_bias_prev_cc'
        load_name=sprintf('%s_Reset_no_bias_prev_cc_VBMC',monkey);
        load(fullfile(load_path,load_name),'Control_model')
        Parameters=Control_model.Mean;
        Parameters.Reset=Parameters.Control_model; %error when fitting it
        clear Control_model
end

clear Reset NoReset

%% get the session list
switch monkey
    case 'Beaker'
        N_tt=8;
        date_list={'181019','181025','181029','181030','181106','181109','181110','181115'};
    case 'Scooter'
        N_tt=9;
        date_list={'181102','181103','181107','181108','181112','181113','181116','181117','181119'};
end

%% load data to analyse

data=load_monkey_data_continuous(fsroot,monkey,N_channels,100);

y = data.choices;

u = [ 0, data.chosen_color(:,1:end-1) ;  % previous choice
    nan, data.payoff(:,1:end-1) ;
    data.options_color;
    data.locations ;
    data.PopLoc ;
    data.PopSize
    data.IsPrev] ; % options offered (no last one)

for i=2:size(u,2)
    if u(11,i)==0
        u(1,i)=0;
        u(2,i)=nan;
    end
end

cd('/Volumes/buschman/Users/Caroline/Learning_attentional_templates/Behavioral_Analysis/Figure1')

%% each session

for n_tt=1:N_tt
    
    ind_this_sess=(data.Sess_number==n_tt);
    this_y=y(:,ind_this_sess);
    this_u=u(:,ind_this_sess);
    this_block=data.X_data(ind_this_sess,5);
    switch switch_type
        case 'Reset'
            Model_predictions = get_QlearningFuncApproxReset_model_RW_VBMC(Parameters,this_u,data,this_block);
        case 'No_reset'
            Model_predictions = get_QlearningFuncApproxNoReset_model_RW_VBMC(Parameters,this_u,data);
        case 'Reset_no_bias_prev_cc'
             Model_predictions = get_Control_no_bias_prev_cc_VBMC(Parameters,this_u,data,this_block);
    end
    Model_predictions.behavior.N_trials=sum(ind_this_sess);
    Model_predictions.behavior.Monkey_accuracy=sum(this_y(1,:)==1)/sum(ind_this_sess);
    Model_predictions.behavior.Best_chosen=this_y(1,:);
    Model_predictions.behavior.SecondBest_chosen=this_y(2,:);
    Model_predictions.behavior.Worst_chosen=this_y(3,:);
    Model_predictions.model_input.u=this_u;
    Model_predictions.model_input.block=this_block;
    Model_predictions.model_input.in=data.X_data;
    
    save_name=sprintf('%s_%s_%d_channels_VBMC',switch_type,model_type,N_channels);
    save_this_session_path=fullfile(save_path,monkey,date_list{n_tt});
    mkdir(save_this_session_path);
    save(fullfile(save_this_session_path,save_name),'Model_predictions');
    
    clear Model_predictions
end

%% all the session together
switch switch_type
    case 'Reset'
        Model_predictions = get_QlearningFuncApproxReset_model_RW_VBMC(Parameters,u,data,data.X_data(:,5));
    case 'No_reset'
        Model_predictions = get_QlearningFuncApproxNoReset_model_RW_VBMC(Parameters,u,data);
    case 'Reset_no_bias_prev_cc'
        Model_predictions = get_Control_no_bias_prev_cc_VBMC(Parameters,u,data,data.X_data(:,5));
end
Model_predictions.behavior.N_trials=length(data.Sess_number);
Model_predictions.behavior.Monkey_accuracy=sum(y(1,:)==1)/length(data.Sess_number);
Model_predictions.behavior.Best_chosen=y(1,:);
Model_predictions.behavior.SecondBest_chosen=y(2,:);
Model_predictions.behavior.Worst_chosen=y(3,:);
Model_predictions.model_input.u=u;
Model_predictions.model_input.block=data.X_data(:,5);
Model_predictions.model_input.in=data.X_data;

save_name=sprintf('All_sessions_%s_%s_%d_channels_VBMC',switch_type,model_type,N_channels);
save_this_session_path=fullfile(save_path,monkey);
mkdir(save_this_session_path);
save(fullfile(save_this_session_path,save_name),'Model_predictions');


