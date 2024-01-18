function []=doExecuteMasterScript_WSLS_VBMC(fsroot,monkey,analysis)
%% Set parameters

PLOT=1;
N_bins=100;

task='Learning_Attentional_Templates';

save_path=fullfile(fsroot,'Projects', task, 'Analysed_data',sprintf('%s',monkey));
switch analysis
    case 'WSLS'
        save_name=sprintf('%s_WSLS_VBMC',monkey);
    case 'WSLS_flip'
        save_name=sprintf('%s_WSLS_flip_VBMC',monkey);
    case 'No_value'
        save_name=sprintf('%s_no_value_VBMC',monkey);
    case 'No_bias_previous_color'
        save_name=sprintf('%s_Reset_no_bias_prev_cc_VBMC',monkey);
    case 'Value_only'
        save_name=sprintf('%s_Reset_value_only_VBMC',monkey);
end

%% load data to analyse

N_channels=6;
data=load_monkey_data_continuous(fsroot,monkey,N_channels,N_bins);

cd(fullfile(fsroot,'/Users/Caroline/Learning_attentional_templates/Behavioral_Analysis/Figure1'))

%% invert the models and save

block_sw_ind=find(data.X_data(:,5)==1);
switch analysis
    case 'WSLS'
        [WSLS.u, WSLS.posterior, WSLS.out] = WSLS_VBMC (data);
    case 'WSLS_flip'
        [WSLS.u, WSLS.posterior, WSLS.out] = WSLS_flip_VBMC (data);
    case 'No_value'
        [Control_model.u, Control_model.posterior, Control_model.out] = No_value_VBMC (data);
    case 'No_bias_previous_color'
        [Control_model.u, Control_model.posterior, Control_model.out] = Reset_no_bias_prev_cc_VBMC (data);
    case 'Value_only'
        [Control_model.u, Control_model.posterior, Control_model.out] = Reset_value_only_VBMC (data);
end

switch analysis
    case {'WSLS', 'WSLS_flip'}
        WSLS.out.y=[data.Best_chosen;data.SecondBest_chosen;data.Worst_chosen];
        for i=1:length(block_sw_ind)
            for l=1:3
                WSLS.Choice_block(:,i,l)=WSLS.out.y(l,block_sw_ind(i):block_sw_ind(i)+80);
                WSLS.Prediction_block(:,i,l)=WSLS.out.suffStat.gx(l,block_sw_ind(i):block_sw_ind(i)+80);
            end
        end
    case {'No_value', 'No_bias_previous_color','Value_only'}
        Control_model.out.y=[data.Best_chosen;data.SecondBest_chosen;data.Worst_chosen];
        for i=1:length(block_sw_ind)
            for l=1:3
                Control_model.Choice_block(:,i,l)=Control_model.out.y(l,block_sw_ind(i):block_sw_ind(i)+80);
                Control_model.Prediction_block(:,i,l)=Control_model.out.suffStat.gx(l,block_sw_ind(i):block_sw_ind(i)+80);
            end
        end
end



if PLOT==1
    
    %% compare the models
    load_name=sprintf('Reset_%s_surprise_RW_%d_channels_VBMC',monkey,6); %load the best model
    load(fullfile(save_path,load_name),'Reset')
    
    switch analysis
        case {'WSLS', 'WSLS_flip'}
            
            figure
            bar(Reset.out.fit.BIC-WSLS.out.fit.BIC)
            title('BIC')
            
            figure
            for l=1:3
                subplot(1,3,l)
                
                plot(mean(Reset.Choice_block(:,:,l),2),'Color','k','LineWidth',2)
                hold on
                plot(mean(Reset.Prediction_block(:,:,l),2),'Color','r','LineWidth',2)
                hold on
                plot(mean(WSLS.Prediction_block(:,:,l),2),'Color','m','LineWidth',2)
                hold on
                
                shadedErrorBar([],mean(Reset.Choice_block(:,:,l),2),std(Reset.Choice_block(:,:,l),0,2)./sqrt(length(block_sw_ind)-1),{'color','k'},2)
                hold on
                shadedErrorBar([],mean(Reset.Prediction_block(:,:,l),2),std(Reset.Prediction_block(:,:,l),0,2)./sqrt(length(block_sw_ind)-1),{'color','r'},2)
                hold on
                shadedErrorBar([],mean(WSLS.Prediction_block(:,:,l),2),std(WSLS.Prediction_block(:,:,l),0,2)./sqrt(length(block_sw_ind)-1),{'color','m'},2)
                legend({'data','reset','WSLS'})
                
            end
            
        case {'No_value', 'No_bias_previous_color', 'Value_only'}
            
            figure
            bar(Reset.out.fit.BIC-Control_model.out.fit.BIC)
            title('BIC')
            
            figure
            for l=1:3
                subplot(1,3,l)
                
                plot(mean(Reset.Choice_block(:,:,l),2),'Color','k','LineWidth',2)
                hold on
                plot(mean(Reset.Prediction_block(:,:,l),2),'Color','r','LineWidth',2)
                hold on
                plot(mean(Control_model.Prediction_block(:,:,l),2),'Color','m','LineWidth',2)
                hold on
                
                shadedErrorBar([],mean(Reset.Choice_block(:,:,l),2),std(Reset.Choice_block(:,:,l),0,2)./sqrt(length(block_sw_ind)-1),{'color','k'},2)
                hold on
                shadedErrorBar([],mean(Reset.Prediction_block(:,:,l),2),std(Reset.Prediction_block(:,:,l),0,2)./sqrt(length(block_sw_ind)-1),{'color','r'},2)
                hold on
                shadedErrorBar([],mean(Control_model.Prediction_block(:,:,l),2),std(Control_model.Prediction_block(:,:,l),0,2)./sqrt(length(block_sw_ind)-1),{'color','m'},2)
                legend({'data','reset','Control'})
                
            end
    end
end

%% rename the prameters

switch analysis
    case {'WSLS', 'WSLS_flip'}
        
        % observation parameters
        WSLS.Mean.B_location(1)=0;
        WSLS.Var.B_location(1)=0;
        for l=1:3
            WSLS.Mean.B_location(l+1)=WSLS.posterior.mean(:,l+1);
            WSLS.Var.B_location(l+1)=WSLS.posterior.var(:,l+1);
        end
        for l=1:2
            WSLS.Mean.B_popout(l)=WSLS.posterior.mean(:,l+4);
            WSLS.Var.B_popout(l)=WSLS.posterior.var(:,l+4);
        end
        WSLS.Mean.color_pref=WSLS.posterior.mean(:,7);
        WSLS.Var.color_pref=WSLS.posterior.var(:,7);
        WSLS.Mean.B_color_pref=WSLS.posterior.mean(:,8);
        WSLS.Var.B_color_pref=WSLS.posterior.var(:,8);
        
        %Evolution parameters
        WSLS.Mean.reward_threshold=WSLS.posterior.mean(:,1);
        WSLS.Var.reward_threshold=WSLS.posterior.var(:,1);
        
        save(fullfile(save_path,save_name),'WSLS','-v7.3');
        
    case 'No_bias_previous_color'
        
        % observation parameters
        Control_model.Mean.B_location(1)=0;
        Control_model.Var.B_location(1)=0;
        for l=1:3
            Control_model.Mean.B_location(l+1)=Control_model.posterior.mean(:,l+2);
            Control_model.Var.B_location(l+1)=Control_model.posterior.var(:,l+2);
        end
        for l=1:2
            Control_model.Mean.B_popout(l)=Control_model.posterior.mean(:,l+5);
            Control_model.Var.B_popout(l)=Control_model.posterior.var(:,l+5);
        end
        Control_model.Mean.color_pref=Control_model.posterior.mean(:,10);
        Control_model.Var.color_pref=Control_model.posterior.var(:,10);
        Control_model.Mean.B_color_pref=Control_model.posterior.mean(:,11);
        Control_model.Var.B_color_pref=Control_model.posterior.var(:,11);
        
        %Evolution parameters
        Control_model.Mean.BF_kappa=Control_model.posterior.mean(:,1);
        Control_model.Var.BF_kappa=Control_model.posterior.var(:,1);
        Control_model.Mean.learning_rate=Control_model.posterior.mean(:,2);
        Control_model.Var.learning_rate=Control_model.posterior.var(:,2);
        
        Control_model.Mean.Reset=Control_model.posterior.mean(:,8);
        Control_model.Var.Reset_model=Control_model.posterior.var(:,8);
        Control_model.Mean.Volatility=Control_model.posterior.mean(:,9);
        Control_model.Var.Volatility=Control_model.posterior.var(:,9);
        
        save(fullfile(save_path,save_name),'Control_model','-v7.3');
        
    case 'No_value' 
        
        % observation parameters
        Control_model.Mean.B_location(1)=0;
        Control_model.Var.B_location(1)=0;
        for l=1:3
            Control_model.Mean.B_location(l+1)=Control_model.posterior.mean(:,l+1);
            Control_model.Var.B_location(l+1)=Control_model.posterior.var(:,l+1);
        end
        for l=1:2
            Control_model.Mean.B_popout(l)=Control_model.posterior.mean(:,l+4);
            Control_model.Var.B_popout(l)=Control_model.posterior.var(:,l+4);
        end
        Control_model.Mean.color_pref=Control_model.posterior.mean(:,7);
        Control_model.Var.color_pref=Control_model.posterior.var(:,7);
        Control_model.Mean.B_color_pref=Control_model.posterior.mean(:,8);
        Control_model.Var.B_color_pref=Control_model.posterior.var(:,8);
        
        Control_model.Mean.B_previous_color=Control_model.posterior.mean(:,1); %this is fixed to 0
        Control_model.Var.B_previous_color=Control_model.posterior.var(:,1); %this is fixed to 0
        
        save(fullfile(save_path,save_name),'Control_model','-v7.3');
    
    case 'Value_only' 
        
        %Evolution parameters
        Control_model.Mean.BF_kappa=Control_model.posterior.mean(:,1);
        Control_model.Var.BF_kappa=Control_model.posterior.var(:,1);
        Control_model.Mean.learning_rate=Control_model.posterior.mean(:,2);
        Control_model.Var.learning_rate=Control_model.posterior.var(:,2);
        
        Control_model.Mean.Reset=Control_model.posterior.mean(:,3);
        Control_model.Var.Reset=Control_model.posterior.var(:,3);
        Control_model.Mean.Volatility=Control_model.posterior.mean(:,4);
        Control_model.Var.Volatility=Control_model.posterior.var(:,4);
        
        save(fullfile(save_path,save_name),'Control_model','-v7.3');
        
end
