function []=doExecuteMasterScriptResetRWVBMC(fsroot,monkey,N_channels,analysis)
%% Set parameters

PLOT=0;
N_bins=100;

task='Learning_Attentional_Templates';

save_path=fullfile(fsroot,'Projects', task, 'Analysed_data',sprintf('%s',monkey));
save_name=sprintf('%s_%s_surprise_RW_%d_channels_VBMC',analysis,monkey,N_channels);

%% load data to analyse

data=load_monkey_data_continuous(fsroot,monkey,N_channels,N_bins);

%% invert the models and save

block_sw_ind=find(data.X_data(:,5)==1);

switch analysis
    case 'Reset'
        [Reset.u, Reset.posterior, Reset.out] = QlearningFuncApproxReset_RW_VBMC (data);
        Reset.out.y=[data.Best_chosen;data.SecondBest_chosen;data.Worst_chosen];
        for i=1:length(block_sw_ind)
            for l=1:3
                Reset.Choice_block(:,i,l)=Reset.out.y(l,block_sw_ind(i):block_sw_ind(i)+80);
                Reset.Prediction_block(:,i,l)=Reset.out.suffStat.gx(l,block_sw_ind(i):block_sw_ind(i)+80);
            end
        end
    case 'No_reset'
        [NoReset.u, NoReset.posterior, NoReset.out] = QlearningFuncApproxNoReset_RW_VBMC (data);
        NoReset.out.y=[data.Best_chosen;data.SecondBest_chosen;data.Worst_chosen];
        for i=1:length(block_sw_ind)
            for l=1:3
                NoReset.Choice_block(:,i,l)=NoReset.out.y(l,block_sw_ind(i):block_sw_ind(i)+80);
                NoReset.Prediction_block(:,i,l)=NoReset.out.suffStat.gx(l,block_sw_ind(i):block_sw_ind(i)+80);
            end
        end
        
end


if PLOT==1
    
    %% compare the models
    
    figure
    subplot(1,3,1)
    bar(Reset.out.F-NoReset.out.F)
    title('LME')
    subplot(1,3,2)
    bar(Reset.out.fit.BIC-NoReset.out.fit.BIC)
    title('BIC')
    subplot(1,3,3)
    bar(Reset.out.fit.AIC-NoReset.out.fit.AIC)
    title('AIC')
    
    figure
    for l=1:3
        subplot(1,3,l)
        
        plot(mean(Reset.Choice_block(:,:,l),2),'Color','k','LineWidth',2)
        hold on
        plot(mean(Reset.Prediction_block(:,:,l),2),'Color','r','LineWidth',2)
        hold on
        plot(mean(NoReset.Prediction_block(:,:,l),2),'Color','m','LineWidth',2)
        hold on
        
        shadedErrorBar([],mean(Reset.Choice_block(:,:,l),2),std(Reset.Choice_block(:,:,l),0,2)./sqrt(length(block_sw_ind)-1),{'color','k'},2)
        hold on
        shadedErrorBar([],mean(Reset.Prediction_block(:,:,l),2),std(Reset.Prediction_block(:,:,l),0,2)./sqrt(length(block_sw_ind)-1),{'color','r'},2)
        hold on
        shadedErrorBar([],mean(NoReset.Prediction_block(:,:,l),2),std(NoReset.Prediction_block(:,:,l),0,2)./sqrt(length(block_sw_ind)-1),{'color','m'},2)
        legend({'data','reset','no reset'})
        
    end
    
end

%% rename the prameters

switch analysis
    case 'Reset'
        
        % observation parameters
        Reset.Mean.B_location(1)=0;
        Reset.Var.B_location(1)=0;
        for l=1:3
            Reset.Mean.B_location(l+1)=Reset.posterior.mean(:,l+2);
            Reset.Var.B_location(l+1)=Reset.posterior.var(:,l+2);
        end
        for l=1:2
            Reset.Mean.B_popout(l)=Reset.posterior.mean(:,l+5);
            Reset.Var.B_popout(l)=Reset.posterior.var(:,l+5);
        end
        Reset.Mean.B_previous_color=Reset.posterior.mean(:,8);
        Reset.Var.B_previous_color=Reset.posterior.var(:,8);
        Reset.Mean.color_pref=Reset.posterior.mean(:,11);
        Reset.Var.color_pref=Reset.posterior.var(:,11);
        Reset.Mean.B_color_pref=Reset.posterior.mean(:,12);
        Reset.Var.B_color_pref=Reset.posterior.var(:,12);
        
        %Evolution parameters
        Reset.Mean.BF_kappa=Reset.posterior.mean(:,1);
        Reset.Var.BF_kappa=Reset.posterior.var(:,1);
        Reset.Mean.learning_rate=Reset.posterior.mean(:,2);
        Reset.Var.learning_rate=Reset.posterior.var(:,2);
        
        Reset.Mean.Reset=Reset.posterior.mean(:,9);
        Reset.Var.Reset=Reset.posterior.var(:,9);
        Reset.Mean.Volatility=Reset.posterior.mean(:,10);
        Reset.Var.Volatility=Reset.posterior.var(:,10);
        
        save(fullfile(save_path,save_name),'Reset','-v7.3');
        
    case 'No_reset'
        
        % observation parameters
        NoReset.Mean.B_location(1)=0;
        NoReset.Var.B_location(1)=0;
        for l=1:3
            NoReset.Mean.B_location(l+1)=NoReset.posterior.mean(:,l+2);
            NoReset.Var.B_location(l+1)=NoReset.posterior.var(:,l+2);
        end
        for l=1:2
            NoReset.Mean.B_popout(l)=NoReset.posterior.mean(:,l+5);
            NoReset.Var.B_popout(l)=NoReset.posterior.var(:,l+5);
        end
        NoReset.Mean.B_previous_color=NoReset.posterior.mean(:,8);
        NoReset.Var.B_previous_color=NoReset.posterior.var(:,8);
        NoReset.Mean.color_pref=NoReset.posterior.mean(:,9);
        NoReset.Var.color_pref=NoReset.posterior.var(:,9);
        NoReset.Mean.B_color_pref=NoReset.posterior.mean(:,10);
        NoReset.Var.B_color_pref=NoReset.posterior.var(:,10);
        
        %Evolution parameters
        NoReset.Mean.BF_kappa=NoReset.posterior.mean(:,1);
        NoReset.Var.BF_kappa=NoReset.posterior.var(:,1);
        NoReset.Mean.learning_rate=NoReset.posterior.mean(:,2);
        NoReset.Var.learning_rate=NoReset.posterior.var(:,2);
        
        
        save(fullfile(save_path,save_name),'NoReset','-v7.3');
end
