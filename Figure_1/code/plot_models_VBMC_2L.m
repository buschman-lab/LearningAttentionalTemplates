%plot_corss_val_VBA

clear all

fsroot='/Volumes/buschman';
task='Learning_Attentional_Templates';

Channels_list=6;%[3 4 5 6 7 8];

%% Beaker

monkey='Beaker';

for n=1:length(Channels_list)
    
    N_channels=Channels_list(n);
        
    save_path=fullfile(fsroot,'Projects', task, 'Analysed_data',sprintf('%s',monkey));
    save_name=sprintf('Reset_%s_surprise_RW_2L_%d_channels_VBMC',monkey,N_channels);
    load(fullfile(save_path,save_name),'Reset')
    Reset_2l=Reset;
    clear Reset
    save_name=sprintf('Reset_%s_surprise_RW_%d_channels_VBMC',monkey,N_channels);
    load(fullfile(save_path,save_name),'Reset')
       
    Beaker.Reset_RW.out(n).BIC=Reset.out.fit.BIC;
    Beaker.Reset_2l_RW.out(n).BIC=Reset_2l.out.fit.BIC;
    Beaker.Reset_RW.out(n).AIC=Reset.out.fit.AIC;
    Beaker.Reset_2l_RW.out(n).AIC=Reset_2l.out.fit.AIC;
    Beaker.Reset_RW.out(n).F=Reset.posterior.F;
    Beaker.Reset_2l_RW.out(n).F=Reset_2l.posterior.F;
    Beaker.Reset_RW.out(n).F_std=Reset.posterior.F_std;
    Beaker.Reset_2l_RW.out(n).F_std=Reset_2l.posterior.F_std;
    Beaker.Reset_RW.Mean(n)=Reset.Mean;
    Beaker.Reset_2l_RW.Mean(n)=Reset_2l.Mean;
    Beaker.Reset_RW.convergence(n)=Reset.posterior.exitflag;
    Beaker.Reset_2l_RW.convergence(n)=Reset_2l.posterior.exitflag;
    
    clear Reset Reset_2l
    
    
end


% Scooter

monkey='Scooter';

for n=1:length(Channels_list)
    
    N_channels=Channels_list(n);
            
    save_path=fullfile(fsroot,'Projects', task, 'Analysed_data',sprintf('%s',monkey));
    save_name=sprintf('Reset_%s_surprise_RW_2L_%d_channels_VBMC',monkey,N_channels);
    load(fullfile(save_path,save_name),'Reset')
    Reset_2l=Reset;
    clear Reset
    save_name=sprintf('Reset_%s_surprise_RW_%d_channels_VBMC',monkey,N_channels);
    load(fullfile(save_path,save_name),'Reset')
       
    Scooter.Reset_RW.out(n).BIC=Reset.out.fit.BIC;
    Scooter.Reset_2l_RW.out(n).BIC=Reset_2l.out.fit.BIC;
    Scooter.Reset_RW.out(n).AIC=Reset.out.fit.AIC;
    Scooter.Reset_2l_RW.out(n).AIC=Reset_2l.out.fit.AIC;
    Scooter.Reset_RW.out(n).F=Reset.posterior.F;
    Scooter.Reset_2l_RW.out(n).F=Reset_2l.posterior.F;
    Scooter.Reset_RW.out(n).F_std=Reset.posterior.F_std;
    Scooter.Reset_2l_RW.out(n).F_std=Reset_2l.posterior.F_std;
    Scooter.Reset_RW.Mean(n)=Reset.Mean;
    Scooter.Reset_2l_RW.Mean(n)=Reset_2l.Mean;
    Scooter.Reset_RW.convergence(n)=Reset.posterior.exitflag;
    Scooter.Reset_2l_RW.convergence(n)=Reset_2l.posterior.exitflag;
    
    clear Reset Reset_2l
    
end

%%

for n=1:length(Channels_list)
    Beaker_LME_diff(n,1)=Beaker.Reset_2l_RW.out(n).F-Beaker.Reset_2l_RW.out(1).F;
    Beaker_LME_diff(n,2)=Beaker.Reset_RW.out(n).F-Beaker.Reset_2l_RW.out(1).F;
    Beaker_LME_diff_std(n,1)=Beaker.Reset_2l_RW.out(n).F_std;
    Beaker_LME_diff_std(n,2)=Beaker.Reset_RW.out(n).F_std;
    
    Scooter_LME_diff(n,1)=Scooter.Reset_2l_RW.out(n).F-Scooter.Reset_2l_RW.out(1).F;
    Scooter_LME_diff(n,2)=Scooter.Reset_RW.out(n).F-Scooter.Reset_2l_RW.out(1).F;
    Scooter_LME_diff_std(n,1)=Beaker.Reset_2l_RW.out(n).F_std;
    Scooter_LME_diff_std(n,2)=Beaker.Reset_RW.out(n).F_std;
    
end

for n=1:length(Channels_list)
    Beaker_AIC_diff(n,1)=Beaker.Reset_2l_RW.out(n).AIC-Beaker.Reset_2l_RW.out(1).AIC;
    Beaker_AIC_diff(n,2)=Beaker.Reset_RW.out(n).AIC-Beaker.Reset_2l_RW.out(1).AIC;
    
    Scooter_AIC_diff(n,1)=Scooter.Reset_2l_RW.out(n).AIC-Scooter.Reset_2l_RW.out(1).AIC;
    Scooter_AIC_diff(n,2)=Scooter.Reset_RW.out(n).AIC-Scooter.Reset_2l_RW.out(1).AIC;
    
end

for n=1:length(Channels_list)
    Beaker_BIC_diff(n,1)=Beaker.Reset_2l_RW.out(n).BIC-Beaker.Reset_2l_RW.out(1).BIC;
    Beaker_BIC_diff(n,2)=Beaker.Reset_RW.out(n).BIC-Beaker.Reset_2l_RW.out(1).BIC;
    
    Scooter_BIC_diff(n,1)=Scooter.Reset_2l_RW.out(n).BIC-Scooter.Reset_2l_RW.out(1).BIC;
    Scooter_BIC_diff(n,2)=Scooter.Reset_RW.out(n).BIC-Scooter.Reset_2l_RW.out(1).BIC;
end

%%

figure
subplot(1,2,1)
for n=1:length(Channels_list)
%     bar(1,Beaker_LME_diff(:,1),'FaceColor',[0.5 0.5 0.5],'LineWidth',2)
    hold on
    bar(1,Beaker_LME_diff(:,2),'FaceColor',[0.75 0.75 0.75],'LineWidth',2)
    bar(2,Scooter_LME_diff(:,2),'FaceColor',[0.75 0.75 0.75],'LineWidth',2)
%     plot(Beaker_LME_diff(:,1)-Beaker_LME_diff_std(:,1),'--','Color',[0.5 0.5 0.5],'LineWidth',2)
%     hold on
%     plot(Beaker_LME_diff(:,1)+Beaker_LME_diff_std(:,1),'--','Color',[0.5 0.5 0.5],'LineWidth',2)
%     hold on
%     plot(Beaker_LME_diff(:,2)-Beaker_LME_diff_std(:,2),'--','Color',[0.75 0.75 0.75],'LineWidth',2)
%     hold on
%     plot(Beaker_LME_diff(:,2)+Beaker_LME_diff_std(:,2),'--','Color',[0.75 0.75 0.75],'LineWidth',2)
%     set(gca,'XTick',(0:1:length(Channels_list)+1));
%     set(gca,'XTickLabel',{'','3','4','5','6','7','8',''})
%     xlabel('# Basis functions')
    ylabel('∆Log model evidence')
%     title('Beaker')
%     legend({'No reset','Reset'})
    box off
%     xlim([0 length(Channels_list)+1])
end
subplot(1,2,2)
for n=1:length(Channels_list)
        hold on
    bar(1,Beaker_BIC_diff(:,2),'FaceColor',[0.75 0.75 0.75],'LineWidth',2)
    bar(2,Scooter_BIC_diff(:,2),'FaceColor',[0.75 0.75 0.75],'LineWidth',2)

%     plot(Scooter_LME_diff(:,1),'Color',[0.5 0.5 0.5],'LineWidth',2)
%     hold on
%     plot(Scooter_LME_diff(:,2),'Color',[0.75 0.75 0.75],'LineWidth',2)
%     hold on
%     plot(Scooter_LME_diff(:,1)-Scooter_LME_diff_std(:,1),'--','Color',[0.5 0.5 0.5],'LineWidth',2)
%     hold on
%     plot(Scooter_LME_diff(:,1)+Scooter_LME_diff_std(:,1),'--','Color',[0.5 0.5 0.5],'LineWidth',2)
%     hold on
%     plot(Scooter_LME_diff(:,2)-Scooter_LME_diff_std(:,2),'--','Color',[0.75 0.75 0.75],'LineWidth',2)
%     hold on
%     plot(Scooter_LME_diff(:,2)+Scooter_LME_diff_std(:,2),'--','Color',[0.75 0.75 0.75],'LineWidth',2)
%     set(gca,'XTick',(0:1:length(Channels_list)+1));
%     set(gca,'XTickLabel',{'','3','4','5','6','7','8',''})
%     xlim([0 length(Channels_list)+1])
    ylabel('∆BIC')
%     title('Beaker')
%     legend({'No reset','Reset'})
    box off
end




%%
save('/Volumes/buschman/Users/Caroline/Learning_attentional_templates/Behavioral_Analysis/All_models_results_2l','Scooter','Beaker')

