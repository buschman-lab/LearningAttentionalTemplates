%plot_corss_val_VBA

clear all

fsroot='/Volumes/buschman';
task='Learning_Attentional_Templates';

Channels_list=[3 4 5 6 7 8];

%% Beaker

monkey='Beaker';

for n=1:length(Channels_list)
    
    N_channels=Channels_list(n);
    
    
    task='Learning_attentional_templates';
    
    save_path=fullfile(fsroot,'Projects', task, 'Analysed_data',sprintf('%s',monkey));
    save_name=sprintf('No_reset_%s_surprise_precision_%d_channels_VBMC',monkey,N_channels);
    load(fullfile(save_path,save_name),'NoReset')
    save_name=sprintf('Reset_%s_surprise_precision_%d_channels_VBMC',monkey,N_channels);
    load(fullfile(save_path,save_name),'Reset')
    
    Beaker.Reset.out(n).BIC=Reset.out.fit.BIC;
    Beaker.NoReset.out(n).BIC=NoReset.out.fit.BIC;
    Beaker.Reset.out(n).AIC=Reset.out.fit.AIC;
    Beaker.NoReset.out(n).AIC=NoReset.out.fit.AIC;
    Beaker.Reset.out(n).F=Reset.posterior.F;
    Beaker.NoReset.out(n).F=NoReset.posterior.F;
    Beaker.Reset.out(n).F_std=Reset.posterior.F_std;
    Beaker.NoReset.out(n).F_std=NoReset.posterior.F_std;
    Beaker.Reset.Mean(n)=Reset.Mean;
    Beaker.NoReset.Mean(n)=NoReset.Mean;
    Beaker.Reset.convergence(n)=Reset.posterior.exitflag;
    Beaker.NoReset.convergence(n)=NoReset.posterior.exitflag;
    
    clear Reset NoReset
    
    save_path=fullfile(fsroot,'Projects', task, 'Analysed_data',sprintf('%s',monkey));
    save_name=sprintf('No_reset_%s_surprise_RW_%d_channels_VBMC',monkey,N_channels);
    load(fullfile(save_path,save_name),'NoReset')
    save_name=sprintf('Reset_%s_surprise_RW_%d_channels_VBMC',monkey,N_channels);
    load(fullfile(save_path,save_name),'Reset')
       
    Beaker.Reset_RW.out(n).BIC=Reset.out.fit.BIC;
    Beaker.NoReset_RW.out(n).BIC=NoReset.out.fit.BIC;
    Beaker.Reset_RW.out(n).AIC=Reset.out.fit.AIC;
    Beaker.NoReset_RW.out(n).AIC=NoReset.out.fit.AIC;
    Beaker.Reset_RW.out(n).F=Reset.posterior.F;
    Beaker.NoReset_RW.out(n).F=NoReset.posterior.F;
    Beaker.Reset_RW.out(n).F_std=Reset.posterior.F_std;
    Beaker.NoReset_RW.out(n).F_std=NoReset.posterior.F_std;
    Beaker.Reset_RW.Mean(n)=Reset.Mean;
    Beaker.NoReset_RW.Mean(n)=NoReset.Mean;
    Beaker.Reset_RW.convergence(n)=Reset.posterior.exitflag;
    Beaker.NoReset_RW.convergence(n)=NoReset.posterior.exitflag;
    
    clear Reset NoReset
    
    save_path=fullfile(fsroot,'Projects', task, 'Analysed_data',sprintf('%s',monkey));
    save_name=sprintf('No_reset_%s_surprise_PH_%d_channels_VBMC',monkey,N_channels);
    load(fullfile(save_path,save_name),'NoReset')
    save_name=sprintf('Reset_%s_surprise_PH_%d_channels_VBMC',monkey,N_channels);
    load(fullfile(save_path,save_name),'Reset')
    
    Beaker.Reset_PH.out(n).BIC=Reset.out.fit.BIC;
    Beaker.NoReset_PH.out(n).BIC=NoReset.out.fit.BIC;
    Beaker.Reset_PH.out(n).AIC=Reset.out.fit.AIC;
    Beaker.NoReset_PH.out(n).AIC=NoReset.out.fit.AIC;
    Beaker.Reset_PH.out(n).F=Reset.posterior.F;
    Beaker.NoReset_PH.out(n).F=NoReset.posterior.F;
    Beaker.Reset_PH.out(n).F_std=Reset.posterior.F_std;
    Beaker.NoReset_PH.out(n).F_std=NoReset.posterior.F_std;
    Beaker.Reset_PH.Mean(n)=Reset.Mean;
    Beaker.NoReset_PH.Mean(n)=NoReset.Mean;
    Beaker.Reset_PH.convergence(n)=Reset.posterior.exitflag;
    Beaker.NoReset_PH.convergence(n)=NoReset.posterior.exitflag;
    
    clear Reset NoReset
    
    
end


%% Scooter

monkey='Scooter';

for n=1:length(Channels_list)
    
    N_channels=Channels_list(n);
    
    
    task='Learning_attentional_templates';
    
    save_path=fullfile(fsroot,'Projects', task, 'Analysed_data',sprintf('%s',monkey));
    save_name=sprintf('No_reset_%s_surprise_precision_%d_channels_VBMC',monkey,N_channels);
    load(fullfile(save_path,save_name),'NoReset')
    save_name=sprintf('Reset_%s_surprise_precision_%d_channels_VBMC',monkey,N_channels);
    load(fullfile(save_path,save_name),'Reset')
    
    Scooter.Reset.out(n).BIC=Reset.out.fit.BIC;
    Scooter.NoReset.out(n).BIC=NoReset.out.fit.BIC;
    Scooter.Reset.out(n).AIC=Reset.out.fit.AIC;
    Scooter.NoReset.out(n).AIC=NoReset.out.fit.AIC;
    Scooter.Reset.out(n).F=Reset.posterior.F;
    Scooter.Reset.out(n).F_std=Reset.posterior.F_std;
    Scooter.NoReset.out(n).F=NoReset.posterior.F;
    Scooter.NoReset.out(n).F_std=NoReset.posterior.F_std;
    Scooter.Reset.Mean(n)=Reset.Mean;
    Scooter.NoReset.Mean(n)=NoReset.Mean;
    Scooter.Reset.convergence(n)=Reset.posterior.exitflag;
    Scooter.NoReset.convergence(n)=NoReset.posterior.exitflag;
    
    clear Reset NoReset
    
    save_path=fullfile(fsroot,'Projects', task, 'Analysed_data',sprintf('%s',monkey));
    save_name=sprintf('No_reset_%s_surprise_RW_%d_channels_VBMC',monkey,N_channels);
    load(fullfile(save_path,save_name),'NoReset')
    save_name=sprintf('Reset_%s_surprise_RW_%d_channels_VBMC',monkey,N_channels);
    load(fullfile(save_path,save_name),'Reset')
    
    Scooter.Reset_RW.out(n).BIC=Reset.out.fit.BIC;
    Scooter.NoReset_RW.out(n).BIC=NoReset.out.fit.BIC;
    Scooter.Reset_RW.out(n).AIC=Reset.out.fit.AIC;
    Scooter.NoReset_RW.out(n).AIC=NoReset.out.fit.AIC;
    Scooter.Reset_RW.out(n).F=Reset.posterior.F;
    Scooter.Reset_RW.out(n).F_std=Reset.posterior.F_std;
    Scooter.NoReset_RW.out(n).F=NoReset.posterior.F;
    Scooter.NoReset_RW.out(n).F_std=NoReset.posterior.F_std;
    Scooter.Reset_RW.Mean(n)=Reset.Mean;
    Scooter.NoReset_RW.Mean(n)=NoReset.Mean;
    Scooter.Reset_RW.convergence(n)=Reset.posterior.exitflag;
    Scooter.NoReset_RW.convergence(n)=NoReset.posterior.exitflag;
    
    clear Reset NoReset
    
    save_path=fullfile(fsroot,'Projects', task, 'Analysed_data',sprintf('%s',monkey));
    save_name=sprintf('No_reset_%s_surprise_PH_%d_channels_VBMC',monkey,N_channels);
    load(fullfile(save_path,save_name),'NoReset')
    save_name=sprintf('Reset_%s_surprise_PH_%d_channels_VBMC',monkey,N_channels);
    load(fullfile(save_path,save_name),'Reset')
   
    Scooter.Reset_PH.out(n).BIC=Reset.out.fit.BIC;
    Scooter.NoReset_PH.out(n).BIC=NoReset.out.fit.BIC;
    Scooter.Reset_PH.out(n).AIC=Reset.out.fit.AIC;
    Scooter.NoReset_PH.out(n).AIC=NoReset.out.fit.AIC;
    Scooter.Reset_PH.out(n).F=Reset.posterior.F;
    Scooter.Reset_PH.out(n).F_std=Reset.posterior.F_std;
    Scooter.NoReset_PH.out(n).F=NoReset.posterior.F;
    Scooter.NoReset_PH.out(n).F_std=NoReset.posterior.F_std;
    Scooter.Reset_PH.Mean(n)=Reset.Mean;
    Scooter.NoReset_PH.Mean(n)=NoReset.Mean;
    Scooter.Reset_PH.convergence(n)=Reset.posterior.exitflag;
    Scooter.NoReset_PH.convergence(n)=NoReset.posterior.exitflag;
    
    clear Reset NoReset
    
end

%%

for n=1:length(Channels_list)
    Beaker_LME_diff(n,1)=Beaker.NoReset_RW.out(n).F-Beaker.NoReset_RW.out(1).F;
    Beaker_LME_diff(n,2)=Beaker.Reset_RW.out(n).F-Beaker.NoReset_RW.out(1).F;
    Beaker_LME_diff_std(n,1)=Beaker.NoReset_RW.out(n).F_std;
    Beaker_LME_diff_std(n,2)=Beaker.Reset.out(n).F_std;
    
    Scooter_LME_diff(n,1)=Scooter.NoReset_RW.out(n).F-Scooter.NoReset_RW.out(1).F;
    Scooter_LME_diff(n,2)=Scooter.Reset_RW.out(n).F-Scooter.NoReset_RW.out(1).F;
    Scooter_LME_diff_std(n,1)=Beaker.NoReset_RW.out(n).F_std;
    Scooter_LME_diff_std(n,2)=Beaker.Reset_RW.out(n).F_std;
    
    Beaker_LME_diff(n,3)=Beaker.NoReset.out(n).F-Beaker.NoReset_RW.out(1).F;
    Beaker_LME_diff(n,4)=Beaker.Reset.out(n).F-Beaker.NoReset_RW.out(1).F;
    Beaker_LME_diff_std(n,3)=Beaker.NoReset.out(n).F_std;
    Beaker_LME_diff_std(n,4)=Beaker.Reset.out(n).F_std;
    
    Scooter_LME_diff(n,3)=Scooter.NoReset.out(n).F-Scooter.NoReset_RW.out(1).F;
    Scooter_LME_diff(n,4)=Scooter.Reset.out(n).F-Scooter.NoReset_RW.out(1).F;
    Scooter_LME_diff_std(n,3)=Beaker.NoReset.out(n).F_std;
    Scooter_LME_diff_std(n,4)=Beaker.Reset.out(n).F_std;
    
    Beaker_LME_diff(n,5)=Beaker.NoReset_PH.out(n).F-Beaker.NoReset_RW.out(1).F;
    Beaker_LME_diff(n,6)=Beaker.Reset_PH.out(n).F-Beaker.NoReset.out(1).F;
    Beaker_LME_diff_std(n,5)=Beaker.NoReset_PH.out(n).F_std;
    Beaker_LME_diff_std(n,6)=Beaker.Reset_PH.out(n).F_std;
    
    Scooter_LME_diff(n,5)=Scooter.NoReset_PH.out(n).F-Scooter.NoReset_RW.out(1).F;
    Scooter_LME_diff(n,6)=Scooter.Reset_PH.out(n).F-Scooter.NoReset.out(1).F;
    Scooter_LME_diff_std(n,5)=Beaker.NoReset_PH.out(n).F_std;
    Scooter_LME_diff_std(n,6)=Beaker.Reset_PH.out(n).F_std;
end

for n=1:length(Channels_list)
    Beaker_AIC_diff(n,1)=Beaker.NoReset_RW.out(n).AIC-Beaker.NoReset_RW.out(1).AIC;
    Beaker_AIC_diff(n,2)=Beaker.Reset_RW.out(n).AIC-Beaker.NoReset_RW.out(1).AIC;
    
    Scooter_AIC_diff(n,1)=Scooter.NoReset_RW.out(n).AIC-Scooter.NoReset_RW.out(1).AIC;
    Scooter_AIC_diff(n,2)=Scooter.Reset_RW.out(n).AIC-Scooter.NoReset_RW.out(1).AIC;
    
    Beaker_AIC_diff(n,3)=Beaker.NoReset.out(n).AIC-Beaker.NoReset_RW.out(1).AIC;
    Beaker_AIC_diff(n,4)=Beaker.Reset.out(n).AIC-Beaker.NoReset_RW.out(1).AIC;
    
    Scooter_AIC_diff(n,3)=Scooter.NoReset.out(n).AIC-Scooter.NoReset_RW.out(1).AIC;
    Scooter_AIC_diff(n,4)=Scooter.Reset.out(n).AIC-Scooter.NoReset_RW.out(1).AIC;
    
    Beaker_AIC_diff(n,5)=Beaker.NoReset_PH.out(n).AIC-Beaker.NoReset_RW.out(1).AIC;
    Beaker_AIC_diff(n,6)=Beaker.Reset_PH.out(n).AIC-Beaker.NoReset_RW.out(1).AIC;
    
    Scooter_AIC_diff(n,5)=Scooter.NoReset_PH.out(n).AIC-Scooter.NoReset_RW.out(1).AIC;
    Scooter_AIC_diff(n,6)=Scooter.Reset_PH.out(n).AIC-Scooter.NoReset_RW.out(1).AIC;
end

for n=1:length(Channels_list)
    Beaker_BIC_diff(n,1)=Beaker.NoReset_RW.out(n).BIC-Beaker.NoReset_RW.out(1).BIC;
    Beaker_BIC_diff(n,2)=Beaker.Reset_RW.out(n).BIC-Beaker.NoReset_RW.out(1).BIC;
    
    Scooter_BIC_diff(n,1)=Scooter.NoReset_RW.out(n).BIC-Scooter.NoReset_RW.out(1).BIC;
    Scooter_BIC_diff(n,2)=Scooter.Reset_RW.out(n).BIC-Scooter.NoReset_RW.out(1).BIC;
    
    Beaker_BIC_diff(n,3)=Beaker.NoReset.out(n).BIC-Beaker.NoReset_RW.out(1).BIC;
    Beaker_BIC_diff(n,4)=Beaker.Reset_RW.out(n).BIC-Beaker.NoReset_RW.out(1).BIC;
    
    Scooter_BIC_diff(n,3)=Scooter.NoReset.out(n).BIC-Scooter.NoReset_RW.out(1).BIC;
    Scooter_BIC_diff(n,4)=Scooter.Reset.out(n).BIC-Scooter.NoReset_RW.out(1).BIC;

    Beaker_BIC_diff(n,5)=Beaker.NoReset_PH.out(n).BIC-Beaker.NoReset_RW.out(1).BIC;
    Beaker_BIC_diff(n,6)=Beaker.Reset_PH.out(n).BIC-Beaker.NoReset_RW.out(1).BIC;
    
    Scooter_BIC_diff(n,5)=Scooter.NoReset_PH.out(n).BIC-Scooter.NoReset_RW.out(1).BIC;
    Scooter_BIC_diff(n,6)=Scooter.Reset_PH.out(n).BIC-Scooter.NoReset_RW.out(1).BIC;

end

%%

figure
subplot(2,3,1)
for n=1:length(Channels_list)
    plot(Beaker_LME_diff(:,1),'Color',[0.5 0.5 0.5],'LineWidth',2)
    hold on
    plot(Beaker_LME_diff(:,2),'Color',[0.75 0.75 0.75],'LineWidth',2)
    hold on
    plot(Beaker_LME_diff(:,1)-Beaker_LME_diff_std(:,1),'--','Color',[0.5 0.5 0.5],'LineWidth',2)
    hold on
    plot(Beaker_LME_diff(:,1)+Beaker_LME_diff_std(:,1),'--','Color',[0.5 0.5 0.5],'LineWidth',2)
    hold on
    plot(Beaker_LME_diff(:,2)-Beaker_LME_diff_std(:,2),'--','Color',[0.75 0.75 0.75],'LineWidth',2)
    hold on
    plot(Beaker_LME_diff(:,2)+Beaker_LME_diff_std(:,2),'--','Color',[0.75 0.75 0.75],'LineWidth',2)
    set(gca,'XTick',(0:1:length(Channels_list)+1));
    set(gca,'XTickLabel',{'','3','4','5','6','7','8',''})
    xlabel('# Basis functions')
    ylabel('∆Log model evidence')
    title('Beaker')
    legend({'No reset','Reset'})
    box off
    xlim([0 length(Channels_list)+1])
end
subplot(2,3,4)
for n=1:length(Channels_list)
    plot(Scooter_LME_diff(:,1),'Color',[0.5 0.5 0.5],'LineWidth',2)
    hold on
    plot(Scooter_LME_diff(:,2),'Color',[0.75 0.75 0.75],'LineWidth',2)
    hold on
    plot(Scooter_LME_diff(:,1)-Scooter_LME_diff_std(:,1),'--','Color',[0.5 0.5 0.5],'LineWidth',2)
    hold on
    plot(Scooter_LME_diff(:,1)+Scooter_LME_diff_std(:,1),'--','Color',[0.5 0.5 0.5],'LineWidth',2)
    hold on
    plot(Scooter_LME_diff(:,2)-Scooter_LME_diff_std(:,2),'--','Color',[0.75 0.75 0.75],'LineWidth',2)
    hold on
    plot(Scooter_LME_diff(:,2)+Scooter_LME_diff_std(:,2),'--','Color',[0.75 0.75 0.75],'LineWidth',2)
    set(gca,'XTick',(0:1:length(Channels_list)+1));
    set(gca,'XTickLabel',{'','3','4','5','6','7','8',''})
    xlim([0 length(Channels_list)+1])
    xlabel('# Basis functions')
    ylabel('∆Log model evidence')
    title('Scooter')
    legend({'No reset','Reset'})
    box off
end

subplot(2,3,2)
for n=1:length(Channels_list)
    plot(Beaker_AIC_diff(:,1),'Color',[0.5 0.5 0.5],'LineWidth',2)
    hold on
    plot(Beaker_AIC_diff(:,2),'Color',[0.75 0.75 0.75],'LineWidth',2)
    set(gca,'XTick',(0:1:length(Channels_list)+1));
    set(gca,'XTickLabel',{'','3','4','5','6','7','8',''})
    xlabel('# Basis functions')
    ylabel('∆AIC')
    title('Beaker')
    legend({'No reset','Reset'})
    box off
    xlim([0 length(Channels_list)+1])
end
subplot(2,3,5)
for n=1:length(Channels_list)
    plot(Scooter_AIC_diff(:,1),'Color',[0.5 0.5 0.5],'LineWidth',2)
    hold on
    plot(Scooter_AIC_diff(:,2),'Color',[0.75 0.75 0.75],'LineWidth',2)
    set(gca,'XTick',(0:1:length(Channels_list)+1));
    set(gca,'XTickLabel',{'','3','4','5','6','7','8',''})
    xlabel('# Basis functions')
    ylabel('∆AIC')
    title('Scooter')
    legend({'No reset','Reset'})
    box off
    xlim([0 length(Channels_list)+1])
end

subplot(2,3,3)
for n=1:length(Channels_list)
    plot(Beaker_BIC_diff(:,1),'Color',[0.5 0.5 0.5],'LineWidth',2)
    hold on
    plot(Beaker_BIC_diff(:,2),'Color',[0.75 0.75 0.75],'LineWidth',2)
    set(gca,'XTick',(0:1:length(Channels_list)+1));
    set(gca,'XTickLabel',{'','3','4','5','6','7','8',''})
    xlabel('# Basis functions')
    ylabel('∆BIC')
    title('Beaker')
    legend({'No reset','Reset'})
    box off
    xlim([0 length(Channels_list)+1])
end
subplot(2,3,6)
for n=1:length(Channels_list)
    plot(Scooter_BIC_diff(:,1),'Color',[0.5 0.5 0.5],'LineWidth',2)
    hold on
    plot(Scooter_BIC_diff(:,2),'Color',[0.75 0.75 0.75],'LineWidth',2)
    set(gca,'XTick',(0:1:length(Channels_list)+1));
    set(gca,'XTickLabel',{'','3','4','5','6','7','8',''})
    xlabel('# Basis functions')
    ylabel('∆BIC')
    title('Scooter')
    legend({'No reset','Reset'})
    box off
    xlim([0 length(Channels_list)+1])
end

%%

figure
subplot(3,2,1)
for n=1:length(Channels_list)
    plot(Beaker_LME_diff(:,1),'Color',[0.5 0.5 0.5],'LineWidth',2)
    hold on
    plot(Beaker_LME_diff(:,2),'Color',[0.75 0.75 0.75],'LineWidth',2)
    hold on
    plot(Beaker_LME_diff(:,3),'Color',[0.25 0.5 0.5],'LineWidth',2)
    hold on
    plot(Beaker_LME_diff(:,4),'Color',[0.5 0.75 0.75],'LineWidth',2)
    hold on
    plot(Beaker_LME_diff(:,5),'Color',[0.5 0.25 0.5],'LineWidth',2)
    hold on
    plot(Beaker_LME_diff(:,6),'Color',[0.75 0.5 0.75],'LineWidth',2)
    set(gca,'XTick',(0:1:length(Channels_list)+1));
    set(gca,'XTickLabel',{'','3','4','5','6','7','8',''})
    legend({'No reset','Reset','Prec no reset','Prec reset','PH no reset', 'PH reset'})
    xlabel('# Basis functions')
    ylabel('∆Log model evidence')
    title('Beaker')
    box off
    
end
subplot(3,2,2)
for n=1:length(Channels_list)
    plot(Scooter_LME_diff(:,1),'Color',[0.5 0.5 0.5],'LineWidth',2)
    hold on
    plot(Scooter_LME_diff(:,2),'Color',[0.75 0.75 0.75],'LineWidth',2)
    hold on
    plot(Scooter_LME_diff(:,3),'Color',[0.25 0.5 0.5],'LineWidth',2)
    hold on
    plot(Scooter_LME_diff(:,4),'Color',[0.5 0.75 0.75],'LineWidth',2)
    hold on
    plot(Scooter_LME_diff(:,5),'Color',[0.5 0.25 0.5],'LineWidth',2)
    hold on
    plot(Scooter_LME_diff(:,6),'Color',[0.75 0.5 0.75],'LineWidth',2)
    set(gca,'XTick',(0:1:length(Channels_list)+1));
    set(gca,'XTickLabel',{'','3','4','5','6','7','8',''})
    xlabel('# Basis functions')
    ylabel('∆Log model evidence')
    title('Scooter')
    legend({'No reset','Reset','Prec no reset','Prec reset','PH no reset', 'PH reset'})
    box off
end

subplot(3,2,3)
for n=1:length(Channels_list)
    plot(Beaker_AIC_diff(:,1),'Color',[0.5 0.5 0.5],'LineWidth',2)
    hold on
    plot(Beaker_AIC_diff(:,2),'Color',[0.75 0.75 0.75],'LineWidth',2)
    hold on
    plot(Beaker_AIC_diff(:,3),'Color',[0.25 0.5 0.5],'LineWidth',2)
    hold on
    plot(Beaker_AIC_diff(:,4),'Color',[0.5 0.75 0.75],'LineWidth',2)
    hold on
    plot(Beaker_AIC_diff(:,5),'Color',[0.5 0.25 0.5],'LineWidth',2)
    hold on
    plot(Beaker_AIC_diff(:,6),'Color',[0.75 0.5 0.75],'LineWidth',2)
    set(gca,'XTick',(0:1:length(Channels_list)+1));
    set(gca,'XTickLabel',{'','3','4','5','6','7','8',''})
    xlabel('# Basis functions')
    ylabel('∆AIC')
    title('Beaker')
    legend({'No reset','Reset','Prec no reset','Prec reset','PH no reset', 'PH reset'})
    box off
end
subplot(3,2,4)
for n=1:length(Channels_list)
    plot(Scooter_AIC_diff(:,1),'Color',[0.5 0.5 0.5],'LineWidth',2)
    hold on
    plot(Scooter_AIC_diff(:,2),'Color',[0.75 0.75 0.75],'LineWidth',2)
    hold on
    plot(Scooter_AIC_diff(:,3),'Color',[0.25 0.5 0.5],'LineWidth',2)
    hold on
    plot(Scooter_AIC_diff(:,4),'Color',[0.5 0.75 0.75],'LineWidth',2)
    hold on
    plot(Scooter_AIC_diff(:,5),'Color',[0.5 0.25 0.5],'LineWidth',2)
    hold on
    plot(Scooter_AIC_diff(:,6),'Color',[0.75 0.5 0.75],'LineWidth',2)
    set(gca,'XTick',(0:1:length(Channels_list)+1));
    set(gca,'XTickLabel',{'','3','4','5','6','7','8',''})
    xlabel('# Basis functions')
    ylabel('∆AIC')
    title('Scooter')
    legend({'No reset','Reset','Prec no reset','Prec reset','PH no reset', 'PH reset'})
    box off
end

subplot(3,2,5)
for n=1:length(Channels_list)
    plot(Beaker_BIC_diff(:,1),'Color',[0.5 0.5 0.5],'LineWidth',2)
    hold on
    plot(Beaker_BIC_diff(:,2),'Color',[0.75 0.75 0.75],'LineWidth',2)
    hold on
    plot(Beaker_BIC_diff(:,3),'Color',[0.25 0.5 0.5],'LineWidth',2)
    hold on
    plot(Beaker_BIC_diff(:,4),'Color',[0.5 0.75 0.75],'LineWidth',2)
    hold on
    plot(Beaker_BIC_diff(:,5),'Color',[0.5 0.25 0.5],'LineWidth',2)
    hold on
    plot(Beaker_BIC_diff(:,6),'Color',[0.75 0.5 0.75],'LineWidth',2)
    set(gca,'XTick',(0:1:length(Channels_list)+1));
    set(gca,'XTickLabel',{'','3','4','5','6','7','8',''})
    xlabel('# Basis functions')
    ylabel('∆BIC')
    title('Beaker')
    legend({'No reset','Reset','Prec no reset','Prec reset','PH no reset', 'PH reset'})
    box off
end
subplot(3,2,6)
for n=1:length(Channels_list)
    plot(Scooter_BIC_diff(:,1),'Color',[0.5 0.5 0.5],'LineWidth',2)
    hold on
    plot(Scooter_BIC_diff(:,2),'Color',[0.75 0.75 0.75],'LineWidth',2)
    hold on
    plot(Scooter_BIC_diff(:,3),'Color',[0.25 0.5 0.5],'LineWidth',2)
    hold on
    plot(Scooter_BIC_diff(:,4),'Color',[0.5 0.75 0.75],'LineWidth',2)
    hold on
    plot(Scooter_BIC_diff(:,5),'Color',[0.5 0.25 0.5],'LineWidth',2)
    hold on
    plot(Scooter_BIC_diff(:,6),'Color',[0.75 0.5 0.75],'LineWidth',2)
    set(gca,'XTick',(0:1:length(Channels_list)+1));
    set(gca,'XTickLabel',{'','3','4','5','6','7','8',''})
    xlabel('# Basis functions')
    ylabel('∆BIC')
    title('Scooter')
    legend({'No reset','Reset','Prec no reset','Prec reset','PH no reset', 'PH reset'})
    box off
end




%%
save('/Volumes/buschman/Users/Caroline/Learning_attentional_templates/Behavioral_Analysis/All_models_results','Scooter','Beaker')

