clear all

arrayID=38;


load('colors.mat')

figure
theta=0:2*pi/100:2*pi-(2*pi/100);
for i=1:length(theta)
    rho(i)=exp(2.5*(cos(theta(i))-1));
end
polarplot(theta,rho,'LineWidth',2,'Color',colors(1,:))
hold on
% for i=1:length(theta)
% rho2(i)=exp(2.5*(cos(theta(i)-pi)-1));
% end
% polarplot(theta,rho2,'LineWidth',2,'Color',colors(20,:))
% hold on
% for i=1:length(theta)
% rho2(i)=exp(2.5*(cos(theta(i)-rule_thr)-1));
% end
% polarplot(theta,rho2,'LineWidth',2,'Color',colors(round(rule_thr/(2*pi)*length(colors)),:))
box off
ax=gca;
ax.ThetaTickLabel={};
ax.RTickLabel={};

clear rho theta

%%
N_channels=6;
theta_channel=0:2*pi/N_channels:2*pi;

kappa=1.82; %Beaker is 1.64, Scooter is 1.82

theta=0:2*pi/100:2*pi-(2*pi/100);
for j=1:N_channels
    for i=1:length(theta)
        rho(i,j)=vonmisespdf(theta(i),theta_channel(j),kappa);
    end
end

%%
figure
for j=1:N_channels
    polarplot(theta,rho(:,j),'LineWidth',2,'Color',colors(1+(j-1)*round(40/N_channels),:))
    hold on
end
box off
ax=gca;
ax.ThetaTickLabel={};
ax.RTickLabel={};
%%



fsroot='/Volumes/buschman';
event='target';
window_start_list=-600;
window_size=900;
this_time=1;

task='Learning_Attentional_Templates';
subtask='exploreexploit/Restricted_FEF/Reset_RW_model';

subsubtask=sprintf('FR_%s_%d_window_from_%d_to_%d',event,window_size,window_start_list(this_time),window_start_list(this_time)+window_size);

dirstem = fullfile('Projects',task,'Analysis','Electrophy_analysis',subtask,subsubtask);
data_path = fullfile(fsroot,dirstem);


data_name=sprintf('ID_%d.mat',arrayID);
load(fullfile(data_path,data_name),'Conditions')

belief_probability=zeros(size(Conditions.belief,2),size(Conditions.belief,1));
for i=1:size(Conditions.stim,2)
    if sum(Conditions.belief(:,i)==0)<size(Conditions.belief(:,i),1)
        belief_probability(i,:)=(Conditions.belief(:,i)-min(Conditions.belief(:,i)))/sum(Conditions.belief(:,i));
    else
        belief_probability(i,:)=Conditions.belief(:,i);
    end
end
mean_value=mean(belief_probability,1);

N_bins=100;
color_binned = 0:(2*pi/N_bins):(2*pi-(2*pi/N_bins));


%%
% for trial=2:2946
%     if max(belief_probability(trial-1,:))>0 && max(belief_probability(trial,:))>0
%     peak_prev(trial)=find(belief_probability(trial-1,:)==max(belief_probability(trial-1,:)));
%     peak_current(trial)=find(belief_probability(trial,:)==max(belief_probability(trial,:)));
%     end
% end
% figure
% plot(Conditions.RPE,abs(peak_prev-peak_current),'o')
% 
%%
% figure
% 
% for trial=705:710
% 
% subplot(1,5,trial-704)
% 
% polarplot(color_binned,exp(2.5*(cos(color_binned-Conditions.rule(trial))))/650,'-','LineWidth',2,'Color',colors(round(Conditions.rule(trial)/2/pi*40),:))
% hold on
% polarplot([Conditions.rule(trial);Conditions.rule(trial)],[0;0.025],'--','LineWidth',3,'Color',colors(round(Conditions.rule(trial)/2/pi*40),:))
% 
% polarplot([Conditions.chosen_color(trial-1);Conditions.chosen_color(trial-1)],[0;0.025],'LineWidth',3,'Color',colors(round(Conditions.chosen_color(trial-1)/2/pi*40),:))
% peak_prev=find(belief_probability(trial-1,:)==max(belief_probability(trial-1,:)));
% % polarplot(color_binned,belief_probability(trial-1,:),'--','LineWidth',2,'Color',colors(round(peak_prev/100*40),:))
% polarplot(color_binned,belief_probability(trial-1,:),'k--','LineWidth',2)
% polarplot([peak_prev/100*2*pi;peak_prev/100*2*pi],[0;0.025],'k--')
% 
% peak_current=find(belief_probability(trial,:)==max(belief_probability(trial,:)));
% % polarplot(color_binned,belief_probability(trial,:),'-','LineWidth',2,'Color',colors(round(peak_current/100*40),:))
% polarplot(color_binned,belief_probability(trial,:),'k-','LineWidth',2)
% polarplot([peak_current/100*2*pi;peak_current/100*2*pi],[0;0.025],'k-')
% title(sprintf('Reward = %0.2g\nRPE = %0.2g',Conditions.Reward(trial-1),Conditions.RPE(trial)))
% 
% box off
% ax=gca;
% ax.ThetaTickLabel={};
% ax.RTickLabel={};
% 
% end


figure

for trial=4

subplot(1,1,trial-3)

polarplot(color_binned,exp(2.5*(cos(color_binned-Conditions.rule(trial))))/100,'-','LineWidth',2,'Color',colors(round(Conditions.rule(trial)/2/pi*40),:))
hold on
polarplot([Conditions.rule(trial);Conditions.rule(trial)],[0;0.15],'--','LineWidth',3,'Color',colors(round(Conditions.rule(trial)/2/pi*40),:))

polarplot([Conditions.chosen_color(trial-1);Conditions.chosen_color(trial-1)],[0;0.15],'LineWidth',3,'Color',colors(round(Conditions.chosen_color(trial-1)/2/pi*40+1),:))
peak_prev=find(Conditions.belief(:,trial-1)==max(Conditions.belief(:,trial-1)));
% polarplot(color_binned,belief_probability(trial-1,:),'--','LineWidth',2,'Color',colors(round(peak_prev/100*40),:))
polarplot(color_binned,Conditions.belief(:,trial-1),'k--','LineWidth',2)

peak_current=find(Conditions.belief(:,trial)==max(Conditions.belief(:,trial)));
% polarplot(color_binned,belief_probability(trial,:),'-','LineWidth',2,'Color',colors(round(peak_current/100*40),:))
polarplot(color_binned,Conditions.belief(:,trial),'k-','LineWidth',2)

title(sprintf('Reward = %0.2g\nRPE = %0.2g',Conditions.Reward(trial-1),Conditions.RPE(trial)))

box off
ax=gca;
ax.ThetaTickLabel={};
ax.RTickLabel={};

end

%%

figure
subplot(1,2,1)
plot(Conditions.weight(:,3),1:6,'-k','LineWidth',2)
xlim([0 0.2])
box off
subplot(1,2,2)
plot(Conditions.weight(:,4),1:6,'-k','LineWidth',2)
box off
xlim([0 0.2])

%%
figure
plot(rho(round(Conditions.chosen_color(3)/2/pi*100),:),1:6,'-k','LineWidth',2,'Color',colors(round(Conditions.chosen_color(trial-1)/2/pi*40+1),:))
box off


%%

arrayID=40;
data_name=sprintf('ID_%d.mat',arrayID);
load(fullfile(data_path,data_name),'Conditions')


figure
hold on
for trial=1:length(Conditions.rule)
    
    x=[trial-1 trial trial trial-1];
    y=[0 0 1 1];
    this_color=colors(round(Conditions.rule(trial)/2/pi*40),:);
    patch(x,y,this_color,'EdgeColor',this_color)

end



