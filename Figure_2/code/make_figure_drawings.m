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


arrayID=3;
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
trial=10;


figure
polarplot(color_binned,belief_probability(trial,:),'LineWidth',2)
hold on
polarplot(color_binned,mean_value(trial).*ones(N_bins,1),'--','LineWidth',2)

box off
ax=gca;
ax.ThetaTickLabel={};
ax.RTickLabel={};
