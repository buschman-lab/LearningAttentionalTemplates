
clear all

event='target';

window_size=200;
for i=1:25
    initial_window=-400;
    event_list{i}='target';
    window_start_list(i)=initial_window+(i-1)*50;
end

fsroot='/Volumes/buschman';
task='Learning_Attentional_Templates';
subtask='exploreexploit/Restricted_FEF/Reset_RW_model';

subsubtask_save=sprintf('Classifier_%s_%d_to_%s_%d',event,window_start_list(1),event,window_start_list(end));
dirstem2 = fullfile('Projects',task,'Analysis','Electrophy_analysis',subtask,subsubtask_save);
load_path=fullfile(fsroot,dirstem2);

count_ROI=ones(length(window_start_list),1);

arrayID_list=[3,8,12,14,18,20,22,26,28,30,32,34,36,38,40,42,44];

color_for_ROI=[0 0.447 0.741; 0.85 0.325 0.098; 0.929 0.694 0.125]; %LIP FEF PFC


color_each_loc(1,:,1)=color_for_ROI(1,:)+0.25;
color_each_loc(1,:,2)=color_for_ROI(1,:)+0.19;
color_each_loc(1,:,3)=color_for_ROI(1,:)+0.06;
color_each_loc(1,:,4)=color_for_ROI(1,:);

color_each_loc(2,:,1)=color_for_ROI(2,:)+0.15;
color_each_loc(2,:,2)=color_for_ROI(2,:)+0.09;
color_each_loc(2,:,3)=color_for_ROI(2,:)-0.03;
color_each_loc(2,:,4)=color_for_ROI(2,:)-0.09;

color_each_loc(3,:,1)=color_for_ROI(3,:)+0.07;
color_each_loc(3,:,2)=color_for_ROI(3,:)+0.01;
color_each_loc(3,:,3)=color_for_ROI(3,:)-0.06;
color_each_loc(3,:,4)=color_for_ROI(3,:)-0.12;


%% load

%LIP
load_name=sprintf('glm_explained_variance_split_boot_%s','LIP');
load(fullfile(load_path,load_name),'p_chosen_contra_boot',...
    'r_chosen_contra_reward_global_boot','p_contra_chosen_global_boot','p_contra_reward_global_boot',...
    'r_chosen_contra_chosen_global_boot','r_chosen_contra_boot','r_chosen_global_boot');

LIP.r_chosen_contra_boot=r_chosen_contra_boot;
LIP.p_chosen_contra_boot=p_chosen_contra_boot;
LIP.r_chosen_contra_chosen_global_boot=r_chosen_contra_chosen_global_boot;
LIP.r_chosen_contra_reward_global_boot=r_chosen_contra_reward_global_boot;
LIP.p_contra_chosen_global_boot=p_contra_chosen_global_boot;
LIP.p_contra_reward_global_boot=p_contra_reward_global_boot;
LIP.r_chosen_global_boot=r_chosen_global_boot;
for t=1:length(window_start_list)
    y=reshape(mean(r_chosen_global_boot(t,:,:),2),size(r_chosen_global_boot,3),1);
    LIP.p_chosen_global_boot(t)=z_test_function_bootstrap(y,0);
end

for t=1:length(window_start_list)
    for n=1:size(r_chosen_contra_boot,3)
        if mean(r_chosen_contra_boot(t,:,n),2)>0 
            if mean(r_chosen_global_boot(t,:,n),2)>0
            sim_r_contra_chosen(t,n)=mean(r_chosen_contra_chosen_global_boot(t,:,n),2)./(sqrt(mean(r_chosen_contra_boot(t,:,n),2).*mean(r_chosen_global_boot(t,:,n),2)));
            else
                sim_r_contra_chosen(t,n)=NaN;
            end
            diss_r_contra_chosen(t,n)=mean(r_chosen_contra_boot(t,:,n),2)-mean(r_chosen_contra_chosen_global_boot(t,:,n),2);
        else
            sim_r_contra_chosen(t,n)=NaN;
            diss_r_contra_chosen(t,n)=NaN;
        end
    end
    LIP.p_chosen_contra_diss(t) = z_test_function_bootstrap(diss_r_contra_chosen(t,:),0);
    LIP.p_chosen_contra_sim(t) = z_test_function_bootstrap(sim_r_contra_chosen(t,:),0);
end

clear *boot* *chosen* *pref* *contra*

%FEF
load_name=sprintf('glm_explained_variance_split_boot_%s','FEF');
load(fullfile(load_path,load_name),'p_chosen_contra_boot',...
    'r_chosen_contra_reward_global_boot','p_contra_chosen_global_boot','p_contra_reward_global_boot',...
    'r_chosen_contra_chosen_global_boot','r_chosen_contra_boot','r_chosen_global_boot');

FEF.r_chosen_contra_boot=r_chosen_contra_boot;
FEF.p_chosen_contra_boot=p_chosen_contra_boot;
FEF.r_chosen_contra_chosen_global_boot=r_chosen_contra_chosen_global_boot;
FEF.r_chosen_contra_reward_global_boot=r_chosen_contra_reward_global_boot;
FEF.p_contra_chosen_global_boot=p_contra_chosen_global_boot;
FEF.p_contra_reward_global_boot=p_contra_reward_global_boot;
FEF.r_chosen_global_boot=r_chosen_global_boot;
for t=1:length(window_start_list)
    y=reshape(mean(r_chosen_global_boot(t,:,:),2),size(r_chosen_global_boot,3),1);
    FEF.p_chosen_global_boot(t)=z_test_function_bootstrap(y,0);
end

for t=1:length(window_start_list)
    for n=1:size(r_chosen_contra_boot,3)
        if mean(r_chosen_contra_boot(t,:,n),2)>0 
            if mean(r_chosen_global_boot(t,:,n),2)>0
            sim_r_contra_chosen(t,n)=mean(r_chosen_contra_chosen_global_boot(t,:,n),2)./(sqrt(mean(r_chosen_contra_boot(t,:,n),2).*mean(r_chosen_global_boot(t,:,n),2)));
            else
                sim_r_contra_chosen(t,n)=NaN;
            end
            diss_r_contra_chosen(t,n)=mean(r_chosen_contra_boot(t,:,n),2)-mean(r_chosen_contra_chosen_global_boot(t,:,n),2);
        else
            sim_r_contra_chosen(t,n)=NaN;
            diss_r_contra_chosen(t,n)=NaN;
        end
    end
    FEF.p_chosen_contra_diss(t) = z_test_function_bootstrap(diss_r_contra_chosen(t,:),0);
    FEF.p_chosen_contra_sim(t) = z_test_function_bootstrap(sim_r_contra_chosen(t,:),0);
end

clear *boot* *chosen* *pref* *contra*

%PFC
load_name=sprintf('glm_explained_variance_split_boot_%s','PFC');
load(fullfile(load_path,load_name),'p_chosen_contra_boot',...
    'r_chosen_contra_reward_global_boot','p_contra_chosen_global_boot','p_contra_reward_global_boot',...
    'r_chosen_contra_chosen_global_boot','r_chosen_contra_boot','r_chosen_global_boot');

PFC.r_chosen_contra_boot=r_chosen_contra_boot;
PFC.p_chosen_contra_boot=p_chosen_contra_boot;
PFC.r_chosen_contra_chosen_global_boot=r_chosen_contra_chosen_global_boot;
PFC.r_chosen_contra_reward_global_boot=r_chosen_contra_reward_global_boot;
PFC.p_contra_chosen_global_boot=p_contra_chosen_global_boot;
PFC.p_contra_reward_global_boot=p_contra_reward_global_boot;
PFC.r_chosen_global_boot=r_chosen_global_boot;
for t=1:length(window_start_list)
    y=reshape(mean(r_chosen_global_boot(t,:,:),2),size(r_chosen_global_boot,3),1);
    PFC.p_chosen_global_boot(t)=z_test_function_bootstrap(y,0);
end

for t=1:length(window_start_list)
    for n=1:size(r_chosen_contra_boot,3)
        if mean(r_chosen_contra_boot(t,:,n),2)>0 
            if mean(r_chosen_global_boot(t,:,n),2)>0
            sim_r_contra_chosen(t,n)=mean(r_chosen_contra_chosen_global_boot(t,:,n),2)./(sqrt(mean(r_chosen_contra_boot(t,:,n),2).*mean(r_chosen_global_boot(t,:,n),2)));
            else
                sim_r_contra_chosen(t,n)=NaN;
            end
            diss_r_contra_chosen(t,n)=mean(r_chosen_contra_boot(t,:,n),2)-mean(r_chosen_contra_chosen_global_boot(t,:,n),2);
        else
            sim_r_contra_chosen(t,n)=NaN;
            diss_r_contra_chosen(t,n)=NaN;
        end
    end
    PFC.p_chosen_contra_diss(t) = z_test_function_bootstrap(diss_r_contra_chosen(t,:),0);
    PFC.p_chosen_contra_sim(t) = z_test_function_bootstrap(sim_r_contra_chosen(t,:),0);
end

clear *boot* *chosen* *pref* *contra*


%% Contra loc

figure
subplot(3,1,1)
n_color=1;
hold on
shadedErrorBar(window_start_list(4:end),mean(mean(LIP.r_chosen_contra_boot(4:end,:,:),2),3),[prctile(mean(LIP.r_chosen_contra_boot(4:end,:,:),2),95,3) prctile(mean(LIP.r_chosen_contra_boot(4:end,:,:),2),5,3)]',{'-','color',[0.5 0.5 0.5],'LineWidth',2},2)
shadedErrorBar(window_start_list(4:end),mean(mean(LIP.r_chosen_global_boot(4:end,:,:),2),3),[prctile(mean(LIP.r_chosen_global_boot(4:end,:,:),2),95,3) prctile(mean(LIP.r_chosen_global_boot(4:end,:,:),2),5,3)]',{'--','color',[0.5 0.5 0.5],'LineWidth',2},2)
shadedErrorBar(window_start_list(4:end),mean(mean(LIP.r_chosen_contra_chosen_global_boot(4:end,:,:),2),3),[prctile(mean(LIP.r_chosen_contra_chosen_global_boot(4:end,:,:),2),95,3) prctile(mean(LIP.r_chosen_contra_chosen_global_boot(4:end,:,:),2),5,3)]',{'-','color',color_for_ROI(n_color,:),'LineWidth',2},2)

plot_significance_level(window_start_list,LIP.p_chosen_contra_boot,0.8,[0.5 0.5 0.5],[0.01,0.05/(length(window_start_list)-3), 0.01/(length(window_start_list)-3)])
plot_significance_level(window_start_list,LIP.p_chosen_global_boot,0.85,[0.5 0.5 0.5],[0.01,0.05/(length(window_start_list)-3), 0.01/(length(window_start_list)-3)])

plot_significance_level(window_start_list,LIP.p_chosen_contra_sim,0.9,color_for_ROI(n_color,:),[0.01,0.05/(length(window_start_list)-3), 0.01/(length(window_start_list)-3)])
plot_significance_level(window_start_list,LIP.p_chosen_contra_diss,0.95,color_for_ROI(n_color,:),[0.01,0.05/(length(window_start_list)-3), 0.01/(length(window_start_list)-3)])

yl.LabelVerticalAlignment = 'middle';
yl.LabelHorizontalAlignment = 'right';
yl.FontSize=12;

text(-575,0.8,'corr(l,l)')
text(-575,0.85,'corr(g,g)')

text(-575,0.9,'sim(l,g)')
text(-575,0.95,'diss(l,g)')

box off
xlabel('Time to target')
xlim([-600 window_start_list(end)+150])
xticks([-600:200:window_start_list(end)+150])
yline(0,'--')
xline(300,'-','Color',[0.25 0.25 0.75],'LineWidth',1)
hold on
xl=xline(400,'-',{'Reward'},'Color',[0.25 0.25 0.75],'LineWidth',1);
xl.LabelVerticalAlignment = 'bottom';
xl.LabelHorizontalAlignment = 'right';
xl.FontSize=12;
ax=gca;
ax.XAxis.FontSize=12;
ax.YAxis.FontSize=12;
ylabel('Correlation chosen value vectors')

ylim([-0.2 1])

subplot(3,1,2)
n_color=2;
hold on
shadedErrorBar(window_start_list(4:end),mean(mean(FEF.r_chosen_contra_boot(4:end,:,:),2),3),[prctile(mean(FEF.r_chosen_contra_boot(4:end,:,:),2),95,3) prctile(mean(FEF.r_chosen_contra_boot(4:end,:,:),2),5,3)]',{'-','color',[0.5 0.5 0.5],'LineWidth',2},2)
shadedErrorBar(window_start_list(4:end),mean(mean(FEF.r_chosen_global_boot(4:end,:,:),2),3),[prctile(mean(FEF.r_chosen_global_boot(4:end,:,:),2),95,3) prctile(mean(FEF.r_chosen_global_boot(4:end,:,:),2),5,3)]',{'--','color',[0.5 0.5 0.5],'LineWidth',2},2)
shadedErrorBar(window_start_list(4:end),mean(mean(FEF.r_chosen_contra_chosen_global_boot(4:end,:,:),2),3),[prctile(mean(FEF.r_chosen_contra_chosen_global_boot(4:end,:,:),2),95,3) prctile(mean(FEF.r_chosen_contra_chosen_global_boot(4:end,:,:),2),5,3)]',{'-','color',color_for_ROI(n_color,:),'LineWidth',2},2)

plot_significance_level(window_start_list,FEF.p_chosen_contra_boot,0.8,[0.5 0.5 0.5],[0.01,0.05/(length(window_start_list)-3), 0.01/(length(window_start_list)-3)])
plot_significance_level(window_start_list,FEF.p_chosen_global_boot,0.85,[0.5 0.5 0.5],[0.01,0.05/(length(window_start_list)-3), 0.01/(length(window_start_list)-3)])

plot_significance_level(window_start_list,FEF.p_chosen_contra_sim,0.9,color_for_ROI(n_color,:),[0.01,0.05/(length(window_start_list)-3), 0.01/(length(window_start_list)-3)])
plot_significance_level(window_start_list,FEF.p_chosen_contra_diss,0.95,color_for_ROI(n_color,:),[0.01,0.05/(length(window_start_list)-3), 0.01/(length(window_start_list)-3)])


yl.LabelVerticalAlignment = 'middle';
yl.LabelHorizontalAlignment = 'right';
yl.FontSize=12;
text(-575,0.8,'corr(l,l)')
text(-575,0.85,'corr(g,g)')

text(-575,0.9,'sim(l,g)')
text(-575,0.95,'diss(l,g)')
box off
xlabel('Time to target')
xlim([-600 window_start_list(end)+150])
xticks([-600:200:window_start_list(end)+150])
yline(0,'--')
xline(300,'-','Color',[0.25 0.25 0.75],'LineWidth',1)
hold on
xl=xline(400,'-',{'Reward'},'Color',[0.25 0.25 0.75],'LineWidth',1);
xl.LabelVerticalAlignment = 'bottom';
xl.LabelHorizontalAlignment = 'right';
xl.FontSize=12;
ax=gca;
ax.XAxis.FontSize=12;
ax.YAxis.FontSize=12;
ylabel('Correlation chosen value vectors')
ylim([-0.2 1])


subplot(3,1,3)
n_color=3;
hold on
shadedErrorBar(window_start_list(4:end),mean(mean(PFC.r_chosen_contra_boot(4:end,:,:),2),3),[prctile(mean(PFC.r_chosen_contra_boot(4:end,:,:),2),95,3) prctile(mean(PFC.r_chosen_contra_boot(4:end,:,:),2),5,3)]',{'-','color',[0.5 0.5 0.5],'LineWidth',2},2)
shadedErrorBar(window_start_list(4:end),mean(mean(PFC.r_chosen_global_boot(4:end,:,:),2),3),[prctile(mean(PFC.r_chosen_global_boot(4:end,:,:),2),95,3) prctile(mean(PFC.r_chosen_global_boot(4:end,:,:),2),5,3)]',{'--','color',[0.5 0.5 0.5],'LineWidth',2},2)
shadedErrorBar(window_start_list(4:end),mean(mean(PFC.r_chosen_contra_chosen_global_boot(4:end,:,:),2),3),[prctile(mean(PFC.r_chosen_contra_chosen_global_boot(4:end,:,:),2),95,3) prctile(mean(PFC.r_chosen_contra_chosen_global_boot(4:end,:,:),2),5,3)]',{'-','color',color_for_ROI(n_color,:),'LineWidth',2},2)

plot_significance_level(window_start_list,PFC.p_chosen_contra_boot,0.8,[0.5 0.5 0.5],[0.01,0.05/(length(window_start_list)-3), 0.01/(length(window_start_list)-3)])
plot_significance_level(window_start_list,PFC.p_chosen_global_boot,0.85,[0.5 0.5 0.5],[0.01,0.05/(length(window_start_list)-3), 0.01/(length(window_start_list)-3)])

plot_significance_level(window_start_list,PFC.p_chosen_contra_sim,0.9,color_for_ROI(n_color,:),[0.01,0.05/(length(window_start_list)-3), 0.01/(length(window_start_list)-3)])
plot_significance_level(window_start_list,PFC.p_chosen_contra_diss,0.95,color_for_ROI(n_color,:),[0.01,0.05/(length(window_start_list)-3), 0.01/(length(window_start_list)-3)])

yl.LabelVerticalAlignment = 'middle';
yl.LabelHorizontalAlignment = 'right';
yl.FontSize=12;
text(-575,0.8,'corr(l,l)')
text(-575,0.85,'corr(g,g)')

text(-575,0.9,'sim(l,g)')
text(-575,0.95,'diss(l,g)')
box off
xlabel('Time to target')
xlim([-600 window_start_list(end)+150])
xticks([-600:200:window_start_list(end)+150])
yline(0,'--')
xline(300,'-','Color',[0.25 0.25 0.75],'LineWidth',1)
hold on
xl=xline(400,'-',{'Reward'},'Color',[0.25 0.25 0.75],'LineWidth',1);
xl.LabelVerticalAlignment = 'bottom';
xl.LabelHorizontalAlignment = 'right';
xl.FontSize=12;
ax=gca;
ax.XAxis.FontSize=12;
ax.YAxis.FontSize=12;
ylabel('Correlation chosen value vectors')
ylim([-0.2 1])

%%
figure
subplot(3,1,1)
n_color=1;
hold on
shadedErrorBar(window_start_list(4:end),mean(mean(LIP.r_chosen_contra_chosen_global_boot(4:end,:,:),2),3),[prctile(mean(LIP.r_chosen_contra_chosen_global_boot(4:end,:,:),2),95,3) prctile(mean(LIP.r_chosen_contra_chosen_global_boot(4:end,:,:),2),5,3)]',{'-','color',color_for_ROI(1,:),'LineWidth',2},2)
shadedErrorBar(window_start_list(4:end),mean(mean(FEF.r_chosen_contra_chosen_global_boot(4:end,:,:),2),3),[prctile(mean(FEF.r_chosen_contra_chosen_global_boot(4:end,:,:),2),95,3) prctile(mean(FEF.r_chosen_contra_chosen_global_boot(4:end,:,:),2),5,3)]',{'-','color',color_for_ROI(2,:),'LineWidth',2},2)
shadedErrorBar(window_start_list(4:end),mean(mean(PFC.r_chosen_contra_chosen_global_boot(4:end,:,:),2),3),[prctile(mean(PFC.r_chosen_contra_chosen_global_boot(4:end,:,:),2),95,3) prctile(mean(PFC.r_chosen_contra_chosen_global_boot(4:end,:,:),2),5,3)]',{'-','color',color_for_ROI(3,:),'LineWidth',2},2)

plot_significance_level(window_start_list,LIP.p_chosen_contra_sim,0.9,color_for_ROI(1,:),[0.01,0.05/(length(window_start_list)-3), 0.01/(length(window_start_list)-3)])
plot_significance_level(window_start_list,LIP.p_chosen_contra_diss,0.7,color_for_ROI(1,:),[0.01,0.05/(length(window_start_list)-3), 0.01/(length(window_start_list)-3)])

plot_significance_level(window_start_list,FEF.p_chosen_contra_sim,0.85,color_for_ROI(2,:),[0.01,0.05/(length(window_start_list)-3), 0.01/(length(window_start_list)-3)])
plot_significance_level(window_start_list,FEF.p_chosen_contra_diss,0.65,color_for_ROI(2,:),[0.01,0.05/(length(window_start_list)-3), 0.01/(length(window_start_list)-3)])

plot_significance_level(window_start_list,PFC.p_chosen_contra_sim,0.8,color_for_ROI(3,:),[0.01,0.05/(length(window_start_list)-3), 0.01/(length(window_start_list)-3)])
plot_significance_level(window_start_list,PFC.p_chosen_contra_diss,0.6,color_for_ROI(3,:),[0.01,0.05/(length(window_start_list)-3), 0.01/(length(window_start_list)-3)])


yl.LabelVerticalAlignment = 'middle';
yl.LabelHorizontalAlignment = 'right';
yl.FontSize=12;

text(-575,0.9,'sim LIP')
text(-575,0.85,'sim FEF')
text(-575,0.8,'sim lPFC')

text(-575,0.7,'diss LIP')
text(-575,0.65,'diss FEF')
text(-575,0.6,'diss lPFC')

box off
xlabel('Time to target')
xlim([-600 window_start_list(end)+150])
xticks([-600:200:window_start_list(end)+150])
yline(0,'--')
xline(300,'-','Color',[0.25 0.25 0.75],'LineWidth',1)
hold on
xl=xline(400,'-',{'Reward'},'Color',[0.25 0.25 0.75],'LineWidth',1);
xl.LabelVerticalAlignment = 'bottom';
xl.LabelHorizontalAlignment = 'right';
xl.FontSize=12;
ax=gca;
ax.XAxis.FontSize=12;
ax.YAxis.FontSize=12;
ylabel('Correlation chosen value vectors')

ylim([-0.2 1])


%% Functions

function plot_significance_level(x,p,a,c,thr)
hold on
for b=1:length(thr)
    this_thr=thr(b);
    for i=4:length(p)-1
        if p(i)<this_thr && p(i+1)<this_thr
            plot(x(i):x(i+1),a*ones(1,x(i+1)-x(i)+1),'-','Color',c,'LineWidth',b)
        end
        if p(i-1)>=this_thr && p(i)<this_thr && p(i+1)>=this_thr
            plot(x(i),a,'.','Color',c,'LineWidth',b)
        end
    end
end
end

function p = z_test_function_bootstrap(dist,null)

if sum(~isnan(dist))>2
    
    dist=dist(~isnan(dist));
    m = mean(dist);
    s = std(dist);
    z = (m-null)/s;
    p=1-normcdf(z);
    
else
    p=NaN;
end

end

