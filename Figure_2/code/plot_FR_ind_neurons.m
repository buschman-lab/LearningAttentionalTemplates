clear all

fsroot='/Volumes/buschman';
arrayID_list=[3 8 12 14 18 20 22 26 28 30 32 34 36 38 40 42 44];
event='target';
this_time=1;

arrayID=22;

task='Learning_Attentional_Templates';
subtask='exploreexploit/Restricted_FEF/Reset_RW_model';

window_start_list=-600;
window_size=900;

subsubtask=sprintf('FR_%s_%d_window_from_%d_to_%d',event,window_size,window_start_list,window_start_list+window_size);

dirstem = fullfile('Projects',task,'Analysis','Electrophy_analysis',subtask,subsubtask);
data_path = fullfile(fsroot,dirstem);

savename=sprintf('belief_peak_decoder_update_NN_relu2_%d',arrayID);

N_channels=3;

par_param=0;
switch fsroot
    case '/Volumes/buschman'
        par_param=1;
end

circular_distance = @(x, y)abs(mod(x-y+pi,2*pi)-pi);

%% session

%Load data
data_name=sprintf('ID_%d.mat',arrayID);

load(fullfile(data_path,data_name),'Conditions','LIP','FEF','PFC')
if arrayID==34 || arrayID==30
    this_ROI=[FEF; PFC];
    this_region=[2*ones(size(FEF,1),1); 3*ones(size(PFC,1),1)];
elseif arrayID==14
    this_ROI=[LIP; PFC];
    this_region=[ones(size(LIP,1),1); 3*ones(size(PFC,1),1)];
else
    this_ROI=[LIP; FEF; PFC];
    this_region=[ones(size(LIP,1),1); 2*ones(size(FEF,1),1); 3*ones(size(PFC,1),1)];
end
clear LIP FEF PFC

for n=1:size(this_ROI,1)
    if sum(~isnan(this_ROI(n,:)))>0
        this_ROI(n,isnan(this_ROI(n,:)))=0;
    end
end

%Extract peak template
N_bins=100;
color_binned = 0:(2*pi/N_bins):(2*pi-(2*pi/N_bins));

for i=1:size(Conditions.stim,2)
    [~, peak_belief_index(i)]=max(Conditions.belief(:,i));
end
peak_belief=color_binned(peak_belief_index);
for i=1:size(Conditions.stim,2)
    if sum(Conditions.belief(:,i)==0)==size(Conditions.belief,1)
        peak_belief(i)=NaN;
    end
end

%Bin peaked template
edges=0:2*pi/3:2*pi;
binned_peak=discretize (peak_belief, edges);

%% fit model

for t=1:size(this_ROI,1)
    if sum(~isnan(this_ROI(t,:)))>500
        [R2(t),~,~]= regf_stim(peak_belief,this_ROI(t,:));
    else
        R2(t)=NaN;
    end
end

for i=1:3
    [~, ind(i)]=max(R2(this_region==i));
    ind(i)=ind(i)+sum(this_region<i);
end

%% load the time serie

event='target';

for i=1:20
    window_start_list(i)=-600+(i-1)*50;
end

target_on=13;
response_on=17;

window_size=300;

subsubtask_classifier=sprintf('PCA_%s_%d_to_%s_%d',event,window_start_list(1),event,window_start_list(end));
for t=1:20
    window_start=window_start_list(t);
    
    subsubtask_data=sprintf('FR_%s_%d_window_from_%d_to_%d',event,window_size,window_start,window_start+window_size);
    data_name=sprintf('ID_%d.mat',arrayID);
    dirstem = fullfile('Projects',task,'Analysis','Electrophy_analysis',subtask,subsubtask_data);
    data_path_time_serie = fullfile(fsroot,dirstem);
    
    load(fullfile(data_path_time_serie,data_name),'Conditions','LIP','FEF','PFC')
    if arrayID==34 || arrayID==30
        this_time_serie(t,:,:)=[FEF; PFC];
    elseif arrayID==14
        this_time_serie(t,:,:)=[LIP; PFC];
    else
        this_time_serie(t,:,:)=[LIP; FEF; PFC];
    end
end


%% Plot the FR
color_for_ROI=[0 0.447 0.741; 0.85 0.325 0.098; 0.929 0.694 0.125]; %LIP FEF PFC

this_time_serie=nanzscore(this_time_serie,3);

figure
for i=1:3
subplot(1,3,i)
for j=1:3
    hold on
shadedErrorBar(window_start_list,nanmean(this_time_serie(:,ind(j),binned_peak==i),3),nanstd(this_time_serie(:,ind(j),binned_peak==i),0,3)./sqrt(sum(binned_peak==i)-1),{'color',color_for_ROI(j,:),'LineWidth',2},2)
end
ylim([-0.5 1]);
end

%%

function [R2,param,exitflag] = regf_stim(x,y) %vm

x0=zeros(4,1);
ub=[10 10 2*pi 2.5];
lb=[-10 -10 -2*pi -2.5];

options = optimoptions('fmincon','Display','off');

fun_fit=@(P)f_full_model_stim_vm(P,y,x);

[param,~,exitflag] = fmincon(fun_fit,x0,[],[],[],[],lb,ub,[],options);

mse=f_full_model_stim_vm(param,y,x);
mse_null=1/length(y).*(sum((y-mean(y)).^2));
R2=1-mse/mse_null;

end


function z = nanzscore(y,dim)

z=(y-nanmean(y,dim))./nanstd(y,0,dim);

end

