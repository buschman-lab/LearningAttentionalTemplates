arrayID_list=[3 8 12 14 18 20 22 26 28 30 32 34 36 38 40 42 44];

fsroot='/Volumes/buschman';

this_time=1;

event='target';
window_start_list=-600;
window_size=900;

task='Learning_Attentional_Templates';
% subtask='exploreexploit/Reset_RW_model';
subtask='exploreexploit/Restricted_FEF/Reset_RW_model';


subsubtask=sprintf('FR_%s_%d_window_from_%d_to_%d',event,window_size,window_start_list(this_time),window_start_list(this_time)+window_size);

dirstem = fullfile('Projects',task,'Analysis','Electrophy_analysis',subtask,subsubtask);
data_path = fullfile(fsroot,dirstem);

savename=sprintf('MDS_template');

N_channels=6;
N_bins=N_channels;

count_LIP=1;
count_FEF=1;
count_PFC=1;

N_sessions=length(arrayID_list);

%% loop over sessions

for n_sess=1:N_sessions
    
    arrayID=arrayID_list(n_sess);
    
    
    %LIP
    if arrayID~=34 && arrayID~=30
        
        data_name=sprintf('ID_%d.mat',arrayID);
        dirstem = fullfile('Projects',task,'Analysis','Electrophy_analysis',subtask,subsubtask);
        data_path_clasifier = fullfile(fsroot,dirstem);
        
        load(fullfile(data_path_clasifier,data_name),'Conditions','LIP')
        
        %get the progression in trial
        for i=1:length(Conditions.block_nb)
            Conditions.block_nb(i)=Conditions.trial_in_block(i)/sum(Conditions.block_nb==Conditions.block_nb(i));
        end
        %get the peak belief
        N_bins=100;
        color_binned =  0:(2*pi/N_bins):(2*pi-(2*pi/N_bins));
        
        for i=1:size(Conditions.stim,2)
            [~, peak_belief_index(i)]=max(Conditions.belief(:,i));
        end
        Conditions.Peak_belief=color_binned(peak_belief_index);
        
        for i=1:size(Conditions.stim,2)
            if sum(Conditions.belief(:,i)==0)==size(Conditions.belief,1)
                Conditions.Peak_belief(i)=NaN;
            end
        end
        
        N_neurons=size(LIP,1);
        
        for nn=1:N_neurons
            for i=1:length(Conditions.block_nb)
                Raw_data_LIP(nn,i)=LIP(nn,i);
            end
        end
        clear LIP
        
        %remove partially active neurons
        data_LIP=Raw_data_LIP(~isnan(sum(Raw_data_LIP(:,:,1),2)),:,:); %remove partially active neurons
        clear Raw_data_LIP
        for nn=1:size(data_LIP,1)
            
            %first get the conditions
            for i=2:length(Conditions.block_nb)
                for t=1:length(window_start_list)
                    Pseudo_pop_LIP(count_LIP,i-1).Classifier_FR=data_LIP(nn,i);
                end
                Pseudo_pop_LIP(count_LIP,i-1).Sess=n_sess;
                Pseudo_pop_LIP(count_LIP,i-1).Block_nb=Conditions.block_nb(i); %we need to take the updated belief here
                Pseudo_pop_LIP(count_LIP,i-1).Peak_belief=mod(Conditions.Peak_belief(i),2*pi);
                for n_c=1:N_channels
                    if mod(Pseudo_pop_LIP(count_LIP,i-1).Peak_belief,2*pi)>=(n_c-1)*2*pi/N_channels && mod(Pseudo_pop_LIP(count_LIP,i-1).Peak_belief,2*pi)<n_c*2*pi/N_channels
                        Pseudo_pop_LIP(count_LIP,i-1).Belief_color=n_c;
                    end
                end
            end
            count_LIP=count_LIP+1;
        end
        
        clear Conditions
        
    end
    
    if arrayID~=14
        %FEF
        data_name=sprintf('ID_%d.mat',arrayID);
        dirstem = fullfile('Projects',task,'Analysis','Electrophy_analysis',subtask,subsubtask);
        data_path_clasifier = fullfile(fsroot,dirstem);
        
        load(fullfile(data_path_clasifier,data_name),'Conditions','FEF')
        
        %get the progression in trial
        for i=1:length(Conditions.block_nb)
            Conditions.block_nb(i)=Conditions.trial_in_block(i)/sum(Conditions.block_nb==Conditions.block_nb(i));
        end
        %get the peak belief
        N_bins=100;
        color_binned =  0:(2*pi/N_bins):(2*pi-(2*pi/N_bins));
        
        for i=1:size(Conditions.stim,2)
            [~, peak_belief_index(i)]=max(Conditions.belief(:,i));
        end
        Conditions.Peak_belief=color_binned(peak_belief_index);
        
        for i=1:size(Conditions.stim,2)
            if sum(Conditions.belief(:,i)==0)==size(Conditions.belief,1)
                Conditions.Peak_belief(i)=NaN;
            end
        end
        
        N_neurons=size(FEF,1);
        
        for nn=1:N_neurons
            for i=1:length(Conditions.block_nb)
                Raw_data_FEF(nn,i)=FEF(nn,i);
            end
        end
        clear FEF
        
        %remove partially active neurons
        data_FEF=Raw_data_FEF(~isnan(sum(Raw_data_FEF(:,:,1),2)),:,:); %remove partially active neurons
        clear Raw_data_FEF
        for nn=1:size(data_FEF,1)
            
            %first get the conditions
            for i=2:length(Conditions.block_nb)
                for t=1:length(window_start_list)
                    Pseudo_pop_FEF(count_FEF,i-1).Classifier_FR=data_FEF(nn,i);
                end
                Pseudo_pop_FEF(count_FEF,i-1).Sess=n_sess;
                Pseudo_pop_FEF(count_FEF,i-1).Block_nb=Conditions.block_nb(i); %we need to take the updated belief here
                Pseudo_pop_FEF(count_FEF,i-1).Peak_belief=mod(Conditions.Peak_belief(i),2*pi);
                for n_c=1:N_channels
                    if mod(Pseudo_pop_FEF(count_FEF,i-1).Peak_belief,2*pi)>=(n_c-1)*2*pi/N_channels && mod(Pseudo_pop_FEF(count_FEF,i-1).Peak_belief,2*pi)<n_c*2*pi/N_channels
                        Pseudo_pop_FEF(count_FEF,i-1).Belief_color=n_c;
                    end
                end
            end
            count_FEF=count_FEF+1;
        end
        
        clear Conditions
        
    end
    
    %PFC
    data_name=sprintf('ID_%d.mat',arrayID);
    dirstem = fullfile('Projects',task,'Analysis','Electrophy_analysis',subtask,subsubtask);
    data_path_clasifier = fullfile(fsroot,dirstem);
    
    load(fullfile(data_path_clasifier,data_name),'Conditions','PFC')
    
    %get the progression in trial
    for i=1:length(Conditions.block_nb)
        Conditions.block_nb(i)=Conditions.trial_in_block(i)/sum(Conditions.block_nb==Conditions.block_nb(i));
    end
    %get the peak belief
    N_bins=100;
    color_binned =  0:(2*pi/N_bins):(2*pi-(2*pi/N_bins));
    
    for i=1:size(Conditions.stim,2)
        [~, peak_belief_index(i)]=max(Conditions.belief(:,i));
    end
    Conditions.Peak_belief=color_binned(peak_belief_index);
    
    for i=1:size(Conditions.stim,2)
        if sum(Conditions.belief(:,i)==0)==size(Conditions.belief,1)
            Conditions.Peak_belief(i)=NaN;
        end
    end
    
    N_neurons=size(PFC,1);
    
    for nn=1:N_neurons
        for i=1:length(Conditions.block_nb)
            Raw_data_PFC(nn,i)=PFC(nn,i);
        end
    end
    clear PFC
    
    %remove partially active neurons
    data_PFC=Raw_data_PFC(~isnan(sum(Raw_data_PFC(:,:,1),2)),:,:); %remove partially active neurons
    clear Raw_data_PFC
    for nn=1:size(data_PFC,1)
        
        %first get the conditions
        for i=2:length(Conditions.block_nb)
            for t=1:length(window_start_list)
                Pseudo_pop_PFC(count_PFC,i-1).Classifier_FR=data_PFC(nn,i);
            end
            Pseudo_pop_PFC(count_PFC,i-1).Sess=n_sess;
            Pseudo_pop_PFC(count_PFC,i-1).Block_nb=Conditions.block_nb(i); %we need to take the updated belief here
            Pseudo_pop_PFC(count_PFC,i-1).Peak_belief=mod(Conditions.Peak_belief(i),2*pi);
            for n_c=1:N_channels
                if mod(Pseudo_pop_PFC(count_PFC,i-1).Peak_belief,2*pi)>=(n_c-1)*2*pi/N_channels && mod(Pseudo_pop_PFC(count_PFC,i-1).Peak_belief,2*pi)<n_c*2*pi/N_channels
                    Pseudo_pop_PFC(count_PFC,i-1).Belief_color=n_c;
                end
            end
        end
        count_PFC=count_PFC+1;
    end
    
    clear Conditions
    
end

%% select the neurons and trials

for nn=1:size(Pseudo_pop_LIP,1)
    for i=1:size(Pseudo_pop_LIP(nn,:),2)
        if ~isempty(Pseudo_pop_LIP(nn,i).Block_nb) && ~isnan(Pseudo_pop_LIP(nn,i).Peak_belief)
            selected_LIP(nn,i,Pseudo_pop_LIP(nn,i).Belief_color)=1;
        end
    end
    Sess_LIP(nn)= Pseudo_pop_LIP(nn,1).Sess;
end

for nn=1:size(Pseudo_pop_FEF,1)
    for i=1:size(Pseudo_pop_FEF(nn,:),2)
        if ~isempty(Pseudo_pop_FEF(nn,i).Block_nb) && ~isnan(Pseudo_pop_FEF(nn,i).Peak_belief)
            selected_FEF(nn,i,Pseudo_pop_FEF(nn,i).Belief_color)=1;
        end
    end
    Sess_FEF(nn)= Pseudo_pop_FEF(nn,1).Sess;
end

for nn=1:size(Pseudo_pop_PFC,1)
    for i=1:size(Pseudo_pop_PFC(nn,:),2)
        if ~isempty(Pseudo_pop_PFC(nn,i).Block_nb) && ~isnan(Pseudo_pop_PFC(nn,i).Peak_belief)
            
            selected_PFC(nn,i,Pseudo_pop_PFC(nn,i).Belief_color)=1;
        end
    end
    Sess_PFC(nn)= Pseudo_pop_PFC(nn,1).Sess;
end

%%
for nn=1:size(Pseudo_pop_LIP,1)
    for i=1:size(Pseudo_pop_LIP(nn,:),2)
        if ~isempty(Pseudo_pop_LIP(nn,i).Block_nb) && ~isnan(Pseudo_pop_LIP(nn,i).Peak_belief)
            selected_rand_lip(nn,i,randi(N_channels))=1;
        end
    end
end
for nn=1:size(Pseudo_pop_FEF,1)
    for i=1:size(Pseudo_pop_FEF(nn,:),2)
        if ~isempty(Pseudo_pop_FEF(nn,i).Block_nb) && ~isnan(Pseudo_pop_FEF(nn,i).Peak_belief)
            selected_rand_fef(nn,i,randi(N_channels))=1;
        end
    end
end
for nn=1:size(Pseudo_pop_PFC,1)
    for i=1:size(Pseudo_pop_PFC(nn,:),2)
        if ~isempty(Pseudo_pop_PFC(nn,i).Block_nb) && ~isnan(Pseudo_pop_PFC(nn,i).Peak_belief)
            selected_rand_pfc(nn,i,randi(N_channels))=1;
        end
    end
end

%% Now create the matrix N x 12 conditions (using a sliding window for the template)


for n_c=1:N_channels
    %LIP
    for nn=1:size(Pseudo_pop_LIP,1)
        if sum(selected_LIP(nn,:,n_c))>1
                list=find(selected_LIP(nn,:,n_c)==1);% | selected_LIP(nn,:,n_c+1)==1 | selected_LIP(nn,:,n_c+2)==1);
            for i=1:length(list)
                buffer(i)=Pseudo_pop_LIP(nn,list(i)).Classifier_FR;
            end
            LIP_MDS_template(nn,n_c)=mean(buffer);
            clear buffer list
        else
            LIP_MDS_template(nn,n_c)=NaN;
        end
    end
    %FEF
    for nn=1:size(Pseudo_pop_FEF,1)
        if sum(selected_FEF(nn,:,n_c))>1
                list=find(selected_FEF(nn,:,n_c)==1);% | selected_FEF(nn,:,n_c+1)==1 | selected_FEF(nn,:,n_c+2)==1);
            for i=1:length(list)
                buffer(i)=Pseudo_pop_FEF(nn,list(i)).Classifier_FR;
            end
            FEF_MDS_template(nn,n_c)=mean(buffer);
            clear buffer list
        else
            FEF_MDS_template(nn,n_c)=NaN;
        end
    end
    %PFC
    for nn=1:size(Pseudo_pop_PFC,1)
        if sum(selected_PFC(nn,:,n_c))>1
                list=find(selected_PFC(nn,:,n_c)==1);% | selected_PFC(nn,:,n_c+1)==1 | selected_PFC(nn,:,n_c+2)==1);
            for i=1:length(list)
                buffer(i)=Pseudo_pop_PFC(nn,list(i)).Classifier_FR;
            end
            PFC_MDS_template(nn,n_c)=mean(buffer);
            clear buffer list
        else
            PFC_MDS_template(nn,n_c)=NaN;
        end
    end
end


%%
Rand_MDS_template(1,1)=NaN;
    %RAnd
    for n_c=1:N_channels
    %LIP
    for nn=1:size(Pseudo_pop_LIP,1)
        if sum(selected_LIP(nn,:,n_c))>1
                list=find(selected_rand_lip(nn,:,n_c)==1);% | selected_LIP(nn,:,n_c+1)==1 | selected_LIP(nn,:,n_c+2)==1);
            for i=1:length(list)
                buffer(i)=Pseudo_pop_LIP(nn,list(i)).Classifier_FR;
            end
            Rand_MDS_template(end+1,n_c)=mean(buffer);
            clear buffer list
        else
            Rand_MDS_template(end+1,n_c)=NaN;
        end
    end
    %FEF
    for nn=1:size(Pseudo_pop_FEF,1)
        if sum(selected_FEF(nn,:,n_c))>1
                list=find(selected_rand_fef(nn,:,n_c)==1);% | selected_FEF(nn,:,n_c+1)==1 | selected_FEF(nn,:,n_c+2)==1);
            for i=1:length(list)
                buffer(i)=Pseudo_pop_FEF(nn,list(i)).Classifier_FR;
            end
            Rand_MDS_template(end+1,n_c)=mean(buffer);
            clear buffer list
        else
            Rand_MDS_template(end+1,n_c)=NaN;
        end
    end
    for nn=1:size(Pseudo_pop_PFC,1)
        if sum(selected_PFC(nn,:,n_c))>1
                list=find(selected_rand_pfc(nn,:,n_c)==1);% | selected_rand(nn,:,n_c+1)==1 | selected_rand(nn,:,n_c+2)==1);
            for i=1:length(list)
                buffer(i)=Pseudo_pop_PFC(nn,list(i)).Classifier_FR;
            end
            Rand_MDS_template(end+1,n_c)=mean(buffer);
            clear buffer list
        else
            Rand_MDS_template(end+1,n_c)=NaN;
        end
    end
    end

%% remove NaN

LIP_MDS_template=LIP_MDS_template(~isnan(sum(sum(LIP_MDS_template,2),3)),:,:);
FEF_MDS_template=FEF_MDS_template(~isnan(sum(sum(FEF_MDS_template,2),3)),:,:);
PFC_MDS_template=PFC_MDS_template(~isnan(sum(sum(PFC_MDS_template,2),3)),:,:);

%%
Rand_MDS_template=Rand_MDS_template(~isnan(sum(sum(Rand_MDS_template,2),3)),:,:);


%% do the pca

[coeff_LIP(:,:), score_LIP(:,:),~,~,explained_LIP(:),~] = pca(LIP_MDS_template(:,:)');
[coeff_FEF(:,:), score_FEF(:,:),~,~,explained_FEF(:),~] = pca(FEF_MDS_template(:,:)');
[coeff_PFC(:,:), score_PFC(:,:),~,~,explained_PFC(:),~] = pca(PFC_MDS_template(:,:)');

[coeff_all(:,:), score_all(:,:),~,~,explained_all(:),~] = pca([LIP_MDS_template; FEF_MDS_template; PFC_MDS_template]');


%%
[coeff_Rand, score_Rand,~,~,explained_Rand,~] = pca(Rand_MDS_template');

%% do the MDS

D_LIP=dist(LIP_MDS_template);
Y_LIP = cmdscale(D_LIP);

D_FEF=dist(FEF_MDS_template);
Y_FEF = cmdscale(D_FEF);

D_PFC=dist(PFC_MDS_template);
Y_PFC = cmdscale(D_PFC);

D_all=dist([LIP_MDS_template; FEF_MDS_template; PFC_MDS_template]);
Y_all = cmdscale(D_all);

opts = statset('MaxIter',5000);

D_LIP2=pdist(LIP_MDS_template');
Y_LIP2 = mdscale(D_LIP2,2,'Options',opts);

D_FEF2=pdist(FEF_MDS_template');
Y_FEF2 = mdscale(D_FEF2,2,'Options',opts);

D_PFC2=pdist(PFC_MDS_template');
Y_PFC2 = mdscale(D_PFC2,2,'Options',opts);

D_all2=pdist([LIP_MDS_template; FEF_MDS_template; PFC_MDS_template]');
Y_all2 = mdscale(D_all2,2,'Options',opts);


%%

D_Rand=dist(Rand_MDS_template);
Y_Rand = cmdscale(D_Rand);

D_Rand2=pdist(Rand_MDS_template');
Y_Rand2 = mdscale(D_Rand2,2,'Options',opts);


%% now plot

load('colors.mat')
for i=1:N_channels
    color_bin(i,:)=colors(1+ceil((i-1)*(40/N_channels)),:);
end


figure

PC_1=1;
PC_2=2;
PC_3=3;

subplot(3,3,1)
hold on


x=score_LIP(:,PC_1);
x(end+1)=score_LIP(1,PC_1);
y=score_LIP(:,PC_2);
y(end+1)=score_LIP(1,PC_2);
z=score_LIP(:,PC_3);
z(end+1)=score_LIP(1,PC_3);
plot(x,y,'-','Color','k','LineWidth',1)
for n_c=1:N_channels
    plot(score_LIP(n_c,PC_1),score_LIP(n_c,PC_2),'-o','MarkerFaceColor',color_bin(n_c,:),'MarkerEdgeColor',color_bin(n_c,:),'Color',color_bin(n_c,:),'MarkerSize',15)
    
end
grid on
lim(1)=min(min(x),min(y));
lim(2)=max(max(x),max(y));
xlim(lim)
ylim(lim)
clear x y z

title(sprintf('LIP PCA'))
xlabel(sprintf('PC %d',PC_1))
ylabel(sprintf('PC %d',PC_2))

subplot(3,3,2)
hold on


x=score_FEF(:,PC_1);
x(end+1)=score_FEF(1,PC_1);
y=score_FEF(:,PC_2);
y(end+1)=score_FEF(1,PC_2);
z=score_FEF(:,PC_3);
z(end+1)=score_FEF(1,PC_3);
plot(x,y,'-','Color','k','LineWidth',1)
for n_c=1:N_channels
    plot(score_FEF(n_c,PC_1),score_FEF(n_c,PC_2),'-o','MarkerFaceColor',color_bin(n_c,:),'MarkerEdgeColor',color_bin(n_c,:),'Color',color_bin(n_c,:),'MarkerSize',15)
    
end
grid on
lim(1)=min(min(x),min(y));
lim(2)=max(max(x),max(y));
xlim(lim)
ylim(lim)
clear x y z

title(sprintf('FEF PCA'))
xlabel(sprintf('PC %d',PC_1))
ylabel(sprintf('PC %d',PC_2))

subplot(3,3,3)
hold on


x=score_PFC(:,PC_1);
x(end+1)=score_PFC(1,PC_1);
y=score_PFC(:,PC_2);
y(end+1)=score_PFC(1,PC_2);
plot(x,y,'-','Color','k','LineWidth',1)
for n_c=1:N_channels
    plot(score_PFC(n_c,PC_1),score_PFC(n_c,PC_2),'-o','MarkerFaceColor',color_bin(n_c,:),'MarkerEdgeColor',color_bin(n_c,:),'Color',color_bin(n_c,:),'MarkerSize',15)
end
grid on
lim(1)=min(min(x),min(y));
lim(2)=max(max(x),max(y));
xlim(lim)
ylim(lim)
clear x y z

title(sprintf('PFC PCA'))
xlabel(sprintf('PC %d',PC_1))
ylabel(sprintf('PC %d',PC_2))
grid on

subplot(3,3,4)
hold on

x=Y_LIP(:,PC_1);
x(end+1)=Y_LIP(1,PC_1);
y=Y_LIP(:,PC_2);
y(end+1)=Y_LIP(1,PC_2);
plot(x,y,'-','Color','k','LineWidth',1)
for n_c=1:N_channels
    plot(Y_LIP(n_c,PC_1),Y_LIP(n_c,PC_2),'-o','MarkerFaceColor',color_bin(n_c,:),'MarkerEdgeColor',color_bin(n_c,:),'Color',color_bin(n_c,:),'MarkerSize',15)
end
grid on
lim(1)=min(min(x),min(y));
lim(2)=max(max(x),max(y));
xlim(lim)
ylim(lim)
clear x y z

title(sprintf('LIP CMDS'))
xlabel(sprintf('dim %d',PC_1))
ylabel(sprintf('dim %d',PC_2))
grid on

subplot(3,3,5)
hold on

x=Y_FEF(:,PC_1);
x(end+1)=Y_FEF(1,PC_1);
y=Y_FEF(:,PC_2);
y(end+1)=Y_FEF(1,PC_2);
plot(x,y,'-','Color','k','LineWidth',1)
for n_c=1:N_channels
    plot(Y_FEF(n_c,PC_1),Y_FEF(n_c,PC_2),'-o','MarkerFaceColor',color_bin(n_c,:),'MarkerEdgeColor',color_bin(n_c,:),'Color',color_bin(n_c,:),'MarkerSize',15)
end
grid on
lim(1)=min(min(x),min(y));
lim(2)=max(max(x),max(y));
xlim(lim)
ylim(lim)
clear x y z

title(sprintf('FEF CMDS'))
xlabel(sprintf('dim %d',PC_1))
ylabel(sprintf('dim %d',PC_2))
grid on

subplot(3,3,6)
hold on

x=Y_PFC(:,PC_1);
x(end+1)=Y_PFC(1,PC_1);
y=Y_PFC(:,PC_2);
y(end+1)=Y_PFC(1,PC_2);
plot(x,y,'-','Color','k','LineWidth',1)
for n_c=1:N_channels
    plot(Y_PFC(n_c,PC_1),Y_PFC(n_c,PC_2),'-o','MarkerFaceColor',color_bin(n_c,:),'MarkerEdgeColor',color_bin(n_c,:),'Color',color_bin(n_c,:),'MarkerSize',15)
end
grid on
lim(1)=min(min(x),min(y));
lim(2)=max(max(x),max(y));
xlim(lim)
ylim(lim)
clear x y z

title(sprintf('PFC CMDS'))
xlabel(sprintf('dim %d',PC_1))
ylabel(sprintf('dim %d',PC_2))
grid on

subplot(3,3,7)
hold on

x=Y_LIP2(:,PC_1);
x(end+1)=Y_LIP2(1,PC_1);
y=Y_LIP2(:,PC_2);
y(end+1)=Y_LIP2(1,PC_2);
plot(x,y,'-','Color','k','LineWidth',1)
for n_c=1:N_channels
    plot(Y_LIP2(n_c,PC_1),Y_LIP2(n_c,PC_2),'-o','MarkerFaceColor',color_bin(n_c,:),'MarkerEdgeColor',color_bin(n_c,:),'Color',color_bin(n_c,:),'MarkerSize',15)
end
grid on
lim(1)=min(min(x),min(y));
lim(2)=max(max(x),max(y));
xlim(lim)
ylim(lim)
clear x y z

title(sprintf('LIP MDS'))
xlabel(sprintf('dim %d',PC_1))
ylabel(sprintf('dim %d',PC_2))
grid on

subplot(3,3,8)
hold on

x=Y_FEF2(:,PC_1);
x(end+1)=Y_FEF2(1,PC_1);
y=Y_FEF2(:,PC_2);
y(end+1)=Y_FEF2(1,PC_2);
plot(x,y,'-','Color','k','LineWidth',1)
for n_c=1:N_channels
    plot(Y_FEF2(n_c,PC_1),Y_FEF2(n_c,PC_2),'-o','MarkerFaceColor',color_bin(n_c,:),'MarkerEdgeColor',color_bin(n_c,:),'Color',color_bin(n_c,:),'MarkerSize',15)
end
grid on
lim(1)=min(min(x),min(y));
lim(2)=max(max(x),max(y));
xlim(lim)
ylim(lim)
clear x y z

title(sprintf('FEF MDS'))
xlabel(sprintf('dim %d',PC_1))
ylabel(sprintf('dim %d',PC_2))
grid on

subplot(3,3,9)
hold on

x=Y_PFC2(:,PC_1);
x(end+1)=Y_PFC2(1,PC_1);
y=Y_PFC2(:,PC_2);
y(end+1)=Y_PFC2(1,PC_2);
plot(x,y,'-','Color','k','LineWidth',1)
for n_c=1:N_channels
    plot(Y_PFC2(n_c,PC_1),Y_PFC2(n_c,PC_2),'-o','MarkerFaceColor',color_bin(n_c,:),'MarkerEdgeColor',color_bin(n_c,:),'Color',color_bin(n_c,:),'MarkerSize',15)
end
grid on
lim(1)=min(min(x),min(y));
lim(2)=max(max(x),max(y));
xlim(lim)
ylim(lim)
clear x y z

title(sprintf('PFC MDS'))
xlabel(sprintf('dim %d',PC_1))
ylabel(sprintf('dim %d',PC_2))
grid on

%%
figure
subplot(1,3,2)
hold on

x=Y_all(:,PC_1);
x(end+1)=Y_all(1,PC_1);
y=Y_all(:,PC_2);
y(end+1)=Y_all(1,PC_2);
z=Y_all(:,PC_3);
z(end+1)=Y_all(1,PC_3);
plot(x,y,'-','Color','k','LineWidth',1)
for n_c=1:N_channels
    plot(Y_all(n_c,PC_1),Y_all(n_c,PC_2),'-o','MarkerFaceColor',color_bin(n_c,:),'MarkerEdgeColor',color_bin(n_c,:),'Color',color_bin(n_c,:),'MarkerSize',15)
end
grid on
lim(1)=min(min(x),min(y));
lim(2)=max(max(x),max(y));
xlim(lim)
ylim(lim)
clear x y z

title(sprintf('All regions CMDS'))
xlabel(sprintf('dim %d',PC_1))
ylabel(sprintf('dim %d',PC_2))
grid on

subplot(1,3,1)
hold on
x=score_all(:,PC_1);
x(end+1)=score_all(1,PC_1);
y=score_all(:,PC_2);
y(end+1)=score_all(1,PC_2);
z=score_all(:,PC_3);
z(end+1)=score_all(1,PC_3);
plot(x,y,'-','Color','k','LineWidth',1)
for n_c=1:N_channels
    plot(score_all(n_c,PC_1),score_all(n_c,PC_2),'-o','MarkerFaceColor',color_bin(n_c,:),'MarkerEdgeColor',color_bin(n_c,:),'Color',color_bin(n_c,:),'MarkerSize',15)
end
grid on
lim(1)=min(min(x),min(y));
lim(2)=max(max(x),max(y));
xlim(lim)
ylim(lim)
clear x y z

title(sprintf('All regions PCA'))
xlabel(sprintf('PC %d',PC_1))
ylabel(sprintf('PC %d',PC_2))
grid on


subplot(1,3,3)
hold on
x=score_all(:,PC_1);
x(end+1)=score_all(1,PC_1);
y=score_all(:,PC_2);
y(end+1)=score_all(1,PC_2);
z=score_all(:,PC_3);
z(end+1)=score_all(1,PC_3);
plot(x,y,'-','Color','k','LineWidth',1)
for n_c=1:N_channels
    plot(score_all(n_c,PC_1),score_all(n_c,PC_2),'-o','MarkerFaceColor',color_bin(n_c,:),'MarkerEdgeColor',color_bin(n_c,:),'Color',color_bin(n_c,:),'MarkerSize',15)
end
grid on
lim(1)=min(min(x),min(y));
lim(2)=max(max(x),max(y));
xlim(lim)
ylim(lim)
clear x y z

title(sprintf('All regions MDS'))
xlabel(sprintf('PC %d',PC_1))
ylabel(sprintf('PC %d',PC_2))
grid on

%%
figure
subplot(1,3,2)
hold on

x=Y_Rand(:,PC_1);
x(end+1)=Y_Rand(1,PC_1);
y=Y_Rand(:,PC_2);
y(end+1)=Y_Rand(1,PC_2);
z=Y_Rand(:,PC_3);
z(end+1)=Y_Rand(1,PC_3);
plot(x,y,'-','Color','k','LineWidth',1)
for n_c=1:N_channels
    plot(Y_Rand(n_c,PC_1),Y_Rand(n_c,PC_2),'-o','MarkerFaceColor',color_bin(n_c,:),'MarkerEdgeColor',color_bin(n_c,:),'Color',color_bin(n_c,:),'MarkerSize',15)
end
grid on
lim(1)=min(min(x),min(y));
lim(2)=max(max(x),max(y));
xlim(lim)
ylim(lim)
clear x y z

title(sprintf('Rand CMDS'))
xlabel(sprintf('dim %d',PC_1))
ylabel(sprintf('dim %d',PC_2))
grid on

subplot(1,3,1)
hold on
x=score_Rand(:,PC_1);
x(end+1)=score_Rand(1,PC_1);
y=score_Rand(:,PC_2);
y(end+1)=score_Rand(1,PC_2);
z=score_Rand(:,PC_3);
z(end+1)=score_Rand(1,PC_3);
plot(x,y,'-','Color','k','LineWidth',1)
for n_c=1:N_channels
    plot(score_Rand(n_c,PC_1),score_Rand(n_c,PC_2),'-o','MarkerFaceColor',color_bin(n_c,:),'MarkerEdgeColor',color_bin(n_c,:),'Color',color_bin(n_c,:),'MarkerSize',15)
end
grid on
lim(1)=min(min(x),min(y));
lim(2)=max(max(x),max(y));
xlim(lim)
ylim(lim)
clear x y z

title(sprintf('Rand PCA'))
xlabel(sprintf('PC %d',PC_1))
ylabel(sprintf('PC %d',PC_2))
grid on

subplot(1,3,3)
hold on
x=score_Rand(:,PC_1);
x(end+1)=score_Rand(1,PC_1);
y=score_Rand(:,PC_2);
y(end+1)=score_Rand(1,PC_2);
z=score_Rand(:,PC_3);
z(end+1)=score_Rand(1,PC_3);
plot(x,y,'-','Color','k','LineWidth',1)
for n_c=1:N_channels
    plot(score_Rand(n_c,PC_1),score_Rand(n_c,PC_2),'-o','MarkerFaceColor',color_bin(n_c,:),'MarkerEdgeColor',color_bin(n_c,:),'Color',color_bin(n_c,:),'MarkerSize',15)
end
grid on
lim(1)=min(min(x),min(y));
lim(2)=max(max(x),max(y));
xlim(lim)
ylim(lim)
clear x y z

title(sprintf('Rand MDS'))
xlabel(sprintf('PC %d',PC_1))
ylabel(sprintf('PC %d',PC_2))
grid on



%% 

color_bin2=rgb2hsv(color_bin);
color_bin2(:,2)=color_bin2(:,2)-0.2;
color_bin2=hsv2rgb(color_bin2);

figure
hold on
x=score_all(:,PC_1);
x(end+1)=score_all(1,PC_1);
y=score_all(:,PC_2);
y(end+1)=score_all(1,PC_2);
z=score_all(:,PC_3);
z(end+1)=score_all(1,PC_3);
plot(x,y,'-','Color','k','LineWidth',1)
for n_c=1:N_channels
    plot(score_all(n_c,PC_1),score_all(n_c,PC_2),'-o','MarkerFaceColor',color_bin(n_c,:),'MarkerEdgeColor',color_bin(n_c,:),'Color',color_bin(n_c,:),'MarkerSize',15)
end
grid on
lim(1)=min(min(x),min(y))-5;
lim(2)=max(max(x),max(y))+5;
xlim(lim)
ylim(lim)
clear x y z
title(sprintf('All regions MDS'))
xlabel(sprintf('dim %d',PC_1))
ylabel(sprintf('dim %d',PC_2))
grid on


