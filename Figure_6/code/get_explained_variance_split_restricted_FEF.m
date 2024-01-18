function get_explained_variance_split_restricted_FEF(fsroot, arrayID, ROI, this_time)

N_shuffles=250;

% %Test param
% clear all
% ROI='LIP';
% fsroot='/Volumes/buschman';
% arrayID=3;
% N_shuffles=10;
% this_time=25;

event='target';
window_size=200;
for i=1:25
    initial_window=-400;
    event_list{i}='target';
    window_start_list(i)=initial_window+(i-1)*50;
end

task='Learning_Attentional_Templates';
subtask='exploreexploit/Restricted_FEF/Reset_RW_model';

subsubtask_save=sprintf('Classifier_%s_%d_to_%s_%d',event,window_start_list(1),event,window_start_list(end));
dirstem2 = fullfile('Projects',task,'Analysis','Electrophy_analysis',subtask,subsubtask_save);
save_path=fullfile(fsroot,dirstem2);
mkdir(save_path)
savename=sprintf('glm_explained_variance_split_%s_%d_%d',ROI,arrayID,this_time);

count_ROI=1;


%% session

if (strncmp(ROI,'LIP',3) && (arrayID==34 || arrayID==30)) || (strncmp(ROI,'FEF',3) && arrayID==14)
    
else
    
    data_name=sprintf('ID_%d.mat',arrayID);
    
    %this is where the waveform is saved
    subsubtask=sprintf('FR_%s_%d_window_from_%d_to_%d','target',300,-300,-300+300);
    dirstem = fullfile('Projects',task,'Analysis','Electrophy_analysis',subtask,subsubtask);
    data_path = fullfile(fsroot,dirstem);
    
    switch ROI
        case 'LIP'
            load(fullfile(data_path,data_name),'Conditions','LIP_wv')
            this_Wv=LIP_wv;
            clear LIP_wv
        case 'FEF'
            load(fullfile(data_path,data_name),'Conditions','FEF_wv')
            this_Wv=FEF_wv;
            clear FEF_wv
        case 'PFC'
            load(fullfile(data_path,data_name),'Conditions','PFC_wv')
            this_Wv=PFC_wv;
            clear PFC_wv
    end
    
    N_bins=100;
    color_binned =  0:(2*pi/N_bins):(2*pi-(2*pi/N_bins));
    
    Conditions.value=zeros(size(Conditions.stim));
    
    %get the value at each location (reorganize the data)
    N_bins=100;
    color_binned =  0:(2*pi/N_bins):(2*pi-(2*pi/N_bins));
    
    Conditions.value=zeros(size(Conditions.stim));
    Conditions.chosen_value=zeros(size(Conditions.chosen_color));
    %bin the stim color
    for i=1:size(Conditions.stim,2)
        for j=1:size(Conditions.stim,1)
            if isnan(Conditions.stim(j,i))
                Conditions.value(j,i)=NaN;
            else
                if Conditions.stim(j,i)>=color_binned(N_bins)
                    this_bin=N_bins;
                else
                    for k=1:N_bins-1
                        if Conditions.stim(j,i)>=color_binned(k) && Conditions.stim(j,i)<color_binned(k+1)
                            this_bin=k;
                        end
                    end
                end
                Conditions.value(j,i)=Conditions.belief(this_bin,i);
            end
        end
        if Conditions.chosen_color(i)>=color_binned(N_bins)
            this_bin=N_bins;
        else
            for k=1:N_bins-1
                if Conditions.chosen_color(i)>=color_binned(k) && Conditions.chosen_color(i)<color_binned(k+1)
                    this_bin=k;
                end
            end
        end
        Conditions.chosen_value(i)=Conditions.belief(this_bin,i);
    end
    for i=1:size(Conditions.value,2)
        sorted_value=sort(Conditions.value(~isnan(Conditions.value(:,i)),i));
        if sorted_value(3)==Conditions.chosen_value(i)
            Conditions.unchosen_value_best(i)=sorted_value(2);
            Conditions.unchosen_value_worst(i)=sorted_value(1);
        elseif sorted_value(2)==Conditions.chosen_value(i)
            Conditions.unchosen_value_best(i)=sorted_value(3);
            Conditions.unchosen_value_worst(i)=sorted_value(1);
        else
            Conditions.unchosen_value_best(i)=sorted_value(3);
            Conditions.unchosen_value_worst(i)=sorted_value(2);
        end
    end
    Conditions.mean_value=mean(Conditions.belief,1);
    
    %let's get the FR at this time
    subsubtask=sprintf('FR_%s_%d_window_from_%d_to_%d',event_list{this_time},window_size,window_start_list(this_time),window_start_list(this_time)+window_size);
    dirstem = fullfile('Projects',task,'Analysis','Electrophy_analysis',subtask,subsubtask);
    data_path = fullfile(fsroot,dirstem);
    
    switch ROI
        case 'LIP'
            load(fullfile(data_path,data_name),'LIP')
            this_ROI_now=LIP;
            clear LIP
        case 'FEF'
            load(fullfile(data_path,data_name),'FEF')
            this_ROI_now=FEF;
            clear FEF
        case 'PFC'
            load(fullfile(data_path,data_name),'PFC')
            this_ROI_now=PFC;
            clear PFC
    end
    %get the baseline
    subsubtask=sprintf('FR_%s_%d_window_from_%d_to_%d',event_list{1},window_size,window_start_list(1),window_start_list(1)+window_size);
    dirstem = fullfile('Projects',task,'Analysis','Electrophy_analysis',subtask,subsubtask);
    data_path = fullfile(fsroot,dirstem);
    
    switch ROI
        case 'LIP'
            load(fullfile(data_path,data_name),'LIP')
            this_ROI_baseline=LIP;
            clear LIP
        case 'FEF'
            load(fullfile(data_path,data_name),'FEF')
            this_ROI_baseline=FEF;
            clear FEF
        case 'PFC'
            load(fullfile(data_path,data_name),'PFC')
            this_ROI_baseline=PFC;
            clear PFC
    end
    %remove the baseline
    for n=1:size(this_ROI_baseline,1)
        this_ROI(n,:)=this_ROI_now(n,:)-this_ROI_baseline(n,:);
    end
    
    clear this_ROI_now this_ROI_baseline
    
    %remove NaNs when FR was just 0
    for n=1:size(this_ROI,1)
        for t=1:size(this_ROI,2)
            if sum(~isnan(this_ROI(n,t,:)))>0
                this_ROI(n,t,isnan(this_ROI(n,t,:)))=0;
            end
        end
    end
    
    
    for n=1:size(this_ROI,1)
        
        this_value_loc1=Conditions.value(1,:);
        this_choice_loc1=(Conditions.choice==1);
        this_value_loc2=Conditions.value(2,:);
        this_choice_loc2=(Conditions.choice==2);
        this_value_loc3=Conditions.value(3,:);
        this_choice_loc3=(Conditions.choice==3);
        this_value_loc4=Conditions.value(4,:);
        this_choice_loc4=(Conditions.choice==4);
        
        this_value_loc=vertcat(this_value_loc1,this_value_loc2,this_value_loc3,this_value_loc4);
        this_choice_loc=vertcat(this_choice_loc1,this_choice_loc2,this_choice_loc3,this_choice_loc4);
        
        this_reward=Conditions.Reward;
        this_chosen_value=Conditions.chosen_value;
        this_unchosen_value_best=Conditions.unchosen_value_best;
        this_mean_value=Conditions.mean_value;
        
        this_prev_reward(1)=NaN;
        this_prev_chosen_value(1)=NaN;
        this_prev_unchosen_value_best(1)=NaN;
        for i=2:length(Conditions.stim)
            this_prev_reward(i)=Conditions.Reward(i-1);
            this_prev_chosen_value(i)=Conditions.chosen_value(i-1);
            this_prev_unchosen_value_best(i)=Conditions.unchosen_value_best(i-1);
        end
        
        if sum(~isnan(this_ROI(n,:)),2)>500 && sum((this_ROI(n,:))==0,2)<size(this_ROI,2)
            %z-score
            y=(this_ROI(n,:)-nanmean(this_ROI(n,:)))/nanstd(this_ROI(n,:));
            
            %reps to make sure there isn't anything special abou the split
            for i=1:10
                
                %split the data in 2
                c=cvpartition(length(y),'Holdout',0.5);
                partition(count_ROI).c(i,1).test=c.training;
                partition(count_ROI).c(i,2).test=c.test;
                
                for p=1:2 %split
                    %How much varaince was there to explain in each half?
                    for j=1:4
                        Explainable_variance_at_loc(count_ROI,j,p,i) = sum((y(partition(count_ROI).c(i,p).test' & ~isnan(this_value_loc(j,:))) - mean(y(partition(count_ROI).c(i,p).test' & ~isnan(this_value_loc(j,:))))).^2);
                    end
                    %Run the glm for each loc
                    for j=1:4
                        mdl=fitglm([this_choice_loc(j,partition(count_ROI).c(i,p).test); (this_choice_loc(j,partition(count_ROI).c(i,p).test)==1).*this_value_loc(j,partition(count_ROI).c(i,p).test); (this_choice_loc(j,partition(count_ROI).c(i,p).test)==0).*this_value_loc(j,partition(count_ROI).c(i,p).test)]',y(partition(count_ROI).c(i,p).test)');
                        beta_split(count_ROI,:,j,p,i)=mdl.Coefficients.Estimate;
                        clear mdl
                    end
                    %Run the glm for all the locs with confounds
                    mdl=fitglm([this_prev_reward(partition(count_ROI).c(i,p).test); this_reward(partition(count_ROI).c(i,p).test); this_prev_chosen_value(partition(count_ROI).c(i,p).test); this_prev_unchosen_value_best(partition(count_ROI).c(i,p).test); this_chosen_value(partition(count_ROI).c(i,p).test); this_unchosen_value_best(partition(count_ROI).c(i,p).test); this_mean_value(partition(count_ROI).c(i,p).test)]',y(partition(count_ROI).c(i,p).test)'); %
                    beta_all_values_split(count_ROI,:,p,i)=mdl.Coefficients.Estimate;
                    clear mdl
                end
                
                %What is the variance explained by the chosen value
                %regressor (within loc and split and across)?
                for p=1:2
                    for q=1:2
                        for j=1:4
                            for k=1:4
                                % Here we replace the chosen value regressor by
                                % that at location k but for each split
                                Explained_variance_across_loc(count_ROI,j,k,p,q,i) = nansum((y(partition(count_ROI).c(i,p).test) - (beta_split(count_ROI,1,j,p,i) + beta_split(count_ROI,2,j,p,i).*this_choice_loc(j,partition(count_ROI).c(i,p).test)  + beta_split(count_ROI,3,k,q,i).*(this_choice_loc(j,partition(count_ROI).c(i,p).test)==1).*this_value_loc(j,partition(count_ROI).c(i,p).test) + beta_split(count_ROI,4,j,p,i).*(this_choice_loc(j,partition(count_ROI).c(i,p).test)==0).*this_value_loc(j,partition(count_ROI).c(i,p).test))).^2);
                                
                                %We shuffle to look at the specific explained
                                %variance by the chosen value regressor
                                for np=1:N_shuffles
                                    y_shuffle=y(partition(count_ROI).c(i,p).test);
                                    choice_shuffle=this_choice_loc(j,partition(count_ROI).c(i,p).test);
                                    value_shuffle=this_value_loc(j,partition(count_ROI).c(i,p).test);
                                    chosen_value_pre_shuffle=(choice_shuffle==1).*value_shuffle;
                                    unchosen_value_shuffle=(choice_shuffle==0).*value_shuffle;
                                    % we want to shuffle only if choice=1
                                    % and isnan(value)
                                    ind_c=find(choice_shuffle==1 & ~isnan(chosen_value_pre_shuffle));
                                    ind_s=ind_c(randperm(length(ind_c)));
                                    chosen_value_shuffle=chosen_value_pre_shuffle;
                                    for ind=1:length(ind_c)
                                        chosen_value_shuffle(ind_s(ind))=chosen_value_pre_shuffle(ind_c(ind));
                                    end
                                    Explained_variance_across_loc_shuffled(count_ROI,j,k,p,q,i,np) = nansum((y_shuffle-(beta_split(count_ROI,1,j,p,i) + beta_split(count_ROI,2,j,p,i).*choice_shuffle  + beta_split(count_ROI,3,k,q,i).*chosen_value_shuffle + beta_split(count_ROI,4,j,p,i).*unchosen_value_shuffle)).^2);
                                end
                            end
                        end
                    end
                end
            end
                                
            Session_id(count_ROI)=arrayID;
            Neuron_id(count_ROI)=n;
            
            clear y x FR *_this_neuron
            
        else
            
            Explainable_variance_at_loc(count_ROI,1:4,1:2,1:10)=NaN;
            beta_split(count_ROI,1:4,1:4,1:2,1:10)=NaN;
            beta_all_values_split(count_ROI,1:8,1:2,1:10)=NaN;
            Explained_variance_across_loc(count_ROI,1:4,1:4,1:2,1:2,1:10)=NaN;
            Explained_variance_across_loc_shuffled(count_ROI,1:4,1:4,1:2,1:2,1:10,1:N_shuffles)=NaN;
            
            Session_id(count_ROI)=arrayID;
            Neuron_id(count_ROI)=n;
            
        end
        
        count_ROI=count_ROI+1;
    end
    clear *_blocks id* Conditions this*
    
    
    
    %% Get the r square
    
    for n=1:size(Explained_variance_across_loc,1) %neuron
        for i=1:10 %rep
            for j=1:4 %loc
                for k=1:4 %loc
                    for s=1:2 %split
                        for p=1:2 %split
                            if Explainable_variance_at_loc(n,j,s,i)<1 %not enough expalinable variance
                                R_square_across_loc(n,j,k,s,p,i)=NaN;
                                z_R_square_across_loc(n,j,k,s,p,i)=NaN;
                                R_square_across_loc_shuffle(n,j,k,s,p,i,1:N_shuffles)=NaN;
                            else
                                R_square_across_loc(n,j,k,s,p,i)=1-Explained_variance_across_loc(n,j,k,s,p,i)/Explainable_variance_at_loc(n,j,s,i);
                                R_square_across_loc_shuffle(n,j,k,s,p,i,:)=1-Explained_variance_across_loc_shuffled(n,j,k,s,p,i,:)/Explainable_variance_at_loc(n,j,s,i);
                                z_R_square_across_loc(n,j,k,s,p,i)=(R_square_across_loc(n,j,k,s,p,i)-mean(R_square_across_loc_shuffle(n,j,k,s,p,i,:)))/std(R_square_across_loc_shuffle(n,j,k,s,p,i,:));
                            end
                        end
                    end
                end
            end
        end
        
        %Mean across repetitions
        %within loc
        for j=1:4
            for s=1:2 %split
                for p=1:2 %split
                    Mean_z_R_square_within_loc(n,j,s,p)=mean(z_R_square_across_loc(n,j,j,s,p,:));
                end
            end
            Mean_z_R_square_within_loc_across_split(n,j)=(Mean_z_R_square_within_loc(n,j,1,2)+Mean_z_R_square_within_loc(n,j,2,1))/2;
        end
        
        %across locs
        for s=1:2 %split
            for p=1:2 %split
                for j=1:4
                    for k=1:4
                        Mean_z_R_square_across_loc(n,j,k,s,p)=mean(z_R_square_across_loc(n,j,k,s,p,:));
                    end
                end
            end
        end
        
        %Compare to how much explained variance could have been explained
        %by this regressor = same split and same loc
        for j=1:4
            for k=1:4
                Mean_z_R_square_explainable_across_loc(n,j,k)=(Mean_z_R_square_across_loc(n,j,k,1,2)/Mean_z_R_square_across_loc(n,j,j,1,1) + Mean_z_R_square_across_loc(n,j,k,2,1)/Mean_z_R_square_across_loc(n,j,j,2,2))/2;
            end
        end
        
        
    end
    
    if exist('partition','var')==0
        partition=NaN;
    end
    
    
    %% save
    
    save(fullfile(save_path,savename),'beta*','partition*','Explain*','*R_square*','*_id','Mean*')
    
end


end

