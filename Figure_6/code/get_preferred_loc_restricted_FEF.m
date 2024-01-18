function get_preferred_loc_restricted_FEF(fsroot, arrayID, ROI)

this_time=10; %[50 250ms] window

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
savename=sprintf('get_pref_loc_%s_%d',ROI,arrayID);

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
                
        if sum(~isnan(this_ROI(n,:)),2)>500 && sum((this_ROI(n,:))==0,2)<size(this_ROI,2)
            %z-score
            y=(this_ROI(n,:)-nanmean(this_ROI(n,:)))/nanstd(this_ROI(n,:));
            
            %get the preferred loc
            loc1=~isnan(Conditions.value(1,:));
            loc2=~isnan(Conditions.value(2,:));
            loc3=~isnan(Conditions.value(3,:));
            loc4=~isnan(Conditions.value(4,:));
            
            m=fitglm(loc1,y);
            T_loc(1,count_ROI)=abs(m.Coefficients.tStat(2));
            p_loc(1,count_ROI)=m.Coefficients.pValue(2);
            m=fitglm(loc2,y);
            T_loc(2,count_ROI)=abs(m.Coefficients.tStat(2));
            p_loc(2,count_ROI)=m.Coefficients.pValue(2);
            m=fitglm(loc3,y);
            T_loc(3,count_ROI)=abs(m.Coefficients.tStat(2));
            p_loc(3,count_ROI)=m.Coefficients.pValue(2);
            m=fitglm(loc4,y);
            T_loc(4,count_ROI)=abs(m.Coefficients.tStat(2));
            p_loc(4,count_ROI)=m.Coefficients.pValue(2);
            [~, Preferred_loc(count_ROI)]=max(T_loc(:,count_ROI));
            if p_loc(Preferred_loc(count_ROI),count_ROI)>=0.05
                Preferred_loc(count_ROI)=NaN;
            end
            
            Session_id(count_ROI)=arrayID;
            Neuron_id(count_ROI)=n;
            
            clear y x FR *_this_neuron
            
        else
            
            Preferred_loc(count_ROI)=NaN;
            T_loc(1:4,count_ROI)=NaN;
            p_loc(1:4,count_ROI)=NaN;
            Session_id(count_ROI)=arrayID;
            Neuron_id(count_ROI)=n;
            
        end
        
        count_ROI=count_ROI+1;
    end
    clear *_blocks id* Conditions this*
    
    
    
    
    %% save
    
    save(fullfile(save_path,savename),'*_loc','*_id')
    
end


end

