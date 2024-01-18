function models_FR_RR_value_function(fsroot, ROI, arrayID, this_time)

% for i=1:7
%     initial_window=0;
%     event_list{i}='reward_end';
%     window_start_list(i)=initial_window+(i-1)*50;
% end
% 
% for i=1:20
%     initial_window=-600;
%     event_list{i+7}='target';
%     window_start_list(i+7)=initial_window+(i-1)*50;
% end
% 
% event=event_list{this_time};
% window_size=300;

event='target';
window_start_list=-600;
window_size=900;

task='Learning_Attentional_Templates';
% subtask='exploreexploit/Reset_RW_model';
subtask='exploreexploit/Restricted_FEF/Reset_RW_model';


subsubtask=sprintf('FR_%s_%d_window_from_%d_to_%d',event,window_size,window_start_list(this_time),window_start_list(this_time)+window_size);

dirstem = fullfile('Projects',task,'Analysis','Electrophy_analysis',subtask,subsubtask);
data_path = fullfile(fsroot,dirstem);

savename=sprintf('models_FR_RR_value_function_%s_%d',ROI,arrayID);

N_param=4;
N_fold=10;
N_models=5;

count_ROI=1;

%% session

if strncmp(ROI,'LIP',3) && arrayID==34
    return
else
    
    data_name=sprintf('ID_%d.mat',arrayID);
    
    switch ROI
        case 'LIP'
            load(fullfile(data_path,data_name),'Conditions','LIP')
            this_ROI=LIP;
            clear LIP
        case 'FEF'
            load(fullfile(data_path,data_name),'Conditions','FEF')
            this_ROI=FEF;
            clear FEF
        case 'PFC'
            load(fullfile(data_path,data_name),'Conditions','PFC')
            this_ROI=PFC;
            clear PFC
    end
    
    for n=1:size(this_ROI,1)
        for t=1:size(this_ROI,2)
            if sum(~isnan(this_ROI(n,t,:)))>0
                this_ROI(n,t,isnan(this_ROI(n,t,:)))=0;
            end
        end
    end
    
    belief_probability=zeros(size(Conditions.belief,2),size(Conditions.belief,1));
    for i=1:size(Conditions.stim,2)
        if sum(Conditions.belief(:,i)==0)<size(Conditions.belief(:,i),1)
            belief_probability(i,:)=(Conditions.belief(:,i)-min(Conditions.belief(:,i)))/sum(Conditions.belief(:,i));
        else
            belief_probability(i,:)=Conditions.belief(:,i);
        end
    end
    belief=Conditions.belief';
    belief_precision=Conditions.belief_precision;
    belief_precision(isnan(belief_precision))=0;
    mean_value=mean(Conditions.belief,1);
    
    N_bins=100;
    color_binned =  0:(2*pi/N_bins):(2*pi-(2*pi/N_bins));
    belief_angle = color_binned+2*pi/N_bins/2; %center of the bin
    
    for i=1:size(Conditions.stim,2)
        [~, peak_belief_index(i)]=max(Conditions.belief(:,i));
    end
    peak_belief=belief_angle(peak_belief_index);
    for i=1:size(Conditions.stim,2)
        if sum(Conditions.belief(:,i)==0)==size(Conditions.belief,1)
            peak_belief(i)=NaN;
        end
    end
    
    value=zeros(size(Conditions.stim));
    Conditions.chosen_value=zeros(size(Conditions.chosen_color));
    %bin the stim color
    for i=1:size(Conditions.stim,2)
        for j=1:size(Conditions.stim,1)
            if isnan(Conditions.stim(j,i))
                value(j,i)=NaN;
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
                value(j,i)=Conditions.belief(this_bin,i);
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
    reward=Conditions.Reward;
    mean_belief=Conditions.mean_belief;
    chosen_color=Conditions.chosen_color;
    
    chosen_value=Conditions.chosen_value;
    
    switch event
        case 'target'
            prev_chosen_color(1)=NaN;
            prev_rpe(1)=NaN;
            prev_reward(1)=NaN;
            prev_chosen_value(1)=NaN;
            for i=2:length(Conditions.stim)
                prev_chosen_color(i)=Conditions.chosen_color(i-1);
                prev_rpe(i)=Conditions.RPE(i); %index misaligned
                prev_reward(i)=Conditions.Reward(i-1);
                prev_chosen_value(i)=chosen_value(i-1);
            end
            for i=1:length(Conditions.stim)-1 %index misaligned
                rpe(i)=Conditions.RPE(i+1);
            end
            rpe(length(Conditions.stim))=NaN;
        case 'reward_end'
            for i=1:length(Conditions.stim)-1
                belief_probability(i,:)=belief_probability(i+1,:);
                belief(i,:)=belief(i+1,:);
                mean_belief(i)=mean_belief(i+1);
                chosen_color(i)=Conditions.chosen_color(i+1);
                reward(i)=Conditions.Reward(i+1);
                chosen_value(i)=Conditions.chosen_value(i+1);
                mean_value(i)=mean_value(i+1);
            end
            for i=1:length(Conditions.stim)-2
                rpe(i)=Conditions.RPE(i+2);
            end
            belief_probability(length(Conditions.stim))=NaN;
            belief(length(Conditions.stim))=NaN;
            mean_belief(length(Conditions.stim))=NaN;
            mean_value(length(Conditions.stim))=NaN;
            chosen_color(length(Conditions.stim))=NaN;
            reward(length(Conditions.stim))=NaN;
            chosen_value(length(Conditions.stim))=NaN;
            rpe(length(Conditions.stim)-1:length(Conditions.stim))=NaN;
            for i=1:length(Conditions.stim)
                prev_chosen_color(i)=Conditions.chosen_color(i);
                
                prev_reward(i)=Conditions.Reward(i);
                prev_chosen_value(i)=Conditions.chosen_value(i);
            end
            for i=1:length(Conditions.stim)-1
                prev_rpe(i)=Conditions.RPE(i+1);
            end
            prev_rpe(length(Conditions.stim))=NaN;
    end
    
    for n=1:size(this_ROI,1)
        
        if sum(~isnan(this_ROI(n,:,1)),2)>500
            y=this_ROI(n,:);
            for i=1:size(y,1)
                if sum(~isnan(y(i)))>0
                    y(i,isnan(y(i)))=0;
                end
            end
            y(1)=NaN; %remove first trial as no prev
            y(end-1:end)=NaN; %remove last 2 trials as no next
            
            fr=y(:)';
            
            selected=logical(~isnan(fr));
%             Belief_this_neuron=belief_probability(selected,:);
            Belief_this_neuron=belief(selected,:);
            Peak_belief_this_neuron=peak_belief(selected)';
            Mean_belief_this_neuron=mean_belief(selected)';
            chosen_color_this_neuron=chosen_color(selected)';
            chosen_value_this_neuron=chosen_value(selected)';
            reward_this_neuron=reward(selected)';
            rpe_this_neuron=rpe(selected)';
            surprise_this_neuron=abs(rpe_this_neuron);
            prev_reward_this_neuron=prev_reward(selected)';
            prev_rpe_this_neuron=prev_rpe(selected)';
            prev_surprise_this_neuron=abs(prev_rpe_this_neuron);
            prev_chosen_color_this_neuron=prev_chosen_color(selected)';
            prev_chosen_value_this_neuron=prev_chosen_value(selected)';
            belief_precision_this_neuron=belief_precision(selected)';
            fr=fr(selected);
            Mean_value_this_neuron=mean_value(selected)';
            
            fr=zscore(fr);
            switch event
                case 'target'
                    [R2(count_ROI,:,:),R2_baseline(count_ROI,:),param(count_ROI,:,:,:),partition_saved(count_ROI)] = compute_log_model_evidence(prev_chosen_color_this_neuron, ...
                        Belief_this_neuron, Peak_belief_this_neuron, Mean_value_this_neuron, fr',N_param, N_fold, N_models, prev_reward_this_neuron);
                case 'reward_end'
                    [R2(count_ROI,:,:),R2_baseline(count_ROI,:),param(count_ROI,:,:,:),partition_saved(count_ROI)] = compute_log_model_evidence(prev_chosen_color_this_neuron, ...
                        Belief_this_neuron, Peak_belief_this_neuron, Mean_value_this_neuron, fr',N_param, N_fold, N_models, prev_reward_this_neuron);
            end
            for j=1:size(R2,2)
                for k=1:size(R2,3)
                    R2_diff(count_ROI,j,k)=R2(count_ROI,j,k)-R2_baseline(count_ROI,k);
                end
            end
            
            Session_id(count_ROI)=arrayID;
            Neuron_id(count_ROI)=n;
            
            clear y x fr *_this_neuron
            
            count_ROI=count_ROI+1;
            
        end
        
    end
    
    if exist('R2')>0
        save(fullfile(data_path,savename),'R2','R2_baseline','R2_diff','param','partition_saved','Neuron_id','Session_id','-v7.3')
    end
    
end

end


function [R2,R2_baseline,param,c,exitflag] = compute_log_model_evidence(Chosen_color, Belief, Peak_belief, Mean_value, y, Nparam, N_fold, N_models, prev_rew)

c = cvpartition(length(y),'KFold',N_fold);

for np=1:N_fold
    
    [R2(:,np),R2_baseline(np),param(:,:,np),exitflag(:,np)]=regf(Chosen_color(c.training(np)), Belief(c.training(np),:), Peak_belief(c.training(np)), Mean_value(c.training(np)), y(c.training(np)),...
        Chosen_color(c.test(np)), Belief(c.test(np),:), Peak_belief(c.test(np)), Mean_value(c.test(np)), y(c.test(np)), Nparam,N_models, prev_rew(c.training(np)), prev_rew(c.test(np)));
    
end

    function [R2,R2_baseline,param,exitflag] = regf(Chosen_color, Belief, Peak_belief, Mean_value, y,...
            Chosen_color_test,  Belief_test, Peak_belief_test, Mean_value_test, ytest, Nparam,Nmodels, prev_rew, prev_rew_test)
        
        R2=zeros(Nmodels,1);
        param=zeros(Nparam,Nmodels);
        
        %Peak belief
        [R2(1),param(1:4,1),exitflag(1)] = regf_stim(Peak_belief,y,Peak_belief_test,ytest);
        %Expected value model
        [R2(2),param(1:4,2),exitflag(2)] = regf_belief(Belief,y,Belief_test,ytest);
        %Mean value
        [R2(3),param(1:2,3),exitflag(3)] = regf_value(Mean_value,y,Mean_value_test,ytest);
        %Chosen color
        [R2(4),param(1:4,4),exitflag(4)] = regf_stim(Chosen_color,y,Chosen_color_test,ytest);
        %Prev chosen_color if rew > 0
        [R2(5),param(1:4,5),exitflag(4)] = regf_stim_if(Chosen_color, prev_rew,y,Chosen_color_test, prev_rew_test,ytest);
        
        mse_baseline_model=1/length(ytest).*(sum((ytest-mean(y)).^2));
        mse_null_model=1/length(ytest).*(sum((ytest-mean(ytest)).^2));
        R2_baseline=1-(mse_baseline_model/mse_null_model);
        
        
        function [R2,param,exitflag] = regf_stim(stim,y,stim_test,ytest) %vm
            
            x0=zeros(4,1);
            ub=[10 10 2*pi 2.5];
            lb=[-10 -10 -2*pi -2.5];
            
            options = optimoptions('fmincon','Display','off');
            
            fun_fit=@(P)f_full_model_stim_vm(P,y,stim);
            
            [param,~,exitflag] = fmincon(fun_fit,x0,[],[],[],[],lb,ub,[],options);
            
            mse=f_full_model_stim_vm(param,ytest,stim_test);
            mse_null=1/length(ytest).*(sum((ytest-mean(ytest)).^2));
            R2=1-mse/mse_null;
            
        end
        
        function [R2,param,exitflag] = regf_belief(Belief,y,Belief_test,ytest) %sum vm weighted by P(b)
            
            x0=zeros(4,1);
            ub=[10 10 2*pi 2.5];
            lb=[-10 -10 -2*pi -2.5];
            
            options = optimoptions('fmincon','Display','off');
            
            fun_fit=@(P)f_full_model_belief_vm(P,y,Belief);
            
            [param,~,exitflag,~] = fmincon(fun_fit,x0,[],[],[],[],lb,ub,[],options);
            
            mse=f_full_model_belief_vm(param,ytest,Belief_test);
            mse_null=1/length(ytest).*(sum((ytest-mean(ytest)).^2));
            R2=1-mse/mse_null;
            
        end
        
        function [R2,param,exitflag] = regf_value(value,y,value_test,ytest)  %linear
            
            x0=zeros(2,1);
            ub=[10 10];
            lb=[-10 -10];
            
            options = optimoptions('fmincon','Display','off');
            
            fun_fit=@(P)f_full_model_value(P,y,value);
            
            [param,~,exitflag,~] = fmincon(fun_fit,x0,[],[],[],[],lb,ub,[],options);
            
            mse=f_full_model_value(param,ytest,value_test);
            mse_null=1/length(ytest).*(sum((ytest-mean(ytest)).^2));
            R2=1-mse/mse_null;
            
        end
        
          function [R2,param,exitflag] = regf_stim_if(stim,rew,y,stim_test, rew_test,ytest) %vm
            
            x0=zeros(4,1);
            ub=[10 10 2*pi 2.5];
            lb=[-10 -10 -2*pi -2.5];
            
            options = optimoptions('fmincon','Display','off');
            
            fun_fit=@(P)f_full_model_stim_vm_if(P,y,stim, rew);
            
            [param,~,exitflag] = fmincon(fun_fit,x0,[],[],[],[],lb,ub,[],options);
            
            mse=f_full_model_stim_vm_if(param,ytest,stim_test,rew_test);
            mse_null=1/length(ytest).*(sum((ytest-mean(ytest)).^2));
            R2=1-mse/mse_null;
            
        end

        
    end

end

