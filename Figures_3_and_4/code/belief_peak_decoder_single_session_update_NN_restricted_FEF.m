function belief_peak_decoder_single_session_update_NN_restricted_FEF(fsroot, arrayID, event, this_time)
% fsroot='/Volumes/buschman';
% arrayID=20;
% event='target';
% this_time=1;

task='Learning_Attentional_Templates';
subtask='exploreexploit/Restricted_FEF/Reset_RW_model';

window_start_list=-600;
window_size=900;

subsubtask=sprintf('FR_%s_%d_window_from_%d_to_%d',event,window_size,window_start_list(this_time),window_start_list(this_time)+window_size);

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
elseif arrayID==14
    this_ROI=[LIP; PFC];
else
    this_ROI=[LIP; FEF; PFC];
end
clear LIP FEF PFC

for n=1:size(this_ROI,1)
    for t=1:size(this_ROI,2)
        if sum(~isnan(this_ROI(n,t,:)))>0
            this_ROI(n,t,isnan(this_ROI(n,t,:)))=0;
        end
    end
end
    
    %Extract peak belief
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
    %Bin the peak belief
    for i=1:size(Conditions.stim,2)
        for n_c=1:N_channels
            if mod(peak_belief(i),2*pi)>=(n_c-1)*2*pi/N_channels && mod(peak_belief(i),2*pi)<n_c*2*pi/N_channels
                peak_belief_bin(i)=n_c;
            end
        end
    end    
    
    for n_c=1:N_channels
        angle_mean_belief_bin(n_c)=((n_c-1)*2*pi/N_channels + n_c*2*pi/N_channels)/2;
    end
    
    %Adapt to event
    switch event
        case 'target'
            prev_chosen_color(1)=NaN;
            prev_rpe(1)=NaN;
            prev_reward(1)=NaN;
            for i=2:length(Conditions.stim)
                prev_chosen_color(i)=Conditions.chosen_color(i-1);
                prev_rpe(i)=Conditions.RPE(i); %index misaligned
                prev_reward(i)=Conditions.Reward(i-1);
            end
            for i=1:length(Conditions.stim)-1 %index misaligned
                rpe(i)=Conditions.RPE(i+1);
            end
            rpe(length(Conditions.stim))=NaN;
        case 'reward_end'
            for i=1:length(Conditions.stim)-1
                peak_belief(i)=peak_belief(i+1);
                chosen_color(i)=Conditions.chosen_color(i+1);
                reward(i)=Conditions.Reward(i+1);
                
            end
            for i=1:length(Conditions.stim)-2
                rpe(i)=Conditions.RPE(i+2);
            end
            peak_belief(length(Conditions.stim))=NaN;
            chosen_color(length(Conditions.stim))=NaN;
            reward(length(Conditions.stim))=NaN;
            rpe(length(Conditions.stim)-1:length(Conditions.stim))=NaN;
            for i=1:length(Conditions.stim)
                prev_chosen_color(i)=Conditions.chosen_color(i);
                prev_reward(i)=Conditions.Reward(i);
            end
            for i=1:length(Conditions.stim)-1
                prev_rpe(i)=Conditions.RPE(i+1);
            end
            prev_rpe(length(Conditions.stim))=NaN;
    end
    
    %FR
    Y=this_ROI;
    Y=Y(~isnan(sum(Y,2)),:); %remove NaN
    
    %Create test/train
    index_start=find(Conditions.trial_in_block==1);
    data_val=zeros(length(Conditions.trial_in_block),length(index_start));
    for i=1:length(index_start)
        data_val(index_start(i):index_start(i)+34,i)=1;
    end
    data_val=logical(data_val);
    
    data_tt=~data_val;
    % remove no peak from training
    for i=1:size(data_tt,1)
        for j=1:size(data_tt,2)
            if data_tt(i,j)==1 && isnan(peak_belief(i))
                data_tt(i,j)=0;
            end
        end
    end

    Length_tt=[0 0 0];
    for i=1:length(index_start)
        for c=1:N_channels
            Length_tt(c)=sum(peak_belief_bin(data_tt(:,i))==c);
        end
        N_tt=min(Length_tt);
        if mod(N_tt,2)~=0
            N_tt=N_tt-1;
        end
        This_TT=zeros(N_tt,3);
        for c=1:N_channels
            buffer=find(peak_belief_bin(data_tt(:,i))==c);
            buffer=buffer(randperm(length(buffer)));
            This_TT(:,c)=buffer(1:N_tt);
            clear buffer
        end
        
        for c=1:N_channels
            data_this_belief=Y(:,This_TT(:,c));
            increment=0;
            for d=1:N_channels
                if d~=c
                    data_other_belief(:,1+increment*N_tt/2:N_tt/2+N_tt/2*increment)=Y(:,This_TT(1:N_tt/2,d));
                    increment=increment+1;
                end
            end
            groups_classifier(1:N_tt,1)=1;
            groups_classifier(N_tt+1:2*N_tt,1)=0;
            data_classifier=[data_this_belief, data_other_belief]';
            
            cv = cvpartition(size(groups_classifier,1),'KFold',10);
            
            %optimize classfier
            
            opts = struct('Optimizer','bayesopt','ShowPlots',false,'CVPartition',cv,...
                'AcquisitionFunctionName','expected-improvement-plus','Verbose',0,'UseParallel',par_param);
            
            svmmod = fitcsvm(data_classifier,groups_classifier,'KernelFunction','linear',...
                'OptimizeHyperparameters','auto','HyperparameterOptimizationOptions',opts);
            
            if ~isempty(svmmod.Beta)
                
                mdl=fitPosterior(svmmod,data_classifier,groups_classifier);
                %
                W(:,c)=mdl.Beta;
                Intercept(:,c)=mdl.Bias;
                
                SVMModel(c).mdl=mdl;
                
            else
                
                W(:,c)=zeros(N_neurons,1);
                Intercept(:,c)=svmmod.Bias;
                
                SVMModel(c).mdl=svmmod;
                
            end
            
            clear cv *_this_belief *_classifier mdl *_other_belief*
        end
        %Now use a NN to determine the relationship between classifiers
        %scores and angles
        Results(i).True_peak_belief_tt=peak_belief(data_tt(:,i));
        Results(i).Coordinates_peak_belief_tt=[cos(peak_belief(data_tt(:,i))); sin(peak_belief(data_tt(:,i)))];
        for c=1:N_channels
            [~, probas,~] = predict(SVMModel(c).mdl,Y(:,data_tt(:,i))');
            Results(i).Classification_tt(c,:)=probas(:,2);
            clear probas
        end
        net=feedforwardnet(10,'trainlm');
                net.layers{1}.transferFcn='poslin';
        %         net.layers{2}.transferFcn='tansig';
        [net,tr]=train(net,Results(i).Classification_tt,Results(i).Coordinates_peak_belief_tt);
        
        %test
        Results(i).Test_net = net(Results(i).Classification_tt);
        Results(i).Decoded_angle_net_tt = angle(Results(i).Test_net(1,:) + 1i.*Results(i).Test_net(2,:));
        Results(i).e_net = circular_distance(Results(i).True_peak_belief_tt,Results(i).Decoded_angle_net_tt);
        Results(i).performance=perform(net,Results(i).Coordinates_peak_belief_tt,Results(i).Test_net);
        Results(i).net=net;
        Results(i).tr=tr;
        %         view(net)
        
        Results(i).Decoded_angle_tt=mod(angle(mean(Results(i).Classification_tt.*exp(1i*angle_mean_belief_bin'))),2*pi);
        Results(i).e_angle = circular_distance(Results(i).True_peak_belief_tt, Results(i).Decoded_angle_tt);
        
        %Now decode the val data
        Results(i).True_peak_belief=peak_belief(data_val(:,i));
        Results(i).Coordinates_peak_belief=[cos(peak_belief(data_val(:,i))); sin(peak_belief(data_val(:,i)))];
        for c=1:N_channels
            [~, probas,~] = predict(SVMModel(c).mdl,Y(:,data_val(:,i))');
            Results(i).Classification(c,:)=probas(:,2);
            clear probas
        end
        
        Results(i).Val_net=net(Results(i).Classification);
        Results(i).Decoded_angle_net=angle(Results(i).Val_net(1,:) + 1i.*Results(i).Val_net(2,:));
        
        Results(i).Decoded_angle=angle(mean(Results(i).Classification.*exp(1i*angle_mean_belief_bin')));
        
        Results(i).e_net_val = circular_distance(Results(i).True_peak_belief,Results(i).Decoded_angle_net);
        Results(i).e_angle_val = circular_distance(Results(i).True_peak_belief,Results(i).Decoded_angle);
%         figure
%         subplot(1,2,1)
%         plot(mod(Decoded_angle(:,i),2*pi))
%         hold on
%         plot(mod(True_peak_belief(:,i),2*pi))
%         title('Decoded angle')
%         subplot(1,2,2)
%         plot(Decoded_norm(:,i))
    end
    
    for i=1:length(index_start)
         Results(i).prev_rpe_val=prev_rpe(data_val(:,i)');
         Results(i).prev_chosen_color_val=prev_chosen_color(data_val(:,i)');
    end
    
    decoded_angle_net_reshaped=Results(1).Decoded_angle_net;
    decoded_angle_reshaped=Results(1).Decoded_angle;
    prev_rpe_val_reshaped=Results(1).prev_rpe_val;
    prev_chosen_color_val_reshaped=Results(1).prev_chosen_color_val;
    true_peak_reshaped=Results(1).True_peak_belief;
    e_net_val_reshaped=Results(1).e_net_val;
    e_angle_val_reshaped=Results(1).e_angle_val;
    block_val_reshaped=1.*ones(1,length(Results(1).Decoded_angle_net));
    
    for i=2:length(index_start)
        decoded_angle_net_reshaped=cat(2,decoded_angle_net_reshaped,Results(i).Decoded_angle_net);
        decoded_angle_reshaped=cat(2,decoded_angle_reshaped,Results(i).Decoded_angle);
        prev_rpe_val_reshaped=cat(2,prev_rpe_val_reshaped,Results(i).prev_rpe_val);
        prev_chosen_color_val_reshaped=cat(2,prev_chosen_color_val_reshaped,Results(i).prev_chosen_color_val);
        true_peak_reshaped=cat(2,true_peak_reshaped,Results(i).True_peak_belief);
        e_net_val_reshaped=cat(2,e_net_val_reshaped,Results(i).e_net_val);
        e_angle_val_reshaped=cat(2,e_angle_val_reshaped,Results(i).e_angle_val);
        block_val_reshaped=cat(2,block_val_reshaped,i.*ones(1,length(Results(1).Decoded_angle_net)));
    end
    
    
%     f=figure;
%     histogram(dist_decoded_angle_chosen_color_rpe_positive)
%     hold on
%     histogram(dist_decoded_angle_chosen_color_rpe_negative)
%     legend({'RPE>0','RPE<0'})
%     xlabel('dist decoded angle to prev chosen color')
%     title(sprintf('belief decoder update %s %d rpe',ROI,arrayID))
% %     txt=sprintf('T=%d \n p=%d',t.tstat, p);
% %     text(0,100,txt)
%     saveas(f,sprintf('belief_decoder_update_%s_%d_rpe',ROI,arrayID),'jpg')

save(fullfile(data_path,savename))

% end


%%
% figure
% subplot(1,2,1)
% hold on
% scatter(mod(true_peak_reshaped,2*pi),mod(decoded_angle_net_reshaped,2*pi))
% scatter(mod(true_peak_reshaped,2*pi),mod(decoded_angle_reshaped,2*pi))
% title(sprintf('%d',arrayID'))
% 
% subplot(1,2,2)
% hold on
% histogram(e_net_val_reshaped)
% histogram(e_angle_val_reshaped)
% xline(pi/2,'LineWidth',2)
% xline(nanmean(e_net_val_reshaped),'b','LineWidth',2)
% xline(nanmean(e_angle_val_reshaped),'r','LineWidth',2)
% title(sprintf('%d',arrayID'))
end


    