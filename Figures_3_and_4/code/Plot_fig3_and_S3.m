clear all

fsroot='/Volumes/buschman';
arrayID_list=[3 8 12 14 18 20 22 26 28 30 32 34 36 38 40 42 44];  % 34  38 %% ADD 38

window_start=-600;
event='target';
task='Learning_Attentional_Templates';
subtask='exploreexploit/Restricted_FEF/Reset_RW_model';

window_size=900;
N_perm=1000;

sig_thr=1.1;

N_bins=100;
color_binned = 0:(2*pi/N_bins):(2*pi-(2*pi/N_bins));

addpath(genpath('CircStat2012a'))
addpath(genpath('Violinplot-Matlab-master'))

%% Functions

circular_distance = @(x, y)mod(x-y+pi,2*pi)-pi;
update_index_distance = @(p_rpe,p_cc,p_t,t)sign(p_rpe).*sign(mod(p_t-p_cc+pi,2*pi)-pi).*(mod(p_t-t+pi,2*pi)-pi);
reshape_across_blocks = @(x)reshape(x,size(x,1)*size(x,2),1);
nanzscore = @(x)(x-nanmean(x))./nanstd(x);

%% Initiate

Decoded_angle=NaN(35,200);
True_peak=NaN(35,200);
Prev_RPE=NaN(35,200);
Prev_chosen_color=NaN(35,200);
p_val_block=NaN(1,200);
p_dist_val_this_block=NaN(1,200);

N_neurons_block=NaN(1,200);
N_sess_block=NaN(1,200);
N_sess_all_block=NaN(35,200);

All_angular_change=NaN(35,200);
All_true_angular_change=NaN(35,200);

All_update=NaN(35,200);
All_true_update=NaN(35,200);

All_dist_prev_decoded_angle_cc=NaN(35,200);
All_dist_prev_true_peak_cc=NaN(35,200);

All_dist_decoded_angle_cc=NaN(35,200);
All_dist_true_peak_cc=NaN(35,200);

All_error=NaN(35,200);
All_prev_error=NaN(35,200);

Precision=NaN(35,200);

count=0;

%% Load and organize

subsubtask=sprintf('FR_%s_%d_window_from_%d_to_%d',event,window_size,window_start,window_start+window_size);
dirstem = fullfile('Projects',task,'Analysis','Electrophy_analysis',subtask,subsubtask);
data_path = fullfile(fsroot,dirstem);

for j=1:length(arrayID_list)
    
    arrayID=arrayID_list(j);
    
    loadname=sprintf('belief_peak_decoder_update_NN_relu_%d.mat',arrayID);
    
    if exist(fullfile(data_path,loadname))==2
        load(fullfile(data_path,loadname),'decoded_angle_net_reshaped','prev_rpe_val_reshaped','prev_chosen_color_val_reshaped','e_net_val_reshaped','e_angle_val_reshaped','true_peak_reshaped','block_val_reshaped','W','Conditions','data_val','Results')
        
        %extract the precision and entropy
        for i=1:size(data_val,2)
            belief_this_val=Conditions.belief(:,data_val(:,i));
            rt(:,i)=Conditions.RT(data_val(:,i));
            reset(i)=(sum(Conditions.switch(data_val(:,i)==1))>0);
            for k=1:size(belief_this_val,2)
                if sum(belief_this_val(:,k)==0)<size(belief_this_val,1)
                    proba_belief_this_val(:,k)=(belief_this_val(:,k)-min(belief_this_val(:,k))+0.01)./sum(belief_this_val(:,k)-min(belief_this_val(:,k))+0.01);
                else
                    proba_belief_this_val(1:size(belief_this_val,1),k)=NaN;
                end
            end
            for k=1:size(belief_this_val,2)
                precision(k,i)=norm(mean(belief_this_val(:,k).*exp(1i*color_binned)));
                entropy(k,i)=-sum(proba_belief_this_val(:,k).*log(proba_belief_this_val(:,k)));
            end
            clear proba_this_belief
        end
        a(j)=nanmean(e_angle_val_reshaped-e_net_val_reshaped);
        [z_dist_val(j), p_dist_val(j)]=compute_shuffled_circ_dist(decoded_angle_net_reshaped,true_peak_reshaped,N_perm);
        [~, pp_ttest_val(j)]=ttest(e_net_val_reshaped,pi/2,'Tail','left');
        this_dist_val(j)=nanmean(abs(mod(decoded_angle_net_reshaped-true_peak_reshaped+pi,2*pi)-pi));
        N_neurons(j)=size(W,1);
        N_template(j)=max(block_val_reshaped);
        
        for n_blocks=1:max(block_val_reshaped)
            [~, p_dist_val_this_block(n_blocks+count)]=ttest(e_net_val_reshaped(block_val_reshaped==n_blocks),pi/2,'Tail','left');
            p_val_block(n_blocks+count)=p_dist_val(j);
            N_neurons_block(n_blocks+count)=N_neurons(j);
            N_sess_block(n_blocks+count)=j;
            Decoded_angle(:,n_blocks+count)=decoded_angle_net_reshaped(block_val_reshaped==n_blocks);
            True_peak(:,n_blocks+count)=true_peak_reshaped(block_val_reshaped==n_blocks);
            Prev_RPE(:,n_blocks+count)=prev_rpe_val_reshaped(block_val_reshaped==n_blocks);
            Prev_chosen_color(:,n_blocks+count)=prev_chosen_color_val_reshaped(block_val_reshaped==n_blocks);
            Precision(1,n_blocks+count)=NaN;
            Prev_entropy(1,n_blocks+count)=NaN;
            Entropy(1,n_blocks+count)=entropy(1,n_blocks);
            RT(1,n_blocks+count)=rt(1,n_blocks);
            Reset(n_blocks+count)=reset(n_blocks);
            if arrayID<=8
                Monkey(1,n_blocks+count)=1;
            else
                Monkey(1,n_blocks+count)=2;
            end

            Decoded_angle(isnan(True_peak(:,n_blocks+count)),n_blocks+count)=NaN; %remove when no peak (flat EV)
            
            All_error(1,n_blocks+count)=abs(circular_distance(True_peak(1,n_blocks+count),Decoded_angle(1,n_blocks+count)));
            
            for i=2:35
                All_angular_change(i,n_blocks+count)=circular_distance(Decoded_angle(i-1,n_blocks+count),Decoded_angle(i,n_blocks+count));
                All_true_angular_change(i,n_blocks+count)=circular_distance(True_peak(i-1,n_blocks+count),True_peak(i,n_blocks+count));
                
                All_update(i,n_blocks+count)=update_index_distance(Prev_RPE(i,n_blocks+count),Prev_chosen_color(i,n_blocks+count),Decoded_angle(i-1,n_blocks+count),Decoded_angle(i,n_blocks+count));
                All_true_update(i,n_blocks+count)=update_index_distance(Prev_RPE(i,n_blocks+count),Prev_chosen_color(i,n_blocks+count),True_peak(i-1,n_blocks+count),True_peak(i,n_blocks+count));
                
                All_dist_prev_decoded_angle_cc(i,n_blocks+count)=circular_distance(Decoded_angle(i-1,n_blocks+count),Prev_chosen_color(i,n_blocks+count));
                All_dist_prev_true_peak_cc(i,n_blocks+count)=circular_distance(True_peak(i-1,n_blocks+count),Prev_chosen_color(i,n_blocks+count));
                
                All_dist_decoded_angle_cc(i,n_blocks+count)=circular_distance(Decoded_angle(i,n_blocks+count),Prev_chosen_color(i,n_blocks+count));
                All_dist_true_peak_cc(i,n_blocks+count)=circular_distance(True_peak(i,n_blocks+count),Prev_chosen_color(i,n_blocks+count));
                
                All_error(i,n_blocks+count)=abs(circular_distance(True_peak(i,n_blocks+count),Decoded_angle(i,n_blocks+count)));
                All_prev_error(i,n_blocks+count)=abs(circular_distance(True_peak(i-1,n_blocks+count),Decoded_angle(i-1,n_blocks+count)));
                
                Precision(i,n_blocks+count)=precision(i-1,n_blocks);
                Prev_entropy(i,n_blocks+count)=entropy(i-1,n_blocks);
                Entropy(i,n_blocks+count)=entropy(i,n_blocks);
                RT(i,n_blocks+count)=rt(i,n_blocks);
                if arrayID<=8
                    Monkey(i,n_blocks+count)=1;
                else
                    Monkey(i,n_blocks+count)=2;
                end

                
                N_sess_all_block(i,n_blocks+count)=j;
            end
        end
        clear precision entropy rt
        count=count+max(block_val_reshaped);
    end
end

%% Reshape

All_angular_change_r=reshape_across_blocks(All_angular_change(:,p_val_block<sig_thr));
All_true_angular_change_r=reshape_across_blocks(All_true_angular_change(:,p_val_block<sig_thr));

All_update_r=reshape_across_blocks(All_update(:,p_val_block<sig_thr));
All_true_update_r=reshape_across_blocks(All_true_update(:,p_val_block<sig_thr));

All_prev_rpe_r=reshape_across_blocks(Prev_RPE(:,p_val_block<sig_thr));
All_prev_chosen_color_r=reshape_across_blocks(Prev_chosen_color(:,p_val_block<sig_thr));

All_dist_prev_decoded_angle_cc_r=reshape_across_blocks(All_dist_prev_decoded_angle_cc(:,p_val_block<sig_thr));
All_dist_prev_true_peak_cc_r=reshape_across_blocks(All_dist_prev_true_peak_cc(:,p_val_block<sig_thr));

All_dist_decoded_angle_cc_r=reshape_across_blocks(All_dist_decoded_angle_cc(:,p_val_block<sig_thr));
All_dist_true_peak_cc_r=reshape_across_blocks(All_dist_true_peak_cc(:,p_val_block<sig_thr));

All_error_r=reshape_across_blocks(All_error(:,p_val_block<sig_thr));
All_prev_error_r=reshape_across_blocks(All_prev_error(:,p_val_block<sig_thr));

N_sess_all_block_r=reshape_across_blocks(N_sess_all_block(:,p_val_block<sig_thr));

All_precision_r=reshape_across_blocks(Precision(:,p_val_block<sig_thr));
All_prev_entropy_r=reshape_across_blocks(Prev_entropy(:,p_val_block<sig_thr));
All_entropy_r=reshape_across_blocks(Entropy(:,p_val_block<sig_thr));
All_rt_r=reshape_across_blocks(RT(:,p_val_block<sig_thr));
All_monkey_r=reshape_across_blocks(Monkey(:,p_val_block<sig_thr));
All_trial_r=reshape_across_blocks(repmat([1:35]',1,size(Monkey,2)));


All_decoded_angle_r=reshape_across_blocks(Decoded_angle(:,p_val_block<sig_thr));
All_true_peak_r=reshape_across_blocks(True_peak(:,p_val_block<sig_thr));


%% fig S3B

figure
subplot(1,2,1)
histogram(All_entropy_r,'FaceColor',[0.5 0.5 0.5])
box off
xlabel('Entropy')
ylabel('# trials')

subplot(1,2,2)
hold on
plot(Entropy(:,Reset==0),'-','Color',[0.5 0.5 0.5])
plot(Entropy(:,Reset>0),'k-')
plot(nanmean(Entropy(:,Reset==0),2),'-','Color',[0.5 0.5 0.5],'LineWidth',2)
plot(nanmean(Entropy(:,Reset>0),2),'k-','LineWidth',2)

xlim([0 35])
ylabel('Entropy')
xlabel('Trials from tmpplate switch')
box off

%% Stats over all trials



disp('%%% Display stats all trials %%%')

%3A
disp('decoding')
[z_dist_decoding, p_dist_decoding] = compute_shuffled_circ_dist(All_decoded_angle_r,All_true_peak_r,N_perm);

%%
[rho_ind, pval_ind] = circ_corrcc(All_decoded_angle_r(~isnan(All_decoded_angle_r)),All_true_peak_r(~isnan(All_decoded_angle_r)));

%3B
disp('improved decoding entropy')
[r_pearson_decoding_entropy_corr, p_pearson_decoding_entropy_corr]=corr(abs(All_error_r(~isnan(All_error_r) & ~isnan(All_entropy_r))), All_entropy_r(~isnan(All_error_r) & ~isnan(All_entropy_r)),'Type','Pearson','Tail','right')

%3D
dist_update_metrics=abs(circular_distance(All_update_r,All_angular_change_r));
disp('update and update in correct direction')
[z_dist_update, p_dist_update] = compute_shuffled_circ_dist(All_update_r,All_angular_change_r,N_perm)

dist_true_update_metrics_pos=abs(circular_distance(All_true_update_r(All_prev_rpe_r>0),All_update_r(All_prev_rpe_r>0)));
disp('update rpe>0')
[z_dist_true_update_pos, p_dist_true_update_pos] = compute_shuffled_circ_dist(All_true_update_r(All_prev_rpe_r>0),All_update_r(All_prev_rpe_r>0),N_perm)

dist_true_update_metrics_neg=abs(circular_distance(All_true_update_r(All_prev_rpe_r<0),All_update_r(All_prev_rpe_r<0)));
disp('update rpe<0')
[z_dist_true_update_neg, p_dist_true_update_neg] = compute_shuffled_circ_dist(All_true_update_r(All_prev_rpe_r<0),All_update_r(All_prev_rpe_r<0),N_perm)

%3E
disp('rpe corr rpe>0')
[r_pearson_corr_pos, p_pearson_corr_pos]=corr(abs(All_angular_change_r(~isnan(All_angular_change_r) & All_prev_rpe_r>0)), abs(All_prev_rpe_r(~isnan(All_update_r) & All_prev_rpe_r>0)),'Type','Pearson','Tail','right')

disp('rpe corr rpe<0')
[r_pearson_corr_neg, p_pearson_corr_neg]=corr(abs(All_angular_change_r(~isnan(All_angular_change_r) & All_prev_rpe_r<0)), abs(All_prev_rpe_r(~isnan(All_update_r) & All_prev_rpe_r<0)),'Type','Pearson','Tail','right')

disp('rpe model corr rpe>0')
[r_pearson_model_corr_pos, p_pearson_model_corr_pos]=corr(abs(All_true_angular_change_r(~isnan(All_angular_change_r) & All_prev_rpe_r>0)), abs(All_prev_rpe_r(~isnan(All_update_r) & All_prev_rpe_r>0)),'Type','Pearson','Tail','right')

disp('rpe model corr rpe<0')
[r_pearson_model_corr_neg, p_pearson_model_corr_neg]=corr(abs(All_true_angular_change_r(~isnan(All_angular_change_r) & All_prev_rpe_r<0)), abs(All_prev_rpe_r(~isnan(All_update_r) & All_prev_rpe_r<0)),'Type','Pearson','Tail','right')

disp('cc corr rpe>0')
[r_pearson_cc_corr_pos, p_pearson_cc_corr_pos]=corr(abs(All_angular_change_r(~isnan(All_angular_change_r) & All_prev_rpe_r>0)), abs(All_dist_prev_decoded_angle_cc_r(~isnan(All_update_r) & All_prev_rpe_r>0)),'Type','Pearson','Tail','right')

disp('cc corr rpe<0')
[r_pearson_cc_corr_neg, p_pearson_cc_corr_neg]=corr(abs(All_angular_change_r(~isnan(All_angular_change_r) & All_prev_rpe_r<0)), abs(All_dist_prev_decoded_angle_cc_r(~isnan(All_update_r) & All_prev_rpe_r<0)),'Type','Pearson','Tail','left')

disp('cc model corr rpe>0')
[r_pearson_cc_update_corr_pos, p_pearson_update_corr_pos]=corr(abs(All_true_angular_change_r(~isnan(All_angular_change_r) & All_prev_rpe_r>0)), abs(All_dist_prev_true_peak_cc_r(~isnan(All_update_r) & All_prev_rpe_r>0)),'Type','Pearson','Tail','right')

disp('cc model corr rpe<0')
[r_pearson_cc_model_corr_neg, p_pearson_cc_model_corr_neg]=corr(abs(All_true_angular_change_r(~isnan(All_angular_change_r) & All_prev_rpe_r<0)), abs(All_dist_prev_true_peak_cc_r(~isnan(All_update_r) & All_prev_rpe_r<0)),'Type','Pearson','Tail','left')



%% 3D

disp('RPE>0')

[h p s t]=ttest(All_true_update_r(All_prev_rpe_r>0),0,'Tail','right')
[h p s t]=ttest(All_update_r(All_prev_rpe_r>0),0,'Tail','right')

disp('RPE<0')

[h p s t]=ttest(All_true_update_r(All_prev_rpe_r<0),0,'Tail','right')
[h p s t]=ttest(All_update_r(All_prev_rpe_r<0),0,'Tail','right')


%% Fig 3D

All_true_update_pos=All_true_update_r(All_prev_rpe_r>0 & All_true_update_r>0);
All_update_pos=All_update_r(All_prev_rpe_r>0 & All_true_update_r>0);

figure
hold on
bar(1,angle(nanmean(exp(1i*All_update_pos(All_true_update_pos<0.1)))),'Facecolor',[0.5 0.5 0.5])
bar(2,angle(nanmean(exp(1i*All_update_pos(All_true_update_pos>=0.1 & All_true_update_pos<0.2)))),'Facecolor',[0.5 0.5 0.5])
bar(3,angle(nanmean(exp(1i*All_update_pos(All_true_update_pos>=0.2)))),'Facecolor',[0.5 0.5 0.5])
errorbar(1,angle(nanmean(exp(1i*All_update_pos(All_true_update_pos<=0.1)))),circ_std(All_update_pos(All_true_update_pos<0.1))/sqrt(sum(All_true_update_pos<=0.1)-1),'k','LineWidth',1)
errorbar(2,angle(nanmean(exp(1i*All_update_pos(All_true_update_pos>=0.1 & All_true_update_pos<=0.2)))),circ_std(All_update_pos(All_true_update_pos>=0.1 & All_true_update_pos<=0.2))/sqrt(sum(All_true_update_pos>=0.1 & All_true_update_pos<=0.2)-1),'k','LineWidth',1)
errorbar(3,angle(nanmean(exp(1i*All_update_pos(All_true_update_pos>=0.2)))),circ_std(All_update_pos(All_true_update_pos>=0.2))/sqrt(sum(All_true_update_pos>=0.2)-1),'k','LineWidth',1)
box off
ylabel('Mean behavioral update toward chosen color')
xlabel('Neural update toward chosen color')
xticks([1 2 3])
xticklabels({'≤0.1','[0.1 0.2]','≥0.2'})

%% Session average

for j=1:length(arrayID_list)
    if p_dist_val(j)<sig_thr
        mean_update_sess(j)=nanmean((All_update_r(N_sess_all_block_r==j)));
        mean_update_pos_sess(j)=nanmean((All_update_r(N_sess_all_block_r==j & All_prev_rpe_r>0)));
        mean_update_neg_sess(j)=nanmean((All_update_r(N_sess_all_block_r==j & All_prev_rpe_r<0)));

        circ_mean_update_sess(j)=angle(nanmean(exp(1i*All_update_r(N_sess_all_block_r==j))));
        circ_mean_update_pos_sess(j)=angle(nanmean(exp(1i*All_update_r(N_sess_all_block_r==j & All_prev_rpe_r>0))));
        circ_mean_update_neg_sess(j)=angle(nanmean(exp(1i*All_update_r(N_sess_all_block_r==j & All_prev_rpe_r<0))));
        
        mean_update_model_sess(j)=nanmean((All_true_update_r(N_sess_all_block_r==j)));
        mean_update_model_pos_sess(j)=nanmean((All_true_update_r(N_sess_all_block_r==j & All_prev_rpe_r>0)));
        mean_update_model_neg_sess(j)=nanmean((All_true_update_r(N_sess_all_block_r==j & All_prev_rpe_r<0)));
        
        circ_mean_update_model_sess(j)=angle(nanmean(exp(1i*All_true_update_r(N_sess_all_block_r==j))));
        circ_mean_update_model_pos_sess(j)=angle(nanmean(exp(1i*All_true_update_r(N_sess_all_block_r==j & All_prev_rpe_r>0))));
        circ_mean_update_model_neg_sess(j)=angle(nanmean(exp(1i*All_true_update_r(N_sess_all_block_r==j & All_prev_rpe_r<0))));
        
        
        pearson_corr_sess(j)=corr(abs(All_angular_change_r(N_sess_all_block_r==j & ~isnan(All_angular_change_r))), abs(All_prev_rpe_r(N_sess_all_block_r==j & ~isnan(All_update_r))),'Type','Pearson');
        pearson_corr_pos_sess(j)=corr(abs(All_angular_change_r(N_sess_all_block_r==j & All_prev_rpe_r>0 & ~isnan(All_angular_change_r))), abs(All_prev_rpe_r(N_sess_all_block_r==j & All_prev_rpe_r>0 & ~isnan(All_update_r))),'Type','Pearson');
        pearson_corr_neg_sess(j)=corr(abs(All_angular_change_r(N_sess_all_block_r==j & All_prev_rpe_r<0 & ~isnan(All_angular_change_r))), abs(All_prev_rpe_r(N_sess_all_block_r==j & All_prev_rpe_r<0 & ~isnan(All_update_r))),'Type','Pearson');
        
        pearson_corr_model_sess(j)=corr(abs(All_true_angular_change_r(N_sess_all_block_r==j & ~isnan(All_angular_change_r))), abs(All_prev_rpe_r(N_sess_all_block_r==j & ~isnan(All_update_r))),'Type','Pearson');
        pearson_corr_model_pos_sess(j)=corr(abs(All_true_angular_change_r(N_sess_all_block_r==j & All_prev_rpe_r>0 & ~isnan(All_angular_change_r))), abs(All_prev_rpe_r(N_sess_all_block_r==j & All_prev_rpe_r>0 & ~isnan(All_update_r))),'Type','Pearson');
        pearson_corr_model_neg_sess(j)=corr(abs(All_true_angular_change_r(N_sess_all_block_r==j & All_prev_rpe_r<0 & ~isnan(All_angular_change_r))), abs(All_prev_rpe_r(N_sess_all_block_r==j & All_prev_rpe_r<0 & ~isnan(All_update_r))),'Type','Pearson');
        
        [z_dist_update_sess(j), p_dist_update_sess(j)] = compute_shuffled_circ_dist(All_true_update_r(N_sess_all_block_r==j),All_update_r(N_sess_all_block_r==j),N_perm);
        [z_dist_update_pos_sess(j), p_dist_update_pos_sess(j)] = compute_shuffled_circ_dist(All_true_update_r(N_sess_all_block_r==j & All_prev_rpe_r>0),All_update_r(N_sess_all_block_r==j & All_prev_rpe_r>0),N_perm);
        [z_dist_update_neg_sess(j), p_dist_update_neg_sess(j)] = compute_shuffled_circ_dist(All_true_update_r(N_sess_all_block_r==j & All_prev_rpe_r<0),All_update_r(N_sess_all_block_r==j & All_prev_rpe_r<0),N_perm);
        
        b_error_prec_sess(j)=corr(All_error_r(N_sess_all_block_r==j & ~isnan(All_error_r)),All_precision_r(N_sess_all_block_r==j & ~isnan(All_error_r)),'Type','Pearson');
        b_error_entropy_sess(j)=corr(All_error_r(N_sess_all_block_r==j & ~isnan(All_error_r) & ~isnan(All_entropy_r)),All_entropy_r(N_sess_all_block_r==j & ~isnan(All_error_r) & ~isnan(All_entropy_r)),'Type','Pearson');
        
        [z_dist_angular_change_sess(j), p_dist_angular_change_sess(j)] = compute_shuffled_circ_dist(All_true_angular_change_r(N_sess_all_block_r==j),All_angular_change_r(N_sess_all_block_r==j),N_perm);
        [z_dist_angular_change_pos_sess(j), p_dist_angular_change_pos_sess(j)] = compute_shuffled_circ_dist(All_true_angular_change_r(N_sess_all_block_r==j & All_prev_rpe_r>0),All_angular_change_r(N_sess_all_block_r==j & All_prev_rpe_r>0),N_perm);
        [z_dist_angular_change_neg_sess(j), p_dist_angular_change_neg_sess(j)] = compute_shuffled_circ_dist(All_true_angular_change_r(N_sess_all_block_r==j & All_prev_rpe_r<0),All_angular_change_r(N_sess_all_block_r==j & All_prev_rpe_r<0),N_perm);
        
        [z_dist_angular_change_prev_to_cc_sess(j), p_dist_angular_change_prev_to_cc_sess(j)] = compute_shuffled_circ_dist(sign(All_prev_rpe_r(N_sess_all_block_r==j)).*All_dist_prev_decoded_angle_cc_r(N_sess_all_block_r==j),All_angular_change_r(N_sess_all_block_r==j),N_perm);
        [z_dist_angular_change_prev_to_cc_pos_sess(j), p_dist_angular_change_prev_to_cc_pos_sess(j)] = compute_shuffled_circ_dist(All_dist_prev_decoded_angle_cc_r(N_sess_all_block_r==j & All_prev_rpe_r>0),All_angular_change_r(N_sess_all_block_r==j & All_prev_rpe_r>0),N_perm);
        [z_dist_angular_change_prev_to_cc_neg_sess(j), p_dist_angular_change_prev_to_cc_neg_sess(j)] = compute_shuffled_circ_dist(All_dist_prev_decoded_angle_cc_r(N_sess_all_block_r==j & All_prev_rpe_r<0),All_angular_change_r(N_sess_all_block_r==j & All_prev_rpe_r<0),N_perm);
        
        [z_dist_true_angular_change_prev_to_cc_sess(j), p_dist_true_angular_change_prev_to_cc_sess(j)] = compute_shuffled_circ_dist(All_dist_prev_true_peak_cc_r(N_sess_all_block_r==j),All_true_angular_change_r(N_sess_all_block_r==j),N_perm);
        [z_dist_true_angular_change_prev_to_cc_pos_sess(j), p_dist_true_angular_change_prev_to_cc_pos_sess(j)] = compute_shuffled_circ_dist(All_dist_prev_true_peak_cc_r(N_sess_all_block_r==j & All_prev_rpe_r>0),All_true_angular_change_r(N_sess_all_block_r==j & All_prev_rpe_r>0),N_perm);
        [z_dist_true_angular_change_prev_to_cc_neg_sess(j), p_dist_true_angular_change_prev_to_cc_neg_sess(j)] = compute_shuffled_circ_dist(All_dist_prev_true_peak_cc_r(N_sess_all_block_r==j & All_prev_rpe_r<0),All_true_angular_change_r(N_sess_all_block_r==j & All_prev_rpe_r<0),N_perm);
        
        [r_dist_abs_angular_change_prev_to_cc_sess(j), p_dist_abs_angular_change_prev_to_cc_sess(j)] = corr(abs(All_dist_prev_decoded_angle_cc_r(N_sess_all_block_r==j & ~isnan(All_angular_change_r))),abs(All_angular_change_r(N_sess_all_block_r==j & ~isnan(All_angular_change_r))));
        [r_dist_abs_angular_change_prev_to_cc_pos_sess(j), p_dist_abs_angular_change_prev_to_cc_pos_sess(j)] = corr(abs(All_dist_prev_decoded_angle_cc_r(N_sess_all_block_r==j & All_prev_rpe_r>0 & ~isnan(All_angular_change_r))),abs(All_angular_change_r(N_sess_all_block_r==j & ~isnan(All_angular_change_r) & All_prev_rpe_r>0)));
        [r_dist_abs_angular_change_prev_to_cc_neg_sess(j), p_dist_abs_angular_change_prev_to_cc_neg_sess(j)] = corr(abs(All_dist_prev_decoded_angle_cc_r(N_sess_all_block_r==j & All_prev_rpe_r<0 & ~isnan(All_angular_change_r))),abs(All_angular_change_r(N_sess_all_block_r==j & ~isnan(All_angular_change_r) & All_prev_rpe_r<0)));
        
        [r_dist_abs_true_angular_change_prev_to_cc_sess(j), p_dist_abs_true_angular_change_prev_to_cc_sess(j)] = corr(abs(All_dist_prev_true_peak_cc_r(N_sess_all_block_r==j & ~isnan(All_true_angular_change_r))),abs(All_true_angular_change_r(N_sess_all_block_r==j & ~isnan(All_true_angular_change_r))));
        [r_dist_abs_true_angular_change_prev_to_cc_pos_sess(j), p_dist_abs_true_angular_change_prev_to_cc_pos_sess(j)] = corr(abs(All_dist_prev_true_peak_cc_r(N_sess_all_block_r==j & All_prev_rpe_r>0 & ~isnan(All_true_angular_change_r))),abs(All_true_angular_change_r(N_sess_all_block_r==j & ~isnan(All_true_angular_change_r) & All_prev_rpe_r>0)));
        [r_dist_abs_true_angular_change_prev_to_cc_neg_sess(j), p_dist_abs_true_angular_change_prev_to_cc_neg_sess(j)] = corr(abs(All_dist_prev_true_peak_cc_r(N_sess_all_block_r==j & All_prev_rpe_r<0 & ~isnan(All_true_angular_change_r))),abs(All_true_angular_change_r(N_sess_all_block_r==j & ~isnan(All_true_angular_change_r) & All_prev_rpe_r<0)));
    end
end

%% Fig 3A and C

figure
subplot(2,1,1)
vs = violinplot(z_dist_val',{'All'},'ViolinColor',[0 0 0],'ShowData',true,'ShowMean',true,'ShowMedian',true);
if ttest(z_dist_val,0,'Tail','left')
    plot(1,2.5,'k*')
end
yline(0,'--')
box off
ylabel(sprintf('z-score decosing distance'))

subplot(2,1,2)
vs = violinplot(b_error_entropy_sess',{'All'},'ViolinColor',[0 0 0],'ShowData',true,'ShowMean',true,'ShowMedian',true);
if ttest(b_error_entropy_sess,0,'Tail','right')
    plot(1,0.5,'k*')
end

yline(0,'--')
box off
ylabel(sprintf('Correlation with precision of template'))



%% figure 3D

figure

subplot(2,1,1)
vs = violinplot([circ_mean_update_model_pos_sess', circ_mean_update_pos_sess'],{'Behavioral','Neural'},'ViolinColor',[0 0 0],'ShowData',true,'ShowMean',true,'ShowMedian',true);
hold on
if ttest(circ_mean_update_model_pos_sess,0,'Tail','right')
    plot(1,0.5,'k*')
end
if ttest(circ_mean_update_pos_sess,0,'Tail','right')
    plot(2,0.5,'k*')
end
yline(0,'--')
box off
ylabel(sprintf('Update toward the chosen color'))

subplot(2,1,2)
vs = violinplot([z_dist_angular_change_pos_sess'],{'Neural ~ beahvioral'},'ViolinColor',[0 0 0],'ShowData',true,'ShowMean',true,'ShowMedian',true);
if ttest(z_dist_angular_change_pos_sess,0,'Tail','left')
    plot(1,1.5,'k*')
end
yline(0,'--')
box off
ylabel(sprintf('Z-score distance between updates'))

%% Figure 3E

figure

subplot(1,2,1)
vs = violinplot([pearson_corr_model_pos_sess', pearson_corr_pos_sess'],{'Behavioral','Neural'},'ViolinColor',[0 0 0],'ShowData',true,'ShowMean',true,'ShowMedian',true);
if ttest(pearson_corr_model_pos_sess,0,'Tail','right')
    plot(1,0.8,'k*')
end
if ttest(pearson_corr_pos_sess,0,'Tail','right')
    plot(2,0.8,'k*')
end
yline(0,'--')
box off
ylabel(sprintf('Pearson correlation with |RPE|'))


subplot(1,2,2)
vs = violinplot([r_dist_abs_true_angular_change_prev_to_cc_pos_sess', r_dist_abs_angular_change_prev_to_cc_pos_sess'],{'Behavioral','Neural'},'ViolinColor',[0 0 0],'ShowData',true,'ShowMean',true,'ShowMedian',true);
hold on
if ttest(r_dist_abs_true_angular_change_prev_to_cc_pos_sess,0,'Tail','right')
    plot(1,0.8,'k*')
end
if ttest(r_dist_abs_angular_change_prev_to_cc_pos_sess,0,'Tail','right')
    plot(2,0.8,'k*')
end
yline(0,'--')
box off
ylabel(sprintf('Pearson correlation with |d(ET(t-1),CC(t-1))|'))

%% Figure 3D and S3C

figure

subplot(2,2,1)
vs = violinplot([circ_mean_update_model_neg_sess'],{'Behavioral'},'ViolinColor',[0 0 0],'ShowData',true,'ShowMean',true,'ShowMedian',true);
hold on
if ttest(circ_mean_update_model_neg_sess,0,'Tail','right')
    plot(1,0.06,'k*')
end
yline(0,'--')
box off
ylabel(sprintf('Update away from the chosen color'))

subplot(2,2,2)
vs = violinplot([circ_mean_update_neg_sess'],{'Neural'},'ViolinColor',[0 0 0],'ShowData',true,'ShowMean',true,'ShowMedian',true);
if ttest(circ_mean_update_neg_sess,0,'Tail','right')
    plot(2,0.4,'k*')
end
yline(0,'--')
box off
ylabel(sprintf('Update away from  chosen color'))




subplot(2,2,[3 4])
vs = violinplot([z_dist_angular_change_neg_sess'],{'Neural ~ beahvioral'},'ViolinColor',[0 0 0],'ShowData',true,'ShowMean',true,'ShowMedian',true);
if ttest(z_dist_angular_change_neg_sess,0,'Tail','left')
    plot(1,1.5,'k*')
end
yline(0,'--')
box off
ylabel(sprintf('Z-score distance between updates'))

%%

figure

subplot(1,2,1)
vs = violinplot([pearson_corr_model_neg_sess', pearson_corr_neg_sess'],{'Behavioral','Neural'},'ViolinColor',[0 0 0],'ShowData',true,'ShowMean',true,'ShowMedian',true);
if ttest(pearson_corr_model_neg_sess,0,'Tail','right')
    plot(1,0.5,'k*')
end
if ttest(pearson_corr_neg_sess,0,'Tail','right')
    plot(2,0.5,'k*')
end
yline(0,'--')
box off
ylabel(sprintf('Pearson correlation with |RPE|'))


subplot(1,2,2)
vs = violinplot([r_dist_abs_true_angular_change_prev_to_cc_neg_sess', r_dist_abs_angular_change_prev_to_cc_neg_sess'],{'Behavioral','Neural'},'ViolinColor',[0 0 0],'ShowData',true,'ShowMean',true,'ShowMedian',true);
hold on
if ttest(r_dist_abs_true_angular_change_prev_to_cc_neg_sess,0,'Tail','left')
    plot(1,0.4,'k*')
end
if ttest(r_dist_abs_angular_change_prev_to_cc_neg_sess,0,'Tail','left')
    plot(2,0.4,'k*')
end
yline(0,'--')
box off
ylabel(sprintf('Pearson correlation with |d(ET(t-1),CC(t-1))|'))


%% Figure 3C


All_prc=prctile(All_entropy_r,0:10:100);

for i=1:length(All_prc)-4
    mean_All_error(i)=nanmean(All_error_r(All_entropy_r>=All_prc(i) & All_entropy_r<=All_prc(i+4)));
    sem_All_error(i)=nanstd(All_error_r(All_entropy_r>=All_prc(i) & All_entropy_r<=All_prc(i+4)))/sqrt(sum(All_entropy_r>=All_prc(i) & All_entropy_r<=All_prc(i+4))-1);
    
    
end

figure
subplot(1,2,2)
shadedErrorBar(All_prc(3:end-2),mean_All_error,sem_All_error,{'-','color',[0 0 0],'LineWidth',2},2)
xlabel('entropy')
ylabel(sprintf('Decoding error'))
box off
title('All')


%% Figure 3E

All_prc=prctile(All_prev_rpe_r(All_prev_rpe_r>0),0:10:100);


for i=1:length(All_prc)-4
    mean_All_update(i)=nanmean(abs(All_angular_change_r(All_prev_rpe_r>=All_prc(i) & All_prev_rpe_r<=All_prc(i+4) & All_prev_rpe_r>0)));
    sem_All_update(i)=nanstd(abs(All_angular_change_r(All_prev_rpe_r>=All_prc(i) & All_prev_rpe_r<=All_prc(i+4) & All_prev_rpe_r>0)))/sqrt(sum(All_prev_rpe_r>=All_prc(i) & All_prev_rpe_r<=All_prc(i+4) & All_prev_rpe_r>0)-1);
    
    
    mean_All_true_update(i)=nanmean(abs(All_true_angular_change_r(All_prev_rpe_r>=All_prc(i) & All_prev_rpe_r<=All_prc(i+4) & All_prev_rpe_r>0)));
    sem_All_true_update(i)=nanstd(abs(All_true_angular_change_r(All_prev_rpe_r>=All_prc(i) & All_prev_rpe_r<=All_prc(i+4) & All_prev_rpe_r>0)))/sqrt(sum(All_prev_rpe_r>=All_prc(i) & All_prev_rpe_r<=All_prc(i+4) & All_prev_rpe_r>0)-1);
    
end


figure
subplot(2,2,1)
shadedErrorBar(All_prc(3:end-2),mean_All_update,sem_All_update,{'-','color',[0 0 0],'LineWidth',2},2)
xlabel('RPE')
ylabel(sprintf('Neural update (RPE > 0)'))
box off
title('All')

subplot(2,2,3)
shadedErrorBar(All_prc(3:end-2),mean_All_true_update,sem_All_true_update,{'-','color',[0 0 0],'LineWidth',2},2)
xlabel('RPE')
ylabel(sprintf('Behavioral update (RPE > 0)'))
box off
title('All')

All_angle=prctile(abs(All_dist_prev_decoded_angle_cc_r(All_prev_rpe_r>0)),0:10:100);
All_angle2=prctile(abs(All_dist_prev_true_peak_cc_r(All_prev_rpe_r>0)),0:10:100);

for i=1:length(All_angle)-4
    mean_All_update_angle(i)=nanmean(abs(All_angular_change_r(abs(All_dist_prev_decoded_angle_cc_r)>=All_angle(i) & abs(All_dist_prev_decoded_angle_cc_r)<=All_angle(i+4) & All_prev_rpe_r>0)));
    sem_All_update_angle(i)=nanstd(abs(All_angular_change_r(abs(All_dist_prev_decoded_angle_cc_r)>=All_angle(i) & abs(All_dist_prev_decoded_angle_cc_r)<=All_angle(i+4) & All_prev_rpe_r>0)))/sqrt(sum(abs(All_dist_prev_decoded_angle_cc_r)>=All_angle(i) & abs(All_dist_prev_decoded_angle_cc_r)<=All_angle(i+4) & All_prev_rpe_r>0)-1);
    
    mean_All_true_update_angle(i)=nanmean(abs(All_true_angular_change_r(abs(All_dist_prev_true_peak_cc_r)>=All_angle2(i) & abs(All_dist_prev_true_peak_cc_r)<=All_angle2(i+4) & All_prev_rpe_r>0)));
    sem_All_true_update_angle(i)=nanstd(abs(All_true_angular_change_r(abs(All_dist_prev_true_peak_cc_r)>=All_angle2(i) & abs(All_dist_prev_true_peak_cc_r)<=All_angle2(i+4) & All_prev_rpe_r>0)))/sqrt(sum(abs(All_dist_prev_true_peak_cc_r)>=All_angle2(i) & abs(All_dist_prev_true_peak_cc_r)<=All_angle2(i+4) & All_prev_rpe_r>0)-1);
    
end

subplot(2,2,2)
shadedErrorBar(All_angle(3:end-2),mean_All_update_angle,sem_All_update_angle,{'-','color',[0 0 0],'LineWidth',2},2)
xlabel(sprintf('Neural |d(ET(t-1), CC(t-1))|'))
ylabel(sprintf('Neural update magnitude'))
box off
xticks([0 pi/4 pi/2 3*pi/4])
xticklabels({'0','π/4','π/2','3π/4'})
xlim([0 pi-pi/8])

subplot(2,2,4)
shadedErrorBar(All_angle(3:end-2),mean_All_true_update_angle,sem_All_true_update_angle,{'-','color',[0 0 0],'LineWidth',2},2)
xlabel(sprintf('Model |d(ET(t-1), CC(t-1))|'))
ylabel(sprintf('Model update magnitude'))
box off
xticks([0 pi/4 pi/2 3*pi/4])
xticklabels({'0','π/4','π/2','3π/4'})
xlim([0 pi-pi/8])


%%

All_prc=prctile(All_prev_rpe_r(All_prev_rpe_r<0),0:10:100);


for i=1:length(All_prc)-4
    mean_All_update(i)=nanmean(abs(All_angular_change_r(All_prev_rpe_r>=All_prc(i) & All_prev_rpe_r<=All_prc(i+4) & All_prev_rpe_r<0)));
    sem_All_update(i)=nanstd(abs(All_angular_change_r(All_prev_rpe_r>=All_prc(i) & All_prev_rpe_r<=All_prc(i+4) & All_prev_rpe_r<0)))/sqrt(sum(All_prev_rpe_r>=All_prc(i) & All_prev_rpe_r<=All_prc(i+4) & All_prev_rpe_r<0)-1);
    
    
    mean_All_true_update(i)=nanmean(abs(All_true_angular_change_r(All_prev_rpe_r>=All_prc(i) & All_prev_rpe_r<=All_prc(i+4) & All_prev_rpe_r<0)));
    sem_All_true_update(i)=nanstd(abs(All_true_angular_change_r(All_prev_rpe_r>=All_prc(i) & All_prev_rpe_r<=All_prc(i+4) & All_prev_rpe_r<0)))/sqrt(sum(All_prev_rpe_r>=All_prc(i) & All_prev_rpe_r<=All_prc(i+4) & All_prev_rpe_r<0)-1);
    
end


figure
subplot(2,2,1)
shadedErrorBar(All_prc(3:end-2),mean_All_update,sem_All_update,{'-','color',[0 0 0],'LineWidth',2},2)
xlabel('RPE')
ylabel(sprintf('Neural update (RPE < 0)'))
box off
title('All')

subplot(2,2,3)
shadedErrorBar(All_prc(3:end-2),mean_All_true_update,sem_All_true_update,{'-','color',[0 0 0],'LineWidth',2},2)
xlabel('RPE')
ylabel(sprintf('Behavioral update (RPE < 0)'))
box off
title('All')

All_angle=prctile(abs(All_dist_prev_decoded_angle_cc_r(All_prev_rpe_r<0)),0:10:100);
All_angle2=prctile(abs(All_dist_prev_true_peak_cc_r(All_prev_rpe_r<0)),0:10:100);

for i=1:length(All_angle)-4
    mean_All_update_angle(i)=nanmean(abs(All_angular_change_r(abs(All_dist_prev_decoded_angle_cc_r)>=All_angle(i) & abs(All_dist_prev_decoded_angle_cc_r)<=All_angle(i+4) & All_prev_rpe_r<0)));
    sem_All_update_angle(i)=nanstd(abs(All_angular_change_r(abs(All_dist_prev_decoded_angle_cc_r)>=All_angle(i) & abs(All_dist_prev_decoded_angle_cc_r)<=All_angle(i+4) & All_prev_rpe_r<0)))/sqrt(sum(abs(All_dist_prev_decoded_angle_cc_r)>=All_angle(i) & abs(All_dist_prev_decoded_angle_cc_r)<=All_angle(i+4) & All_prev_rpe_r<0)-1);
    
    mean_All_true_update_angle(i)=nanmean(abs(All_true_angular_change_r(abs(All_dist_prev_true_peak_cc_r)>=All_angle2(i) & abs(All_dist_prev_true_peak_cc_r)<=All_angle2(i+4) & All_prev_rpe_r<0)));
    sem_All_true_update_angle(i)=nanstd(abs(All_true_angular_change_r(abs(All_dist_prev_true_peak_cc_r)>=All_angle2(i) & abs(All_dist_prev_true_peak_cc_r)<=All_angle2(i+4) & All_prev_rpe_r<0)))/sqrt(sum(abs(All_dist_prev_true_peak_cc_r)>=All_angle2(i) & abs(All_dist_prev_true_peak_cc_r)<=All_angle2(i+4) & All_prev_rpe_r<0)-1);
    
end

subplot(2,2,2)
shadedErrorBar(All_angle(3:end-2),mean_All_update_angle,sem_All_update_angle,{'-','color',[0 0 0],'LineWidth',2},2)
xlabel(sprintf('Neural |d(ET(t-1), CC(t-1))|'))
ylabel(sprintf('Neural update magnitude (RPE < 0)'))
box off
xticks([0 pi/4 pi/2 3*pi/4])
xticklabels({'0','π/4','π/2','3π/4'})
xlim([0 pi-pi/8])

subplot(2,2,4)
shadedErrorBar(All_angle(3:end-2),mean_All_true_update_angle,sem_All_true_update_angle,{'-','color',[0 0 0],'LineWidth',2},2)
xlabel(sprintf('Model |d(ET(t-1), CC(t-1))|'))
ylabel(sprintf('Model update magnitude (RPE < 0)'))
box off
xticks([0 pi/4 pi/2 3*pi/4])
xticklabels({'0','π/4','π/2','3π/4'})
xlim([0 pi-pi/8])


%% figure 3B

load('colors.mat')
mean_decoded_angle=mod(angle(nanmean(exp(1i*Decoded_angle),1)),2*pi);
mean_mean_decoded_angle=mod(angle(nanmean(exp(1i*mean_decoded_angle),2)),2*pi);
mean_true_angle=mod(angle(nanmean(exp(1i*True_peak),1)),2*pi);

%%
disp('mean block color')
[z_dist_mean, p_dist_mean] = compute_shuffled_circ_dist(mean_decoded_angle,mean_true_angle,1000)

mean_dist=nanmean(abs(circular_distance(mean_decoded_angle,mean_true_angle)));

%%

[rho_mean, pval_mean] = circ_corrcc(mean_decoded_angle(~isnan(mean_decoded_angle)), mean_true_angle(~isnan(mean_decoded_angle)));
nanmean(abs(circular_distance(mean_decoded_angle,mean_true_angle)))

% binned_angle=0:pi/8:2*pi;
% mean_decoded_angle_binned=NaN(1,length(binned_angle));
% for i=3:length(binned_angle)-2
%     mean_decoded_angle_binned(i)=mod(angle(nanmean(exp(1i*Decoded_angle(mod(True_peak,2*pi)>=binned_angle(i-2) & mod(True_peak,2*pi)<binned_angle(i+2))),1)),2*pi);
% end

figure
subplot(1,2,1)
for i=1:length(mean_true_angle)
    if ~isnan(mean_true_angle(i))
        hold on        
        plot(mean_true_angle(i),mean_decoded_angle(i),'o','Color',colors(round(mod(mean_true_angle(i),2*pi)/2/pi*40),:),'MarkerFaceColor',colors(round(mod(mean_true_angle(i),2*pi)/2/pi*40),:))
    end
end
% plot(binned_angle,mean_decoded_angle_binned,'Color',[0.5 0.5 0.5],'LineWidth',2)
predicted=rho_mean*([0:1:8]-pi)+mean_mean_decoded_angle;

plot(0:8,predicted,'k--','LineWidth',2)
plot(0:8,0:8,'k-')

xticks([0 pi/2 pi 3*pi/2 2*pi])
xticklabels({'0','π/2','π','3π/2','2π'})
yticks([0 pi/2 pi 3*pi/2 2*pi])
yticklabels({'0','π/2','π','3π/2','2π'})
ylim([0 2*pi])
xlim([0 2*pi])
xlabel('Mean behavioral ET')
ylabel('Mean neural ET')
box off


%% figure 3B

count_distance=1;

for i=1:length(mean_true_angle)-1
    for j=i+1:length(mean_true_angle)
        if ~isnan(mean_true_angle(i)) && ~isnan(mean_true_angle(j))
            hold on
            dist_true_angles_on_wheel(count_distance)=circular_distance(mean_true_angle(i),mean_true_angle(j));
            dist_decoded_angles_on_wheel(count_distance)=circular_distance(mean_decoded_angle(i),mean_decoded_angle(j));
            count_distance=count_distance+1;
        end
    end
end

disp('distance between block colors')

%%
[z_dist_on_wheel, p_dist_on_wheel] = compute_shuffled_circ_dist(dist_true_angles_on_wheel,dist_decoded_angles_on_wheel,1000)

mean_dist_on_wheel=nanmean(abs(circular_distance(dist_true_angles_on_wheel,dist_decoded_angles_on_wheel)));

[rho_dist, pval_dist] = circ_corrcc(dist_decoded_angles_on_wheel(~isnan(dist_decoded_angles_on_wheel)), dist_true_angles_on_wheel(~isnan(dist_decoded_angles_on_wheel)));

%%

values=hist3([dist_decoded_angles_on_wheel' dist_true_angles_on_wheel' ],[9 9]);
figure
subplot(1,2,1)
imagesc(values)
hold on
plot(0:10,0:10,'r-','Linewidth',2)
set(gca,'Ydir','normal')
set(gca,'Xdir','normal')
colorbar
xticks([1 3 5 7 9])
xticklabels({'-π','-π/2','0','π/2','π'})
yticks([1 3 5 7 9])
yticklabels({'-π','-π/2','0','π/2','π'})
xlabel('d(ET(i),ET(j)) behavioral')
ylabel('d(ET(i),ET(j)) neural')
box off
colormap hot

subplot(1,2,2)
plot(dist_true_angles_on_wheel,dist_decoded_angles_on_wheel,'k.')
hold on
plot(-5:5,-5:5,'r-','Linewidth',2)
xticks([-pi -pi/2 0 pi/2 pi])
xticklabels({'-π','-π/2','0','π/2','π'})
yticks([-pi -pi/2 0 pi/2 pi])
yticklabels({'-π','-π/2','0','π/2','π'})
ylim([-pi pi])
xlim([-pi pi])
xlabel('d(ET(i),ET(j)) behavioral')
ylabel('d(ET(i),ET(j)) neural')
box off




%% Functions

function [z_circ_dist, p_circ_dist] = compute_shuffled_circ_dist(dist1,dist2,N_perm)

dist_diff=nanmean(abs(mod(dist1-dist2+pi,2*pi)-pi));
shuffled_difference=zeros(N_perm,1);
for n=1:N_perm
    dist2=dist2(randperm(length(dist2)));
    shuffled_difference(n)=nanmean(abs(mod(dist1-dist2+pi,2*pi)-pi));
end

z_circ_dist=(dist_diff-mean(shuffled_difference))/std(shuffled_difference);
p_circ_dist=sum(dist_diff>shuffled_difference)/N_perm;

% figure
% histogram(shuffled_difference,'FaceColor',[0.5 0.5 0.5])
% hold on
% xline(dist_diff,'r','LineWidth',2)
% box off
% 
end


