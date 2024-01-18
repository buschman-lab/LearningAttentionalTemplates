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
Reset=NaN(35,200);

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
        
        %extract the reset
        for i=1:size(data_val,2)
            reset(:,i)=Conditions.switch(:,data_val(:,i));
        end
        
        for n_blocks=1:max(block_val_reshaped)
            N_sess_block(n_blocks+count)=j;
            Decoded_angle(:,n_blocks+count)=decoded_angle_net_reshaped(block_val_reshaped==n_blocks);
            True_peak(:,n_blocks+count)=true_peak_reshaped(block_val_reshaped==n_blocks);
            Prev_RPE(:,n_blocks+count)=prev_rpe_val_reshaped(block_val_reshaped==n_blocks);
            Prev_chosen_color(:,n_blocks+count)=prev_chosen_color_val_reshaped(block_val_reshaped==n_blocks);
            Reset(1,n_blocks+count)=reset(1,n_blocks);
            
%             Decoded_angle(isnan(True_peak(:,n_blocks+count)),n_blocks+count)=NaN; %remove when no peak (flat EV)
            
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
                
                Reset(i,n_blocks+count)=reset(i,n_blocks);
                
                N_sess_all_block(i,n_blocks+count)=j;
            end
        end
        clear  reset
        count=count+max(block_val_reshaped);
    end
end


%%

All_angular_change_around_reset=NaN(10,size(All_angular_change,2));
All_true_angular_change_around_reset=NaN(10,size(All_angular_change,2));

id_reset=NaN(1,size(All_angular_change,2));
for n=1:size(All_angular_change_around_reset,2) %blocks
    if ~isempty(find(Reset(:,n)==1,1))
        id_reset(n)=find(Reset(:,n)==1,1,'first');
        All_angular_change_around_reset(:,n)=abs(All_angular_change(id_reset(n):id_reset(n)+9,n));
         All_true_angular_change_around_reset(:,n)=abs(All_true_angular_change(id_reset(n):id_reset(n)+9,n));
         
         All_angular_change_around_reset(:,n)=All_angular_change_around_reset(:,n)-nanmean(All_angular_change_around_reset(:,n));
         All_true_angular_change_around_reset(:,n)=All_true_angular_change_around_reset(:,n)-nanmean(All_true_angular_change_around_reset(:,n));
         
%         if id_reset(n)==2
%             All_angular_change_around_reset(2:7,n)=All_angular_change(id_reset(n)-1:id_reset(n)+4,n);
%             All_true_angular_change_around_reset(2:7,n)=All_true_angular_change(id_reset(n)-1:id_reset(n)+4,n);
%         elseif id_reset(n)>=3
%             All_angular_change_around_reset(:,n)=All_angular_change(id_reset(n)-2:id_reset(n)+4,n);
%             All_true_angular_change_around_reset(:,n)=All_true_angular_change(id_reset(n)-2:id_reset(n)+4,n);
%         end
    end
end

%%
    figure
    subplot(1,2,1)
shadedErrorBar(0:9,nanmean(All_angular_change_around_reset,2),nanstd(All_angular_change_around_reset,0,2)./sqrt(sum(~isnan(All_angular_change_around_reset),2)-1),{'-','color',[0 0 0],'LineWidth',2},2)
xlabel(sprintf('Trials from reset'))
ylabel(sprintf('Angular change in ET'))
title('Neural update')
box off
    subplot(1,2,2)
shadedErrorBar(0:9,nanmean(All_true_angular_change_around_reset,2),nanstd(All_true_angular_change_around_reset,0,2)./sqrt(sum(~isnan(All_true_angular_change_around_reset),2)-1),{'-','color',[0 0 0],'LineWidth',2},2)
xlabel(sprintf('Trials from reset'))
ylabel(sprintf('Angular change in ET'))
box off
title('Behavioral update')

%%
y=reshape(All_angular_change_around_reset,1,size(All_angular_change_around_reset,1)*size(All_angular_change_around_reset,2));
basis=1:10;
x=reshape(basis'.*ones(size(All_angular_change_around_reset)),1,size(All_angular_change_around_reset,1)*size(All_angular_change_around_reset,2));

m=fitglm(x,y)

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

end


