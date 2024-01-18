%check number of blocks in pseudo-pop
fsroot='/Volumes/buschman';


for i=1:7
    initial_window=0;
    event_list{i}='reward_end';
    window_start_list(i)=initial_window+(i-1)*50;
end

for i=1:20
    initial_window=-600;
    event_list{i+7}='target';
    window_start_list(i+7)=initial_window+(i-1)*50;
end

task='Learning_Attentional_Templates';
% subtask='exploreexploit/Reset_RW_model';
subtask='exploreexploit/Restricted_FEF/Reset_RW_model';

window_size=300;

subsubtask_classifier=sprintf('Classifier_%s_%d_to_%s_%d',event_list{1},window_start_list(1),event_list{end},window_start_list(end));
% subsubtask_classifier=sprintf('Classifier_%s_%d_to_%s_%d',event,window_start_list(1),event,window_start_list(end));
dirstem = fullfile('Projects',task,'Analysis','Electrophy_analysis',subtask,subsubtask_classifier);
data_path_clasifier = fullfile(fsroot,dirstem);
% boot_name='Pseudo_pop_peak_belief_prog_classifier_bootstrap_more_trials'; %added peak

boot_name='Pseudo_pop_peak_belief_prog_classifier_bootstrap_4bins';

load(fullfile(data_path_clasifier,boot_name),'Pseudo_pop_bootstrap')

color_for_ROI=[0 0.447 0.741; 0.85 0.325 0.098; 0.929 0.694 0.125]; %LIP FEF PFC


%% Blocks the trials were taken from

N_prog=3;
N_bins=size(Pseudo_pop_bootstrap.LIP.progression.TT_trials,5);

for n=1:100
    clear Block*
    N_neurons=size(Pseudo_pop_bootstrap.LIP.progression.TT_trials,2);
    
    for k=1%:N_prog
        for nn=1:N_neurons
            for n_c=1:N_bins
                for i=1:size(Pseudo_pop_bootstrap.LIP.progression.TT_trials,3)
                    Block_TT_LIP(i,nn,n_c,k)=Pseudo_pop_bootstrap.LIP.Pseudo_pop(Pseudo_pop_bootstrap.LIP.Neurons(n,nn),Pseudo_pop_bootstrap.LIP.progression.TT_trials(n,nn,i,3,n_c)).Block_nb;
                end
                
            end
        end
    end
%     %combine across progression levels
%     Block_TT_LIP=vertcat(Block_TT_LIP(:,:,:,1), Block_TT_LIP(:,:,:,2), Block_TT_LIP(:,:,:,3));
    
    for nn=1:N_neurons
        for n_c=1:N_bins
            Unique_block_LIP(nn,n_c,n)=length(unique(Block_TT_LIP(:,nn,n_c)));
            
        end
    end
    
    N_neurons=size(Pseudo_pop_bootstrap.FEF.progression.TT_trials,2);
    
    for k=1%:N_prog
        for nn=1:N_neurons
            for n_c=1:N_bins
                for i=1:size(Pseudo_pop_bootstrap.FEF.progression.TT_trials,3)
                    Block_TT_FEF(i,nn,n_c,k)=Pseudo_pop_bootstrap.FEF.Pseudo_pop(Pseudo_pop_bootstrap.FEF.Neurons(n,nn),Pseudo_pop_bootstrap.FEF.progression.TT_trials(n,nn,i,3,n_c)).Block_nb;
                end
                
            end
        end
    end
%     %combine across progression levels
%     Block_TT_FEF=vertcat(Block_TT_FEF(:,:,:,1), Block_TT_FEF(:,:,:,2), Block_TT_FEF(:,:,:,3));
    
    for nn=1:N_neurons
        for n_c=1:N_bins
            Unique_block_FEF(nn,n_c,n)=length(unique(Block_TT_FEF(:,nn,n_c)));
            
        end
    end
    
    
    N_neurons=size(Pseudo_pop_bootstrap.PFC.progression.TT_trials,2);
    
    for k=1%:N_prog
        for nn=1:N_neurons
            for n_c=1:N_bins
                for i=1:size(Pseudo_pop_bootstrap.PFC.progression.TT_trials,3)
                    Block_TT_PFC(i,nn,n_c,k)=Pseudo_pop_bootstrap.PFC.Pseudo_pop(Pseudo_pop_bootstrap.PFC.Neurons(n,nn),Pseudo_pop_bootstrap.PFC.progression.TT_trials(n,nn,i,3,n_c)).Block_nb;
                end
                
            end
        end
    end
%     %combine across progression levels
%     Block_TT_PFC=vertcat(Block_TT_PFC(:,:,:,1), Block_TT_PFC(:,:,:,2), Block_TT_PFC(:,:,:,3));
    
    for nn=1:N_neurons
        for n_c=1:N_bins
            Unique_block_PFC(nn,n_c,n)=length(unique(Block_TT_PFC(:,nn,n_c)));
            
        end
    end
    
end

%%
figure
for i=1:4
    subplot(3,4,i)
    histogram(Unique_block_LIP(:,i,3),'Facecolor',color_for_ROI(1,:))
    box off
    ylabel('# neurons')
    xlabel('# blocks in bin')
    
    subplot(3,4,i+4)
    histogram(Unique_block_FEF(:,i,3),'Facecolor',color_for_ROI(2,:))
    box off
    ylabel('# neurons')
    xlabel('# blocks in bin')
    
    subplot(3,4,i+8)
    histogram(Unique_block_PFC(:,i,3),'Facecolor',color_for_ROI(3,:))
    box off
    ylabel('# neurons')
    xlabel('# blocks in bin')
    
    
end

for nc=1:N_bins
    for n=1:100
        one_block_LIP(nc,n)=length(find(Unique_block_LIP(:,nc,n)==1));
        one_block_FEF(nc,n)=length(find(Unique_block_FEF(:,nc,n)==1));
        one_block_PFC(nc,n)=length(find(Unique_block_PFC(:,nc,n)==1));
    end
end


%%
figure
for i=1:3
    subplot(3,3,i)
    histogram(mean(Unique_block_LIP(:,i,:),1),'Facecolor',color_for_ROI(1,:))
    box off
    ylabel('# neurons')
    xlabel('# blocks in bin')
    
    subplot(3,3,i+3)
    histogram(mean(Unique_block_FEF(:,i,:),1),'Facecolor',color_for_ROI(2,:))
    box off
    ylabel('# neurons')
    xlabel('# blocks in bin')
    
    subplot(3,3,i+6)
    histogram(mean(Unique_block_PFC(:,i,:),1),'Facecolor',color_for_ROI(3,:))
    box off
    ylabel('# neurons')
    xlabel('# blocks in bin')
    
    
end



