function get_boot_glm_all_values_split_correct(ROI)

N_boot=5000;

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
save_path=fullfile(fsroot,dirstem2);

save_name=sprintf('glm_explained_variance_split_boot_%s',ROI);

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

switch ROI
    case 'LIP'
        color_for_plot=color_for_ROI(1,:);
        n_color=1;
    case 'FEF'
        color_for_plot=color_for_ROI(2,:);
        n_color=2;
    case 'PFC'
        color_for_plot=color_for_ROI(3,:);
        n_color=3;
end



%%

for i=1:length(arrayID_list)
    arrayID=arrayID_list(i);
    
    if (arrayID==34 && strcmp(ROI,'LIP')) || (arrayID==30 && strcmp(ROI,'LIP')) || (arrayID==14 && strcmp(ROI,'FEF'))
        
    else
        for this_time=1:length(window_start_list)
            load_name=sprintf('glm_explained_variance_split_%s_%d_%d',ROI,arrayID,this_time);
            load(fullfile(save_path,load_name),'beta_split','beta_all_values_split')
            
            for n=1:size(beta_split,1)
                
                All_beta_split(count_ROI(this_time),:,:,:,:,this_time)=beta_split(n,:,:,:,:);
                All_beta_all_values_split(count_ROI(this_time),:,:,:,this_time)=beta_all_values_split(n,:,:,:);
                                
                count_ROI(this_time)=count_ROI(this_time)+1;
            end
            
            clear  beta_split beta_all_values_split
        end
    end
end



%%
All_beta_split_contra=vertcat(All_beta_split(:,:,:,:,:,:),All_beta_split(:,:,:,:,:,:));
All_beta_all_values_split_contra=vertcat(All_beta_all_values_split(:,:,:,:,:),All_beta_all_values_split(:,:,:,:,:));

contra_beta_split=vertcat(ones(size(All_beta_split,1),1).*1,ones(size(All_beta_split,1),1).*2);

for i=1:size(contra_beta_split,1)
    if contra_beta_split(i)==1
        Ipsi_contra_beta_split(i)=4;
        contra_diag_contra_beta_split(i)=3;
        contra_contra_beta_split(i)=2;
    elseif contra_beta_split(i)==2
        Ipsi_contra_beta_split(i)=3;
        contra_diag_contra_beta_split(i)=4;
        contra_contra_beta_split(i)=1;
    end
end

%% Contra

for i=1:size(All_beta_split,5) %rep
    for t=1:size(All_beta_split,6) %time
        r_chosen_contra(i,t)=calculate_chosen_corr(All_beta_split_contra,contra_beta_split,contra_beta_split,1,2,i,t);
        r_chosen_contra_ipsi(i,t)=calculate_chosen_corr(All_beta_split_contra,contra_beta_split,Ipsi_contra_beta_split,1,2,i,t);
        r_chosen_contra_contra(i,t)=calculate_chosen_corr(All_beta_split_contra,contra_beta_split,contra_contra_beta_split,1,2,i,t);
        r_chosen_contra_contra_diag(i,t)=calculate_chosen_corr(All_beta_split_contra,contra_beta_split,contra_diag_contra_beta_split,1,2,i,t);
        
        r_chosen_contra_ipsi(i+size(All_beta_split,5), t)=calculate_chosen_corr(All_beta_split_contra,contra_beta_split,Ipsi_contra_beta_split,2,1,i,t);
        r_chosen_contra_contra(i+size(All_beta_split,5), t)=calculate_chosen_corr(All_beta_split_contra,contra_beta_split,contra_contra_beta_split,2,1,i,t);
        r_chosen_contra_contra_diag(i+size(All_beta_split,5), t)=calculate_chosen_corr(All_beta_split_contra,contra_beta_split,contra_diag_contra_beta_split,2,1,i,t);
        
        r_chosen_contra_chosen_global(i,t)=calculate_chosen_corr_global(All_beta_split_contra,All_beta_all_values_split_contra,contra_beta_split,1,2,i,t,6);
        r_chosen_contra_reward_global(i,t)=calculate_chosen_corr_global(All_beta_split_contra,All_beta_all_values_split_contra,contra_beta_split,1,2,i,t,3);
        
        r_chosen_contra_chosen_global(i+size(All_beta_split,5), t)=calculate_chosen_corr_global(All_beta_split_contra,All_beta_all_values_split_contra,contra_beta_split,2,1,i,t,6);
        r_chosen_contra_reward_global(i+size(All_beta_split,5), t)=calculate_chosen_corr_global(All_beta_split_contra,All_beta_all_values_split_contra,contra_beta_split,2,1,i,t,3);
        
        r_unchosen_contra(i,t)=calculate_unchosen_corr(All_beta_split_contra,contra_beta_split,contra_beta_split,1,2,i,t);
        r_chosen_unchosen_contra(i,t)=calculate_chosen_unchosen_corr(All_beta_split_contra,contra_beta_split,contra_beta_split,1,2,i,t);
        r_chosen_unchosen_contra(i+size(All_beta_split,5), t)=calculate_chosen_unchosen_corr(All_beta_split_contra,contra_beta_split,contra_beta_split,2,1,i,t);
        
        r_chosen_choice_contra(i,t)=calculate_chosen_choice_corr(All_beta_split_contra,contra_beta_split,contra_beta_split,1,2,i,t);
        r_chosen_choice_contra(i+size(All_beta_split,5),t)=calculate_chosen_choice_corr(All_beta_split_contra,contra_beta_split,contra_beta_split,2,1,i,t);
       
        r_global_chosen_choice_contra(i,t)=calculate_chosen_choice_corr_global(All_beta_split_contra,All_beta_all_values_split_contra,contra_beta_split,2,1,i,t);
        r_global_chosen_choice_contra(i+size(All_beta_split,5),t)=calculate_chosen_choice_corr_global(All_beta_split_contra,All_beta_all_values_split_contra,contra_beta_split,1,2,i,t);
        
    end
end

%%
%cross corr
for i=1:size(All_beta_split,5) %rep
    for t=1:size(All_beta_split,6) %time
        for p=1:size(All_beta_split,6) %time
            r_chosen_contra_cross_corr(t,p,i)=nancorr(All_beta_split_contra(:,3,contra_beta_split(i),1,i,t),All_beta_split_contra(:,3,contra_beta_split(i),2,i,p));
            r_chosen_contra_cross_corr(t,p,i+size(All_beta_split,5))=nancorr(All_beta_split_contra(:,3,contra_beta_split(i),2,i,t),All_beta_split_contra(:,3,contra_beta_split(i),1,i,p));
            
            r_chosen_contra_global_cross_corr(t,p,i)=nancorr(All_beta_split_contra(:,3,contra_beta_split(i),1,i,t),All_beta_all_values_split_contra(:,6,2,i,p));
            r_chosen_contra_global_cross_corr(t,p,i+size(All_beta_split,5))=nancorr(All_beta_split_contra(:,3,contra_beta_split(i),2,i,t),All_beta_all_values_split_contra(:,6,1,i,p));
            
            r_chosen_choice_contra_cross_corr(t,p,i)=nancorr(All_beta_split_contra(:,3,contra_beta_split(i),1,i,t),All_beta_split_contra(:,2,contra_beta_split(i),2,i,p));
            r_chosen_choice_contra_cross_corr(t,p,i+size(All_beta_split,5))=nancorr(All_beta_split_contra(:,3,contra_beta_split(i),2,i,t),All_beta_split_contra(:,2,contra_beta_split(i),1,i,p));
            
            r_chosen_global_cross_corr(t,p,i)=nancorr(All_beta_all_values_split_contra(:,6,1,i,t),All_beta_all_values_split_contra(:,6,2,i,p));
            r_chosen_global_cross_corr(t,p,i+size(All_beta_split,5))=nancorr(All_beta_all_values_split_contra(:,6,2,i,t),All_beta_all_values_split_contra(:,6,1,i,p));
        end
    end
end

for t=1:size(All_beta_split,6) %time
    for p=1:size(All_beta_split,6) %time
        if mean(r_chosen_contra_cross_corr(t,p,:))>0 && mean(r_chosen_global_cross_corr(t,p,:))>0
            diss_chosen_contra_global_cross_corr(t,p)=mean(r_chosen_contra_global_cross_corr(t,p,:))/sqrt(mean(r_chosen_contra_cross_corr(t,p,:))*mean(r_chosen_global_cross_corr(t,p,:)));
        end
    end
end


%%
ipsi_beta_split=vertcat(ones(size(All_beta_split,1),1).*3,ones(size(All_beta_split,1),1).*4);

for i=1:size(ipsi_beta_split,1)
    if ipsi_beta_split(i)==3
        ipsi_ipsi_beta_split(i)=4;
        ipsi_diag_contra_beta_split(i)=1;
        ipsi_contra_beta_split(i)=2;
    elseif ipsi_beta_split(i)==4
        ipsi_ipsi_beta_split(i)=3;
        ipsi_diag_contra_beta_split(i)=2;
        ipsi_contra_beta_split(i)=1;
    end
end

for i=1:size(All_beta_split,5) %rep
    for t=1:size(All_beta_split,6) %time
        r_chosen_ipsi(i,t)=calculate_chosen_corr(All_beta_split_contra,ipsi_beta_split,ipsi_beta_split,1,2,i,t);
        r_chosen_ipsi_ipsi(i,t)=calculate_chosen_corr(All_beta_split_contra,ipsi_beta_split,ipsi_ipsi_beta_split,1,2,i,t);
        r_chosen_ipsi_contra(i,t)=calculate_chosen_corr(All_beta_split_contra,ipsi_beta_split,ipsi_contra_beta_split,1,2,i,t);
        r_chosen_ipsi_contra_diag(i,t)=calculate_chosen_corr(All_beta_split_contra,ipsi_beta_split,ipsi_diag_contra_beta_split,1,2,i,t);
        
        r_chosen_ipsi_ipsi(i+size(All_beta_split,5), t)=calculate_chosen_corr(All_beta_split_contra,ipsi_beta_split,ipsi_ipsi_beta_split,2,1,i,t);
        r_chosen_ipsi_contra(i+size(All_beta_split,5), t)=calculate_chosen_corr(All_beta_split_contra,ipsi_beta_split,ipsi_contra_beta_split,2,1,i,t);
        r_chosen_ipsi_contra_diag(i+size(All_beta_split,5), t)=calculate_chosen_corr(All_beta_split_contra,ipsi_beta_split,ipsi_diag_contra_beta_split,2,1,i,t);
        
        r_unchosen_ipsi(i,t)=calculate_unchosen_corr(All_beta_split_contra,ipsi_beta_split,ipsi_beta_split,1,2,i,t);
        r_chosen_unchosen_ipsi(i,t)=calculate_chosen_unchosen_corr(All_beta_split_contra,ipsi_beta_split,ipsi_ipsi_beta_split,1,2,i,t);
        r_chosen_unchosen_ipsi(i+size(All_beta_split,5), t)=calculate_chosen_unchosen_corr(All_beta_split_contra,ipsi_beta_split,ipsi_ipsi_beta_split,2,1,i,t);
        
        r_chosen_ipsi_chosen_global(i,t)=calculate_chosen_corr_global(All_beta_split_contra,All_beta_all_values_split_contra,ipsi_beta_split,1,2,i,t,6);
        r_chosen_ipsi_reward_global(i,t)=calculate_chosen_corr_global(All_beta_split_contra,All_beta_all_values_split_contra,ipsi_beta_split,1,2,i,t,3);
        
        r_chosen_ipsi_chosen_global(i+size(All_beta_split,5), t)=calculate_chosen_corr_global(All_beta_split_contra,All_beta_all_values_split_contra,ipsi_beta_split,2,1,i,t,6);
        r_chosen_ipsi_reward_global(i+size(All_beta_split,5), t)=calculate_chosen_corr_global(All_beta_split_contra,All_beta_all_values_split_contra,ipsi_beta_split,2,1,i,t,3);
    end
end


%%
for i=1:size(All_beta_split,5) %rep
    for l=1:size(All_beta_split,3) %loc
        for k=1:size(All_beta_split,3) %loc
            for t=1:size(All_beta_split,6) %time
                r_chosen(t,l,k,i)=nancorr(All_beta_split(:,3,l,1,i,t),All_beta_split(:,3,k,2,i,t));
            end
        end
    end
end

for i=1:size(All_beta_split,5) %rep
    for t=1:size(All_beta_split,6) %time
        r_chosen_global(i,t)=nancorr(All_beta_all_values_split(:,6,1,i,t),All_beta_all_values_split(:,6,2,i,t));
        r_unchosen_global(i,t)=nancorr(All_beta_all_values_split(:,7,1,i,t),All_beta_all_values_split(:,7,2,i,t));
        r_chosen_unchosen_global(i,t)=nancorr(All_beta_all_values_split(:,6,1,i,t),All_beta_all_values_split(:,7,2,i,t));
        r_chosen_unchosen_global(i+size(All_beta_split,5),t)=nancorr(All_beta_all_values_split(:,7,1,i,t),All_beta_all_values_split(:,6,2,i,t));
        r_reward_global(i,t)=nancorr(All_beta_all_values_split(:,3,1,i,t),All_beta_all_values_split(:,3,2,i,t));
    end
end

for i=1:size(All_beta_split,5) %rep
    for l=1:size(All_beta_split,3) %loc
        for t=1:size(All_beta_split,6) %time
            r_chosen_global_across(t,l,i)=nancorr(All_beta_split(:,3,l,1,i,t),All_beta_all_values_split(:,6,2,i,t));
            r_chosen_global_across(t,l,i+size(All_beta_split,5))=nancorr(All_beta_split(:,3,l,2,i,t),All_beta_all_values_split(:,6,1,i,t));
            
            r_reward_global_across(t,l,i)=nancorr(All_beta_split(:,3,l,1,i,t),All_beta_all_values_split(:,3,2,i,t));
            r_reward_global_across(t,l,i+size(All_beta_split,5))=nancorr(All_beta_split(:,3,l,2,i,t),All_beta_all_values_split(:,3,1,i,t));
        end
    end
end

%% Bootstrap for stability

for n=1:N_boot
    ind=randsample(size(All_beta_split,1),size(All_beta_split,1),'true');
    
    for i=1:size(All_beta_split,5) %rep
        for l=1:size(All_beta_split,3) %loc
            for k=1:size(All_beta_split,3) %loc
                for t=1:size(All_beta_split,6) %time
                    r_chosen_boot(t,l,k,i,n)=nancorr(All_beta_split(ind,3,l,1,i,t),All_beta_split(ind,3,k,2,i,t));
                end
            end
        end
    end
    
    %global non concatenated
    for i=1:size(All_beta_split,5) %rep
        for t=1:size(All_beta_split,6) %time
            r_chosen_global_global_boot(t,i,n)=nancorr(All_beta_all_values_split(ind,6,1,i,t),All_beta_all_values_split(ind,6,2,i,t));
            r_reward_global_global_boot(t,i,n)=nancorr(All_beta_all_values_split(ind,3,1,i,t),All_beta_all_values_split(ind,3,2,i,t));
            r_unchosen_global_global_boot(t,i,n)=nancorr(All_beta_all_values_split(ind,7,1,i,t),All_beta_all_values_split(ind,7,2,i,t));
            r_chosen_unchosen_global_global_boot(t,i,n)=nancorr(All_beta_all_values_split(ind,6,1,i,t),All_beta_all_values_split(ind,7,2,i,t));
            r_chosen_unchosen_global_global_boot(t,i+size(All_beta_split,5),n)=nancorr(All_beta_all_values_split(ind,7,1,i,t),All_beta_all_values_split(ind,6,2,i,t));
        end
    end
    
    for i=1:size(All_beta_split,5) %rep
        for l=1:size(All_beta_split,3) %loc
            for t=1:size(All_beta_split,6) %time
                r_chosen_global_across_boot(t,l,i,n)=nancorr(All_beta_split(ind,3,l,1,i,t),All_beta_all_values_split(ind,6,2,i,t));
                r_chosen_global_across_boot(t,l,i+size(All_beta_split,5),n)=nancorr(All_beta_split(ind,3,l,2,i,t),All_beta_all_values_split(ind,6,1,i,t));
                
                r_reward_global_across_boot(t,l,i,n)=nancorr(All_beta_split(ind,3,l,1,i,t),All_beta_all_values_split(ind,3,2,i,t));
                r_reward_global_across_boot(t,l,i+size(All_beta_split,5),n)=nancorr(All_beta_split(ind,3,l,2,i,t),All_beta_all_values_split(ind,3,1,i,t));
            end
        end
    end
     
    ind=randsample(size(All_beta_split_contra,1),size(All_beta_split_contra,1),'true');
    for i=1:size(All_beta_split,5) %rep
        for t=1:size(All_beta_split,6) %time
            
            %contra
            r_chosen_contra_boot(t,i,n)=calculate_chosen_corr_boot(All_beta_split_contra,contra_beta_split,contra_beta_split,1,2,i,t,ind);
            r_chosen_contra_ipsi_boot(t,i,n)=calculate_chosen_corr_boot(All_beta_split_contra,contra_beta_split,Ipsi_contra_beta_split,1,2,i,t,ind);
            r_chosen_contra_contra_boot(t,i,n)=calculate_chosen_corr_boot(All_beta_split_contra,contra_beta_split,contra_contra_beta_split,1,2,i,t,ind);
            r_chosen_contra_contra_diag_boot(t,i,n)=calculate_chosen_corr_boot(All_beta_split_contra,contra_beta_split,contra_diag_contra_beta_split,1,2,i,t,ind);
            
            r_chosen_contra_ipsi_boot(t,i+size(All_beta_split,5),n)=calculate_chosen_corr_boot(All_beta_split_contra,contra_beta_split,Ipsi_contra_beta_split,2,1,i,t,ind);
            r_chosen_contra_contra_boot(t,i+size(All_beta_split,5),n)=calculate_chosen_corr_boot(All_beta_split_contra,contra_beta_split,contra_contra_beta_split,2,1,i,t,ind);
            r_chosen_contra_contra_diag_boot(t,i+size(All_beta_split,5),n)=calculate_chosen_corr_boot(All_beta_split_contra,contra_beta_split,contra_diag_contra_beta_split,2,1,i,t,ind);
            
            r_unchosen_contra_boot(t,i,n)=calculate_unchosen_corr_boot(All_beta_split_contra,contra_beta_split,contra_beta_split,1,2,i,t,ind);
            r_chosen_unchosen_contra_boot(t,i,n)=calculate_chosen_unchosen_corr_boot(All_beta_split_contra,contra_beta_split,contra_beta_split,1,2,i,t,ind);
            r_chosen_unchosen_contra_boot(t,i+size(All_beta_split,5),n)=calculate_chosen_unchosen_corr_boot(All_beta_split_contra,contra_beta_split,contra_beta_split,2,1,i,t,ind);
            
            r_chosen_choice_contra_boot(t,i,n)=calculate_chosen_choice_corr_boot(All_beta_split_contra,contra_beta_split,contra_beta_split,1,2,i,t,ind);
            r_chosen_choice_contra_boot(t,i+size(All_beta_split,5),n)=calculate_chosen_choice_corr_boot(All_beta_split_contra,contra_beta_split,contra_beta_split,2,1,i,t,ind);
                        
            r_chosen_contra_chosen_global_boot(t,i,n)=calculate_chosen_corr_global_boot(All_beta_split_contra,All_beta_all_values_split_contra,contra_beta_split,1,2,i,t,6,ind);
            r_chosen_contra_reward_global_boot(t,i,n)=calculate_chosen_corr_global_boot(All_beta_split_contra,All_beta_all_values_split_contra,contra_beta_split,1,2,i,t,3,ind);
            
            r_chosen_contra_chosen_global_boot(t,i+size(All_beta_split,5),n)=calculate_chosen_corr_global_boot(All_beta_split_contra,All_beta_all_values_split_contra,contra_beta_split,2,1,i,t,6,ind);
            r_chosen_contra_reward_global_boot(t,i+size(All_beta_split,5),n)=calculate_chosen_corr_global_boot(All_beta_split_contra,All_beta_all_values_split_contra,contra_beta_split,2,1,i,t,3,ind);
            
            r_global_chosen_choice_contra_boot(t,i,n)=calculate_chosen_choice_corr_global_boot(All_beta_split_contra,All_beta_all_values_split_contra,contra_beta_split,2,1,i,t,ind);
            r_global_chosen_choice_contra_boot(t,i+size(All_beta_split,5),n)=calculate_chosen_choice_corr_global_boot(All_beta_split_contra,All_beta_all_values_split_contra,contra_beta_split,1,2,i,t,ind);
        
            %global concatenated to match the dim of contra
            r_chosen_global_boot(t,i,n)=calculate_corr_global_boot(All_beta_all_values_split_contra,1,2,i,t,6,ind);
            r_chosen_global_boot(t,i+size(All_beta_split,5),n)=calculate_corr_global_boot(All_beta_all_values_split_contra,2,1,i,t,6,ind);
            
            r_reward_global_boot(t,i,n)=calculate_corr_global_boot(All_beta_all_values_split_contra,1,2,i,t,3,ind);
            r_reward_global_boot(t,i+size(All_beta_split,5),n)=calculate_corr_global_boot(All_beta_all_values_split_contra,2,1,i,t,3,ind);
            
            %ipsi
            r_chosen_ipsi_boot(t,i,n)=calculate_chosen_corr_boot(All_beta_split_contra,ipsi_beta_split,ipsi_beta_split,1,2,i,t,ind);
            r_chosen_ipsi_ipsi_boot(t,i,n)=calculate_chosen_corr_boot(All_beta_split_contra,ipsi_beta_split,ipsi_ipsi_beta_split,1,2,i,t,ind);
            r_chosen_ipsi_contra_boot(t,i,n)=calculate_chosen_corr_boot(All_beta_split_contra,ipsi_beta_split,ipsi_contra_beta_split,1,2,i,t,ind);
            r_chosen_ipsi_contra_diag_boot(t,i,n)=calculate_chosen_corr_boot(All_beta_split_contra,ipsi_beta_split,ipsi_diag_contra_beta_split,1,2,i,t,ind);
            
            r_chosen_ipsi_ipsi_boot(t,i+size(All_beta_split,5),n)=calculate_chosen_corr_boot(All_beta_split_contra,ipsi_beta_split,ipsi_ipsi_beta_split,2,1,i,t,ind);
            r_chosen_ipsi_contra_boot(t,i+size(All_beta_split,5),n)=calculate_chosen_corr_boot(All_beta_split_contra,ipsi_beta_split,ipsi_contra_beta_split,2,1,i,t,ind);
            r_chosen_ipsi_contra_diag_boot(t,i+size(All_beta_split,5),n)=calculate_chosen_corr_boot(All_beta_split_contra,ipsi_beta_split,ipsi_diag_contra_beta_split,2,1,i,t,ind);
            
            r_unchosen_ipsi_boot(t,i,n)=calculate_unchosen_corr_boot(All_beta_split_contra,ipsi_beta_split,ipsi_ipsi_beta_split,1,2,i,t,ind);
            r_chosen_unchosen_ipsi_boot(t,i,n)=calculate_chosen_unchosen_corr_boot(All_beta_split_contra,ipsi_beta_split,ipsi_ipsi_beta_split,1,2,i,t,ind);
            r_chosen_unchosen_ipsi_boot(t,i+size(All_beta_split,5),n)=calculate_chosen_unchosen_corr_boot(All_beta_split_contra,ipsi_beta_split,ipsi_ipsi_beta_split,2,1,i,t,ind);

            r_chosen_ipsi_chosen_global_boot(t,i,n)=calculate_chosen_corr_global_boot(All_beta_split_contra,All_beta_all_values_split_contra,ipsi_beta_split,1,2,i,t,6,ind);
            r_chosen_ipsi_reward_global_boot(t,i,n)=calculate_chosen_corr_global_boot(All_beta_split_contra,All_beta_all_values_split_contra,ipsi_beta_split,1,2,i,t,3,ind);
            
            r_chosen_ipsi_chosen_global_boot(t,i+size(All_beta_split,5),n)=calculate_chosen_corr_global_boot(All_beta_split_contra,All_beta_all_values_split_contra,ipsi_beta_split,2,1,i,t,6,ind);
            r_chosen_ipsi_reward_global_boot(t,i+size(All_beta_split,5),n)=calculate_chosen_corr_global_boot(All_beta_split_contra,All_beta_all_values_split_contra,ipsi_beta_split,2,1,i,t,3,ind);
        end
    end
    
    %cross corr
    for i=1:size(All_beta_split,5) %rep
        for t=1:size(All_beta_split,6) %time
            for p=1:size(All_beta_split,6) %time
                r_chosen_contra_cross_corr_boot(t,p,i,n)=nancorr(All_beta_split_contra(ind,3,contra_beta_split(i),1,i,t),All_beta_split_contra(ind,3,contra_beta_split(i),2,i,p));
                r_chosen_contra_cross_corr_boot(t,p,i+size(All_beta_split,5),n)=nancorr(All_beta_split_contra(ind,3,contra_beta_split(i),2,i,t),All_beta_split_contra(ind,3,contra_beta_split(i),1,i,p));
                
                r_chosen_contra_global_cross_corr_boot(t,p,i,n)=nancorr(All_beta_split_contra(ind,3,contra_beta_split(i),1,i,t),All_beta_all_values_split_contra(ind,6,2,i,p));
                r_chosen_contra_global_cross_corr_boot(t,p,i+size(All_beta_split,5),n)=nancorr(All_beta_split_contra(ind,3,contra_beta_split(i),2,i,t),All_beta_all_values_split_contra(ind,6,1,i,p));
                
                r_chosen_choice_contra_cross_corr_boot(t,p,i,n)=nancorr(All_beta_split_contra(ind,3,contra_beta_split(i),1,i,t),All_beta_split_contra(ind,2,contra_beta_split(i),2,i,p));
                r_chosen_choice_contra_cross_corr_boot(t,p,i+size(All_beta_split,5),n)=nancorr(All_beta_split_contra(ind,3,contra_beta_split(i),2,i,t),All_beta_split_contra(ind,2,contra_beta_split(i),1,i,p));
                
                r_chosen_global_cross_corr_boot(t,p,i,n)=nancorr(All_beta_all_values_split_contra(ind,6,1,i,t),All_beta_all_values_split_contra(ind,6,2,i,p));
                r_chosen_global_cross_corr_boot(t,p,i+size(All_beta_split,5),n)=nancorr(All_beta_all_values_split_contra(ind,6,2,i,t),All_beta_all_values_split_contra(ind,6,1,i,p));
            end
        end
    end
    

end


%%
for n=1:N_boot
    for t=1:size(All_beta_split,6) %time
        for p=1:size(All_beta_split,6) %time
            if mean(r_chosen_contra_cross_corr_boot(t,p,:,n))>0 && mean(r_chosen_global_cross_corr_boot(t,p,:,n))>0
                diss_chosen_contra_global_cross_corr_boot(t,p,n)=mean(r_chosen_contra_global_cross_corr_boot(t,p,:,n))/sqrt(mean(r_chosen_contra_cross_corr_boot(t,p,:,n))*mean(r_chosen_global_cross_corr_boot(t,p,:,n)));
            else
                diss_chosen_contra_global_cross_corr_boot(t,p,n)=NaN;
            end
        end
    end
end

%%

for t=1:length(window_start_list)
    for p=1:size(All_beta_split,6) %time
        y=reshape(diss_chosen_contra_global_cross_corr_boot(t,p,:),size(r_chosen_contra_boot,3),1);
        p_diss_chosen_contra_global_cross_corr_boot(t,p) = z_test_function_bootstrap(y,0);
        
        y=reshape(mean(r_chosen_contra_cross_corr_boot(t,p,:,:),3),size(r_chosen_contra_boot,3),1);
        p_r_chosen_contra_cross_corr_boot(t,p) = z_test_function_bootstrap(y,0);
        
        y=reshape(mean(r_chosen_global_cross_corr_boot(t,p,:,:),3),size(r_chosen_contra_boot,3),1);
        p_r_chosen_global_cross_corr_boot(t,p) = z_test_function_bootstrap(y,0);
    end    
end



%% Dissatenuation coeff

for t=1:length(window_start_list)
    for n=1:N_boot
        if mean(r_chosen_contra_boot(t,:,n),2)>0 && mean(r_chosen_global_boot(t,:,n),2)>0
            diss_r_contra_chosen(t,n)=mean(r_chosen_contra_chosen_global_boot(t,:,n),2)./(sqrt(mean(r_chosen_contra_boot(t,:,n),2).*mean(r_chosen_global_boot(t,:,n),2)));
        else
            diss_r_contra_chosen(t,n)=NaN;
        end
        if mean(r_chosen_contra_boot(t,:,n),2)>0 && mean(r_reward_global_boot(t,:,n),2)>0
            diss_r_contra_reward(t,n)=mean(r_chosen_contra_chosen_global_boot(t,:,n),2)./(sqrt(mean(r_chosen_contra_boot(t,:,n),2).*mean(r_reward_global_boot(t,:,n),2)));
        else
            diss_r_contra_reward(t,n)=NaN;
        end
        if mean(r_chosen_contra_boot(t,:,n),2)>0 && mean(r_unchosen_contra_boot(t,:,n),2)>0
            diss_r_contra_chosen_unchosen(t,n)=mean(r_chosen_unchosen_contra_boot(t,:,n),2)./(sqrt(mean(r_chosen_contra_boot(t,:,n),2).*mean(r_unchosen_contra_boot(t,:,n),2)));
        else
            diss_r_contra_chosen_unchosen(t,n)=NaN;
        end
        if mean(r_chosen_global_global_boot(t,:,n),2)>0 && mean(r_unchosen_global_global_boot(t,:,n),2)>0
            diss_r_global_chosen_unchosen(t,n)=mean(r_chosen_unchosen_global_global_boot(t,:,n),2)./(sqrt(mean(r_chosen_global_global_boot(t,:,n),2).*mean(r_unchosen_global_global_boot(t,:,n),2)));
        else
            diss_r_global_chosen_unchosen(t,n)=NaN;
        end
    end
    p_diss_r_contra_chosen(t) = z_test_function_bootstrap(diss_r_contra_chosen(t,:),0);
    p_diss_r_contra_reward(t) = z_test_function_bootstrap(diss_r_contra_reward(t,:),0);
    
    p_diss_r_contra_chosen_unchosen(t) = z_test_function_bootstrap(diss_r_contra_chosen_unchosen(t,:),0);
    p_diss_r_global_chosen_unchosen(t) = z_test_function_bootstrap(diss_r_global_chosen_unchosen(t,:),0);
    
    for n=1:N_boot
        if mean(r_chosen_ipsi_boot(t,:,n),2)>0 && mean(r_chosen_global_boot(t,:,n),2)>0
            diss_r_ipsi_chosen(t,n)=mean(r_chosen_ipsi_chosen_global_boot(t,:,n),2)./(sqrt(mean(r_chosen_ipsi_boot(t,:,n),2).*mean(r_chosen_global_boot(t,:,n),2)));
        else
            diss_r_ipsi_chosen(t,n)=NaN;
        end
        if mean(r_chosen_ipsi_boot(t,:,n),2)>0 && mean(r_reward_global_boot(t,:,n),2)>0
            diss_r_ipsi_reward(t,n)=mean(r_chosen_ipsi_chosen_global_boot(t,:,n),2)./(sqrt(mean(r_chosen_ipsi_boot(t,:,n),2).*mean(r_reward_global_boot(t,:,n),2)));
        else
            diss_r_ipsi_reward(t,n)=NaN;
        end
    end
    p_diss_r_ipsi_chosen(t) = z_test_function_bootstrap(diss_r_ipsi_chosen(t,:),0);
    p_diss_r_ipsi_reward(t) = z_test_function_bootstrap(diss_r_ipsi_reward(t,:),0);
end

%%

for t=1:length(window_start_list)
    y=reshape(mean(r_chosen_contra_boot(t,:,:),2),size(r_chosen_contra_boot,3),1);
    p_chosen_contra_boot(t) = z_test_function_bootstrap(y,0);
    y=reshape(mean(r_chosen_contra_ipsi_boot(t,:,:),2),size(r_chosen_contra_boot,3),1);
    p_chosen_contra_ipsi_boot(t) = z_test_function_bootstrap(y,0);
    y=reshape(mean(r_chosen_contra_contra_boot(t,:,:),2),size(r_chosen_contra_boot,3),1);
    p_chosen_contra_contra_boot(t) = z_test_function_bootstrap(y,0);
    y=reshape(mean(r_chosen_contra_contra_diag_boot(t,:,:),2),size(r_chosen_contra_boot,3),1);
    p_chosen_contra_contra_diag_boot(t) = z_test_function_bootstrap(y,0);
    
    y=reshape(mean(r_unchosen_contra_boot(t,:,:),2),size(r_chosen_contra_boot,3),1);
    p_unchosen_contra_boot(t) = z_test_function_bootstrap(y,0);
    
    y=reshape(mean(r_chosen_contra_chosen_global_boot(t,:,:),2),size(r_chosen_contra_boot,3),1);
    p_contra_chosen_global_boot(t) = z_test_function_bootstrap(y,0);
    y=reshape(mean(r_chosen_contra_reward_global_boot(t,:,:),2),size(r_chosen_contra_boot,3),1);
    p_contra_reward_global_boot(t) = z_test_function_bootstrap(y,0);
    
    y=reshape(mean(r_chosen_contra_boot(t,:,:),2)-mean(r_chosen_contra_chosen_global_boot(t,:,:),2),size(r_chosen_contra_boot,3),1);
    p_chosen_vs_glaobal_contra_boot(t) = z_test_function_bootstrap(y,0);
    
    y=reshape(mean(r_chosen_global_global_boot(t,:,:),2),size(r_chosen_contra_boot,3),1);
    p_chosen_global_global_boot(t) = z_test_function_bootstrap(y,0);
    y=reshape(mean(r_unchosen_global_global_boot(t,:,:),2),size(r_chosen_contra_boot,3),1);
    p_unchosen_global_global_boot(t) = z_test_function_bootstrap(y,0);
    y=reshape(mean(r_reward_global_boot(t,:,:),2),size(r_chosen_contra_boot,3),1);
    p_reward_global_global_boot(t) = z_test_function_bootstrap(y,0);
end


for t=1:length(window_start_list)
    y=reshape(mean(r_chosen_ipsi_boot(t,:,:),2),size(r_chosen_contra_boot,3),1);
    p_chosen_ipsi_boot(t) = z_test_function_bootstrap(y,0);
    y=reshape(mean(r_chosen_ipsi_ipsi_boot(t,:,:),2),size(r_chosen_contra_boot,3),1);
    p_chosen_ipsi_ipsi_boot(t) = z_test_function_bootstrap(y,0);
    y=reshape(mean(r_chosen_ipsi_contra_boot(t,:,:),2),size(r_chosen_contra_boot,3),1);
    p_chosen_ipsi_contra_boot(t) = z_test_function_bootstrap(y,0);
    y=reshape(mean(r_chosen_ipsi_contra_diag_boot(t,:,:),2),size(r_chosen_contra_boot,3),1);
    p_chosen_ipsi_contra_diag_boot(t) = z_test_function_bootstrap(y,0);
    
    y=reshape(mean(r_unchosen_ipsi_boot(t,:,:),2),size(r_chosen_ipsi_boot,3),1);
    p_unchosen_ipsi_boot(t) = z_test_function_bootstrap(y,0);
    
    y=reshape(mean(r_chosen_ipsi_chosen_global_boot(t,:,:),2),size(r_chosen_contra_boot,3),1);
    p_ipsi_chosen_global_boot(t) = z_test_function_bootstrap(y,0);
    y=reshape(mean(r_chosen_ipsi_reward_global_boot(t,:,:),2),size(r_chosen_contra_boot,3),1);
    p_ipsi_reward_global_boot(t) = z_test_function_bootstrap(y,0);
    
    y=reshape(mean(r_chosen_ipsi_boot(t,:,:),2)-mean(r_chosen_ipsi_chosen_global_boot(t,:,:),2),size(r_chosen_contra_boot,3),1);
    p_chosen_vs_global_ipsi_boot(t) = z_test_function_bootstrap(y,0);
end


%% save

save(fullfile(save_path, save_name))

end



%% functions


function r=calculate_chosen_corr(beta,pref,other,ind1,ind2,this_rep,this_time)

for i=1:size(beta,1)
    x(i)=beta(i,3,pref(i),ind1,this_rep,this_time);
    y(i)=beta(i,3,other(i),ind2,this_rep,this_time);
end

if sum(~isnan(x) & ~isnan(y))>2
    x2=x(~isnan(x) & ~isnan(y));
    y2=y(~isnan(x) & ~isnan(y));
    
    r=corr(x2',y2');
else
    r=NaN;
end
end

function r=calculate_unchosen_corr(beta,pref,other,ind1,ind2,this_rep,this_time)
for i=1:size(beta,1)
    x(i)=beta(i,4,pref(i),ind1,this_rep,this_time);
    y(i)=beta(i,4,other(i),ind2,this_rep,this_time);
end

if sum(~isnan(x) & ~isnan(y))>2
    x2=x(~isnan(x) & ~isnan(y));
    y2=y(~isnan(x) & ~isnan(y));
    
    r=corr(x2',y2');
else
    r=NaN;
end
end

function r=calculate_chosen_unchosen_corr(beta,pref,other,ind1,ind2,this_rep,this_time)
for i=1:size(beta,1)
    x(i)=beta(i,3,pref(i),ind1,this_rep,this_time);
    y(i)=beta(i,4,other(i),ind2,this_rep,this_time);
end

if sum(~isnan(x) & ~isnan(y))>2
    x2=x(~isnan(x) & ~isnan(y));
    y2=y(~isnan(x) & ~isnan(y));
    
    r=corr(x2',y2');
else
    r=NaN;
end
end

function r=calculate_chosen_choice_corr(beta,pref,other,ind1,ind2,this_rep,this_time)
for i=1:size(beta,1)
    x(i)=beta(i,2,pref(i),ind1,this_rep,this_time);
    y(i)=beta(i,3,other(i),ind2,this_rep,this_time);
end

if sum(~isnan(x) & ~isnan(y))>2
    x2=x(~isnan(x) & ~isnan(y));
    y2=y(~isnan(x) & ~isnan(y));
    
    r=corr(x2',y2');
else
    r=NaN;
end
end

function r = calculate_chosen_corr_global(beta,beta_global,pref,ind1,ind2,this_rep,this_time,this_reg)
for i=1:size(beta,1)
    x(i)=beta(i,3,pref(i),ind1,this_rep,this_time); %chosen
    y(i)=beta_global(i,this_reg,ind2,this_rep,this_time); %chosen or reward global
end
if sum(~isnan(x) & ~isnan(y))>2
    
    x2=x(~isnan(x) & ~isnan(y));
    y2=y(~isnan(x) & ~isnan(y));
    
    r=corr(x2',y2');
else
    r=NaN;
end
end

function r = calculate_chosen_choice_corr_global(beta,beta_global,pref,ind1,ind2,this_rep,this_time)
for i=1:size(beta,1)
    x(i)=beta(i,2,pref(i),ind1,this_rep,this_time); %choice
    y(i)=beta_global(i,6,ind2,this_rep,this_time); %chosen global
end
if sum(~isnan(x) & ~isnan(y))>2
    
    x2=x(~isnan(x) & ~isnan(y));
    y2=y(~isnan(x) & ~isnan(y));
    
    r=corr(x2',y2');
else
    r=NaN;
end
end

function r=calculate_chosen_corr_boot(beta,pref,other,ind1,ind2,this_rep,this_time,ind_boot)

for i=1:size(beta,1)
    x(i)=beta(ind_boot(i),3,pref(ind_boot(i)),ind1,this_rep,this_time);
    y(i)=beta(ind_boot(i),3,other(ind_boot(i)),ind2,this_rep,this_time);
end
if sum(~isnan(x) & ~isnan(y))>2
    
    x2=x(~isnan(x) & ~isnan(y));
    y2=y(~isnan(x) & ~isnan(y));
    
    r=corr(x2',y2');
else
    r=NaN;
end
end

function r=calculate_unchosen_corr_boot(beta,pref,other,ind1,ind2,this_rep,this_time,ind_boot)

for i=1:size(beta,1)
    x(i)=beta(ind_boot(i),4,pref(ind_boot(i)),ind1,this_rep,this_time);
    y(i)=beta(ind_boot(i),4,other(ind_boot(i)),ind2,this_rep,this_time);
end
if sum(~isnan(x) & ~isnan(y))>2
    
    x2=x(~isnan(x) & ~isnan(y));
    y2=y(~isnan(x) & ~isnan(y));
    
    r=corr(x2',y2');
else
    r=NaN;
end
end

function r=calculate_chosen_unchosen_corr_boot(beta,pref,other,ind1,ind2,this_rep,this_time,ind_boot)

for i=1:size(beta,1)
    x(i)=beta(ind_boot(i),3,pref(ind_boot(i)),ind1,this_rep,this_time);
    y(i)=beta(ind_boot(i),4,other(ind_boot(i)),ind2,this_rep,this_time);
end
if sum(~isnan(x) & ~isnan(y))>2
    
    x2=x(~isnan(x) & ~isnan(y));
    y2=y(~isnan(x) & ~isnan(y));
    
    r=corr(x2',y2');
else
    r=NaN;
end
end

function r=calculate_chosen_choice_corr_boot(beta,pref,other,ind1,ind2,this_rep,this_time,ind_boot)

for i=1:size(beta,1)
    x(i)=beta(ind_boot(i),2,pref(ind_boot(i)),ind1,this_rep,this_time);
    y(i)=beta(ind_boot(i),3,other(ind_boot(i)),ind2,this_rep,this_time);
end
if sum(~isnan(x) & ~isnan(y))>2
    
    x2=x(~isnan(x) & ~isnan(y));
    y2=y(~isnan(x) & ~isnan(y));
    
    r=corr(x2',y2');
else
    r=NaN;
end
end


function r = calculate_chosen_corr_global_boot(beta,beta_global,pref,ind1,ind2,this_rep,this_time,this_reg,ind_boot)
for i=1:size(beta,1)
    x(i)=beta(ind_boot(i),3,pref(ind_boot(i)),ind1,this_rep,this_time);
    y(i)=beta_global(ind_boot(i),this_reg,ind2,this_rep,this_time);
end
if sum(~isnan(x) & ~isnan(y))>2
    
    x2=x(~isnan(x) & ~isnan(y));
    y2=y(~isnan(x) & ~isnan(y));
    
    r=corr(x2',y2');
else
    r=NaN;
end
end

function r = calculate_chosen_choice_corr_global_boot(beta,beta_global,pref,ind1,ind2,this_rep,this_time,ind_boot)
for i=1:size(beta,1)
    x(i)=beta(ind_boot(i),2,pref(ind_boot(i)),ind1,this_rep,this_time); %choice
    y(i)=beta_global(ind_boot(i),6,ind2,this_rep,this_time); %global chosen
end
if sum(~isnan(x) & ~isnan(y))>2
    
    x2=x(~isnan(x) & ~isnan(y));
    y2=y(~isnan(x) & ~isnan(y));
    
    r=corr(x2',y2');
else
    r=NaN;
end
end


function r = calculate_corr_global_boot(beta_global,ind1,ind2,this_rep,this_time,this_reg,ind_boot)
for i=1:size(beta_global,1)
    x(i)=beta_global(ind_boot(i),this_reg,ind1,this_rep,this_time);
    y(i)=beta_global(ind_boot(i),this_reg,ind2,this_rep,this_time);
end
if sum(~isnan(x) & ~isnan(y))>2
    
    x2=x(~isnan(x) & ~isnan(y));
    y2=y(~isnan(x) & ~isnan(y));
    
    r=corr(x2',y2');
else
    r=NaN;
end
end


function r = nancorr(x,y)
if sum(~isnan(x) & ~isnan(y))>2
    
    x2=x(~isnan(x) & ~isnan(y));
    y2=y(~isnan(x) & ~isnan(y));
    
    r=corr(x2,y2);
else
    r=NaN;
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


