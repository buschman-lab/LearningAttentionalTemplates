function get_boot_glm_all_values_split_pref(ROI)

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

save_name=sprintf('glm_explained_variance_split_boot_pref_%s',ROI);

count_ROI=ones(length(window_start_list),1);
count_ROI_pref=1;

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



%% Load everything

for i=1:length(arrayID_list)
    arrayID=arrayID_list(i);
    
    if (arrayID==34 && strcmp(ROI,'LIP')) || (arrayID==30 && strcmp(ROI,'LIP')) || (arrayID==14 && strcmp(ROI,'FEF'))
        
    else
        for this_time=1:length(window_start_list)
            load_name=sprintf('glm_explained_variance_split_%s_%d_%d',ROI,arrayID,this_time);
            load(fullfile(save_path,load_name),'beta_split')
            
            for n=1:size(beta_split,1)
                
                All_beta_split(count_ROI(this_time),:,:,:,:,this_time)=beta_split(n,:,:,:,:);
                
                count_ROI(this_time)=count_ROI(this_time)+1;
            end
            
            clear  beta_split beta_all_values_split
        end
        
        load_name=sprintf('get_pref_loc_%s_%d',ROI,arrayID);
        load(fullfile(save_path,load_name),'Preferred_loc','T_loc','p_loc')
        
        for n=1:size(Preferred_loc,2)
            
            All_Preferred_loc(count_ROI_pref)=Preferred_loc(n);
            All_T_loc(count_ROI_pref,:)=T_loc(:,n);
            All_p_loc(count_ROI_pref,:)=p_loc(:,n);
            count_ROI_pref=count_ROI_pref+1;
        end
        clear  Preferred_loc T_loc p_loc
    end
end



%% select based on pref

All_beta_split_all=vertcat(All_beta_split(:,:,:,:,:,:),All_beta_split(:,:,:,:,:,:),All_beta_split(:,:,:,:,:,:),All_beta_split(:,:,:,:,:,:)); %concatenate all 4 loc 
this_loc=vertcat(ones(size(All_beta_split,1),1).*1,ones(size(All_beta_split,1),1).*2,ones(size(All_beta_split,1),1).*3,ones(size(All_beta_split,1),1).*4);

select_pref=reshape(All_p_loc,size(All_p_loc,1)*size(All_p_loc,2),1);
select_pref=select_pref<0.05; %selct significant locations
select_unpref=~select_pref;

All_beta_split_pref=All_beta_split_all(select_pref,:,:,:,:,:);
All_beta_split_unpref=All_beta_split_all(select_unpref,:,:,:,:,:);

this_loc_pref=this_loc(select_pref);
this_loc_unpref=this_loc(select_unpref);


%% Bootstrap for stability

for n=1:N_boot         
    ind_p=randsample(size(All_beta_split_pref,1),size(All_beta_split_pref,1),'true');
    ind_u=randsample(size(All_beta_split_unpref,1),size(All_beta_split_unpref,1),'true');
    ind_u=ind_u(1:length(ind_p)); %subsample for equal N
    
    for i=1:size(All_beta_split,5) %rep
        for t=1:size(All_beta_split,6) %time
            
            %pref
            r_chosen_pref_boot(t,i,n)=calculate_chosen_corr_boot(All_beta_split_pref,this_loc_pref,1,2,i,t,ind_p);            
            r_unchosen_pref_boot(t,i,n)=calculate_unchosen_corr_boot(All_beta_split_pref,this_loc_pref,1,2,i,t,ind_p);
            
            %unpref
            r_chosen_unpref_boot(t,i,n)=calculate_chosen_corr_boot(All_beta_split_unpref,this_loc_unpref,1,2,i,t,ind_u);
            r_unchosen_unpref_boot(t,i,n)=calculate_unchosen_corr_boot(All_beta_split_unpref,this_loc_unpref,1,2,i,t,ind_u);

        end
    end
    
end


%%

for t=1:length(window_start_list)
    y=reshape(mean(r_chosen_pref_boot(t,:,:),2),size(r_chosen_pref_boot,3),1);
    p_chosen_pref_boot(t) = z_test_function_bootstrap(y,0);
    y=reshape(mean(r_unchosen_pref_boot(t,:,:),2),size(r_chosen_pref_boot,3),1);
    p_unchosen_ref_boot(t) = z_test_function_bootstrap(y,0);
    
    
    y=reshape(mean(r_chosen_unpref_boot(t,:,:),2),size(r_chosen_pref_boot,3),1);
    p_chosen_unpref_boot(t) = z_test_function_bootstrap(y,0);
    y=reshape(mean(r_unchosen_unpref_boot(t,:,:),2),size(r_chosen_pref_boot,3),1);
    p_unchosen_unpref_boot(t) = z_test_function_bootstrap(y,0);
end



%% save

save(fullfile(save_path, save_name))

end



%% functions


function r=calculate_chosen_corr_boot(beta,loc,ind1,ind2,this_rep,this_time,ind_boot)

for i=1:size(ind_boot,1)
    x(i)=beta(ind_boot(i),3,loc(ind_boot(i)),ind1,this_rep,this_time);
    y(i)=beta(ind_boot(i),3,loc(ind_boot(i)),ind2,this_rep,this_time);
end
if sum(~isnan(x) & ~isnan(y))>2
    
    x2=x(~isnan(x) & ~isnan(y));
    y2=y(~isnan(x) & ~isnan(y));
    
    r=corr(x2',y2');
else
    r=NaN;
end
end

function r=calculate_unchosen_corr_boot(beta,loc,ind1,ind2,this_rep,this_time,ind_boot)

for i=1:size(ind_boot,1)
    x(i)=beta(ind_boot(i),4,loc(ind_boot(i)),ind1,this_rep,this_time);
    y(i)=beta(ind_boot(i),4,loc(ind_boot(i)),ind2,this_rep,this_time);
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


