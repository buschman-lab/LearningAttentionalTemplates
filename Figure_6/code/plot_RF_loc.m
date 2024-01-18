clear all

ROI='PFC';

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

%% Plot

figure
hold on
for i=1:4
    bar(i,sum(All_p_loc(:,i)<0.05)/sum(~isnan(All_p_loc(:,i))),'FaceColor',color_for_plot)
end
ylim([0 0.3])



