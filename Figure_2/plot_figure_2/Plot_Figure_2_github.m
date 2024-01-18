%Plot Figure 2

load('Figure_2.mat')
load('Figure_2_PCA.mat')
load('Figure_2_decoding_accuracy.mat')

addpath(genpath('Violinplot-Matlab-master'))

%% Figure 2C

figure
subplot(3,1,1)
hold on
for j=1:3
    [X, Y, Z]=ellipsoid(Mean_LIP_proj(j,1),Mean_LIP_proj(j,3),Mean_LIP_proj(j,2),Std_LIP_proj(j,1),Std_LIP_proj(j,3),Std_LIP_proj(j,2));%
    line([0.5 Mean_LIP_proj(j,1)], [0.5, Mean_LIP_proj(j,3)], [0.5 Mean_LIP_proj(j,2)],'Color',color_for_ROI(1,:),'LineWidth',2)
    plot3(Mean_LIP_proj(j,1),Mean_LIP_proj(j,3),Mean_LIP_proj(j,2),'o','color',color_template(j,:),'MarkerFacecolor',color_template(j,:),'MarkerSize',10)
    surf(X,Y,Z,'Facecolor',color_template(j,:),'FaceAlpha',0.35,'EdgeColor','none')
end
xlabel('Pink template')
ylabel('Blue template')
zlabel('Brown template')
ylim([0 1])
xlim([0 1])
zlim([0 1])
grid on
xticks(0:0.2:1)
yticks(0:0.2:1)
zticks(0:0.2:1)
view([-0.5 3 1])

subplot(3,1,2)
hold on
for j=1:3
    [X, Y, Z]=ellipsoid(Mean_FEF_proj(j,1),Mean_FEF_proj(j,3),Mean_FEF_proj(j,2),Std_FEF_proj(j,1),Std_FEF_proj(j,3),Std_FEF_proj(j,2));%
    line([0.5 Mean_FEF_proj(j,1)], [0.5, Mean_FEF_proj(j,3)], [0.5 Mean_FEF_proj(j,2)],'Color',color_for_ROI(2,:),'LineWidth',2)
    plot3(Mean_FEF_proj(j,1),Mean_FEF_proj(j,3),Mean_FEF_proj(j,2),'o','color',color_template(j,:),'MarkerFacecolor',color_template(j,:),'MarkerSize',10)
    surf(X,Y,Z,'Facecolor',color_template(j,:),'FaceAlpha',0.35,'EdgeColor','none')
end
xlabel('Pink template')
ylabel('Blue template')
zlabel('Brown template')
ylim([0 1])
xlim([0 1])
zlim([0 1])
grid on
xticks(0:0.2:1)
yticks(0:0.2:1)
zticks(0:0.2:1)
view([-0.5 3 1])

subplot(3,1,3)
hold on
for j=1:3
    [X, Y, Z]=ellipsoid(Mean_PFC_proj(j,1),Mean_PFC_proj(j,3),Mean_PFC_proj(j,2),Std_PFC_proj(j,1),Std_PFC_proj(j,3),Std_PFC_proj(j,2));%
    line([0.5 Mean_PFC_proj(j,1)], [0.5, Mean_PFC_proj(j,3)], [0.5 Mean_PFC_proj(j,2)],'Color',color_for_ROI(3,:),'LineWidth',2)
    plot3(Mean_PFC_proj(j,1),Mean_PFC_proj(j,3),Mean_PFC_proj(j,2),'o','color',color_template(j,:),'MarkerFacecolor',color_template(j,:),'MarkerSize',10)
    surf(X,Y,Z,'Facecolor',color_template(j,:),'FaceAlpha',0.35,'EdgeColor','none')
end
xticks(0:0.2:1)
yticks(0:0.2:1)
zticks(0:0.2:1)

xlabel('Pink template')
ylabel('Blue template')
zlabel('Brown template')
ylim([0 1])
xlim([0 1])
zlim([0 1])
grid on

view([-0.5 3 1])



%% Figure 2D

figure
subplot(3,1,1)
vs = violinplot(LIP_mean_classification_net',{'1/3','2/3','3/3'},'ViolinColor',color_for_ROI_2(1,:,2),'ShowData',false,'ViolinAlpha',0.5,'HalfViolin','left','BoxColor',[0 0 0],'MedianColor',[0.75 0.75 0.75]);
hold on
vs = violinplot(LIP_mean_classification_net_each_prog',{'1/3','2/3','3/3'},'ViolinColor',color_for_ROI_2(1,:,1),'ShowData',false,'ViolinAlpha',0.5,'HalfViolin','right','BoxColor',[0.5 0.5 0.5],'MedianColor',[1 1 1]);

yline(1/3,'--',{'Chance level'},'Color',[0.5 0.5 0.5],'LabelVerticalAlignment','middle','LabelHorizontalAlignment','center')
box off
% title('LIP','FontSize',12)
ylabel(sprintf('Expected template\nclassification accuracy'),'FontSize',12)
xlabel('Progression in block','FontSize',12)
ylim([0 1])

subplot(3,1,2)
vs = violinplot(FEF_mean_classification_net',{'1/3','2/3','3/3'},'ViolinColor',color_for_ROI_2(2,:,2),'ShowData',false,'ViolinAlpha',0.5,'HalfViolin','left','BoxColor',[0 0 0],'MedianColor',[0.75 0.75 0.75]);
hold on
vs = violinplot(FEF_mean_classification_net_each_prog',{'1/3','2/3','3/3'},'ViolinColor',color_for_ROI_2(2,:,1),'ShowData',false,'ViolinAlpha',0.5,'HalfViolin','right','BoxColor',[0.5 0.5 0.5],'MedianColor',[1 1 1]);

yline(1/3,'--',{'Chance level'},'Color',[0.5 0.5 0.5],'LabelVerticalAlignment','middle','LabelHorizontalAlignment','center')
box off
% title('FEF','FontSize',12)
ylabel(sprintf('Expected template\nclassification accuracy'),'FontSize',12)
xlabel('Progression in block','FontSize',12)
ylim([0 1])

subplot(3,1,3)
vs = violinplot(PFC_mean_classification_net',{'1/3','2/3','3/3'},'ViolinColor',color_for_ROI_2(3,:,2),'ShowData',false,'ViolinAlpha',0.5,'HalfViolin','left','BoxColor',[0 0 0],'MedianColor',[0.75 0.75 0.75]);
hold on
vs = violinplot(PFC_mean_classification_net_each_prog',{'1/3','2/3','3/3'},'ViolinColor',color_for_ROI_2(3,:,1),'ShowData',false,'ViolinAlpha',0.5,'HalfViolin','right','BoxColor',[0.5 0.5 0.5],'MedianColor',[1 1 1]);

yline(1/3,'--',{'Chance level'},'Color',[0.5 0.5 0.5],'LabelVerticalAlignment','middle','LabelHorizontalAlignment','center')
box off
% title('PFC','FontSize',12)
ylabel(sprintf('Expected template\nclassification accuracy'),'FontSize',12)
xlabel('Progression in block','FontSize',12)
ylim([0 1])

%% Figure 2E - PCA
color_bin(1,:)=colors(1,:);
color_bin(2,:)=colors(14,:);
color_bin(3,:)=colors(27,:);

figure
%Figure S2D
% PC_1=1;
% PC_2=2;
% PC_3=3;

PC_1=2;
PC_2=3;
PC_3=4;


subplot(3,1,1)
hold on
for n_c=1:3
    
    x(:)=PC_LIP(n_c,:,PC_1);
    y(:)=PC_LIP(n_c,:,PC_2);
    z(:)=PC_LIP(n_c,:,PC_3);
    if n_c==1
        plot3(x,y,z,'-','MarkerFaceColor',color_bin(n_c,:),'MarkerEdgeColor',color_bin(n_c,:),'Color',color_bin(n_c,:),'LineWidth',2)
    elseif n_c==2
        plot3(x,y,z,'-','MarkerFaceColor',color_bin(n_c,:),'MarkerEdgeColor',color_bin(n_c,:),'Color',color_bin(n_c,:),'LineWidth',2)
    else
        plot3(x,y,z,'-','MarkerFaceColor',color_bin(n_c,:),'MarkerEdgeColor',color_bin(n_c,:),'Color',color_bin(n_c,:),'LineWidth',2)
    end
    plot3(PC_LIP(n_c,target_on,PC_1),PC_LIP(n_c,target_on,PC_2),PC_LIP(n_c,target_on,PC_3),'^','MarkerSize',7,'MarkerFaceColor',[0.75 0.75 0.75],'MarkerEdgeColor',[0.75 0.75 0.75])
    hold on
    plot3(PC_LIP(n_c,response_on,PC_1),PC_LIP(n_c,response_on,PC_2),PC_LIP(n_c,response_on,PC_3),'^','MarkerSize',7,'MarkerFaceColor',[0.25 0.75 0.25],'MarkerEdgeColor',[0.25 0.75 0.25])
    grid on
    clear x y z
    
end
title(sprintf('LIP '))
xlabel(sprintf('PC %d',PC_1))
ylabel(sprintf('PC %d',PC_2))
zlabel(sprintf('PC %d',PC_3))
% view([-1 -1 1])
view([-1 -1 1])

subplot(3,1,2)
hold on
for n_c=1:3
    x(:)=PC_FEF(n_c,:,PC_1);
    y(:)=PC_FEF(n_c,:,PC_2);
    z(:)=PC_FEF(n_c,:,PC_3);
    if n_c==1
        plot3(x,y,z,'-','MarkerFaceColor',color_bin(n_c,:),'MarkerEdgeColor',color_bin(n_c,:),'Color',color_bin(n_c,:),'LineWidth',2)
    elseif n_c==2
        plot3(x,y,z,'-','MarkerFaceColor',color_bin(n_c,:),'MarkerEdgeColor',color_bin(n_c,:),'Color',color_bin(n_c,:),'LineWidth',2)
    else
        plot3(x,y,z,'-','MarkerFaceColor',color_bin(n_c,:),'MarkerEdgeColor',color_bin(n_c,:),'Color',color_bin(n_c,:),'LineWidth',2)
    end
    plot3(PC_FEF(n_c,target_on,PC_1),PC_FEF(n_c,target_on,PC_2),PC_FEF(n_c,target_on,PC_3),'^','MarkerSize',7,'MarkerFaceColor',[0.75 0.75 0.75],'MarkerEdgeColor',[0.75 0.75 0.75])
    hold on
    plot3(PC_FEF(n_c,response_on,PC_1),PC_FEF(n_c,response_on,PC_2),PC_FEF(n_c,response_on,PC_3),'^','MarkerSize',7,'MarkerFaceColor',[0.25 0.75 0.25],'MarkerEdgeColor',[0.25 0.75 0.25])
    
    clear x y z
    grid on
end
title(sprintf('FEF '))
xlabel(sprintf('PC %d',PC_1))
ylabel(sprintf('PC %d',PC_2))
zlabel(sprintf('PC %d',PC_3))
view([-1 1 1])



subplot(3,1,3)
hold on

for n_c=1:3
    x(:)=PC_PFC(n_c,:,PC_1);
    y(:)=PC_PFC(n_c,:,PC_2);
    z(:)=PC_PFC(n_c,:,PC_3);
    if n_c==1
        plot3(x,y,z,'-','MarkerFaceColor',color_bin(n_c,:),'MarkerEdgeColor',color_bin(n_c,:),'Color',color_bin(n_c,:),'LineWidth',2)
    elseif n_c==2
        plot3(x,y,z,'-','MarkerFaceColor',color_bin(n_c,:),'MarkerEdgeColor',color_bin(n_c,:),'Color',color_bin(n_c,:),'LineWidth',2)
    else
        plot3(x,y,z,'-','MarkerFaceColor',color_bin(n_c,:),'MarkerEdgeColor',color_bin(n_c,:),'Color',color_bin(n_c,:),'LineWidth',2)
    end
    plot3(PC_PFC(n_c,target_on,PC_1),PC_PFC(n_c,target_on,PC_2),PC_PFC(n_c,target_on,PC_3),'^','MarkerSize',7,'MarkerFaceColor',[0.75 0.75 0.75],'MarkerEdgeColor',[0.75 0.75 0.75])
    hold on
    plot3(PC_PFC(n_c,response_on,PC_1),PC_PFC(n_c,response_on,PC_2),PC_PFC(n_c,response_on,PC_3),'^','MarkerSize',7,'MarkerFaceColor',[0.25 0.75 0.25],'MarkerEdgeColor',[0.25 0.75 0.25])
    
    clear x y z
    
end
title(sprintf('PFC '))
xlabel(sprintf('PC %d',PC_1))
ylabel(sprintf('PC %d',PC_2))
zlabel(sprintf('PC %d',PC_3))
% view([-1 1 1])
view([-1 1 1])

grid on

%% Figure 2E - decoding accuracy

figure
subplot(3,1,1)
hold on
shadedErrorBar(window_start_list(8:end),mean(LIP_across_time_diag(8:end,:),2)',[prctile(LIP_across_time_diag(8:end,:),95,2),prctile(LIP_across_time_diag(8:end,:),5,2)]' ,{'color',color_for_ROI_early(1,:),'LineWidth',2},2)
shadedErrorBar(window_start_list(8:end),mean(LIP_across_time_dyn(8:end,:),2)',[prctile(LIP_across_time_dyn(8:end,:),95,2),prctile(LIP_across_time_dyn(8:end,:),5,2)]' ,{'color',color_for_ROI_late(1,:),'LineWidth',2},2)
plot_significance_level(window_start_list(8:end),P_LIP_diag(8:end),0.8,color_for_ROI_early(1,:),[0.01,0.05/length(window_start_list), 0.01/length(window_start_list)])
plot_significance_level(window_start_list(8:end),P_LIP_dyn(8:end),0.85,color_for_ROI_late(1,:),[0.01,0.05/length(window_start_list), 0.01/length(window_start_list)])
yline(1/3,'--')
ylabel({'Decoding accuracy'})
box off
xlabel({'Time to targets onset'})

subplot(3,1,2)
hold on
shadedErrorBar(window_start_list(8:end),mean(FEF_across_time_diag(8:end,:),2)',[prctile(FEF_across_time_diag(8:end,:),95,2),prctile(FEF_across_time_diag(8:end,:),5,2)]' ,{'color',color_for_ROI_early(2,:),'LineWidth',2},2)
shadedErrorBar(window_start_list(8:end),mean(FEF_across_time_dyn(8:end,:),2)',[prctile(FEF_across_time_dyn(8:end,:),95,2),prctile(FEF_across_time_dyn(8:end,:),5,2)]' ,{'color',color_for_ROI_late(2,:),'LineWidth',2},2)
plot_significance_level(window_start_list(8:end),P_FEF_diag(8:end),0.8,color_for_ROI_early(2,:),[0.01,0.05/length(window_start_list), 0.01/length(window_start_list)])
plot_significance_level(window_start_list(8:end),P_FEF_dyn(8:end),0.85,color_for_ROI_late(2,:),[0.01,0.05/length(window_start_list), 0.01/length(window_start_list)])
yline(1/3,'--')
ylabel({'Decoding accuracy'})
box off
xlabel({'Time to targets onset'})

subplot(3,1,3)
hold on
shadedErrorBar(window_start_list(8:end),mean(PFC_across_time_diag(8:end,:),2)',[prctile(PFC_across_time_diag(8:end,:),95,2),prctile(PFC_across_time_diag(8:end,:),5,2)]' ,{'color',color_for_ROI_early(3,:),'LineWidth',2},2)
shadedErrorBar(window_start_list(8:end),mean(PFC_across_time_dyn(8:end,:),2)',[prctile(PFC_across_time_dyn(8:end,:),95,2),prctile(PFC_across_time_dyn(8:end,:),5,2)]' ,{'color',color_for_ROI_late(3,:),'LineWidth',2},2)
plot_significance_level(window_start_list(8:end),P_PFC_diag(8:end),0.8,color_for_ROI_early(3,:),[0.01,0.05/length(window_start_list), 0.01/length(window_start_list)])
plot_significance_level(window_start_list(8:end),P_PFC_dyn(8:end),0.85,color_for_ROI_late(3,:),[0.01,0.05/length(window_start_list), 0.01/length(window_start_list)])
yline(1/3,'--')
ylabel({'Decoding accuracy'})
box off
xlabel({'Time to targets onset'})

%% Functions


function plot_significance_level(x,p,a,c,thr)
hold on
for b=1:length(thr)
    this_thr=thr(b);
    for i=1:length(p)-1
        if p(i)<=this_thr && p(i+1)<=this_thr
            plot(x(i):x(i+1),a*ones(1,x(i+1)-x(i)+1),'-','Color',c,'LineWidth',b)
        end
        if i>1
            if p(i-1)>=this_thr && p(i)<=this_thr && p(i+1)>=this_thr
                plot(x(i),a,'.','Color',c,'LineWidth',b)
            end
        end
    end
end
end


