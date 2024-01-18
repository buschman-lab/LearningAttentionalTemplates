%Plot Figure 1

%load data
load('Monkey_B')
load('Monkey_S')

%set colors
cfp=[0.5 0.5 0.5];
load('colors')
color_current=colors(17,:); %green
color_old=colors(1,:); %pink

color_current_data=rgb2hsv(color_current);
color_current_data=hsv2rgb(color_current_data(1),color_current_data(2),0.65);

color_old_data=rgb2hsv(color_old);
color_old_data=hsv2rgb(color_old_data(1),color_old_data(2),0.9);

color_current_model=rgb2hsv(color_current);
color_current_model=hsv2rgb(color_current_model(1),color_current_model(2),0.35);

color_old_model=rgb2hsv(color_old);
color_old_model=hsv2rgb(color_old_model(1),color_old_model(2),0.5);

N_plot=80;

%% Figure 1F

figure
subplot(2,1,1)
shadedErrorBar([],Monkey_B.Figure_1_F_mean,Monkey_B.Figure_1_F_sem,{'color',cfp,'LineWidth',2},2)
box off
ylabel(sprintf('Distance between template estimate \n and the true template (in rad)'),'FontSize',20)
xlabel('Progression in block','FontSize',20)
set(gca,'XTick',(0:10:30));
set(gca,'XTickLabel',{'0','1/3','2/3','1'})
ylim([0 pi/2+pi/8])
set(gca,'YTick',(0:pi/4:pi/2));
set(gca,'YTickLabel',{'0','π/4','π/2'})
ax=gca;
ax.XAxis.FontSize=16;
ax.YAxis.FontSize=16;

subplot(2,1,2)
shadedErrorBar([],Monkey_S.Figure_1_F_mean,Monkey_S.Figure_1_F_sem,{'color',cfp,'LineWidth',2},2)
box off
ylabel(sprintf('Distance between template estimate \n and the true template (in rad)'),'FontSize',20)
xlabel('Progression in block','FontSize',20)
set(gca,'XTick',(0:10:30));
set(gca,'XTickLabel',{'0','1/3','2/3','1'})
ylim([0 pi/2+pi/8])
set(gca,'YTick',(0:pi/4:pi/2));
set(gca,'YTickLabel',{'0','π/4','π/2'})
ax=gca;
ax.XAxis.FontSize=16;
ax.YAxis.FontSize=16;


%% Figure 1 E

figure
subplot(2,1,1)
plot(nanmean(Monkey_B.Figure_1_E_current_data,2),'Color',color_current_data,'Linewidth',2)
hold on
plot(nanmean(Monkey_B.Figure_1_E_previous_data,2),'Color',color_old_data,'Linewidth',2)
hold on

plot(nanmean(Monkey_B.Figure_1_E_current_model,2),'--','Linewidth',2,'Color',color_current_model)
hold on
plot(nanmean(Monkey_B.Figure_1_E_previous_model,2),'--','Linewidth',2,'Color',color_old_model)
hold on

shadedErrorBar([],nanmean(Monkey_B.Figure_1_E_current_data,2),nanstd(Monkey_B.Figure_1_E_current_data,0,2)./sqrt(size(Monkey_B.Figure_1_E_current_data,2)-1),{'color',color_current_data},0.5)
hold on
shadedErrorBar([],nanmean(Monkey_B.Figure_1_E_previous_data,2),nanstd(Monkey_B.Figure_1_E_previous_data,0,2)./sqrt(size(Monkey_B.Figure_1_E_previous_data,2)-1),{'color',color_old_data},0.5)
hold on
shadedErrorBar([],nanmean(Monkey_B.Figure_1_E_current_model,2),nanstd(Monkey_B.Figure_1_E_current_model,0,2)./sqrt(size(Monkey_B.Figure_1_E_current_model,2)-1),{'--','LineWidth',2,'color',color_current_model},0.5)
hold on
shadedErrorBar([],nanmean(Monkey_B.Figure_1_E_previous_model,2),nanstd(Monkey_B.Figure_1_E_previous_model,0,2)./sqrt(size(Monkey_B.Figure_1_E_previous_model,2)-1),{'--','LineWidth',2,'color',color_old_model},0.5)

yline(0.33,'--')
legend({'Data: current template','Data: previous template', 'Model: current template','Model: previous template'})
ylim([0 1])
xlim([0 N_plot+1])
ylabel(sprintf('Proportion of trials in which\nthe best target was chosen'),'FontSize',20)
xlabel('Trial from template switch','FontSize',20)
box off
ax=gca;
ax.XAxis.FontSize=16;
ax.YAxis.FontSize=16;

subplot(2,1,2)
plot(nanmean(Monkey_S.Figure_1_E_current_data,2),'Color',color_current_data,'Linewidth',2)
hold on
plot(nanmean(Monkey_S.Figure_1_E_previous_data,2),'Color',color_old_data,'Linewidth',2)
hold on

plot(nanmean(Monkey_S.Figure_1_E_current_model,2),'--','Linewidth',2,'Color',color_current_model)
hold on
plot(nanmean(Monkey_S.Figure_1_E_previous_model,2),'--','Linewidth',2,'Color',color_old_model)
hold on

shadedErrorBar([],nanmean(Monkey_S.Figure_1_E_current_data,2),nanstd(Monkey_S.Figure_1_E_current_data,0,2)./sqrt(size(Monkey_S.Figure_1_E_current_data,2)-1),{'color',color_current_data},0.5)
hold on
shadedErrorBar([],nanmean(Monkey_S.Figure_1_E_previous_data,2),nanstd(Monkey_S.Figure_1_E_previous_data,0,2)./sqrt(size(Monkey_S.Figure_1_E_previous_data,2)-1),{'color',color_old_data},0.5)
hold on
shadedErrorBar([],nanmean(Monkey_S.Figure_1_E_current_model,2),nanstd(Monkey_S.Figure_1_E_current_model,0,2)./sqrt(size(Monkey_S.Figure_1_E_current_model,2)-1),{'--','LineWidth',2,'color',color_current_model},0.5)
hold on
shadedErrorBar([],nanmean(Monkey_S.Figure_1_E_previous_model,2),nanstd(Monkey_S.Figure_1_E_previous_model,0,2)./sqrt(size(Monkey_S.Figure_1_E_previous_model,2)-1),{'--','LineWidth',2,'color',color_old_model},0.5)

yline(0.33,'--')
legend({'Data: current template','Data: previous template', 'Model: current template','Model: previous template'})
ylim([0 1])
xlim([0 N_plot+1])
ylabel(sprintf('Proportion of trials in which\nthe best target was chosen'),'FontSize',20)
xlabel('Trial from template switch','FontSize',20)
box off
ax=gca;
ax.XAxis.FontSize=16;
ax.YAxis.FontSize=16;
