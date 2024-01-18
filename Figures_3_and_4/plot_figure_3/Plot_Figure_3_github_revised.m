%Plot Figure 3

load('Figure_3.mat')

addpath(genpath('Violinplot-Matlab-master'))

%% Figure 3A

N_channels=6;
for i=1:N_channels
    color_bin(i,:)=colors(1+ceil((i-1)*(40/N_channels)),:);
end

PC_1=1;
PC_2=2;
PC_3=3;

color_bin2=rgb2hsv(color_bin);
color_bin2(:,2)=color_bin2(:,2)-0.2;
color_bin2=hsv2rgb(color_bin2);

figure
hold on
x=score_all(:,PC_1);
x(end+1)=score_all(1,PC_1);
y=score_all(:,PC_2);
y(end+1)=score_all(1,PC_2);
z=score_all(:,PC_3);
z(end+1)=score_all(1,PC_3);
plot(x,y,'-','Color','k','LineWidth',1)
for n_c=1:N_channels
    plot(score_all(n_c,PC_1),score_all(n_c,PC_2),'-o','MarkerFaceColor',color_bin(n_c,:),'MarkerEdgeColor',color_bin(n_c,:),'Color',color_bin(n_c,:),'MarkerSize',15)
end
grid on
lim(1)=min(min(x),min(y))-5;
lim(2)=max(max(x),max(y))+5;
xlim(lim)
ylim(lim)
clear x y z
title(sprintf('All regions MDS'))
xlabel(sprintf('dim %d',PC_1))
ylabel(sprintf('dim %d',PC_2))
grid on



%% Figure 3B

figure
subplot(2,1,1)
vs = violinplot(z_dist_val',{'All'},'ViolinColor',[0 0 0],'ShowData',true,'ShowMean',true,'ShowMedian',true);
if ttest(z_dist_val,0,'Tail','left')
    plot(1,2.5,'k*')
end
yline(0,'--')
box off
ylabel(sprintf('z-score(distance neural and\nbehavioral estimated templates'))

%% Figure 3D

mean_mean_decoded_angle=mod(angle(nanmean(exp(1i*mean_decoded_angle),2)),2*pi);

figure
subplot(1,2,1)
for i=1:length(mean_true_angle)
    if ~isnan(mean_true_angle(i))
        hold on
        plot(mean_true_angle(i),mean_decoded_angle(i),'o','Color',colors(round(mod(mean_true_angle(i),2*pi)/2/pi*40),:),'MarkerFaceColor',colors(round(mod(mean_true_angle(i),2*pi)/2/pi*40),:))
    end
end
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

%% Figure 3E

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

%figure S3C

% subplot(1,2,2)
% plot(dist_true_angles_on_wheel,dist_decoded_angles_on_wheel,'k.')
% hold on
% plot(-5:5,-5:5,'r-','Linewidth',2)
% xticks([-pi -pi/2 0 pi/2 pi])
% xticklabels({'-π','-π/2','0','π/2','π'})
% yticks([-pi -pi/2 0 pi/2 pi])
% yticklabels({'-π','-π/2','0','π/2','π'})
% ylim([-pi pi])
% xlim([-pi pi])
% xlabel('d(ET(i),ET(j)) behavioral')
% ylabel('d(ET(i),ET(j)) neural')
% box off

