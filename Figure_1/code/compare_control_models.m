%Plot the results of the control models

clear all;
monkey='Scooter';

fsroot='/Volumes/buschman';
task='Learning_Attentional_Templates';
save_path=fullfile(fsroot,'Projects', task, 'Analysed_data',sprintf('%s',monkey));

%% compare the models
load_name=sprintf('Reset_%s_surprise_RW_%d_channels_VBMC',monkey,6); %load the best model
load(fullfile(save_path,load_name),'Reset')
load_name=sprintf('%s_WSLS_flip_VBMC',monkey);
load(fullfile(save_path,load_name),'WSLS')
WSLS_flip=WSLS; %rename
clear WSLS
load_name=sprintf('%s_WSLS_VBMC',monkey);
load(fullfile(save_path,load_name),'WSLS')
load_name=sprintf('%s_no_value_VBMC',monkey);
load(fullfile(save_path,load_name),'Control_model')
Control_no_value=Control_model;
clear Control_model
load_name=sprintf('%s_Reset_value_only_VBMC',monkey);
load(fullfile(save_path,load_name),'Control_model')
Control_value_only=Control_model;
clear Control_model
load_name=sprintf('%s_Reset_no_bias_prev_cc_VBMC',monkey);
load(fullfile(save_path,load_name),'Control_model')
Control_no_bias_prev_cc=Control_model;
clear Control_model

block_nb=size(Reset.Prediction_block,2);
N_plot=size(Reset.Prediction_block,1);
%%

figure
subplot(1,2,1)
hold on
bar(1,WSLS_flip.out.fit.BIC-WSLS.out.fit.BIC,'FaceColor',[0.75 0 0])
bar(2,Control_no_value.out.fit.BIC-WSLS.out.fit.BIC,'FaceColor',[0 0 0.75])
bar(3,Control_value_only.out.fit.BIC-WSLS.out.fit.BIC,'FaceColor',[0 0.25 0])
bar(4,Control_no_bias_prev_cc.out.fit.BIC-WSLS.out.fit.BIC,'FaceColor',[0 0.5 0])
bar(5,Reset.out.fit.BIC-WSLS.out.fit.BIC,'FaceColor',[0 0.75 0])

set(gca,'XTick',(1:1:5));
set(gca,'XTickLabel',{'WSLS flip','No value','No biases','No prev color bias','Reset'})

title('BIC')

subplot(1,2,2)
hold on
l=1;    %Best option only
plot(mean(Reset.Choice_block(:,:,l),2),'Color','k','LineWidth',2) %data

plot(mean(WSLS.Prediction_block(:,:,l),2),'Color',[0.5 0 0],'LineWidth',2)
plot(mean(WSLS_flip.Prediction_block(:,:,l),2),'Color',[0.75 0 0],'LineWidth',2)
plot(mean(Control_no_value.Prediction_block(:,:,l),2),'Color',[0 0 0.75],'LineWidth',2)
plot(mean(Control_value_only.Prediction_block(:,:,l),2),'Color',[0 0.25 0],'LineWidth',2)
plot(mean(Control_no_bias_prev_cc.Prediction_block(:,:,l),2),'Color',[0 0.5 0],'LineWidth',2)

plot(mean(Reset.Prediction_block(:,:,l),2),'Color',[0 0.75 0],'LineWidth',2) %best model

shadedErrorBar([],mean(Reset.Choice_block(:,:,l),2),std(Reset.Choice_block(:,:,l),0,2)./sqrt(block_nb-1),{'color','k'},2) %data

shadedErrorBar([],mean(WSLS.Prediction_block(:,:,l),2),std(WSLS.Prediction_block(:,:,l),0,2)./sqrt(block_nb-1),{'color',[0.5 0 0]},2)
shadedErrorBar([],mean(WSLS_flip.Prediction_block(:,:,l),2),std(WSLS_flip.Prediction_block(:,:,l),0,2)./sqrt(block_nb-1),{'color',[0.75 0 0]},2)
shadedErrorBar([],mean(Control_no_value.Prediction_block(:,:,l),2),std(Control_no_value.Prediction_block(:,:,l),0,2)./sqrt(block_nb-1),{'color',[0 0 0.75]},2)
shadedErrorBar([],mean(Control_value_only.Prediction_block(:,:,l),2),std(Control_value_only.Prediction_block(:,:,l),0,2)./sqrt(block_nb-1),{'color',[0 0.25 0]},2)
shadedErrorBar([],mean(Control_no_bias_prev_cc.Prediction_block(:,:,l),2),std(Control_no_bias_prev_cc.Prediction_block(:,:,l),0,2)./sqrt(block_nb-1),{'color',[0 0.5 0]},2)

shadedErrorBar([],mean(Reset.Prediction_block(:,:,l),2),std(Reset.Prediction_block(:,:,l),0,2)./sqrt(block_nb-1),{'color',[0 0.75 0]},2) %best model

yline(0.33,'--')
ylim([0 1])
xlim([0 N_plot+1])
ylabel(sprintf('Proportion of trials in which\nthe best target was chosen'),'FontSize',20)
xlabel('Trial from template switch','FontSize',20)
box off
ax=gca;
ax.XAxis.FontSize=16;
ax.YAxis.FontSize=16;

legend({'data','WSLS','WSLS flip','No value','Value only','No bias prev color','Reset'},'Location','Southeast')


%%
figure
subplot(1,2,1)
hold on
bar(1,WSLS_flip.out.fit.BIC-WSLS.out.fit.BIC,'FaceColor',[0.75 0 0])
bar(2,Control_no_value.out.fit.BIC-WSLS.out.fit.BIC,'FaceColor',[0 0 0.75])
bar(3,Control_value_only.out.fit.BIC-WSLS.out.fit.BIC,'FaceColor',[0 0.25 0])
bar(4,Control_no_bias_prev_cc.out.fit.BIC-WSLS.out.fit.BIC,'FaceColor',[0 0.5 0])
bar(5,Reset.out.fit.BIC-WSLS.out.fit.BIC,'FaceColor',[0 0.75 0])

set(gca,'XTick',(1:1:5));
set(gca,'XTickLabel',{'WSLS flip','No value','No biases','No prev color bias','Reset'})

title('BIC')

subplot(1,2,2)
hold on
l=1;    %Best option only
plot(mean(Reset.Choice_block(:,:,l),2),'Color','k','LineWidth',2) %data

plot(mean(WSLS.Prediction_block(:,:,l),2),'Color',[0.5 0 0],'LineWidth',2)
plot(mean(WSLS_flip.Prediction_block(:,:,l),2),'Color',[0.75 0 0],'LineWidth',2)
plot(mean(Control_no_value.Prediction_block(:,:,l),2),'Color',[0 0 0.75],'LineWidth',2)
plot(mean(Control_value_only.Prediction_block(:,:,l),2),'Color',[0 0.25 0],'LineWidth',2)
plot(mean(Control_no_bias_prev_cc.Prediction_block(:,:,l),2),'Color',[0 0.5 0],'LineWidth',2)

plot(mean(Reset.Prediction_block(:,:,l),2),'Color',[0 0.75 0],'LineWidth',2) %best model

% shadedErrorBar([],mean(Reset.Choice_block(:,:,l),2),std(Reset.Choice_block(:,:,l),0,2)./sqrt(block_nb-1),{'color','k'},2) %data
% 
% shadedErrorBar([],mean(WSLS.Prediction_block(:,:,l),2),std(WSLS.Prediction_block(:,:,l),0,2)./sqrt(block_nb-1),{'color',[0.5 0 0]},2)
% shadedErrorBar([],mean(WSLS_flip.Prediction_block(:,:,l),2),std(WSLS_flip.Prediction_block(:,:,l),0,2)./sqrt(block_nb-1),{'color',[0.75 0 0]},2)
% shadedErrorBar([],mean(Control_no_value.Prediction_block(:,:,l),2),std(Control_no_value.Prediction_block(:,:,l),0,2)./sqrt(block_nb-1),{'color',[0 0 0.75]},2)
% shadedErrorBar([],mean(Control_value_only.Prediction_block(:,:,l),2),std(Control_value_only.Prediction_block(:,:,l),0,2)./sqrt(block_nb-1),{'color',[0 0.25 0]},2)
% shadedErrorBar([],mean(Control_no_bias_prev_cc.Prediction_block(:,:,l),2),std(Control_no_bias_prev_cc.Prediction_block(:,:,l),0,2)./sqrt(block_nb-1),{'color',[0 0.5 0]},2)
% 
% shadedErrorBar([],mean(Reset.Prediction_block(:,:,l),2),std(Reset.Prediction_block(:,:,l),0,2)./sqrt(block_nb-1),{'color',[0 0.75 0]},2) %best model

yline(0.33,'--')
ylim([0 1])
xlim([0 N_plot+1])
ylabel(sprintf('Proportion of trials in which\nthe best target was chosen'),'FontSize',20)
xlabel('Trial from template switch','FontSize',20)
box off
ax=gca;
ax.XAxis.FontSize=16;
ax.YAxis.FontSize=16;

legend({'data','WSLS','WSLS flip','No value','Value only','No bias prev color','Reset'},'Location','Southeast')


%% compare ateentiuonal templates

load(fullfile(save_path,'All_sessions_Reset_no_bias_prev_cc_RW_6_channels_VBMC.mat'));
Model_predictions_control=Model_predictions;
load(fullfile(save_path,'All_sessions_Reset_RW_6_channels_VBMC.mat'));

for i=1:size(Model_predictions_control.model_outputs.Value_for_choice,2)
[~, AT_control(i)]=max(Model_predictions_control.model_outputs.Value_for_choice(:,i));
end
for i=1:size(Model_predictions_control.model_outputs.Value_for_choice,2)
[~, AT(i)]=max(Model_predictions.model_outputs.Value_for_choice(:,i));
end

%in rad
AT=AT/100*2*pi;
AT_control=AT_control/100*2*pi;

%%
values=hist3([AT' AT_control'],[100 100]);
x=0:2*pi/100:2*pi-2*pi/100;

figure
hold on
imagesc(x,x,values./sum(values,2))
colormap(flipud(gray))
ylim([0 2*pi])
set(gca,'YTick',(0:pi/2:pi*2));
set(gca,'YTickLabel',{'0','π/2','π','3π/2'})
xlim([0 2*pi])
set(gca,'XTick',(0:pi/2:pi*2));
set(gca,'XTickLabel',{'0','π/2','π','3π/2'})
xlabel({'Attentional template'})
ylabel({'Control attentional template'})

ax=gca;
ax.XAxis.FontSize=16;
ax.YAxis.FontSize=16;


[r p]=circ_corrcc(AT,AT_control)
