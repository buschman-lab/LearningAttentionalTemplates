
clear all

fsroot='/Volumes/buschman';
task='Learning_attentional_templates';

%set monkey
monkey='Beaker';

N_channels=6;
model_type='RW';
switch_type='Reset';
N_bins=100;

% load data
load_path = sprintf('/Volumes/buschman/Projects/Learning_Attentional_Templates/Analysed_data/%s',monkey);
load_name=sprintf('All_sessions_%s_%s_%d_channels_VBMC',switch_type,model_type,N_channels);
load(fullfile(load_path,load_name));
data=load_monkey_data_continuous(fsroot,monkey,N_channels,100);

load('colors')
color_current=colors(17,:); %green
color_old=colors(1,:); %pink

switch monkey
    case 'Beaker'
        cfp=[0.5 0.5 0.5];
    case 'Scooter'
        cfp=[0.5 0.5 0.5];
end

cd('/Volumes/buschman/Users/Caroline/Learning_attentional_templates/Behavioral_Analysis/Figure1')

%% Look for difficulty effects

positive_values=Model_predictions.model_outputs.Option_value-min(Model_predictions.model_outputs.Option_value,[],'all');

for i=1:length(data.Chosen_ind)
    V_normalized(i)=positive_values(1,i)/sum(positive_values(:,i));
end
V_normalized(sum(positive_values,2)<0.1)=NaN;


figure
histogram(V_normalized)

%% Look for entropy effects

for i=1:length(data.Chosen_ind)
    value_this_trial=Model_predictions.model_outputs.Value_for_choice(:,i);
    if sum(value_this_trial==0)<size(value_this_trial,1)
        proba_value_this_val(:,i)=(value_this_trial-min(value_this_trial)+0.01)./sum(value_this_trial-min(value_this_trial)+0.01);
    else
        proba_value_this_val(1:size(value_this_trial,1),i)=NaN;
    end
    entropy(i)=-sum(proba_value_this_val(:,i).*log(proba_value_this_val(:,i)));
    
end

figure
histogram(entropy)

%% AT location

for i=1:length(data.Chosen_ind)
    value_this_trial=Model_predictions.model_outputs.Value_for_choice(:,i);
    if sum(value_this_trial==0)<size(value_this_trial,1)
        [~, expected_AT(i)]=max(value_this_trial);
    else
        expected_AT(i)=NaN;
    end
end

expected_AT=expected_AT/100*2*pi;
% 
% figure
% histogram(expected_AT)


%% Plot

clear std* accuracy*

bin_V=prctile(V_normalized,0:20:100);
for j=1:length(bin_V)-1
    accuracy_binned_choice_difficulty(j)=mean(data.Accuracy(V_normalized>=bin_V(j) & V_normalized<=bin_V(j+1)));
    std_accuracy_binned_choice_difficulty(j)=std(data.Accuracy(V_normalized>=bin_V(j) & V_normalized<=bin_V(j+1)))./sqrt(sum(V_normalized>=bin_V(j) & V_normalized<=bin_V(j+1))-1);
end

figure
subplot(2,3,1)
shadedErrorBar(bin_V(1:end-1),accuracy_binned_choice_difficulty,std_accuracy_binned_choice_difficulty,{'color',cfp,'LineWidth',2},2)
box off
ylabel(sprintf('Accuracy'),'FontSize',20)
xlabel('Vbest / sum all values','FontSize',20)
ax=gca;

bin_E=prctile(entropy,0:20:100);
for j=1:length(bin_E)-1
    accuracy_binned_entropy(j)=mean(data.Accuracy(entropy>=bin_E(j) & entropy<=bin_E(j+1)));
    std_accuracy_binned_entropy(j)=std(data.Accuracy(entropy>=bin_E(j) & entropy<=bin_E(j+1)))./sqrt(sum(entropy>=bin_E(j) & entropy<=bin_E(j+1))-1);
end

subplot(2,3,2)
shadedErrorBar(bin_E(1:end-1),accuracy_binned_entropy,std_accuracy_binned_entropy,{'color',cfp,'LineWidth',2},2)
box off
ylabel(sprintf('Accuracy'),'FontSize',20)
xlabel('Entropy','FontSize',20)
ax=gca;

bin_AT=0:pi/2:2*pi;
for j=1:length(bin_AT)-1
    accuracy_binned_expected_AT(j)=mean(data.Accuracy(expected_AT>=bin_AT(j) & expected_AT<=bin_AT(j+1)));
    std_accuracy_binned_expected_AT(j)=std(data.Accuracy(expected_AT>=bin_AT(j) & expected_AT<=bin_AT(j+1)))./sqrt(sum(expected_AT>=bin_AT(j) & expected_AT<=bin_AT(j+1))-1);
    sum(expected_AT>=bin_AT(j) & expected_AT<=bin_AT(j+1))
end

subplot(2,3,3)
polarplot([bin_AT(1:end-1) bin_AT(1)],[accuracy_binned_expected_AT accuracy_binned_expected_AT(1)],'Color',cfp,'LineWidth',2)
hold on
polarplot([bin_AT(1:end-1) bin_AT(1)],[accuracy_binned_expected_AT+std_accuracy_binned_expected_AT accuracy_binned_expected_AT(1)+std_accuracy_binned_expected_AT(1)],'Color',cfp,'LineWidth',1)
polarplot([bin_AT(1:end-1) bin_AT(1)],[accuracy_binned_expected_AT-std_accuracy_binned_expected_AT accuracy_binned_expected_AT(1)-std_accuracy_binned_expected_AT(1)],'Color',cfp,'LineWidth',1)

box off
pax = gca;
pax.ThetaAxisUnits = 'radians';
title('Accuracy','FontSize',20,'FontWeight','Normal')

reset_ind=find(Model_predictions.model_outputs.Switch==1);
reset_ind=reset_ind+1; %look at trial after the reset

for j=1:length(reset_ind)
    Accuracy_from_reset(:,j)=data.Accuracy(reset_ind(j):reset_ind(j)+50);
    Entropy_from_reset(:,j)=entropy(reset_ind(j):reset_ind(j)+50);
end

subplot(2,3,4)
shadedErrorBar([],nanmean(Accuracy_from_reset,2),nanstd(Accuracy_from_reset,0,2)./sqrt(length(reset_ind)-1),{'color',cfp,'LineWidth',2},2)
box off
ylabel(sprintf('Accuracy'),'FontSize',20)
xlabel('Trial from reset','FontSize',20)
ax=gca;

subplot(2,3,5)
shadedErrorBar([],nanmean(Entropy_from_reset,2),nanstd(Entropy_from_reset,0,2)./sqrt(length(reset_ind)-1),{'color',cfp,'LineWidth',2},2)
box off
ylabel(sprintf('Entropy'),'FontSize',20)
xlabel('Trial from reset','FontSize',20)
ax=gca;

for j=1:length(bin_AT)-1
    [beta_AT(:,j),~,s]=glmfit(zscore(V_normalized(expected_AT>=bin_AT(j) & expected_AT<=bin_AT(j+1))),data.Accuracy(expected_AT>=bin_AT(j) & expected_AT<=bin_AT(j+1))','binomial');
    se_beta_AT(:,j)=s.se;
end

subplot(2,3,6)
polarplot([bin_AT(1:end-1) bin_AT(1)],[beta_AT(2,:) beta_AT(2,1)],'Color',cfp,'LineWidth',2)
hold on
polarplot([bin_AT(1:end-1) bin_AT(1)],[beta_AT(2,:)+se_beta_AT(2,:) beta_AT(2,1)+se_beta_AT(2,1)],'Color',cfp,'LineWidth',1)
polarplot([bin_AT(1:end-1) bin_AT(1)],[beta_AT(2,:)-se_beta_AT(2,:) beta_AT(2,1)-se_beta_AT(2,1)],'Color',cfp,'LineWidth',1)

box off
pax = gca;
pax.ThetaAxisUnits = 'radians';
title('Vbest / sum all values','FontSize',20,'FontWeight','Normal')

%% Stats
m=fitglm(zscore(positive_values,[],2)',data.Accuracy,'Distribution','binomial')

m1=fitglm(V_normalized,data.Accuracy,'Distribution','binomial')

m3=fitglm(entropy,data.Accuracy,'Distribution','binomial')

distances_to_template=1-abs(mod(expected_AT-data.options_color+pi,2*pi)-pi)/pi;

m4=fitglm(nanzscore([distances_to_template(1:2,:); entropy]',1),data.Accuracy,'interactions','Distribution','binomial')

%%

m_model=fitglm(zscore(positive_values,[],2)',zscore(Model_predictions.model_outputs.Choice_prediction(1,:)))

m1_model=fitglm(V_normalized,zscore(Model_predictions.model_outputs.Choice_prediction(1,:)))

m3_model=fitglm(entropy,zscore(Model_predictions.model_outputs.Choice_prediction(1,:)))

m4_model=fitglm(nanzscore([distances_to_template(1:2,:); entropy]',1),zscore(Model_predictions.model_outputs.Choice_prediction(1,:)),'interactions')




%%
figure
subplot(1,2,1)
hold on
bar(m.Coefficients.Estimate(2:end),'FaceColor',[0.5 0.5 0.5]);
errorbar(m.Coefficients.Estimate(2:end), m.Coefficients.SE(2:end),'.k','LineWidth',2);
title('data')

subplot(1,2,2)
hold on
bar(m_model.Coefficients.Estimate(2:end),'FaceColor',[0.5 0.5 0.5]);
errorbar(m_model.Coefficients.Estimate(2:end), m_model.Coefficients.SE(2:end),'.k','LineWidth',2);
title('model')

figure
subplot(1,2,1)
hold on
bar(m4.Coefficients.Estimate(2:end),'FaceColor',[0.5 0.5 0.5]);
errorbar(m4.Coefficients.Estimate(2:end), m4.Coefficients.SE(2:end),'.k','LineWidth',2);
title('data')

subplot(1,2,2)
hold on
bar(m4_model.Coefficients.Estimate(2:end),'FaceColor',[0.5 0.5 0.5]);
errorbar(m4_model.Coefficients.Estimate(2:end), m4_model.Coefficients.SE(2:end),'.k','LineWidth',2);
title('model')


%% how often does the model make the riht prediciton?
    y=[];
    y_model=[];
    y_model_binary=[];
    
    for j=1:3
        for i=1:size(data.choices,2)
            y(end+1)=data.choices(j,i);
            y_model(end+1)=Model_predictions.model_outputs.Choice_prediction(j,i);
            if Model_predictions.model_outputs.Choice_prediction(j,i)==max(Model_predictions.model_outputs.Choice_prediction(:,i))
                y_model_binary(end+1)=1;
            else
                y_model_binary(end+1)=0;
            end
        end
    end
correct_prediction=sum(y_model_binary==1 & y==1)/sum(y==1);    
    




%%
function z = nanzscore(x,dim)

z=(x-nanmean(x,dim))./nanstd(x,0,dim);

end
