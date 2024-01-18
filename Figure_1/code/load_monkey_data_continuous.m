function [ inG ] = load_monkey_data_continuous(fsroot,monkey,N_channels,N_bins)

raw_data_path=fullfile(fsroot,'Projects/Learning_Attentional_Templates/Data/');

addpath(fullfile(fsroot,'Users/Caroline/Learning_attentional_templates/Behavioral_Analysis'));


switch monkey
    case 'Beaker'
        date_list={'181019','181025','181029','181030','181106','181109','181110','181115'}; %'181017',
        data_path=fullfile(raw_data_path, monkey, '/');
    case 'Scooter'
        date_list={'181102','181103','181107','181108','181112','181113','181116','181117','181119'};  %,'181126','181130'
        data_path=fullfile(raw_data_path, monkey, '/');
    case 'Beaker_dummy'
        date_list={'181017'};
        data_path=fullfile(raw_data_path, monkey, '/');
    case 'Both'
        date_list={'181019','181025','181029','181030','181106','181109','181110','181115','181102','181103','181107','181108','181112','181113','181116','181117','181119'};
        monkey_list={'Beaker','Beaker','Beaker','Beaker','Beaker','Beaker','Beaker','Beaker','Scooter','Scooter','Scooter','Scooter','Scooter','Scooter','Scooter','Scooter','Scooter'};
end

file_list={};
file_path_list={};

for i=1:size(date_list,2)
    
    if strmatch(monkey,'Both')
        monkey_for_file=monkey_list{1,i};
    else
        monkey_for_file=monkey;
    end
    
    data_path=fullfile(raw_data_path, monkey, '/');
    cd(data_path)
    
    file_to_add_path=[date_list{1,i}  '/bhv'];
    cd(file_to_add_path)
    
    if exist([monkey_for_file '_' date_list{1,i} '_01_bhv.mat'])==2
        file_to_add=[monkey_for_file '_' date_list{1,i} '_01_bhv.mat'];
    elseif exist([monkey_for_file '_' date_list{1,i} '_00_bhv.mat'])==2
        file_to_add=[monkey_for_file '_' date_list{1,i} '_00_bhv.mat'];
    end
    
    file_list{end+1}=file_to_add;
    file_path_list{end+1}=file_to_add_path;
    
    cd ../..
    
end

cd(data_path)
data = ParseExploreExploitData(file_list, file_path_list);

data_all_trials=data;
data_all_trials.TrialNumber=ones(size(data.StopCond,1),1);
for i=1:size(data.StopCond,1)
    data_all_trials.TrialNumber(i)=i;
end

%First remove non 1/-1 stop cond trials
data.Block=data.Block(abs(data.StopCond)==1,:);
data.BestColor=data.BestColor(abs(data.StopCond)==1,:);
data.TargetColors=data.TargetColors(abs(data.StopCond)==1,:);
data.TargetDist=data.TargetDist(abs(data.StopCond)==1,:);
data.TargetLoc=data.TargetLoc(abs(data.StopCond)==1,:);
data.PopoutInd=data.PopoutInd(abs(data.StopCond)==1,:);
data.PopoutSize=data.PopoutSize(abs(data.StopCond)==1,:);
data.Resp=data.Resp(abs(data.StopCond)==1,:);
data.RespErr=data.RespErr(abs(data.StopCond)==1,:);
data.RT=data.RT(abs(data.StopCond)==1,:);
data.ChosenColor=data.ChosenColor(abs(data.StopCond)==1,:);
data.Reward=data.Reward(abs(data.StopCond)==1,:);
data.TrialNumber=data_all_trials.TrialNumber(abs(data.StopCond)==1,:);
data.IdleTime=data.IdleTime(abs(data.StopCond)==1,:);

data.StopCond=data.StopCond(abs(data.StopCond)==1,:);

%% Find transitions
day_sw = find(diff([0; data.Block(:, 1)]));
block_sw = find(diff([0; data.Block(:, 2)]));

%Calculate the length of blocks
block_len = NaN*ones(size(block_sw));
for cur_block_ind = 1:(length(block_sw)-1)
    %Is this the last block in the day?
    if data.Block(block_sw(cur_block_ind+1), 1) ~= data.Block(block_sw(cur_block_ind), 1)
        continue;
    end
    block_len(cur_block_ind) = block_sw(cur_block_ind+1) - block_sw(cur_block_ind);
end

%% Find performance around a block switch

block_best_col = data.BestColor(block_sw);
prev_best_col = NaN*ones(size(block_best_col));
for cur_block_ind = 2:length(block_sw)
    %Is this the first block in the day?
    if data.Block(block_sw(cur_block_ind), 1) ~= data.Block(block_sw(cur_block_ind-1), 1)
        continue;
    end
    prev_best_col(cur_block_ind) = data.BestColor(block_sw(cur_block_ind - 1));
end

%Order the colors
rel_err = data.RespErr;
[~, best_color_target] = min(data.TargetDist,[],2);
[~, worst_color_target] = max(data.TargetDist,[],2);
second_best_color_target=abs((best_color_target+worst_color_target)-6);

for i=1:length(best_color_target)
    best_color_target_distance(i,1)=data.TargetDist(i,best_color_target(i));
    second_best_color_target_distance(i,1)=data.TargetDist(i,second_best_color_target(i));
    worst_color_target_distance(i,1)=data.TargetDist(i,worst_color_target(i));
    best_color_target_loc(i,1)=data.TargetLoc(i,best_color_target(i));
    second_best_color_target_loc(i,1)=data.TargetLoc(i,second_best_color_target(i));
    worst_color_target_loc(i,1)=data.TargetLoc(i,worst_color_target(i));
    
    if ~isnan(data.Resp(i))
        chosen_distance(i,1)=data.TargetDist(i,data.Resp(i));
        chosen_loc(i,1)=data.TargetLoc(i,data.Resp(i));
    else
        chosen_distance(i,1)=NaN;
        chosen_loc(i,1)=NaN;
    end
end

%% Get the data

Choices = data.ChosenColor;

OptionsColors = data.TargetColors;

Best = data.BestColor;

X_pre(:,1)=Best(:); % Best color in the block

for i=1:length(best_color_target)
    BestOption(i,1)=OptionsColors(i,best_color_target(i,1));
    SecondBestOption(i,1)=OptionsColors(i,second_best_color_target(i,1));
    WorstOption(i,1)=OptionsColors(i,worst_color_target(i,1));
end

X_pre(:,2)=BestOption(:); % best colored offered: if chosen accuracy =1
X_pre(:,3)=SecondBestOption(:); % second best
X_pre(:,4)=WorstOption(:); %worst

Accuracy=zeros(size(Choices,1),1);
for i=1:size(Choices,1)
    if Choices(i)==BestOption(i)
        Accuracy(i)=1;
    end
end

Small_popout_chosen=NaN(size(Choices,1),1);
Big_popout_chosen=NaN(size(Choices,1),1);
for i=1:length(data.PopoutInd)
    if data.PopoutInd(i)>0 && data.PopoutSize(i)==0.75
        if data.TargetLoc(i,data.PopoutInd(i))==chosen_loc(i)
            Small_popout_chosen(i,1)=1;
        else
            Small_popout_chosen(i,1)=0;
        end
    elseif data.PopoutInd(i)>0 && data.PopoutSize(i)==2.25
        if data.TargetLoc(i,data.PopoutInd(i))==chosen_loc(i)
            Big_popout_chosen(i,1)=1;
        else
            Big_popout_chosen(i,1)=0;
        end
    end
end

Y_continuous(:,1)=Choices(:);

TrialInBlock=zeros(size(Choices,1),1);
BlockLength=zeros(size(Choices,1),1);
IsPrev=zeros(size(Choices,1),1);
PrevBest=zeros(size(Choices,1),1);

for j=1:(size(block_sw,1)-1) %block_switch_list
    for i=block_sw(j):block_sw(j+1)-1
        if ~isnan(block_len(j)) %remove the last block of the day (unfinished)
            TrialInBlock(i)=i-block_sw(j)+1;
            BlockLength(i)=block_sw(j+1)-block_sw(j);
            if isnan(prev_best_col(j))
                IsPrev(i)=0;
            else
                IsPrev(i)=1;
                PrevBest(i)=prev_best_col(j);
            end
        else
            TrialInBlock(i)=NaN;
            BlockLength(i)=NaN;
            IsPrev(i)=NaN;
            PrevBest(i)=NaN;
        end
    end
end
for i=block_sw(end):size(Choices,1)
    TrialInBlock(i)=NaN; %i-block_sw(end)+1;
    BlockLength(i)=NaN; %size(Choices,1)-block_sw(end);
    if isnan(prev_best_col(end))
        IsPrev(i)=NaN; %0;
    else
        IsPrev(i)=NaN; %1;
        PrevBest(i)=NaN; %prev_best_col(end);
    end
end

X_pre(:,5)=TrialInBlock(:);
X_pre(:,6)=BlockLength(:);

X_pre(:,7)=IsPrev(:);
X_pre(:,8)=PrevBest(:);

Reward=data.Reward;

X_pre(:,9)=Reward(:);

Rmax=zeros(size(Choices,1),1);
for i=1:size(block_sw,1)-1
    Rmax(block_sw(i):block_sw(i+1))=max(Reward(block_sw(i):block_sw(i+1)));
end
Rmax(block_sw(end):size(Choices,1))=max(Reward(block_sw(end):size(Choices,1)));

X_pre(:,10)=Rmax(:);
X_pre(:,11)=0;
X_pre(:,12)=Best(:);

%location of targets
X_pre(:,13)=best_color_target_loc(:);
X_pre(:,14)=second_best_color_target_loc(:);
X_pre(:,15)=worst_color_target_loc(:);

%popout
for i=1:length(data.PopoutInd)
    if data.PopoutInd(i)>0
        X_pre(i,16)=data.TargetLoc(i,data.PopoutInd(i));
    else
        X_pre(i,16)=0;
    end
end
X_pre(:,17)=data.PopoutSize(:);

X_pre(:,18)=data.RT(:);

for i=1:size(X_pre,1)
    if isnan(X_pre(i,16))
        X_pre(i,16)=0;
        X_pre(i,17)=0;
    end
end

X_pre(1,19)=0;
for i=2:size(X_pre,1)
    if IsPrev(i)==1
        X_pre(i,19)=chosen_loc(i-1);
    else
        X_pre(i,19)=NaN;
    end
end

%find out what the attempt rate in the X past trials
AttemptedHistory=5;

data_all_trials.AttemptedTrial=zeros(size(data_all_trials.StopCond,1),1);
data_all_trials.SuccessTrial=zeros(size(data_all_trials.StopCond,1),1);

for i=1:size(data_all_trials.StopCond,1)
    if data_all_trials.StopCond(i)==1 || data_all_trials.StopCond(i)==-1 || data_all_trials.StopCond(i)==-5 || data_all_trials.StopCond(i)==-6
        data_all_trials.AttemptedTrial(i)=1;
    else
        data_all_trials.AttemptedTrial(i)=0;
    end
    if data_all_trials.StopCond(i)==1 || data_all_trials.StopCond(i)==-1
        data_all_trials.SuccessTrial(i)=1;
    else
        data_all_trials.SuccessTrial(i)=0;
    end
end

for i=1:size(data.StopCond,1)
    flag=100;
    for j=1:length(day_sw)
        if data_all_trials.TrialNumber(data.TrialNumber(i))-day_sw(j)>=0 && data_all_trials.TrialNumber(data.TrialNumber(i))-day_sw(j)<=AttemptedHistory
            flag=data_all_trials.TrialNumber(data.TrialNumber(i))-day_sw(j);
        end
    end
    if flag==100
        X_pre(i,20)=mean(data_all_trials.AttemptedTrial(data_all_trials.TrialNumber(data.TrialNumber(i))-AttemptedHistory:data_all_trials.TrialNumber(data.TrialNumber(i))));
        X_pre(i,21)=mean(data_all_trials.SuccessTrial(data_all_trials.TrialNumber(data.TrialNumber(i))-AttemptedHistory:data_all_trials.TrialNumber(data.TrialNumber(i))));
    else
        X_pre(i,20)=mean(data_all_trials.AttemptedTrial(data_all_trials.TrialNumber(data.TrialNumber(i))-flag:data_all_trials.TrialNumber(data.TrialNumber(i))));
        X_pre(i,21)=mean(data_all_trials.SuccessTrial(data_all_trials.TrialNumber(data.TrialNumber(i))-flag:data_all_trials.TrialNumber(data.TrialNumber(i))));
    end
end
X_pre(:,22)=data.IdleTime;
X_pre(:,23)=chosen_loc;

%%
%remove nan values
Index = ~isnan(X_pre(:,6));
y_continuous = Y_continuous(Index);
X_data=X_pre(Index,:);
Accuracy=Accuracy(Index);

X_data(:,11)=y_continuous(:); %chosen color

%find corresponding block switch values when you remove NaN
J = cumsum(Index);
Block_switch=zeros(size(block_sw));
Day_switch=zeros(size(day_sw));
for i=1:size(block_sw)
    if Index(block_sw(i))==1
        Block_switch(i)=J(block_sw(i));
    end
end
Block_switch=Block_switch(Block_switch>0);
for i=1:size(day_sw)
    Day_switch(i)=J(day_sw(i));
end

Block_day=zeros(size(Block_switch,1),1);
for i=1:size(Block_switch)
    for j=1:size(Day_switch)-1
        if Block_switch(i)>=Day_switch(j) && Block_switch(i)<Day_switch(j+1)
            Block_day(i)=j;
        end
    end
    if Block_switch(i)>=Day_switch(end)
        Block_day(i)=size(Day_switch,1);
    end
end

X_data(Block_switch,24)=1;

N_forward_trials=size(X_data,1);

Choices_mat=zeros(3,N_bins,N_forward_trials);

Chosen=zeros(1,N_forward_trials);

Best_chosen=NaN(1,N_forward_trials);
SecondBest_chosen=NaN(1,N_forward_trials);
Worst_chosen=NaN(1,N_forward_trials);

Reward=zeros(1,N_forward_trials);
Reward_max=zeros(1,N_forward_trials);

Chosen_ind=NaN(N_forward_trials,1);


%% Bin the colors
%Oustanding questions = how many bins? where do we start? for now we set it
%to  0
color_binned =  0:(2*pi/N_bins):(2*pi-(2*pi/N_bins));

for i=1:N_forward_trials
    if ~isnan(X_data(i,11))
        %Reward nat for update rule
        if X_data(i,11)>=color_binned(N_bins)
            Reward(1,i)=X_data(i,9)/X_data(i,10);
            Reward_max(1,i)=X_data(i,10);
            Chosen(1,i)=N_bins;
        else
            for k=1:N_bins-1
                if X_data(i,11)>=color_binned(k) && X_data(i,11)<color_binned(k+1)
                    Reward(1,i)=X_data(i,9)/X_data(i,10);
                    Reward_max(1,i)=X_data(i,10);
                    Chosen(1,i)=k;
                end
            end
        end
        %Choice offered
        for l=1:3
            if X_data(i,l+1)>=color_binned(N_bins)
                Choices_mat(l,i)=N_bins;
            else
                for k=1:N_bins-1
                    if X_data(i,l+1)>=color_binned(k) && X_data(i,l+1)<color_binned(k+1)
                        Choices_mat(l,i)=k;
                    end
                end
            end
        end
        for k=1:N_bins
            if Chosen(1,i)==k
                if Choices_mat(1,i)==k
                    Best_chosen(1,i)=1;
                    SecondBest_chosen(1,i)=0;
                    Worst_chosen(1,i)=0;
                elseif Choices_mat(2,i)==k
                    SecondBest_chosen(1,i)=1;
                    Best_chosen(1,i)=0;
                    Worst_chosen(1,i)=0;
                else
                    Worst_chosen(1,i)=1;
                    Best_chosen(1,i)=0;
                    SecondBest_chosen(1,i)=0;
                end
            end
        end
        if Best_chosen(1,i)==1
            Chosen_ind(i)=1;
        elseif SecondBest_chosen(1,i)==1
            Chosen_ind(i)=2;
        elseif Worst_chosen(1,i)==1
            Chosen_ind(i)=3;
        end
    end
end

%Is not first of the day
IsNotFirst=ones(N_forward_trials,1);
IsNotFirst(Day_switch)=0;

%find the days
Sess_number=ones(N_forward_trials,1).*length(Day_switch);
for j=2:length(Day_switch)
    Sess_number(Day_switch(j-1):Day_switch(j)-1)=j-1;
end

%inupt to VBA
inG.nOptions=3;
inG.nTrials=N_forward_trials;
inG.N_channels=N_channels;

inG.chosen_color=y_continuous';
inG.payoff=Reward;

inG.options_color=X_data(:,2:4)';
inG.locations=X_data(:,13:15)';
inG.isPrev=X_data(:,7)';
inG.PopLoc=X_data(:,16)';
inG.PopSize=X_data(:,17)';
inG.IsPrev=IsNotFirst';

inG.choices=[Best_chosen;SecondBest_chosen;Worst_chosen];

inG.Reward=Reward;
inG.Reward_max=Reward_max;
inG.color_binned=color_binned;
inG.N_forward_trials=N_forward_trials;
inG.N_learning=length(Block_switch);
inG.X_data=X_data;
inG.N_bins=N_bins;
inG.Best_chosen=Best_chosen;
inG.SecondBest_chosen=SecondBest_chosen;
inG.Worst_chosen=Worst_chosen;
inG.N_channels=N_channels;
inG.Chosen_ind=Chosen_ind';
inG.Accuracy=Accuracy';

inG.Sess_number=Sess_number';
end

