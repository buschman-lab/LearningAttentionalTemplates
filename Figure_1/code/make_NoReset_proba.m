function [Model_predictions, u, y] = make_NoReset_proba(Parameters,data,feedback_each_option)
%Here we genrate the data from the model rather tahn fitting the monkeys'
%choices

% High beta_reset => categorical
beta_reset=10000;

Location(1:4)=1:4;
Popout(1)=0.75;
Popout(2)=2.25;

Parameters.Precision_effect=0;

%% rename entries

% task dimension
N_channels=6;
N_bins=100;

angle_rule=zeros(N_bins,1);
for i=1:N_bins
    angle_rule(i,1)=(i-1)*2*pi/N_bins;
end

%% Function approximation basis
%each bin s a channel, the update is on the channel, with a width defined
%by kappa
channels_binned =  0:(2*pi/N_channels):(2*pi-(2*pi/N_channels));
color_binned =  0:(2*pi/N_bins):(2*pi-(2*pi/N_bins));

basis_functions=zeros(N_channels,N_bins);

for i=1:N_bins
    for j=1:N_channels
        basis_functions(j,i) = vonmisespdf(color_binned(i),channels_binned(j),Parameters.BF_kappa);
    end
end

%% Initialize

Switch=zeros(1,length(data.IsPrev));
Weights=zeros(N_channels,length(data.IsPrev));
Trials_since_switch=zeros(1,length(data.IsPrev));
Prediction_error=zeros(1,length(data.IsPrev));
Precision=zeros(1,length(data.IsPrev));
Value=zeros(N_bins,length(data.IsPrev));
Option_value=zeros(3,length(data.IsPrev));
Choice_value=zeros(3,length(data.IsPrev));
Choice_prediction=zeros(3,length(data.IsPrev));

previous_chosen_color=zeros(1,length(data.IsPrev));
previous_feedback=zeros(1,length(data.IsPrev));

chosen_id=zeros(1,length(data.IsPrev));
chosen_color=zeros(1,length(data.IsPrev));
feedback=zeros(1,length(data.IsPrev));
y=zeros(3,length(data.IsPrev));

%% Update of belief => evolution function

for t=1:length(data.IsPrev)
    if data.IsPrev(t)==0
        % In the case there was no feedback (first trial), reset/initalize
        % =========================================================================
        Weights(1:N_channels,t)=0;
        Trials_since_switch(t)=1;
        Prediction_error(t)=0;
        Value(1:N_bins,t)=0;
        Precision(t)=0;
        Switch(t)=0;
        %choose an option randomly
        chosen_id(t)=randi(3);
        chosen_color(t)=data.options_color(chosen_id(t),t);
        feedback(t)=feedback_each_option(chosen_id(t),t);
        y(chosen_id(t),t)=1;
        previous_chosen_color(t)=NaN;
        previous_feedback(t)=NaN;
    else
        % Apply the RL
        % =========================================================================
        
        previous_chosen_color(t)=chosen_color(t-1);
        previous_feedback(t)=feedback(t-1);
        
        %bin the previous chosen color color
        if previous_chosen_color(t)>=color_binned(N_bins)
            prev_chosen_color_bin=N_bins;
        else
            for k=1:N_bins-1
                if previous_chosen_color(t)>=color_binned(k) && previous_chosen_color(t)<color_binned(k+1)
                    prev_chosen_color_bin=k;
                end
            end
        end
        
        
        weights=Weights(:,t-1);
        
        value=zeros(N_bins,1);
        for c=1:N_bins
            value(c)=sum(Weights(:,t-1).*basis_functions(:,c));
        end
        
        pe = previous_feedback(t) - value(prev_chosen_color_bin);
                
        Weights(:,t)= weights + Parameters.learning_rate*pe*basis_functions(:,prev_chosen_color_bin);
        
        Prediction_error(t)=pe; %that's the previous trial PE!!!
        
        for c=1:N_bins
            Value(c,t)=sum(Weights(:,t).*basis_functions(:,c));
        end
                
        
        % Extract the value => observation function
        % inputs
        options_color=data.options_color(:,t);
        options_location=data.locations(:,t);
        popout_location=data.PopLoc(t);
        popout_size=data.PopSize(t);
        isPrev=data.IsPrev(t);
        previous_color=previous_chosen_color(t);
        
        options_color_bin=[0 0 0];
        %bin the offered color
        for l=1:3
            if options_color(l)>=color_binned(N_bins)
                options_color_bin(l)=N_bins;
            else
                for k=1:N_bins-1
                    if options_color(l)>=color_binned(k) && options_color(l)<color_binned(k+1)
                        options_color_bin(l)=k;
                    end
                end
            end
        end
        
        value=Value(:,t);
        %get the value and add the biases
        for l=1:3
            Option_value(l,t)=value(options_color_bin(l));
            Choice_value(l,t)=value(options_color_bin(l))+sum(Parameters.B_location(:).*(options_location(l)==Location(:)))+...
                (options_location(l)==popout_location).*sum(Parameters.B_popout(:).*(popout_size==Popout(:)))+...
                isPrev.*Parameters.B_previous_color.*(1-abs(mod(options_color(l)-previous_color+pi,2*pi)-pi)/pi)+...
                Parameters.B_color_pref.*(1-abs(mod(options_color(l)-Parameters.color_pref+pi,2*pi)-pi)/pi);
        end
        
        % inverse temperature
        beta = 0.3;
        
        % apply softmax
        % -------------------------------------------------------------------------
        Choice_prediction(:,t) = exp(Choice_value(:,t)./beta)./(sum(exp(Choice_value(:,t)./beta)));
        
        %Softmax for choice
        x=rand;
        if x<=Choice_prediction(1,t)
            chosen_id(t)=1;
        elseif x<=Choice_prediction(1,t)+Choice_prediction(2,t)
            chosen_id(t)=2;
        else
            chosen_id(t)=3;
        end
        chosen_color(t)=data.options_color(chosen_id(t),t);
        feedback(t)=feedback_each_option(chosen_id(t),t);
        y(chosen_id(t),t)=1;
    end
end
    
    %% get the model predictions in the same format as before
    
    Model_predictions.parameter.BF_kappa=Parameters.BF_kappa;
    Model_predictions.parameter.softmax_reset=beta_reset;
    Model_predictions.parameter.Volatility=Parameters.Volatility;
    Model_predictions.parameter.learning_rate=Parameters.learning_rate;
    Model_predictions.parameter.Precision_effect=Parameters.Precision_effect;
    
    Model_predictions.parameter.Noise_softmax_action=beta;
    Model_predictions.parameter.B_location=Parameters.B_location;
    Model_predictions.parameter.B_popout=Parameters.B_popout;
    Model_predictions.parameter.B_previous_color=Parameters.B_previous_color;
    Model_predictions.parameter.Color_pref=Parameters.color_pref;
    Model_predictions.parameter.B_color_pref=Parameters.B_color_pref;
    
    Model_predictions.model_specificities.N_bins=N_bins;
    Model_predictions.model_specificities.N_channels=N_channels;
    
    Model_predictions.model_outputs.Weight=Weights;
    Model_predictions.model_outputs.Option_value=Option_value;
    Model_predictions.model_outputs.Value_for_choice=Value;
    Model_predictions.model_outputs.Choice_value=Choice_value;
    Model_predictions.model_outputs.Choice_prediction=Choice_prediction;
        
    Model_predictions.model_outputs.RPE=Prediction_error;
    Model_predictions.model_outputs.Precision=Precision;
    Model_predictions.model_outputs.Switch=Switch;
    Model_predictions.model_outputs.Trials_since_switch=Trials_since_switch;
    
    %genrate u for inversion
    
u = [ 0, chosen_color(1:end-1) ;  % previous choice
    nan, feedback(1:end-1) ;
    data.options_color;
    data.locations ;
    data.PopLoc ;
    data.PopSize
    data.IsPrev] ; % options offered (no last one)
    
    
    
    
    
    
    
    
    
    
    
