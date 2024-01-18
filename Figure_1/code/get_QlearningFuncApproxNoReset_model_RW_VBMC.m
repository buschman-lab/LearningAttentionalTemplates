function Model_predictions = get_QlearningFuncApproxNoReset_model_RW_VBMC(Parameters,u,in)

% High beta_reset => categorical
beta_reset=10000;

Location(1:4)=1:4;
Popout(1)=0.75;
Popout(2)=2.25;

Parameters.Precision_effect=0;

%% rename entries

% task dimension
N_channels=in.N_channels;
N_bins=in.N_bins;

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

Weights=zeros(N_channels,size(u,2));
Prediction_error=zeros(1,size(u,2));
Value=zeros(N_bins,size(u,2));
Option_value=zeros(3,size(u,2));
Choice_value=zeros(3,size(u,2));
Choice_prediction=zeros(3,size(u,2));

%% Update of belief => evolution function

for t=1:size(u,2)
    
    chosen_color = u(1,t); %previous trial
    feedback = u(2,t); %previous trial
    
    %bin the offered color
    if chosen_color>=color_binned(N_bins)
        chosen_color_bin=N_bins;
    else
        for k=1:N_bins-1
            if chosen_color>=color_binned(k) && chosen_color<color_binned(k+1)
                chosen_color_bin=k;
            end
        end
    end
    
    if isnan(feedback)
        % In the case there was no feedback (first trial), reset
        % =========================================================================
        Weights(1:N_channels,t)=0;
        Prediction_error(t)=0;
        Value(1:N_bins,t)=0;
    else
        % Apply the RL
        % =========================================================================
        
        weights=Weights(:,t-1);
        
        value=zeros(N_bins,1);
        for c=1:N_bins
            value(c)=sum(Weights(:,t-1).*basis_functions(:,c));
        end
        
        pe = feedback - value(chosen_color_bin);
                
        Weights(:,t)= weights + Parameters.learning_rate*pe*basis_functions(:,chosen_color_bin);
        
        Prediction_error(t)=pe; %that's the previous trial PE!!!
        
        for c=1:N_bins
            Value(c,t)=sum(Weights(:,t).*basis_functions(:,c));
        end
        
    end
    
end


%% Extract the value => observation function

for t=1:size(u,2)
    
    % inputs
    options_color=u(3:5,t);
    options_location=u(6:8,t);
    popout_location=u(9,t);
    popout_size=u(10,t);
    isPrev=u(11,t);
    previous_color=u(1,t);
    
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
    
end

%% get the model predictions in the same format as before

down_sampled_Value_for_choice=NaN(N_channels,size(Value,2));
for i=1:size(Value,2)
    for n=1:N_channels-1
        down_sampled_Value_for_choice(n,i)=mean(Value(1+(n-1)*round(N_bins/N_channels):n*round(N_bins/N_channels),i),1);
    end
    down_sampled_Value_for_choice(N_channels,i)=mean(Value(1+(n-1)*round(N_bins/N_channels):end,i),1);
end

Model_predictions.parameter.BF_kappa=Parameters.BF_kappa;
Model_predictions.parameter.softmax_reset=beta_reset;
Model_predictions.parameter.learning_rate=Parameters.learning_rate;

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

Model_predictions.model_outputs.down_sampled_Value_for_choice=down_sampled_Value_for_choice;

Model_predictions.model_outputs.RPE=Prediction_error;













