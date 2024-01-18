function Model_predictions = get_WSLS_model_VBMC(Parameters,u,in)

%% Task variables
Location(1:4)=1:4;
Popout(1)=0.75;
Popout(2)=2.25;


%% task dimension
N_bins=in.N_bins;
color_binned =  0:(2*pi/N_bins):(2*pi-(2*pi/N_bins));

%% Initialize
Value=zeros(N_bins,size(u,2));
Option_value=zeros(3,size(u,2));
Choice_value=zeros(3,size(u,2));
Choice_prediction=zeros(3,size(u,2));

%% Update of belief => evolution function
%The belief is the previous color if the re

for t=1:size(u,2)
    
    chosen_color = u(1,t); %previous trial
    feedback = u(2,t); %previous trial
    
    if isnan(feedback)
        % In the case there was no feedback (first trial), reset
        % =========================================================================
        Value(1:N_bins,t)=0;
    else
        % Apply the conditional valuation: distance to previous chosen
        % color if the reward is high enough, 0 otherwise
        % =========================================================================
        if feedback > Parameters.reward_threshold
            for c=1:N_bins
                Value(c,t)=abs(mod(color_binned(c)-chosen_color+pi,2*pi)-pi)/pi; %normalized absolute angular distance to previously chosen color
            end
        else
            Value(1:N_bins,t)=0;
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
            Parameters.B_color_pref.*(1-abs(mod(options_color(l)-Parameters.color_pref+pi,2*pi)-pi)/pi);
    end
    
    % inverse temperature
    beta = 0.3;
    
    % apply softmax
    % -------------------------------------------------------------------------
    Choice_prediction(:,t) = exp(Choice_value(:,t)./beta)./(sum(exp(Choice_value(:,t)./beta)));
    
end

%% get the model predictions in the same format as before

Model_predictions.parameter.reward_threshold=Parameters.reward_threshold;

Model_predictions.parameter.Noise_softmax_action=beta;
Model_predictions.parameter.B_location=Parameters.B_location;
Model_predictions.parameter.B_popout=Parameters.B_popout;
Model_predictions.parameter.Color_pref=Parameters.color_pref;
Model_predictions.parameter.B_color_pref=Parameters.B_color_pref;

Model_predictions.model_specificities.N_bins=N_bins;

Model_predictions.model_outputs.Option_value=Option_value;
Model_predictions.model_outputs.Value_for_choice=Value;
Model_predictions.model_outputs.Choice_value=Choice_value;
Model_predictions.model_outputs.Choice_prediction=Choice_prediction;













