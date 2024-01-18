function [ll, gx] = f_NoReset_fit_softmax(P, u, y, in)

% task dimension
N_channels=in.N_channels;
N_bins=in.N_bins;

% inputs
chosen_color = u(1,:); %previous trial
feedback = u(2,:); %previous trial

options_color=u(3:5,:);
options_location=u(6:8,:);
popout_location=u(9,:);
popout_size=u(10,:);
isPrev=u(11,:);
previous_color=u(1,:);

% observation parameters
% =========================================================================

%Location bias
B_location(1)=0;
B_location(2:4)=P(3:5);
Location(1:4)=1:4;

%popout bias
B_popout(1:2)=P(6:7);
Popout(1)=0.75;
Popout(2)=2.25;

%previous chosen color
B_previous_color=P(8);

%preference for color
Color_pref=mod(P(9),2*pi);
B_color_pref=abs(P(10));

% inverse temperature is fitted in this version
beta = P(11);

% evolution parameters
% =========================================================================
%Basis functions (bounded between 0 and 1)
BF_kappa=P(1);

%Learning (bounded between 0 and 1)
learning_rate = abs(P(2));

%precision effect on learning rate
Precision_effect=0;

%angle for calculating the confidence in rule
angle_rule=zeros(N_bins,1);
for i=1:N_bins
    angle_rule(i,1)=(i-1)*2*pi/N_bins;
end

%% Binning
color_binned =  0:(2*pi/N_bins):(2*pi-(2*pi/N_bins));

%% Function approximation basis
channels_binned =  0:(2*pi/N_channels):(2*pi-(2*pi/N_channels));
basis_functions=zeros(N_channels,N_bins);

for i=1:N_bins
    for j=1:N_channels
        basis_functions(j,i) = vonmisespdf(color_binned(i),channels_binned(j),BF_kappa);
    end
end

%% Loop over trials

for t=1:size(chosen_color,2)
    
    % Learning
    
    % In the case there was no feedback (first trial), do nothing
    % =========================================================================
    
    if isnan(feedback(t)) %1st trial of the sesssion
        weights(:,t)=zeros(N_channels,1);
    else
        
        % Apply the RL
        % =========================================================================
        
        %bin the offered color
        if chosen_color(t)>=color_binned(N_bins)
            chosen_color_bin=N_bins;
        else
            for k=1:N_bins-1
                if chosen_color(t)>=color_binned(k) && chosen_color(t)<color_binned(k+1)
                    chosen_color_bin=k;
                end
            end
        end
        
        
        value_old=zeros(N_bins,1);
        for c=1:N_bins
            value_old(c)=sum(weights(:,t-1).*basis_functions(:,c));
        end
        precision=norm(sum(value_old.*exp(1i*angle_rule)));
        
        pe = feedback(t) - value_old(chosen_color_bin);
        
        weights(:,t)= weights(:,t-1) + (learning_rate+Precision_effect/(precision+1))*pe*basis_functions(:,chosen_color_bin);
        
    end
    
    % Choice
    
    % Extract value of each option
    % -------------------------------------------------------------------------
    
    value_new=zeros(N_bins,1);
    for c=1:N_bins
        value_new(c)=sum(weights(:,t).*basis_functions(:,c));
    end
    
    %bin the offered color
    for l=1:3
        if options_color(l,t)>=color_binned(N_bins)
            options_color_bin(l)=N_bins;
        else
            for k=1:N_bins-1
                if options_color(l,t)>=color_binned(k) && options_color(l,t)<color_binned(k+1)
                    options_color_bin(l)=k;
                end
            end
        end
    end
    
    %get the value and add the biases
    for l=1:3
        choice_value(l)=value_new(options_color_bin(l))+...
            sum(B_location(:).*(options_location(l,t)==Location(:)))+...
            (options_location(l,t)==popout_location(t)).*sum(B_popout(:).*(popout_size(t)==Popout(:)))+...
            isPrev(t).*B_previous_color.*(1-abs(mod(options_color(l,t)-previous_color(t)+pi,2*pi)-pi)/pi)+...
            B_color_pref.*(1-abs(mod(options_color(l,t)-Color_pref+pi,2*pi)-pi)/pi);
    end
    
    % apply softmax
    % -------------------------------------------------------------------------
    gx(:,t) =  exp(choice_value./beta)./(sum(exp(choice_value./beta)));
    
    
end

ll=sum(log(gx((y==1))));


end


