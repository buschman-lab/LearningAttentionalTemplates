function [ll, gx] = f_WSLS_flip_VBMC(P, u, y, in)

% task dimension
N_bins=in.N_bins;

% inputs
chosen_color = u(1,:); %previous trial
feedback = u(2,:); %previous trial

options_color=u(3:5,:);
options_location=u(6:8,:);
popout_location=u(9,:);
popout_size=u(10,:);

% observation parameters
% =========================================================================

%Location bias
B_location(1)=0;
B_location(2:4)=P(2:4);
Location(1:4)=1:4;

%popout bias
B_popout(1:2)=P(5:6);
Popout(1)=0.75;
Popout(2)=2.25;

%preference for color
Color_pref=mod(P(7),2*pi);
B_color_pref=abs(P(8));

% inverse temperature
beta = 0.3;

% evolution parameters
% =========================================================================

%reward level to remember the previous color chosen
reward_threshold=P(1);

%% Binning
color_binned =  0:(2*pi/N_bins):(2*pi-(2*pi/N_bins));

%% Loop over trials

for t=1:size(chosen_color,2)
        
    if isnan(feedback)
        % In the case there was no feedback (first trial), reset
        % =========================================================================
        Value(1:N_bins)=0;
    else
        % Apply the conditional valuation: distance to previous chosen
        % color if the reward is high enough, 0 otherwise
        % =========================================================================
        if feedback(t) > reward_threshold
            for c=1:N_bins
                Value(c)=1-abs(mod(color_binned(c)-chosen_color(t)+pi,2*pi)-pi)/pi; %normalized absolute angular distance to previously chosen color
            end
        else
            for c=1:N_bins
                Value(c)=abs(mod(color_binned(c)-chosen_color(t)+pi,2*pi)-pi)/pi; %if the reward is too low, repulsion
            end
        end
    end
    
    % Choice
    % Extract value of each option
    % -------------------------------------------------------------------------
       
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
        choice_value(l,1)=Value(options_color_bin(l))+...
            sum(B_location(:).*(options_location(l,t)==Location(:)))+...
            (options_location(l,t)==popout_location(t)).*sum(B_popout(:).*(popout_size(t)==Popout(:)))+...
            B_color_pref.*(1-abs(mod(options_color(l,t)-Color_pref+pi,2*pi)-pi)/pi);
    end
    
    % apply softmax
    % -------------------------------------------------------------------------
    gx(:,t) =  exp(choice_value./beta)./(sum(exp(choice_value./beta)));
    
end

ll=sum(log(gx((y==1))));


end


