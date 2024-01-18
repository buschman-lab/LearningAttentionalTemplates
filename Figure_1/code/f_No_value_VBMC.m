function [ll, gx] = f_No_value_VBMC(P, u, y, in)

% task dimension
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

%Previous chosen color
B_previous_color=P(1);

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


%% Binning
color_binned =  0:(2*pi/N_bins):(2*pi-(2*pi/N_bins));

%% Loop over trials

for t=1:size(chosen_color,2)
        
    %No value computation
    
    Value(1:N_bins)=0;    
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
            isPrev(t).*B_previous_color.*(1-abs(mod(options_color(l,t)-previous_color(t)+pi,2*pi)-pi)/pi)+...
            B_color_pref.*(1-abs(mod(options_color(l,t)-Color_pref+pi,2*pi)-pi)/pi);
    end
    
    % apply softmax
    % -------------------------------------------------------------------------
    gx(:,t) =  exp(choice_value./beta)./(sum(exp(choice_value./beta)));
    
end

ll=sum(log(gx((y==1))));


end


