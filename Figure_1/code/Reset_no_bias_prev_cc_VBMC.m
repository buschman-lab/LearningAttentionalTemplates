function [u, posterior, out] = Reset_no_bias_prev_cc_VBMC (data)
% // VBMC toolbox //////////////////////////////////////////////////////////
%
% OUT:
%   - data: same strucutre as the input with defaults filled in.
%   - posterior, out: results of the inversion
%
% /////////////////////////////////////////////////////////////////////////

N_param=11;
Nrunsmax=5;

% prepare data for model inversion
% =========================================================================
% observations
y = data.choices;

% inputs

% prev_chosen_color_binned=u(1);
% prev_feedback_binned=u(2);
% options_color=u(3:5);
% options_location=u(6:8);
% popout_location=u(9);
% popout_size=u(10);
% isPrev=u(11);


u = [ 0, data.chosen_color(:,1:end-1) ;  % previous choice
    nan, data.payoff(:,1:end-1) ;
    data.options_color;
    data.locations ;
    data.PopLoc ;
    data.PopSize
    data.IsPrev] ; % options offered (no last one)

%intersect nan for each switch in the days

for i=2:size(u,2)
    if u(11,i)==0
        u(1,i)=0;
        u(2,i)=nan;
    end
end

% note: the first input is a nan as there is no "previous trial" for the
% first trial of the experiment

% specify model
% =========================================================================
% observation and evolution functions
% -------------------------------------------------------------------------

in = struct ('nOptions', data.nOptions, 'N_channels', data.N_channels, 'N_bins', data.N_bins);

lower_bound=-1;
upper_bound=1;

xstart=zeros(1,N_param);
xlb=ones(1,N_param).*lower_bound;
xub=ones(1,N_param).*upper_bound;

xplb=ones(1,N_param).*lower_bound*0.68;
xpub=ones(1,N_param).*upper_bound*0.68;

%Kappa
xstart(1)=2.5;
xlb(1)=0.5;
xub(1)=5;
xplb(1)=1;
xpub(1)=3.5;

%Learing_rate
xstart(2)=0.5;
xlb(2)=0;
xub(2)=6;
xplb(2)=0.1;
xpub(2)=1;

%location
xplb(3:5)=-0.2;
xpub(3:5)=0.2;

%popout
xplb(6)=-0.5;
xpub(6)=0.05;
xplb(7)=-0.05;
xpub(7)=0.5;

%Reset
xstart(8)=0.7;
xlb(8)=0;
xub(8)=1.2;
xplb(8)=0.5;
xpub(8)=0.9;

%Volatility
xstart(9)=0.05;
xlb(9)=0;
xub(9)=1;
xplb(9)=0.01;
xpub(9)=0.5;

%Color preference
xstart(10)=0;
xplb(10)=-1;
xpub(10)=1;
xlb(10)=-pi;
xub(10)=pi;

%priors
priors.x0 = xstart;
priors.xplb = xplb;
priors.xpub = xpub;
priors.xlb = xlb;
priors.xub = xub;


    function ll = wraper_f_Reset_no_bias_prev_cc_VBMC(P, u, y, in)
        [ll,~]=f_Reset_no_bias_prev_cc_VBMC(P, u, y, in);
    end


%% Parameters estimates 1st round

fun_bads=@(x) - wraper_f_Reset_no_bias_prev_cc_VBMC(x, u, y, in);

%Note: bads minimise a function so we use -LL
opt_options = [];
opt_options.Display = 'final';
priors.x0_bads = bads(fun_bads, priors.x0,priors.xlb,priors.xub,priors.xplb, priors.xpub,[],opt_options);

%% Parameters estimates 2nd round
%center the prioirs but keeps a wide search
priors.x0_bads(2)=abs(priors.x0_bads(2));
priors.x0_bads(10)=mod(priors.x0_bads(10)+pi,2*pi)-pi;
priors.x0_bads(11)=abs(priors.x0_bads(11));

priors.xlb_bads=priors.x0_bads-1;
priors.xub_bads=priors.x0_bads+1;
priors.xplb_bads=priors.x0_bads-0.68;
priors.xpub_bads=priors.x0_bads+0.68;

%Kappa must be >0
if priors.xlb_bads(1)<0.25
    priors.xlb_bads(1)=0.25;
    priors.xplb_bads(1)=0.3;
end

%Learnign rate must be >0
if priors.xlb_bads(2)<0
    priors.xlb_bads(2)=0;
    priors.xplb_bads(2)=0.05;
end

%Volatility is tricky, extend the bounds
priors.xlb_bads(9)=0;
priors.xub_bads(9)=1;
if priors.xplb_bads(9)<0.001
    priors.xplb_bads(9)=0.001;
end
    
%Color preference
priors.xplb_bads(10)=priors.x0_bads(10)-pi/2;
priors.xpub_bads(10)=priors.x0_bads(10)+pi/2;
priors.xlb_bads(10)=priors.x0_bads(10)-pi;
priors.xub_bads(10)=priors.x0_bads(10)+pi;
   
fun_vbmc = @(x) wraper_f_Reset_no_bias_prev_cc_VBMC(x, u, y, in);

posterior.exitflag=0;
options = vbmc('defaults');
options.RetryMaxFunEvals = options.MaxFunEvals;

for iRun=1:Nrunsmax
    if posterior.exitflag~=1
        
        fprintf('  VBMC run #%d/%d...\n', iRun, Nrunsmax);
        
        if iRun == 1
            [posterior.vp,posterior.F,posterior.F_std,posterior.exitflag,posterior.output] = vbmc(fun_vbmc,priors.x0_bads,priors.xlb_bads,priors.xub_bads,priors.xplb_bads,priors.xpub_bads,options);
        else
            [posterior.vp,posterior.F,posterior.F_std,posterior.exitflag,posterior.output] = vbmc(fun_vbmc,posterior.vp,[],[],[],[],options);
        end
    end
end

posterior.Xs = vbmc_rnd(posterior.vp,3e5);
for i=1:N_param
    posterior.mean(i) = mean(posterior.Xs(:,i));  % Posterior mean
    posterior.var(i) = std(posterior.Xs(:,i));  % Posterior mean
end
posterior.mean(2) = mean(abs(posterior.Xs(:,2)));
posterior.var(2) = std(abs(posterior.Xs(:,2)));
posterior.mean(10) = angle(mean(exp(1i.*posterior.Xs(:,10))));
posterior.var(10)=0;
posterior.mean(11) = mean(abs(posterior.Xs(:,11)));
posterior.var(11) = std(abs(posterior.Xs(:,11)));

[out.suffStat.LL, out.suffStat.gx]=f_Reset_no_bias_prev_cc_VBMC(posterior.mean,u,y,in);

out.fit.BIC = N_param*log(length(y))-2*out.suffStat.LL; 
out.fit.AIC = 2*N_param-2*out.suffStat.LL; 


end


