function [posterior, out] = QlearningFuncApproxReset_RW_VBMC_recovery (data, u, y)
% // VBMC toolbox //////////////////////////////////////////////////////////
%
% OUT:
%   - data: same strucutre as the input with defaults filled in.
%   - posterior, out: results of the inversion
%
% /////////////////////////////////////////////////////////////////////////

%set random seed
rng('shuffle')

N_param=12;
Nrunsmax=5;

% prepare data for model inversion
% =========================================================================
% we provide the observations

% inputs
% weprovide u

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

%Learning_rate
xstart(2)=0.5;
xlb(2)=0;
xub(2)=6;
xplb(2)=0.1;
xpub(2)=1.5;

%location
xplb(3:5)=-0.2;
xpub(3:5)=0.2;

%popout
xplb(6)=-0.5;
xpub(6)=0.05;
xplb(7)=-0.05;
xpub(7)=0.5;

%previous color
xplb(8)=-0.05;
xpub(8)=0.5;

%Reset
xstart(9)=0.7;
xlb(9)=0;
xub(9)=1.2;
xplb(9)=0.3;
xpub(9)=0.9;

%Volatility
xstart(10)=0.05;
xlb(10)=0;
xub(10)=1;
xplb(10)=0.01;
xpub(10)=0.5;

%Color preference
xstart(11)=0;
xplb(11)=-1;
xpub(11)=1;
xlb(11)=-pi;
xub(11)=pi;

%priors
priors.x0 = xstart;
priors.xplb = xplb;
priors.xpub = xpub;
priors.xlb = xlb;
priors.xub = xub;


    function ll = wraper_f_QleanringFuncApproxReset_VBMC(P, u, y, in)
        [ll,~]=f_QleanringFuncApproxReset_RW_VBMC(P, u, y, in);
    end


%% Parameters estimates 1st round

fun_bads=@(x) - wraper_f_QleanringFuncApproxReset_VBMC(x, u, y, in);

%Note: bads minimise a function so we use -LL
opt_options = [];
opt_options.Display = 'all'; %'final'
priors.x0_bads = bads(fun_bads, priors.x0,priors.xlb,priors.xub,priors.xplb, priors.xpub,[],opt_options);

%% Parameters estimates 2nd round
%center the prioirs but keeps a wide search
priors.x0_bads(2)=abs(priors.x0_bads(2));
priors.x0_bads(11)=mod(priors.x0_bads(11)+pi,2*pi)-pi;
priors.x0_bads(12)=abs(priors.x0_bads(12));

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
priors.xlb_bads(10)=0;
priors.xub_bads(10)=1;
priors.xpub_bads(10)=0.9;
if priors.xplb_bads(10)<0.001
    priors.xplb_bads(10)=0.001;
end
if priors.x0_bads(10)<0.001
    priors.x0_bads(10)=0.01;
end

%Color preference
priors.xplb_bads(11)=priors.x0_bads(11)-pi/2;
priors.xpub_bads(11)=priors.x0_bads(11)+pi/2;
priors.xlb_bads(11)=priors.x0_bads(11)-pi;
priors.xub_bads(11)=priors.x0_bads(11)+pi;
   
fun_vbmc = @(x) wraper_f_QleanringFuncApproxReset_VBMC(x, u, y, in);

posterior.exitflag=0;
options = vbmc('defaults');
options.RetryMaxFunEvals = options.MaxFunEvals;
options.Display='all'; 

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
posterior.mean(11) = angle(mean(exp(1i.*posterior.Xs(:,11))));
posterior.var(11)=0;
posterior.mean(12) = mean(abs(posterior.Xs(:,12)));
posterior.var(12) = std(abs(posterior.Xs(:,12)));

[out.suffStat.LL, out.suffStat.gx]=f_QleanringFuncApproxReset_RW_VBMC(posterior.mean,u,y,in);

out.fit.BIC = N_param*log(length(y))-2*out.suffStat.LL;
out.fit.AIC = 2*N_param-2*out.suffStat.LL;


end


