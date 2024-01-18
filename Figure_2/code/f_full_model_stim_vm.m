function [gx] = f_full_model_stim_vm(P,y,stim)

stim_on=~isnan(stim);
stim(isnan(stim))=0;

kappa=10^(P(4));

gx = 1/length(y).*sum((y - ( P(1) + P(2).*stim_on.*vonmisespdf(stim, mod(P(3),2*pi), kappa) )).^2) ;

end

