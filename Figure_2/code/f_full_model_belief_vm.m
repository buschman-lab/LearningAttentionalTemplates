function [gx] = f_full_model_belief_vm(P,y,X)

b_prob=X;
b_angles=0:2*pi/size(X,2):2*pi-2*pi/size(X,2);

kappa=10^(P(4));

gx = 1/length(y).*sum((y - (P(1) + P(2).*sum(b_prob.*vonmisespdf (b_angles, mod(P(3),2*pi), kappa),2))).^2) ;

end

