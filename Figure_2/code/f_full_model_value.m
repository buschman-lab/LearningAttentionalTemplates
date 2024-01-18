function [gx] = f_full_model_value(P,y,value)

value(isnan(value))=0;

gx = 1/length(y).*sum((y - (P(1) + P(2).*value)).^2) ;

end
