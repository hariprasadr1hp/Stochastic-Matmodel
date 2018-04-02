function [L] = function_L(a,b)
w1 = (a-b)/(a+b);
L = pi*(a+b)*(1+ (3*w1*w1)/(10+sqrt(4-3*w1*w1)));
end