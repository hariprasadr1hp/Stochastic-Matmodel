function fac = factor(a,b)

w= (a-b)/(a+b);

fac = pi*(a+b)*(1+(3*w*w/(10+sqrt(4-(3*w*w)))));
end