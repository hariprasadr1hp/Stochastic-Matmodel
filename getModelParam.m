% The function estimates the model parameters 'lamda' and 'p' of the 
% given outcome


function param = getModelParam(modPar,minFunc,noOfPoints, spacing,a1,a2,b)

intsity = modPar(1);
p = modPar(2);
L1 = factor(a1,b);
L2 = factor(a2,b);
L3 =pi*a1*b;
L4 = pi*a2*b;
 minContrast =0;

for i= 0: noOfPoints-1
    
    r = i* spacing;
    
    fac1 = 1- exp(-intsity*(r*r + (2*r/pi)*(p*L1 + (1-p)*L2) + p*L3 +(1-p)*L4));
     %fac1 = 1- exp(-lam(r*r ));
    
    minContrast = minContrast + ((fac1 - minFunc(i+1,2))*(fac1 - minFunc(i+1,2)));
    
    
end

param = minContrast;

end


function fac = factor(a,b)

w= (a-b)/(a+b);

fac = pi*(a+b)*(1+(3*w*w/(10+sqrt(4-(3*w*w)))));
end


