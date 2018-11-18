function [D1E] = Delta_fun(X)

epsilon = 1; 

D1E =(epsilon/pi)./(epsilon^2.+X.^2);

%D2E =  Epsilon ./(pi .* (Epsilon.^2 + X.^2) );
end