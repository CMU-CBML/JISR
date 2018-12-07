function [H1E] = regular_Heiviside_fun3D(X)
 
epsilon = 1; % the number is not fixed.
%Epsilon = 1.5; 
%%%%%%%%%%%% Create the first type of regulation of Heivide function %%%%%%%%%%%%  
H1E = 0.5*(1+(2/pi)*atan(X./epsilon));

end


