
function [H1E] = regular_Heiviside_fun(X)
 
epsilon = 1; % the number is not fixed.
%Epsilon = 1.5; 
%%%%%%%%%%%% Create the first type of regulation of Heivide function %%%%%%%%%%%% 

     
H1E = 0.5*(1+(2/pi)*atan(X./epsilon));

%H2E = 1/2 .* (1 + 2./pi .* atan(X./Epsilon));
% H2E = 1/2 .* (1 + 2./pi .* atan(Z./Epsilon)); % calculate the second kind of regulation of Heiviside function

% figure  
% axis([-3 3 -0.1 1.1])
% plot(Z,H1E,'b-','LineWidth',2)
% 
% hold on 
% 
% plot(Z,H2E,'r--','LineWidth',2)
% 
% axis([-2 2 -0.1 1.1])
% xlabel('z')
% ylabel('H_{1,\epsilon}(z) or H_{2,\epsilon}(z)')
% title('Regularizations of H_{1,\epsilon}(z) and H_{2,\epsilon}(z)')
% legend('H_{1,\epsilon}(z)','H_{2,\epsilon}(z)',2)
% 
% hold off


end


