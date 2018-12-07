% function [D1E] = Delta_fun3D(X)
% 
% Epsilon = (max(X(:))-min(X(:)))*0.5; % the number is not fixed.
% %Epsilon = 0.1;
% 
% N = size(X,1);
% M = size(X,2);
% O = size(X,3);
% 
% D1E = zeros(N,M,O);
% 
% for i = 1:N
%     
%     for j = 1:M
%         
%         for k = 1:O
%             if (abs(X(i,j,k)) <= Epsilon)
%                 
%                 D1E(i,j,k) = 1/(2*Epsilon) * (1 + cos((pi * X(i,j,k))/Epsilon));
%                 
%             end
%         end
%     end
%     
% end
% 
% %D2E =  Epsilon ./(pi .* (Epsilon.^2 + X.^2) );
% 
% 
% end

function [D1E] = Delta_fun3D(X)

epsilon = 1; 
D11 = zeros(size(X));

for k = 1:size(X,3)
    for j = 1:size(X,2)
        for i = 1:size(X,1)
D11(i,j,k) =(epsilon/pi)./(epsilon^2.+X(i,j,k).^2);
        end
    end
end

D1E = D11;

%D2E =  Epsilon ./(pi .* (Epsilon.^2 + X.^2) );


end