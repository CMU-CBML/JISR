function [ norm ] = norm3d( X )
% square norm of a 3D matrix array

norm = sqrt(sum(sum(sum(X.^2))));

end