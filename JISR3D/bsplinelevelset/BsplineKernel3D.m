function [ kernel ] = BsplineKernel3D( )
%create 4x4x4 Bspline kernel 

b = [1/24, 11/24, 11/24, 1/24];

kernel=zeros(4,4,4);

for i=1:4
    for j=1:4
        for k=1:4
            kernel(j,i,k)=b(i)*b(j)*b(k);
        end
    end
end
            
end
