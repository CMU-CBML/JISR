function [ Ilabel ] = ConvertIntensity3D_binary( I0 )
%ConverIntensity convert intensity 0-255 to label 0-11
% input I0= 3 dimensional image
% Output Ilabel = Image with label 0-11
nx=size(I0,2);
ny=size(I0,1);
nz = size(I0,3);
mean1 = mean(I0(:));
max1 = max(I0(:));
Ilabel=zeros(size(I0));
for i=1:nx
    for j=1:ny
        for k=1:nz
            if (I0(j,i,k)>=0 && I0(j,i,k)<=0.5)
                Ilabel(j,i,k)=0;
            end
            if (I0(j,i,k)>0.5 && I0(j,i,k)<=1)
                Ilabel(j,i,k)=1;
            end
        end
    end
end
end

