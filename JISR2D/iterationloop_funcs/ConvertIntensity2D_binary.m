function [ Ilabel ] = ConvertIntensity2D_binary( I0 )
%ConverIntensity convert intensity 0-255 to label 0-11
% input I0= 3 dimensional image
% Output Ilabel = Image with label 0-11
nx=size(I0,2);
ny=size(I0,1);

Ilabel=zeros(size(I0));
for i=1:nx
    for j=1:ny
        if I0(j,i)>=0 && I0(j,i)<=mean(I0(:))
            Ilabel(j,i)=0;
        end
        if I0(j,i)>mean(I0(:)) && I0(j,i)<=max(I0(:))
            Ilabel(j,i)=1;
        end
    end
end
end

