function coef = img2coef3D(I, nx, ny, nz)
% img2coef_3D compute coeffient of B-spline level set representation for a
% given image I 

%input: I = 3 dimensional image
%      nx, ny, nz = size of image I in x, y, and z direction
% output: coef = coefficient of B-spline level set representation. 


I = reshape(I, ny, nx, nz);
coef = zeros(ny+3, nx+3, nz+3);


%filter in the x-direction
parfor  k=1:nz
        for i=1:ny
            coef(i,:,k) = bspline1dfilt(I(i,:,k));
        end
end
 
% filter in the y-direction
for k=1:nz
    for i=1:nx+3
    coef(:,i,k) = bspline1dfilt(coef(1:end-3,i,k)')';
    end
end

%filter in the z-direction
for j=1:ny+3
    for i=1:nx+3
        temp_coef=zeros(1,nz);
        for k=1:nz
            temp_coef(k)=coef(j,i,k);
        end
        coef(j,i,1:nz+3)=bspline1dfilt(temp_coef);
       
    end
end