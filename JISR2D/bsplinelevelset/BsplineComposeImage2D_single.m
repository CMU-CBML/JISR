function [ I0_new, I0dx_new, I0dy_new] = BsplineComposeImage2D_single( Vx, Vy, coef, nx, ny)
%BsplineComposeImage2D compose the transformation field Vx, Vy, with the image (intensity) given in terms of a B-Spline level-set by coef
%input:  Vx, Vy = transformation field
%        coef   = coef of Bspline level set of the image
%        nx, ny = size of the image
%Output: I0_new = image after composition with transfomation field.
%Note: Vx, Vy should be between 0...nx, and 0..ny respectively

% p=3;
% q=3;

%number of knotspans in each direction
Vx = reshape(Vx, nx, ny);
Vy = reshape(Vy, nx, ny);

%coef = reshape(coef, ny+q, nx+p);
I0_new = zeros(size(Vx));
I0dx_new = zeros(size(Vx));
I0dy_new = zeros(size(Vx));

parfor j=1:ny
    for i=1:nx
        %[i, j]
        
        x_disp = Vx(i,j);
        y_disp = Vy(i,j);
        
        
        %keep the displacement within the bounds of the image
        x_disp = max(0, x_disp);
        y_disp = max(0, y_disp);
        
        x_disp = min(x_disp, size(coef,1)-3-eps);
        y_disp = min(y_disp, size(coef,2)-3-eps);
        
        
        k1 = ceil(x_disp-2+eps);
        l1 = ceil(y_disp-2+eps);
        
        
        x = x_disp-k1;
        y = y_disp-l1;
        
        
        basis_x = [(2-x)^3/6, 2/3-(x-1)^2+(x-1)^3/2, 2/3-(x-2)^2+(2-x)^3/2, (x-1)^3/6];
        basis_y = [(2-y)^3/6, 2/3-(y-1)^2+(y-1)^3/2, 2/3-(y-2)^2+(2-y)^3/2, (y-1)^3/6];
        
        
        basis_dx = [-0.5*(2-x)^2, -2*(x-1)+1.5*(x-1)^2, -2*(x-2)-1.5*(2-x)^2, 0.5*(x-1)^2];
        basis_dy = [-0.5*(2-y)^2, -2*(y-1)+1.5*(y-1)^2, -2*(y-2)-1.5*(2-y)^2, 0.5*(y-1)^2];
        
        
        coef_loc= coef(l1+2:l1+5,k1+2:k1+5);
        
        I0_new(i,j) = basis_y*coef_loc*basis_x';
        I0dx_new(i,j) = basis_y*coef_loc*basis_dx';
        I0dy_new(i,j) = basis_dy*coef_loc*basis_x';
    end
end

