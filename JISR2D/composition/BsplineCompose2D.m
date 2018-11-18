function [ Vx_new, Vy_new ] = BsplineCompose2D( Vx, Vy, coef_x, coef_y, nx, ny)
%BsplineCompose2D compose the transformation field Vx, Vy, with the update field given by
%coef_x, coef_y.

% Input: Vx, Vy = transformation field from previous iteration 
%        coef_x, coef_y = update transformation field 
%        nx, ny = image size
% Output: Vx_new, Vy_new = composed transformation field

%Note: Vx, Vy should be between 0...nx, and 0..ny respectively


p=3;
q=3;

Vx = reshape(Vx, ny, nx);
Vy = reshape(Vy, ny, nx);
coef_x = reshape(coef_x, ny+q, nx+p);
coef_y = reshape(coef_y, ny+q, nx+p);

Vx_new = zeros(size(Vx));
Vy_new = zeros(size(Vy));

[Vx_ident, Vy_ident] = meshgrid((0.5:nx-0.5), (0.5:ny-0.5));

parfor j=1:ny
    for i=1:nx
        %[i, j]
        x_disp = Vx(j,i);
        y_disp = Vy(j,i);
        
        %keep the displacement within the bounds of the image
        x_disp = max(0, x_disp);
        y_disp = max(0, y_disp);
        
        x_disp = min(x_disp, nx-eps);
        y_disp = min(y_disp, ny-eps);
        
        k1 = ceil(x_disp-2+eps);
        l1 = ceil(y_disp-2+eps);
        
        x = x_disp-k1;
        y = y_disp-l1;
        
        basis_x = [(2-x)^3/6, 2/3-(x-1)^2+(x-1)^3/2, 2/3-(x-2)^2+(2-x)^3/2, (x-1)^3/6];
        basis_y = [(2-y)^3/6, 2/3-(y-1)^2+(y-1)^3/2, 2/3-(y-2)^2+(2-y)^3/2, (y-1)^3/6];

        %basis_y = basis_y(end:-1:1);
        coef_loc_x = coef_x(l1+2:l1+5,k1+2:k1+5);
        coef_loc_y = coef_y(l1+2:l1+5,k1+2:k1+5);
        
        Vx_new(j,i) = basis_y*coef_loc_x*basis_x';
        Vy_new(j,i) = basis_y*coef_loc_y*basis_x';
  
        
    end
end

Vx_new(:,1:3)=Vx_ident(:,1:3);
Vx_new(:,nx-2:nx)=Vx_ident(:,nx-2:nx);
Vy_new(1:3,:)=Vy_ident(1:3,:);
Vy_new(ny-2:ny,:)=Vy_ident(ny-2:ny,:);
