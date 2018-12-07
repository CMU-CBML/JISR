function [ I0_new, I0dx_new, I0dy_new, I0dz_new] = BsplineComposeImage3D_single( Vx, Vy, Vz, coef, nx, ny, nz)
%BsplineComposeImage2D compose the transformation field Vx, Vy, with the image (intensity) given in terms of a B-Spline level-set by coef
%input:  Vx, Vy = transformation field
%        coef   = coef of Bspline level set of the image
%        nx, ny = size of the image
%Output: I0_new = image after composition with transfomation field.
%Note: Vx, Vy should be between 0...nx, and 0..ny respectively

pU=3;
pV=3;
pW=3;

%number of knotspans in each direction
Vx = reshape(Vx, nx, ny, nz);
Vy = reshape(Vy, nx, ny, nz);
Vz = reshape(Vz, nx, ny, nz);
%coef = reshape(coef, ny+q, nx+p);

I0_new = zeros(size(Vx));
I0dx_new = zeros(size(Vx));
I0dy_new = zeros(size(Vx));
I0dz_new = zeros(size(Vx));

parfor k=1:nz
    for j = 1:ny
        for i=1:nx
            %[i, j]
            
            x_disp = Vx(i,j,k);
            y_disp = Vy(i,j,k);
            z_disp = Vz(i,j,k);
            
            %keep the displacement within the bounds of the image
            x_disp = max(0, x_disp);
            y_disp = max(0, y_disp);
            z_disp = max(0, z_disp);
            
            x_disp = min(x_disp, size(coef,1)-3-eps);
            y_disp = min(y_disp, size(coef,2)-3-eps);
            z_disp = min(z_disp, size(coef,3)-3-eps);
            
            k1 = ceil(x_disp-2+eps);
            l1 = ceil(y_disp-2+eps);
            m1 = ceil(z_disp-2+eps);
            
            x = x_disp-k1;
            y = y_disp-l1;
            z = z_disp-m1;

            basis_x = [(2-x)^3/6, 2/3-(x-1)^2+(x-1)^3/2, 2/3-(x-2)^2+(2-x)^3/2, (x-1)^3/6];
            basis_y = [(2-y)^3/6, 2/3-(y-1)^2+(y-1)^3/2, 2/3-(y-2)^2+(2-y)^3/2, (y-1)^3/6];
            basis_z = [(2-z)^3/6, 2/3-(z-1)^2+(z-1)^3/2, 2/3-(z-2)^2+(2-z)^3/2, (z-1)^3/6];
            
            basis_dx = [-0.5*(2-x)^2, -2*(x-1)+1.5*(x-1)^2, -2*(x-2)-1.5*(2-x)^2, 0.5*(x-1)^2];
            basis_dy = [-0.5*(2-y)^2, -2*(y-1)+1.5*(y-1)^2, -2*(y-2)-1.5*(2-y)^2, 0.5*(y-1)^2];
            basis_dz = [-0.5*(2-z)^2, -2*(z-1)+1.5*(z-1)^2, -2*(z-2)-1.5*(2-z)^2, 0.5*(z-1)^2];

            coef_loc= coef(l1+2:l1+5,k1+2:k1+5,m1+2:m1+5);
            
            for aa = 1:pU+1
                for bb = 1:pV+1
                    for cc = 1:pW+1
                        I0_new(i,j,k) = I0_new(i,j,k) + basis_y(1,bb)*basis_z(1,cc)*coef_loc(aa,bb,cc)*basis_x(1,aa);
                        I0dx_new(i,j,k) = I0dx_new(i,j,k) + basis_y(1,bb)*basis_z(1,cc)*coef_loc(aa,bb,cc)*basis_dx(1,aa);
                        I0dy_new(i,j,k) = I0dy_new(i,j,k) + basis_dy(1,bb)*basis_z(1,cc)*coef_loc(aa,bb,cc)*basis_x(1,aa);
                        I0dz_new(i,j,k) = I0dz_new(i,j,k) + basis_y(1,bb)*basis_dz(1,cc)*coef_loc(aa,bb,cc)*basis_x(1,aa);
                    end
                end
            end
        end
    end
end
