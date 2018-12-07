function [ I0_new ] = BsplineComposeImage3D(nx, ny, nz,  VXL, VYL, VZL, CI0)
%BsplineComposeImage3D compose the transformation field VXL, VYL, VZL, with the image (intensity) given in terms of a B-Spline level-set by coef
%Note:VXL, VYL, VZL should be between 0...nx, 0..ny, 0..nz respectively

% input:   nx, ny, nz = size of image correspond to coefficent (CI0)
%         VXL,VYL,VZL = transformation field in x, y, and z direction respectively. 
%                 CI0 = coefficient of the image.

% output:      I0_new = new image given in terms of intensity after composed with the given transformation field.  

%degree of Bspline level set
p=3;
q=3;
r=3;

Vx=reshape(VXL, [ny nx nz]); 
Vy=reshape(VYL, [ny nx nz]);
Vz=reshape(VZL, [ny nx nz]);

coef=reshape(CI0,[ny+q nx+p nz+r]);

I0_new = zeros(size(Vx));


parfor j=1:ny
    for i=1:nx
        for k=1:nz
            %[i, j]
            x_disp = Vx(j,i,k);
            y_disp = Vy(j,i,k);
            z_disp = Vz(j,i,k);
            
            %keep the displacement within the bounds of the image
            x_disp = max(0, x_disp);
            y_disp = max(0, y_disp);
            z_disp = max(0, z_disp);
            
            x_disp = min(x_disp, nx-eps);
            y_disp = min(y_disp, ny-eps);
            z_disp = min(z_disp, nz-eps);
            
            k1 = ceil(x_disp-2+eps);
            l1 = ceil(y_disp-2+eps);
            m1 = ceil(z_disp-2+eps);
            
            x = x_disp-k1;
            y = y_disp-l1;
            z = z_disp-m1;
            
            basis_x = [(2-x)^3/6, 2/3-(x-1)^2+(x-1)^3/2, 2/3-(x-2)^2+(2-x)^3/2, (x-1)^3/6];
            basis_y = [(2-y)^3/6, 2/3-(y-1)^2+(y-1)^3/2, 2/3-(y-2)^2+(2-y)^3/2, (y-1)^3/6];
            basis_z = [(2-z)^3/6, 2/3-(z-1)^2+(z-1)^3/2, 2/3-(z-2)^2+(2-z)^3/2, (z-1)^3/6];
            
            coef_loc= coef(l1+2:l1+5,k1+2:k1+5,m1+2:m1+5);

            I0_temp = 0;
            for u=1:4
                for v=1:4
                    for w=1:4
                        I0_temp=I0_temp+coef_loc(v,u,w)*basis_x(u)*basis_y(v)*basis_z(w);
                        
                    end
                end
            end
            I0_new(j,i,k)=I0_temp;
                        

        end
    end
end

