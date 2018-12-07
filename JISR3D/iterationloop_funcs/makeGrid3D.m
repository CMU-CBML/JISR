function [ gridX, gridY, gridZ, pointsX, pointsY,pointsZ ] = makeGrid3D(ngx, ngy,ngz, nx, ny ,nz )
%makes a grid of points

pointsX = round(linspace(1,nx,ngx));
pointsY = round(linspace(1,ny,ngy));
pointsZ = round(linspace(1,nz,ngz));

[gridX, gridY, gridZ] = meshgrid(pointsX, pointsY, pointsZ);

end