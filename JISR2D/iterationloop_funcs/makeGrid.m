function [ gridX, gridY, pointsX, pointsY ] = makeGrid(ngx, ngy, nx, ny  )
%makes a grid of points

pointsX = round(linspace(1,nx,ngx));
pointsY = round(linspace(1,ny,ngy));

[gridX, gridY] = meshgrid(pointsX, pointsY);

end

