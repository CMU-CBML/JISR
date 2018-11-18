function plotGrid(gridX, gridY)
%plots the grid 

ngx = size(gridX,1);
ngy = size(gridX,2);

%plot the vertical lines
for i=1:ngx
    lineVectorX = zeros(1,ngy);
    lineVectorY = zeros(1,ngy);    
    for j=1:ngy
        lineVectorX(j)=gridX(i,j);
        lineVectorY(j)=gridY(i,j);
    end
    line(lineVectorX,lineVectorY);
    hold on
end

%plot the horizontal lines
for i=1:ngy
    lineVectorX = zeros(1,ngx);
    lineVectorY = zeros(1,ngx);    
    for j=1:ngx
        lineVectorX(j)=gridX(j,i);
        lineVectorY(j)=gridY(j,i);
    end
    line(lineVectorX,lineVectorY);
    hold on
end

        
