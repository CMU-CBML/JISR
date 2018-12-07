filename = sprintf('Sphere_Cube_mesh%d.vtk',multil+1);
fileID = fopen(filename,'w');
fprintf(fileID, '# vtk DataFile Version 2.0\n');
fprintf(fileID, 'mesh\n');
fprintf(fileID, 'ASCII\n');
fprintf(fileID, 'DATASET UNSTRUCTURED_GRID\n\n');
fprintf(fileID, 'POINTS %d float\n',size(ActiveNodes,1));

for i = 1:size(ActiveNodes,1),
    fprintf(fileID, '%f %f %f\n', Node(ActiveNodes(i,1),1), Node(ActiveNodes(i,1),2), Node(ActiveNodes(i,1),3));
end

fprintf(fileID, 'CELLS  %d %d\n',ac_ct,9*ac_ct);

for i = 1:ac_ct,
    nodes = Em(ac(i,2)).nodes(ac(i,1),:);
    fprintf(fileID, '8 %d %d %d %d %d %d %d %d\n',Node(nodes(1,1),4)-1,Node(nodes(1,2),4)-1,Node(nodes(1,3),4)-1,Node(nodes(1,4),4)-1,Node(nodes(1,5),4)-1,Node(nodes(1,6),4)-1,Node(nodes(1,7),4)-1,Node(nodes(1,8),4)-1);
end

fprintf(fileID, 'CELL_TYPES %d\n',ac_ct);
for i = 1:ac_ct,
    fprintf(fileID, '12 \n');
end

fclose(fileID);