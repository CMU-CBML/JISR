%Initialize the parametric space

sizeImage = [size(M,1),size(M,2),size(M,3)];
disp('Store the knotvectors...');
[CP, knotvectorU, knotvectorV, knotvectorW, uknotvectorU, uknotvectorV, uknotvectorW] = storeKnotArray(param,sizeImage,maxlev);

%% Compute Elem, Basis, Control point data structure for THBS splines
disp('Set B-spline grid...');

[Dm,Em] = setBsplineGrid_func(knotvectorU, knotvectorV, knotvectorW, uknotvectorU,uknotvectorV,uknotvectorW,param, maxlev);
%[Dm,Em] = setBsplineGrid_func_withoutNode(knotvectorU,knotvectorV,knotvectorW,param,maxlev);

%% Store level 1 control points
disp('Compute control points at level 1 using Greville coordinates...');
knotu = knotvectorU{1,1};
knotv = knotvectorV{1,1};
knotw = knotvectorW{1,1};

pp = zeros(param.nobU(1,1)*param.nobV(1,1)*param.nobW(1,1),3);
for k = 1:param.nobW(1,1)
    for j =1:param.nobV(1,1)
        for i = 1:param.nobU(1,1)
            index =(k-1)*param.nobV(1,1)*param.nobU(1,1) + (j-1)*param.nobU(1,1) + i;
            %Greville Abscissae
            coordx = sum(knotu(i+1:i+param.pU))./param.pU;
            coordy = sum(knotv(j+1:j+param.pV))./param.pV;
            coordz = sum(knotw(k+1:k+param.pW))./param.pW;
            pp(index,:) = [coordx,coordy,coordz];
        end
    end
end

%Control points stored for each refienement level
CP(1).pts = pp;
Pm = CP;
Pm_old = Pm;

