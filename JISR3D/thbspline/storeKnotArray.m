function [ControlPointArray, KnotvectorU, KnotvectorV, KnotvectorW,uknot, vknot, wknot] = storeKnotArray(param,sizeImage,maxlevel)
% This function stores the knot vectors for each refinement level and the
% control points for each refinement level
%Input: param struct array storing all the parameters
% sizeImage: image size 
%Output: knotvectors in u, v, w parametric directions for all the
%refinement levels stored in cell arrays: KnotvectorU, KnotvectorV,
%KnotvectorW
%ControlPointArray stores the control points for each refinement level in
%struct array 

sx = sizeImage(1,1);
sy = sizeImage(1,2);
sz = sizeImage(1,3);
p = param.pU;
q = param.pV;
r = param.pW;
xelem = param.nelemU(1,1);
yelem = param.nelemV(1,1);
zelem = param.nelemW(1,1);

CP = cell(maxlevel,1); %initial control points
knotvectorU = cell(maxlevel,1);
knotvectorV = cell(maxlevel,1);
knotvectorW = cell(maxlevel,1);
uknot1 = cell(maxlevel,1);
vknot1 = cell(maxlevel,1);
wknot1 = cell(maxlevel,1);
for level = 1:maxlevel
    knotvectorU{level,1} = [0.*ones(1,p), 0:0.5^(level-1)*(sx)/xelem:sx, sx.*ones(1,p)];   
    knotvectorV{level,1} = [0.*ones(1,q), 0:0.5^(level-1)*(sy)/yelem:sy, sy.*ones(1,q)];
    knotvectorW{level,1} = [0.*ones(1,r), 0:0.5^(level-1)*(sz)/zelem:sz, sz.*ones(1,r)];
    uknot1{level,1} = unique(knotvectorU{level,1});
    vknot1{level,1} = unique(knotvectorV{level,1});
    wknot1{level,1} = unique(knotvectorW{level,1});
    px = zeros(param.nobU(level,1)*param.nobV(level,1)*param.nobW(level,1),3);
    CP{level,1} = px;
end

ControlPointArray = cell2struct(CP,'pts',2);
KnotvectorU = knotvectorU;
KnotvectorV = knotvectorV;
KnotvectorW = knotvectorW;
uknot = uknot1;
vknot = vknot1;
wknot = wknot1;
end