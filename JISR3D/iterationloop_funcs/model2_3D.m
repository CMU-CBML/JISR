function [SXY,DUSXY,DVSXY,DWSXY] = model2_3D(arf,beta,gma,uspan,vspan,wspan,dersU,dersV,dersW,pU,pV,pW)

arf_length = length(arf);
beta_length = length(beta);
gma_length = length(gma);

PX = zeros(gma_length,beta_length,arf_length);
PY = zeros(gma_length,beta_length,arf_length);
PZ = zeros(gma_length,beta_length,arf_length);

for i = 1:arf_length
    
    term1 = arf(i);
    
    for j = 1:beta_length
        
        term2 = beta(j);
        
        for k = 1:gma_length
            
            term3 = gma(k);
            
            PX(i,j,k) = term1;
            PY(i,j,k) = term2;
            PZ(i,j,k) = term3;
        end
    end
    
end

%uspan = FindSpan(arf_length-1,pU,plotpoints(iter,1),knotvectorU);
%vspan = FindSpan(beta_length-1,pV,plotpoints(iter,2),knotvectorV);

%dersU = Der1BasisFun(uspan,plotpoints(iter,1),pU,knotvectorU);
%dersV = Der1BasisFun(vspan,plotpoints(iter,2),pV,knotvectorV);

%dersU = BasisFun(uspan,plotpoints(iter,1),pU,knotvectorU);
%dersV = BasisFun(vspan,plotpoints(iter,2),pV,knotvectorV);

Nu = dersU(1,:);
Nv = dersV(1,:);
Nw = dersW(1,:);

DNu = dersU(2,:);
DNv = dersV(2,:);
DNw = dersW(2,:);

wind = wspan - pW;

S1 = 0;
S2 = 0;
S3 = 0;
S4 = 0;
S5 = 0;
S6 = 0;
S7 = 0;
S8 = 0;
S9 = 0;
S10 = 0;
S11 = 0;
S12 = 0;

for m = 1:pU+1
    
    temp1 = 0;
    temp2 = 0;
    temp3 = 0;
    temp4 = 0;
    temp5 = 0;
    temp6 = 0;
    temp7 = 0;
    temp8 = 0;
    temp9 = 0;
    temp10 = 0;
    temp11a = 0;
    temp12 = 0;
    
    uind = uspan - pU + m;
    
    for k = 1:pV+1
        
        temp11 = 0;
        temp21 = 0;
        temp31 = 0;
        temp41 = 0;
        temp51 = 0;
        temp61 = 0;
        temp71 = 0;
        temp81 = 0;
        temp91 = 0;
        temp101 = 0;
        temp111 = 0;
        temp121 = 0;
        
        vind = vspan-pV+k;
        
        for l = 1:pW+1
            
            temp11 = temp11 + Nw(l)* PX(uind,vind,wind+l);
            temp21 = temp21 + Nw(l)* PY(uind,vind,wind+l);
            temp31 = temp31 + Nw(l)* PZ(uind,vind,wind+l);
            
            temp41 = temp41 + Nw(l)* PX(uind,vind,wind+l);
            temp51 = temp51 + Nw(l)* PY(uind,vind,wind+l);
            temp61 = temp61 + Nw(l)* PZ(uind,vind,wind+l);
            
            temp71 = temp71 + Nw(l)* PX(uind,vind,wind+l);
            temp81 = temp81 + Nw(l)* PY(uind,vind,wind+l);
            temp91 = temp91 + Nw(l)* PZ(uind,vind,wind+l);
            
            temp101 = temp101 + DNw(l)* PX(uind,vind,wind+l);
            temp111 = temp111 + DNw(l)* PY(uind,vind,wind+l);
            temp121 = temp121 + DNw(l)* PZ(uind,vind,wind+l);
        end
        
        temp1 = temp1 + Nv(k)* temp11;
        temp2 = temp2 + Nv(k)* temp21;
        temp3 = temp3 + Nv(k)* temp31;
        
        temp4 = temp4 + Nv(k)* temp41;
        temp5 = temp5 + Nv(k)* temp51;
        temp6 = temp6 + Nv(k)* temp61;
        
        temp7 = temp7 + DNv(k)* temp71;
        temp8 = temp8 + DNv(k)* temp81;
        temp9 = temp9 + DNv(k)* temp91;
        
        temp10 = temp10 + Nv(k)* temp101;
        temp11a = temp11a + Nv(k)* temp111;
        temp12 = temp12 + Nv(k)* temp121;
    end
    
    S1 = S1 + Nu(m)* temp1;
    S2 = S2 + Nu(m)* temp2;
    S3 = S3 + Nu(m) * temp3;
    
    S4 = S4 + DNu(m)* temp4;
    S5 = S5 + DNu(m)* temp5;
    S6 = S6 + DNu(m)* temp6;
    
    S7 = S7 + Nu(m)* temp7;
    S8 = S8 + Nu(m)* temp8;
    S9 = S9 + Nu(m)* temp9;
    
    S10 = S10 + Nu(m)* temp10;
    S11 = S11 + Nu(m)* temp11a;
    S12 = S12 + Nu(m)* temp12;
end


SXY = [S1,S2,S3];
DUSXY = [S4,S5,S6];
DVSXY = [S7,S8,S9];
DWSXY = [S10,S11,S12];
end