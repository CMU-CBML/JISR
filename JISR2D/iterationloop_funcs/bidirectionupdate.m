function [FBXX,BFXX,PHIF,PHIB,PHIUF,PHIVF,PHIUB,PHIVB]  = bidirectionupdate(ac, parameters, Jm, Coeff, Em, knotvectorU, knotvectorV, nobU, nobV, Pmf1,Pmb1,BIGXXfb,BIGYYfb,BIGXXbf,BIGYYbf)

ac_ct = size(Jm,1);
orderGauss = 6;
xlen = 6;
pU = parameters.pU;
pV = parameters.pV;

PHIF =  cell(ac_ct,1);
PHIB =  cell(ac_ct,1);
PHIUF = cell(ac_ct,1);
PHIVF = cell(ac_ct,1);
PHIUB = cell(ac_ct,1);
PHIVB = cell(ac_ct,1);

FBX = zeros(ac_ct*xlen,xlen);
FBY = zeros(ac_ct*xlen,xlen);
FBXDX = zeros(ac_ct*xlen,xlen);
FBXDY = zeros(ac_ct*xlen,xlen);
FBYDX = zeros(ac_ct*xlen,xlen);
FBYDY = zeros(ac_ct*xlen,xlen);

BFX = zeros(ac_ct*xlen,xlen);
BFY = zeros(ac_ct*xlen,xlen);
BFXDX = zeros(ac_ct*xlen,xlen);
BFXDY = zeros(ac_ct*xlen,xlen);
BFYDX = zeros(ac_ct*xlen,xlen);
BFYDY = zeros(ac_ct*xlen,xlen);

FBXX = zeros(ac_ct*xlen,xlen,6);
BFXX = zeros(ac_ct*xlen,xlen,6);

% % %Starting the image registration part now....
for i = 1:ac_ct,
    
    cell_index = ac(i,1);
    cell_level = ac(i,2);
    
    SBXf = BIGXXfb(1+(i-1)*xlen:i*xlen,1:6);
    SBYf = BIGYYfb(1+(i-1)*xlen:i*xlen,1:6);
    SBXb = BIGXXbf(1+(i-1)*xlen:i*xlen,1:6);
    SBYb = BIGYYbf(1+(i-1)*xlen:i*xlen,1:6);
    
    SB = Jm{i,1};
    supp_size = size(SB,1);
    gg_coeff = Coeff{i,1};
    
    [s1,w1] = ggquad(orderGauss);
    [s2,w2] = ggquad(orderGauss);
    xlen = size(s1,1);
    
    EE = Em{cell_level,1}; %inputing cells at a level
    knotU = knotvectorU{cell_level,1};
    knotV = knotvectorV{cell_level,1};
    unobg = nobU(cell_level,1);
    vnobg = nobV(cell_level,1);
    u = EE{cell_index,2};
    v = EE{cell_index,3};
    
    g_phif = zeros(supp_size,6,6);
    g_phiuf = zeros(supp_size,6,6);
    g_phivf = zeros(supp_size,6,6);
    
    g_phib = zeros(supp_size,6,6);
    g_phiub = zeros(supp_size,6,6);
    g_phivb = zeros(supp_size,6,6);
    
    for gg1 = 1:xlen,
        for gg2 = 1:xlen,
            uuf = SBXf(gg2,gg1);
            vvf = SBYf(gg2,gg1);
            uub = SBXb(gg2,gg1);
            vvb = SBYb(gg2,gg1);
            
            uknotf = FindSpan(unobg-1,pU,uuf,knotU) + 1;
            vknotf = FindSpan(vnobg-1,pV,vvf,knotV) + 1;
            uknotb = FindSpan(unobg-1,pU,uub,knotU) + 1;
            vknotb = FindSpan(vnobg-1,pV,vvb,knotV) + 1;
            
            RRD1f = Der1BasisFun(uknotf-1,uuf,pU,knotU);
            RRD2f = Der1BasisFun(vknotf-1,vvf,pV,knotV);
            RRD1b = Der1BasisFun(uknotb-1,uub,pU,knotU);
            RRD2b = Der1BasisFun(vknotb-1,vvb,pV,knotV);
            
            RRDf = (RRD1f(1,:)')*(RRD2f(1,:));
            RRDUf = (RRD1f(1,:)')*(RRD2f(2,:));
            RRDVf = (RRD1f(2,:)')*(RRD2f(1,:));
            
            RRDb = (RRD1b(1,:)')*(RRD2b(1,:));
            RRDUb = (RRD1b(1,:)')*(RRD2b(2,:));
            RRDVb = (RRD1b(2,:)')*(RRD2b(1,:));
            
            inc=0;
            phiif = zeros(16,1);
            phiiuf = zeros(16,1);
            phiivf = zeros(16,1);
            
            phiib = zeros(16,1);
            phiiub = zeros(16,1);
            phiivb = zeros(16,1);
            
            for m1 = 1:size(RRDf,1),
                for m2=1:size(RRDf,2),
                    phiif(16-inc,1)= RRDf(m2,m1);
                    phiiuf(16-inc,1)= RRDUf(m2,m1);
                    phiivf(16-inc,1)= RRDVf(m2,m1);
                    
                    phiib(16-inc,1)= RRDb(m2,m1);
                    phiiub(16-inc,1)= RRDUb(m2,m1);
                    phiivb(16-inc,1)= RRDVb(m2,m1);
                    inc=inc+1;
                end
            end
            
            phi_pif = gg_coeff*phiif;
            phi_piuf =gg_coeff*phiiuf;
            phi_pivf = gg_coeff*phiivf;
            
            phi_pib = gg_coeff*phiib;
            phi_piub =gg_coeff*phiiub;
            phi_pivb = gg_coeff*phiivb;
            
            g_phif(:,gg2,gg1) = phi_pif;
            g_phiuf(:,gg2,gg1) = phi_piuf;
            g_phivf(:,gg2,gg1) = phi_pivf;
            
            g_phib(:,gg2,gg1) = phi_pib;
            g_phiub(:,gg2,gg1) = phi_piub;
            g_phivb(:,gg2,gg1) = phi_pivb;
        end
    end
    
    PHIF{i,1} = g_phif;
    PHIUF{i,1} = g_phiuf;
    PHIVF{i,1} = g_phivf;
    
    PHIB{i,1} = g_phib;
    PHIUB{i,1} = g_phiub;
    PHIVB{i,1} = g_phivb;
    
    %Now to compute the vectors, fx, fux, fvx
    SBXf = zeros(xlen,xlen);
    SBYf = zeros(xlen,xlen);
    SBUXf = zeros(xlen,xlen);
    SBUYf = zeros(xlen,xlen);
    SBVXf = zeros(xlen,xlen);
    SBVYf = zeros(xlen,xlen);
    
    SBXb = zeros(xlen,xlen);
    SBYb = zeros(xlen,xlen);
    SBUXb = zeros(xlen,xlen);
    SBUYb = zeros(xlen,xlen);
    SBVXb = zeros(xlen,xlen);
    SBVYb = zeros(xlen,xlen);
    
    supp_phif = PHIF{i,1};
    supp_phiuf = PHIUF{i,1};
    supp_phivf =PHIVF{i,1};
    
    supp_phib = PHIB{i,1};
    supp_phiub = PHIUB{i,1};
    supp_phivb =PHIVB{i,1};
    
    for gg1 = 1:xlen
        for gg2 = 1:xlen
            
            sumbxf = 0;
            sumbyf = 0;
            sumbuxf = 0;
            sumbuyf = 0;
            sumbvxf = 0;
            sumbvyf = 0;
            
            sumbxb = 0;
            sumbyb = 0;
            sumbuxb = 0;
            sumbuyb = 0;
            sumbvxb = 0;
            sumbvyb = 0;
            
            for kg = 1:supp_size
                CEbf = Pmf1{SB(kg,2),1};
                CEbb = Pmb1{SB(kg,2),1};
                
                pif = CEbf(SB(kg,1),1);
                pjf = CEbf(SB(kg,1),2);
                
                pib = CEbb(SB(kg,1),1);
                pjb = CEbb(SB(kg,1),2);
                
                sumbxf = sumbxf + pif*supp_phif(kg,gg2,gg1);
                sumbyf = sumbyf + pjf*supp_phif(kg,gg2,gg1);
                sumbuxf = sumbuxf + pif*supp_phiuf(kg,gg2,gg1);
                sumbuyf = sumbuyf + pjf*supp_phiuf(kg,gg2,gg1);
                sumbvxf = sumbvxf + pif*supp_phivf(kg,gg2,gg1);
                sumbvyf = sumbvyf + pjf*supp_phivf(kg,gg2,gg1);
                
                sumbxb = sumbxb + pib*supp_phib(kg,gg2,gg1);
                sumbyb = sumbyb + pjb*supp_phib(kg,gg2,gg1);
                sumbuxb = sumbuxb + pib*supp_phiub(kg,gg2,gg1);
                sumbuyb = sumbuyb + pjb*supp_phiub(kg,gg2,gg1);
                sumbvxb = sumbvxb + pib*supp_phivb(kg,gg2,gg1);
                sumbvyb = sumbvyb + pjb*supp_phivb(kg,gg2,gg1);
                
            end
            
            SBXf(gg2,gg1) = sumbxf;
            SBYf(gg2,gg1)  = sumbyf;
            SBUXf(gg2,gg1)  = sumbuxf;
            SBUYf(gg2,gg1)  = sumbuyf;
            SBVXf(gg2,gg1)  = sumbvxf;
            SBVYf(gg2,gg1)  = sumbvyf;
            
            SBXb(gg2,gg1)  = sumbxb;
            SBYb(gg2,gg1)  = sumbyb;
            SBUXb(gg2,gg1)  = sumbuxb;
            SBUYb(gg2,gg1)  = sumbuyb;
            SBVXb(gg2,gg1)  = sumbvxb;
            SBVYb(gg2,gg1)  = sumbvyb;
        end
    end
    
    FBX(1+(i-1)*xlen:i*xlen,1:6) = SBXf;
    FBY(1+(i-1)*xlen:i*xlen,1:6) = SBYf;
    FBXDX(1+(i-1)*xlen:i*xlen,1:6) = SBUXf;
    FBXDY(1+(i-1)*xlen:i*xlen,1:6) = SBUYf;
    FBYDX(1+(i-1)*xlen:i*xlen,1:6) = SBVXf;
    FBYDY(1+(i-1)*xlen:i*xlen,1:6) = SBVYf;
    
    BFX(1+(i-1)*xlen:i*xlen,1:6) = SBXb;
    BFY(1+(i-1)*xlen:i*xlen,1:6) = SBYb;
    BFXDX(1+(i-1)*xlen:i*xlen,1:6) = SBUXb;
    BFXDY(1+(i-1)*xlen:i*xlen,1:6) = SBUYb;
    BFYDX(1+(i-1)*xlen:i*xlen,1:6) = SBVXb;
    BFYDY(1+(i-1)*xlen:i*xlen,1:6) = SBVYb;
    
end

FBXX(:,:,1) = FBX;
FBXX(:,:,2) = FBY;
FBXX(:,:,3) = FBXDX;
FBXX(:,:,4) = FBXDY;
FBXX(:,:,5) = FBYDX;
FBXX(:,:,6) = FBYDY;

BFXX(:,:,1) = BFX;
BFXX(:,:,2) = BFY;
BFXX(:,:,3) = BFXDX;
BFXX(:,:,4) = BFXDY;
BFXX(:,:,5) = BFYDX;
BFXX(:,:,6) = BFYDY;
end