function Jdet_vector = computeJacobian(parameters,ac,Coeff,Em,Jm,knotvectorU,knotvectorV,Pmf1)
pU = parameters.pU;
pV = parameters.pV;
ac_ct = size(ac,1);
Jdet = zeros(ac_ct,1);
for i = 1:ac_ct
    cell_ind = ac(i,1);
    cell_lev = ac(i,2);
    c1 = Coeff{i,1};
    EE = Em{cell_lev,1}; %inputing cells at a level
    %CE = Pm{cell_lev,1};  %inputing control points at a level
    %SB = EE{cell_ind,6};%inputing support bsplines for an active cell
    %supp_size = size(SB,1);
    SB1 = Jm{i,1};
    supp_size1 = size(SB1,1);
    knotU1 = knotvectorU{cell_lev,1};
    knotV1 = knotvectorV{cell_lev,1};
    
    u = EE{cell_ind,2};
    v = EE{cell_ind,3};
    
    x1 = knotU1(1,u(1,1));
    x2 = knotU1(1,u(2,1));
    y1 = knotV1(1,v(1,1));
    y2 = knotV1(1,v(2,1));
    
    fx(1,1) = x1;
    fx(1,2) = y1;
    fx(2,1) = x2;
    fx(2,2) = y1;
    fx(3,1) = x2;
    fx(3,2) = y2;
    fx(4,1) = x1;
    fx(4,2) = y2;
    
    fmidx = x1+0.5*(-x1+x2);
    fmidy = y1+0.5*(-y1+y2);
    
    
    uu = fmidx;
    vv = fmidy;
    
    sumbxu = 0;
    sumbyu = 0;
    sumbxv = 0;
    sumbyv = 0;
    sumbx = 0;
    sumby = 0;
    
    RRD1 = Der1BasisFun(u(1,1)-1,uu,pU,knotU1);
    RRD2 = Der1BasisFun(v(1,1)-1,vv,pV,knotV1);
    RRDU = (RRD1(1,:)')*(RRD2(2,:));
    RRDV = (RRD1(2,:)')*(RRD2(1,:));
    RR =   (RRD1(1,:)')*(RRD2(1,:));
    inc=0;
    phiiu=zeros(16,1);
    phiiv=zeros(16,1);
    phii = zeros(16,1);
    for m1=1:size(RR,1)
        for m2=1:size(RR,2)
            phiiu(16-inc,1)= RRDU(m2,m1);
            phiiv(16-inc,1)= RRDV(m2,m1);
            phii(16-inc,1) = RR(m2,m1);
            inc = inc+1;
        end
    end
    PHIUC = c1*phiiu;
    PHIVC = c1*phiiv;
    PHIIC = c1*phii;
    for k = 1:supp_size1
        CEb = Pmf1{SB1(k,2),1};
        pi = CEb(SB1(k,1),1);
        pj = CEb(SB1(k,1),2);
        sumbxu = sumbxu + pi*PHIUC(k,1);
        sumbyu = sumbyu + pj*PHIUC(k,1);
        sumbxv = sumbxv + pi*PHIVC(k,1);
        sumbyv = sumbyv + pj*PHIVC(k,1);
        sumbx = sumbx + pi*PHIIC(k,1);
        sumby = sumby + pj*PHIIC(k,1);
    end
    J = [sumbxu,sumbxv;sumbyu,sumbyv];
    Jdet(i,1) = det(J);
end
Jdet_vector = Jdet;
end