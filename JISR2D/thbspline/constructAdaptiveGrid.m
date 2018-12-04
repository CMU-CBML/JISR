function [Jm,Coeff,Pixel,HH,PHI,PHIU,PHIV,BIGX,BIGY] = constructAdaptiveGrid(ac,parameters,Dm,Em,M,knotvectorU,knotvectorV,multilev,nobU,nobV,nelemU)

ac_ct = size(ac,1);
maxlev=size(knotvectorU,1);
Jm = cell(ac_ct,1);
orderGauss = parameters.orderGauss;
pU = parameters.pU;
pV = parameters.pV;
[Vx_ident, Vy_ident] = meshgrid((0.5:size(M,1)-0.5),(0.5:size(M,2)-0.5));

for i = 1:ac_ct,
    
    cell_ind = ac(i,1);
    cell_lev = ac(i,2);
    counter = 0;
    DE = Dm{cell_lev,1};
    EE = Em{cell_lev,1};
    
    %Collect the non-zeros HB splines over each active cell
    local_b = EE{cell_ind,6};
    lb_size = size(local_b,2);
    cell_supp = zeros(1,2);
    
    for j = 1:lb_size,
        if(DE{local_b(1,j),3}==1),
            counter= counter+1;
            cell_supp(counter,1) = local_b(1,j);
            cell_supp(counter,2) = cell_lev;
        end
    end
    curr_cell = cell_ind;
    
    if(cell_lev~=1),
        for m = (cell_lev-1):-1:1,
            EEM = Em{m,1};
            EEMM = Em{m+1,1};
            DEM = Dm{m,1};
            
            curr_cell = EEMM{curr_cell,10};
            local_supp = EEM{curr_cell,6};
            local_supp_size = size(local_supp,2);
            for j = 1:local_supp_size,
                if(DEM{local_supp(1,j),3}==1),
                    counter=counter+1;
                    cell_supp(counter,1) = local_supp(1,j);
                    cell_supp(counter,2) = m;
                end
            end
        end
    end
    Jm{i,1} = cell_supp;
end
%==========================================================================
Coeff = cell(ac_ct,1);
for  i = 1:ac_ct,
    count=0;
    cell_ind = ac(i,1);
    cell_lev = ac(i,2);
    ctemp = cell(maxlev,1);
    EE = Em{cell_lev,1};
    DE = Dm{cell_lev,1};
    local_b = EE{cell_ind,6};
    lb_size = size(local_b,2);
    ct = zeros(16,16);
    ct11 = zeros(16,16);
    coef_arr = zeros(1,16);
    for j = 1:lb_size,
        ct(j,j) = 1;
    end
    ctemp{cell_lev,1} = ct;
    for j =1:lb_size,
        if(DE{local_b(1,j),3}==1),
            count= count+1;
            coef_arr(count,:) = ct(j,:);
        end
    end
    curr_cell = cell_ind;
    curr_cell_level = cell_lev;
    if(cell_lev~=1),
        for m=(cell_lev-1):-1:1,
            EEM = Em{m,1};
            DEM = Dm{m,1};
            DEMM = Dm{m+1,1};
            EEMM = Em{m+1,1};
            local_s1 = EEMM{curr_cell,6};
            curr_cell = EEMM{curr_cell,10};
            local_s = EEM{curr_cell,6};
            ls_size = size(local_s,2);
            ct = zeros(16,16);
            ct11 = zeros(16,16);
            for j=1:ls_size,
                cmat = DEM{local_s(1,j),12};
                for k=1:ls_size,
                    for k1=1:size(cmat,1),
                        if(cmat(k1,2)==local_s1(1,k)),
                            ct(j,k) = cmat(k1,1);
                            if(DEMM{local_s1(1,k),3}==0&&DEMM{local_s1(1,k),7}==0),
                                ct11(j,k) = ct(j,k);
                            end
                        end
                    end
                end
            end
            ct1 = ctemp{m+1,1};
            ct = ct*ct1;
            ctnew = ct11*ct1;
            ctemp{m,1} = ct;
            for j=1:ls_size,
                if(DEM{local_s(1,j),3}==1),
                    count = count+1;
                    coef_arr(count,:) = ctnew(j,:);
                end
            end
        end
    end
    Coeff{i,1} = coef_arr;
end

Pixel = cell(size(M,1)*size(M,2),2);
kuMAX = knotvectorU{multilev+1,1};
kvMAX = knotvectorV{multilev+1,1};
unobMAX = nobU(multilev+1,1);
vnobMAX = nobV(multilev+1,1);
EE = Em{multilev+1,1};
p_ind = 0;
for i = 1:size(M,2),
    for j = 1:size(M,1),
        p_ind = p_ind+1;
        uu = FindSpan(unobMAX-1,pU,Vx_ident(j,i),kuMAX) + 1;
        vv = FindSpan(vnobMAX-1,pV,Vy_ident(j,i),kvMAX) + 1;
        
        cellx = (uu-pU);
        celly = (vv-pV);
        cell_ind = nelemU(multilev+1,1)*(celly-1)+cellx;
        curr_cell = cell_ind;
        if(EE{cell_ind,4} == 1),
            
            act_ind = EE{cell_ind,11};
            SB = Jm{act_ind,1};
            sb_size = size(SB,1);
            pix_coeff = Coeff{act_ind,1};
            SB1 = EE{cell_ind,6};
            sb_size1 = size(SB1,1);
            RR1 = BasisFun((uu-1),Vx_ident(j,i),pU,kuMAX);
            RR2 = BasisFun((vv-1),Vy_ident(j,i),pV,kvMAX);
            RR = (RR1')*(RR2);
            inc=0;
            phii = zeros(16,1);
            for m1 = 1:size(RR,1),
                for m2=1:size(RR,2),
                    phii(16-inc,1)= RR(m2,m1);
                    inc=inc+1;
                end
            end
            phi_pi = pix_coeff*phii;
            Pixel{p_ind,1} = act_ind;
            Pixel{p_ind,2} = phi_pi;
            
        else
            
            for m = (multilev+1):-1:2,
                EEM = Em{m,1};
                EEMM = Em{(m-1),1};
                curr_cell = EEM{curr_cell,10};
                if(EEMM{curr_cell,4} == 1),
                    
                    act_ind = EEMM{curr_cell,11};
                    SB = Jm{act_ind,1};
                    sb_size = size(SB,1);
                    knotuu = knotvectorU{(m-1),1};
                    knotvv = knotvectorV{(m-1),1};
                    unob1 = nobU(m-1,1);
                    vnob1 = nobV(m-1,1);
                    u1 = FindSpan(unob1-1,pU,Vx_ident(j,i),knotuu) + 1;
                    v1 = FindSpan(vnob1-1,pV,Vy_ident(j,i),knotvv) + 1;
                    pix_coeff = Coeff{act_ind,1};
                    SB1 = EE{cell_ind,6};
                    sb_size1 = size(SB1,1);
                    RR1 = BasisFun((u1-1),Vx_ident(j,i),pU,knotuu);
                    RR2 = BasisFun((v1-1),Vy_ident(j,i),pV,knotvv);
                    RR = (RR1')*(RR2);
                    inc=0;
                    phii = zeros(16,1);
                    for m1 = 1:size(RR,1),
                        for m2=1:size(RR,2),
                            phii(16-inc,1)= RR(m2,m1);
                            inc=inc+1;
                        end
                    end
                    phi_pi = pix_coeff*phii;
                    
                    Pixel{p_ind,1} = act_ind;
                    Pixel{p_ind,2} = phi_pi;
                    break;
                end
            end
        end
    end
end

xlen = 6;
HH = zeros(ac_ct,2);
PHI =  cell(ac_ct,1);
PHIU = cell(ac_ct,1);
PHIV = cell(ac_ct,1);
BIGX = zeros(ac_ct*xlen,xlen);
BIGY = zeros(ac_ct*xlen,xlen);
% % %Starting the image registration part now....
for i = 1:ac_ct,
    
    cell_index = ac(i,1);
    cell_level = ac(i,2);
    
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
    
    a1 = knotU(1,u(1,1));
    b1 = knotU(1,u(2,1));
    h1 = (b1-a1)/2;
    m1 = (a1+b1)/2;
    
    a2 = knotV(1,v(1,1));
    b2 = knotV(1,v(2,1));
    h2 = (b2-a2)/2;
    m2 = (a2+b2)/2;
    
    HH(i,:) = [h1,h2];
    s1 = s1*h1 + m1;
    s2 = s2*h2 + m2;
    
    SB1 = EE{cell_index,6};
    supp_size1 = size(SB1,1);
    [GX,GY] = meshgrid(s1,s2);
    BIGX(1+(i-1)*xlen:i*xlen,1:6) = GX;
    BIGY(1+(i-1)*xlen:i*xlen,1:6) = GY;
    
    g_phi = zeros(supp_size,6,6);
    g_phiu = zeros(supp_size,6,6);
    g_phiv = zeros(supp_size,6,6);
    
    
    for gg1 = 1:xlen,
        for gg2 = 1:xlen,
            uu = GX(gg2,gg1);
            vv = GY(gg2,gg1);
            uknot = FindSpan(unobg-1,pU,uu,knotU) + 1;
            vknot = FindSpan(vnobg-1,pV,vv,knotV) + 1;
            RRD1 = Der1BasisFun(uknot-1,uu,pU,knotU);
            RRD2 = Der1BasisFun(vknot-1,vv,pV,knotV);
            RRD = (RRD1(1,:)')*(RRD2(1,:));
            RRDU = (RRD1(1,:)')*(RRD2(2,:));
            RRDV = (RRD1(2,:)')*(RRD2(1,:));
            inc=0;
            phii = zeros(16,1);
            phiiu = zeros(16,1);
            phiiv = zeros(16,1);
            for m1 = 1:size(RR,1),
                for m2=1:size(RR,2),
                    phii(16-inc,1)= RRD(m2,m1);
                    phiiu(16-inc,1)= RRDU(m2,m1);
                    phiiv(16-inc,1)= RRDV(m2,m1);
                    inc=inc+1;
                end
            end
            phi_pi = gg_coeff*phii;
            phi_piu =gg_coeff*phiiu;
            phi_piv = gg_coeff*phiiv;
            g_phi(:,gg2,gg1) = phi_pi;
            g_phiu(:,gg2,gg1) = phi_piu;
            g_phiv(:,gg2,gg1) = phi_piv;
        end
    end
    
    PHI{i,1} = g_phi;
    PHIU{i,1} = g_phiu;
    PHIV{i,1} = g_phiv;
end
end

