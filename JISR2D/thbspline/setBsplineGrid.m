function [Dm,Pm,Em,Bvect,knotvectorU,knotvectorV,nobU,nobV,nelemU] = setBsplineGrid(maxlev,parameters,F)

CP = cell(maxlev,1);

knotvectorU = cell(maxlev,1);
nelemU = zeros(maxlev,1);
uknotvectorU = cell(maxlev,1);
kU = zeros(maxlev,1);
nobU = zeros(maxlev,1);

knotvectorV = cell(maxlev,1);
nelemV = zeros(maxlev,1);
uknotvectorV = cell(maxlev,1);
kV = zeros(maxlev,1);
nobV = zeros(maxlev,1);
pU = parameters.pU;
pV = parameters.pV;
nelemx = parameters.nelemx;
nelemy = parameters.nelemy;

%% Knotvectors
for lev = 1:maxlev,
    knotvectorU{lev,1} = [0.*ones(1,pU), 0:0.5^(lev-1)*size(F,2)/nelemx:size(F,2), size(F,2).*ones(1,pU)];
    uknotvectorU{lev,1} = unique(knotvectorU{lev,1});
    nelemU(lev,1) = size(uknotvectorU{lev,1},2) -1;
    kU(lev,1) = length(knotvectorU{lev,1});
    nobU(lev,1) = kU(lev,1) - pU -1;
    
    knotvectorV{lev,1} = [0.*ones(1,pV), 0:0.5^(lev-1)*size(F,1)/nelemy:size(F,1), size(F,1).*ones(1,pV)];
    uknotvectorV{lev,1} = unique(knotvectorV{lev,1});
    nelemV(lev,1) = size(uknotvectorV{lev,1},2) -1;
    kV(lev,1) = length(knotvectorV{lev,1});
    nobV(lev,1) = kV(lev,1) - pV -1;
    
    px = zeros(nobU(lev,1)*nobV(lev,1),2);
    BBv = zeros(nobU(lev,1)*nobV(lev,1),2);
    CP{lev,1} = px;
    Bvect{lev,1} = BBv;
end

%% Element
Elem = cell(maxlev,1);
AA = cell(maxlev,1);
PC = cell(maxlev,1);
Basis = cell(maxlev,1);
ch = zeros(4,1);
connectU = zeros(2,1);
connectV = zeros(2,1);

for lev = 1:maxlev,
    knotU = knotvectorU{lev,1};
    knotV = knotvectorV{lev,1};
    supp_ct = zeros(nobU(lev,1)*nobV(lev,1),1);
    supp_i = zeros(nobU(lev,1)*nobV(lev,1),1);
    
    if(lev<=maxlev-1)
        Parcell = zeros(nelemU(lev+1,1)*nelemV(lev+1,1),1);
    end
    %elem = cell(nelemU(level,1),1);
    glnumnodes = 0;
    E = 0;
    BBvector = zeros(1,(1+pU)*(1+pV));
    EEN = cell(nelemU(lev,1)*nelemV(lev,1),9);
    for intU = 1:nobU(lev,1),
        for intV = 1:nobV(lev,1),
            glnumnodes = glnumnodes + 1;
            if((intU>=(pU+1))&&(intV>=(pV+1))),
                E = E+1;
                if(lev>1)
                    pcell = PC{lev,1};
                    EEN{E,10} = pcell(E,1);
                end
                yind = ceil(E/nelemU(lev,1));
                xind = mod(E,nelemU(lev,1));
                if(xind==0),
                    xind = nelemU(lev,1);
                end
                i = 1;
                EEN{E,1} = E;
                
                connectU(1,1) = intU;
                connectU(2,1) = connectU(1,1)+1;
                
                
                connectV(1,1) = intV;
                connectV(2,1) = connectV(1,1)+1;
                
                EEN{E,2} = connectV;
                EEN{E,3} = connectU;
                
                if(lev==1),
                    EEN{E,4} = 1;
                else
                    EEN{E,4} = 0;
                end
                if(lev<=maxlev-1)
                    
                    ch(1,1) = nelemU(lev+1,1)*((2*yind-1)-1)+(2*xind-1);
                    ch(2,1) = nelemU(lev+1,1)*((2*yind-1)-1)+(2*xind);
                    ch(3,1) = nelemU(lev+1,1)*((2*yind)-1)+(2*xind);
                    ch(4,1) = nelemU(lev+1,1)*((2*yind)-1)+(2*xind-1);
                    for pct = 1:4,
                        Parcell(ch(pct,1),1) = E;
                    end
                else
                    ch = zeros(4,1);
                end
                EEN{E,5} = ch;
                
                
                for loci = 1:pU+1,
                    for locj = 1:pV+1,
                        B  = glnumnodes - (loci-1)*(kV(lev,1)-pV-1) -locj + 1;
                        BBvector(1,i) = B;
                        supp_ct(B,1) = supp_ct(B,1)+1;
                        supp_i(B,supp_ct(B,1)) = E;
                        i=i+1;
                    end
                end
                EEN{E,6} = BBvector;
                
                EEN{E,7} = lev;
                
                
                cellcx = knotU(intV)+(knotU(intV+1)-knotU(intV))/2;
                cellcy = knotV(intU)+(knotV(intU+1)-knotV(intU))/2;
                
                EEN{E,8} = cellcx;
                EEN{E,9} = cellcy;
                EEN{E,11} = 0;
            end
        end
    end
    Elem{lev,1} = EEN;
    EEM = Elem{lev,1};
    if(lev<=maxlev-1)
        PC{lev+1,1} = Parcell;
    end
    CCX = zeros(sqrt(size(EEM,1)));
    CCY = zeros(sqrt(size(EEM,1)));
    for i = 1:size(EEM,1),
        yind = ceil(i/nelemU(lev,1));
        xind = mod(i,nelemU(lev,1));
        if(xind==0),
            xind = nelemU(lev,1);
        end
        CCX(xind,yind) = EEM{i,8};
        CCY(xind,yind) = EEM{i,9};
    end
    
    Basis_m = cell(nobU(lev,1)*nobV(lev,1),7);
    if(lev<=maxlev-1)
        Anc = zeros(nobU(lev+1)*nobV(lev+1),1);
        counta = zeros(nobU(lev+1)*nobV(lev+1),1);
    end
    knot_cu = knotvectorU{lev,1};
    knot_cv = knotvectorV{lev,1};
    if(lev<=maxlev-1)
        knot_fu = knotvectorU{lev+1,1};
        knot_fv = knotvectorV{lev+1,1};
    end
    
    
    bc = 0;
    for basis_i =1:nobU(lev,1),
        for basis_j = 1:nobV(lev,1),
            bc=bc+1;
            if(lev>1),
                AC = AA{lev,1};
                Basis_m{bc,5} = AC(bc,:);
            end
            Basis_m{bc,6} = supp_i(bc,:);
            BB(1,1) = basis_j;
            BB(1,2) = basis_i;
            Basis_m{bc,1} = BB;
            
            Basis_m{bc,2} = lev;
            Basis_m{bc,8} = 0;
            Basis_m{bc,9} = 0;
            if(lev==1),
                Basis_m{bc,3} = 1;
            else
                Basis_m{bc,3} = 0;
            end
            
            Basis_m{bc,7} = 0;
            for iU = 1:nobU(lev,1),
                dof1 = iU;
                Basis_m{dof1,9} = 1;
                dof2 = iU + nobU(lev,1)*(nobV(lev,1)-1);
                Basis_m{dof2,9} = 1;
            end
            
            for iV = 1:nobV(lev,1),
                dof3 = 1+nobU(lev,1)*(iV-1);
                Basis_m{dof3,9}  =1;
                dof4 = nobU(lev,1)+nobU(lev,1)*(iV-1);
                Basis_m{dof4,9}  =1;
            end
            if(lev<=maxlev-1)
                intc_u1 = knot_cu(1,basis_j);
                intc_u2 = knot_cu(1,basis_j+pU+1);
                
                intc_v1 = knot_cv(1,basis_i);
                intc_v2 = knot_cv(1,basis_i+pV+1);
                
                if(basis_j < pU+1),
                    uindex_start = basis_j;
                else
                    uindex_start = FindSpan(nobU(lev+1,1),pU,intc_u1,knot_fu)+1;
                end
                
                if((nobU(lev,1)-basis_j)<pU)
                    uindex_end = size(knotvectorU{lev+1,1},2)-(nobU(lev,1)-basis_j);
                else
                    uindex_end = FindSpan(nobU(lev+1,1),pU,intc_u2,knot_fu)+1;
                end
                
                if(basis_i < pV+1),
                    vindex_start = basis_i;
                else
                    vindex_start = FindSpan(nobV(lev+1,1),pV,intc_v1,knot_fv)+1;
                end
                if((nobV(lev,1)-basis_i)<pV)
                    vindex_end = size(knotvectorV{lev+1,1},2)-(nobV(lev,1)-basis_i);
                else
                    vindex_end = FindSpan(nobV(lev+1,1),pV,intc_v2,knot_fv)+1;
                end
                
                Knotu = knot_cu(basis_j:(basis_j+pU+1));
                Knotv = knot_cv(basis_i:(basis_i+pV+1));
                
                newKnotu = knot_fu(uindex_start:uindex_end);
                newKnotv = knot_fv(vindex_start:vindex_end);
                
                nob_childu = size(newKnotu,2)-pU-1;
                nob_childv = size(newKnotv,2)-pV-1;
                nu = size(Knotu,2)-pU-1;
                nv = size(Knotv,2)-pV-1;
                TmatU =  Tmatrix(Knotu,newKnotu,pU);
                TmatV =  Tmatrix(Knotv,newKnotv,pV);
                
                TmatU = TmatU(1:nob_childu,1:nu);
                TmatV = TmatV(1:nob_childv,1:nv);
                Tmat = TmatU*TmatV';
                Coeff = zeros(nob_childu*nob_childv,2);
                C = zeros(nob_childu*nob_childv,2);
                cc=0;
                for c_childu = 1:nob_childv,
                    for c_childv = 1:nob_childu,
                        cc= cc+1;
                        C(cc,1) = uindex_start+(c_childv-1);
                        C(cc,2) = vindex_start+(c_childu-1);
                        
                        cbb = nobU(lev+1,1)*(C(cc,2)-1) + (C(cc,1));
                        Coeff(cc,1) = Tmat(c_childv,c_childu);
                        counta(cbb,1) = counta(cbb,1)+1;
                        Coeff(cc,2) = cbb;
                        Anc(cbb,counta(cbb,1)) = bc;
                        
                    end
                end
                
                Basis_m{bc,4} = C;
                Basis_m{bc,12} = Coeff;
                
            else
                Basis_m{bc,4} = 0;
            end
        end
        
    end
    Basis{lev,1} = Basis_m;
    if(lev<=maxlev-1)
        
        AA{lev+1,1} = Anc;
        
    end
end

Dm = Basis;
Em = Elem;

knotu = knotvectorU{1,1};
knotv = knotvectorV{1,1};

pp = zeros(nobU(1,1)*nobV(1,1),2);
for i = 1:nobU(1,1),
    for j =1:nobV(1,1),
        index27=(i-1)*nobV(1,1)+j;
        coordx = sum(knotu(i+1:i+pU))./pU;
        coordy = sum(knotv(j+1:j+pV))./pV;
        pp(index27,:) = [coordx,coordy];
    end
end

CP{1,1} = pp;
Pm = CP;

end