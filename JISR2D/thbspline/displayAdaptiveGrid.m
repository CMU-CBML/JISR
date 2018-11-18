figure;
CF = cell(ac_ct,1);


for i = 1:ac_ct,
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
    cfx = zeros(4,2);
    
    for j =1:4,
        
        uu = fx(j,1);
        vv = fx(j,2);
        
        sumbx = 0;
        sumby = 0;
        RR1 = BasisFun(u(1,1)-1,uu,pU,knotU1);
        RR2 = BasisFun(v(1,1)-1,vv,pV,knotV1);
        RR = RR1'*RR2;
        inc=0;
        phii=zeros(16,1);
        for m1=1:size(RR,1),
            for m2=1:size(RR,2),
                
                phii(16-inc,1) = RR(m2,m1);
                inc = inc+1;
            end
        end
        PHID = c1*phii;
        for k = 1:supp_size1,
            CEb = Pmold{SB1(k,2),1};
            pi = CEb(SB1(k,1),1);
            pj = CEb(SB1(k,1),2);
            sumbx = sumbx + pj*PHID(k,1);
            sumby = sumby + pi*PHID(k,1);

            
        end
        cfx(j,1) = sumbx;
        cfx(j,2) = sumby;
    end
    x1=cfx(1,1);
    x2=cfx(2,1);
    y1=cfx(1,2);
    y2=cfx(3,2);
    x = [x1, x2, x2, x1, x1];
    y = [y1, y1, y2, y2, y1];
    
    plot(x,y,'-','color','black','MarkerEdgeColor','k')
    axis([0 nx 0 ny])
    set(gca,'position',[0 0 1 1],'units','normalized')
    hold on;
    CF{i,1} = cfx;
    % end
end

hold off