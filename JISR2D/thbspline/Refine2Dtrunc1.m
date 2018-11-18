function [Ifinal,Efinal, Pfinal] = Refine2Dtrunc1(indexi,indexj,level,Dem,Elem,P,knotvectorU,knotvectorV,pU,pV)

Em = Elem; %transfer the element hierarchy
Bm = Dem; % transfer the bspline hierarchy
Pm = P; %transfer the control point hierarchy
P1 = Pm{level,1};  %control points for current and next level refinement
P2 = Pm{level+1,1};
El1 = Em{level,1}; %element indices for current and next level refinement
El2 = Em{level+1,1};
B1 = Bm{level,1};  %bspline indices for current and next level refinement
B2 = Bm{level+1,1};

knotcu = knotvectorU{level,1}; %knot vectors for current and next level refinement
knotfu = knotvectorU{level+1,1};
knotcv = knotvectorV{level,1};
knotfv = knotvectorV{level+1,1};
Bindex = sqrt(size(B1,1))*(indexj-1)+indexi; %calculate the global index of a bspline
% global_index = (nobx)*(indexj-1)+ (indexi)
supp_b = B1{Bindex,6};
b_ct = 0;
Btemp = zeros(1,1);
%supp_b is the support cell indices of bspline with global index Bindex
%============================================================================
% Cell refinement
for i=1:size(supp_b,2),%Loop over th support cells of the bsr.
    if(supp_b(1,i)~=0 && El1{supp_b(1,i),4}==1), %I have zeros values appending, that's why only consider non-zero
        El1{supp_b(1,i),4}=0; % Make them inactive
        CIcell = El1{supp_b(1,i),5};
        noci = size(CIcell,1);
        for k = 1:noci,  %Loop over their children cells
            ind = CIcell(k,1);
            El2{ind,4} = 1;  %Make them active
        end
    end
end
%============================================================================
%Bspline refinement
for i =1:size(B1,1),
    supp_b = B1{i,6};
    supp_ct = 0;
    inact = 0;
    for j = 1:size(supp_b,2),
        if(supp_b(1,j)~=0),
            supp_ct = supp_ct+1;
            if(El1{supp_b(1,j),4}==0),
                inact = inact+1;
            end
        end
    end
    if(supp_ct==inact),
        if(B1{i,3}==1),
            B1{i,3} = 0;
            B1{i,7} = 1;
            
            nob_child = size(B1{i,4},1);
            Ci = B1{i,4};
            for j = 1:nob_child,
                Bind1 = Ci(j,1);
                Bind2 = Ci(j,2);
                Bcind = sqrt(size(B2,1))*(Bind2-1)+Bind1;
                B2{Bcind,3} = 1;
            end
            
            if(B1{i,8}==0),
                b_ct = b_ct + 1;
                Btemp(b_ct,1) = i;
            end
            
            idu = B1{i,1}(1,1);
            idv = B1{i,1}(1,2);
            nobc= sqrt(size(B1,1));
            nobf= sqrt(size(B2,1));
            
            if (idu==1),
                ofx1 = 0;
            elseif (idu==2),
                ofx1 = 1;
            else
                ofx1 = 2;
            end
            
            if (idu==nobc),
                ofx2 = 0;
            elseif (idu==(nobc-1)),
                ofx2 = 1;
            else
                ofx2 = 2;
            end
            
            if (idv==1),
                ofy1 = 0;
            elseif (idv==2),
                ofy1 = 1;
            else
                ofy1 = 2;
            end
            
            if (idv==nobc),
                ofy2 = 0;
            elseif (idv==(nobc-1)),
                ofy2 = 1;
            else
                ofy2 = 2;
            end
            for t1 = (idu-ofx1):1:(idu+ofx2),
                for t2 = (idv-ofy1):1:(idv+ofy2),
                    gg = nobc*(t2-1)+t1;
                    if(B1{gg,3}==1 && B1{gg,8}==0),
                        B1{gg,8} = 1;
                        b_ct = b_ct+1;
                        Btemp(b_ct,1) = gg;
                    end
                end
            end
        end
    end
end
%Set truncation flags
for i = 1:size(Btemp,1),
    if(Btemp(1,1)~=0)
        
        itdu = B1{Btemp(i,1),1}(1,1);
        itdv = B1{Btemp(i,1),1}(1,2);
        newkcu = knotcu(itdu:itdu+pU+1);
        newkcv = knotcv(itdv:itdv+pV+1);
        
        if(itdu < pU+1),
            usu = itdu;
        else
            usu = FindSpan(nobf,pU,knotcu(1,itdu),knotfu)+1;
        end
        
        if((nobc-itdu)<pU)
            ueu= size(knotvectorU{level+1,1},2)-(nobc-itdu);
        else
            ueu = FindSpan(nobf,pU,knotcu(1,itdu+pU+1),knotfu)+1;
        end
        
        if(itdv < pV+1),
            usv = itdv;
        else
            usv = FindSpan(nobf,pV,knotcv(1,itdv),knotfv)+1;
        end
        
        if((nobc-itdv)<pV)
            uev = size(knotvectorV{level+1,1},2)-(nobc-itdv);
        else
            uev = FindSpan(nobf,pV,knotcv(1,itdv+pV+1),knotfv)+1;
        end
        
        newkfu = knotfu(usu:ueu);
        newkfv = knotfv(usv:uev);
        
        unb1 = size(newkcu,2)-pU-1;
        unb2 = size(newkfu,2)-pU-1;
        
        vnb1 = size(newkcv,2)-pV-1;
        vnb2 = size(newkfv,2)-pV-1;
        
        TmatU =  Tmatrix(newkcu,newkfu,pU);
        TmatV =  Tmatrix(newkcv,newkfv,pV);
        
        TmatU = TmatU(1:unb2,1:unb1);
        TmatV = TmatV(1:vnb2,1:vnb1);
        
        
        Pnew1 = ((TmatU*P1(Btemp(i,1),1))*TmatV');
        Pnew2 = ((TmatU*P1(Btemp(i,1),2))*TmatV');
        
        for kk1=1:unb2,
            for kk2 =1:vnb2,
                ing = nobf*(usv+(kk2-1)-1)+(usu+(kk1-1));
                
                P2(ing,1) = P2(ing,1)+Pnew1(kk1,kk2);
                P2(ing,2) = P2(ing,2)+Pnew2(kk1,kk2);
                
            end
        end
    end
end


Em{level,1} = El1;
Em{level+1,1} = El2;
Bm{level,1} = B1;
Bm{level+1,1} = B2;
Pm{level+1,1} = P2;
Pfinal = Pm;
Ifinal = Bm;
Efinal = Em;
end





