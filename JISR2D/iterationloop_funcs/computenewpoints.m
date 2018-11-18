function [PXF,PXB] = computenewpoints(M,Pixel,Jm,Pmf,Pmb)

px = 0;
pxxf= zeros(size(M,1),size(M,2));
pyyf= zeros(size(M,1),size(M,2));
pxxb= zeros(size(M,1),size(M,2));
pyyb= zeros(size(M,1),size(M,2));
for i = 1:size(M,1),
    for j = 1:size(M,2),
        px = px+1;
        ac_ind = Pixel{px,1};
        supp = Pixel{px,2};
        SB = Jm{ac_ind,1};
        ss = size(SB,1);
        fxxf = 0;
        fyyf = 0;
        
        fxxb = 0;
        fyyb = 0;
        
        for k = 1:ss,
            
            CEbf = Pmf{SB(k,2),1};
            pif = CEbf(SB(k,1),1);
            pjf = CEbf(SB(k,1),2);
            
            fxxf = fxxf + pif*supp(k,1);
            fyyf = fyyf + pjf*supp(k,1);
            
            CEbb = Pmb{SB(k,2),1};
            pib = CEbb(SB(k,1),1);
            pjb = CEbb(SB(k,1),2);
            
            fxxb = fxxb + pib*supp(k,1);
            fyyb = fyyb + pjb*supp(k,1);
        end
        
        pxxf(i,j) = fxxf;
        pyyf(i,j) = fyyf;
        
        pxxb(i,j) = fxxb;
        pyyb(i,j) = fyyb;
    end
end
PXF(:,:,1) = pxxf;
PXF(:,:,2) = pyyf;
PXB(:,:,1) = pxxb;
PXB(:,:,2) = pyyb;
end