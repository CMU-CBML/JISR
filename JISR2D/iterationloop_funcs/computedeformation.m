function [BIGXf,BIGYf,BIGXb,BIGYb] = computedeformation(Jm,PHI,PHIU,PHIV,Pmf,Pmb)

ac_ct = size(Jm,1);
xlen = 6;

BIGGXf = zeros(ac_ct*xlen,xlen,3);
BIGGYf = zeros(ac_ct*xlen,xlen,3);
BIGGXb = zeros(ac_ct*xlen,xlen,3);
BIGGYb = zeros(ac_ct*xlen,xlen,3);

BIGXXf = zeros(ac_ct*xlen,xlen);
BIGYYf = zeros(ac_ct*xlen,xlen);
BIGMUXf = zeros(ac_ct*xlen,xlen);
BIGMUYf = zeros(ac_ct*xlen,xlen);
BIGMVXf = zeros(ac_ct*xlen,xlen);
BIGMVYf = zeros(ac_ct*xlen,xlen);

BIGXXb = zeros(ac_ct*xlen,xlen);
BIGYYb = zeros(ac_ct*xlen,xlen);
BIGMUXb = zeros(ac_ct*xlen,xlen);
BIGMUYb = zeros(ac_ct*xlen,xlen);
BIGMVXb = zeros(ac_ct*xlen,xlen);
BIGMVYb = zeros(ac_ct*xlen,xlen);

for i = 1:ac_ct
    
    SB = Jm{i,1};
    supp_phi = PHI{i,1};
    supp_phiu = PHIU{i,1};
    supp_phiv = PHIV{i,1};
    supp_size = size(SB,1);
    
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
                CEbf = Pmf{SB(kg,2),1};
                CEbb = Pmb{SB(kg,2),1};
                
                pif = CEbf(SB(kg,1),1);
                pjf = CEbf(SB(kg,1),2);
                
                pib = CEbb(SB(kg,1),1);
                pjb = CEbb(SB(kg,1),2);
                
                sumbxf = sumbxf + pif*supp_phi(kg,gg1,gg2);
                sumbyf = sumbyf + pjf*supp_phi(kg,gg1,gg2);
                sumbuxf = sumbuxf + pif*supp_phiu(kg,gg1,gg2);
                sumbuyf = sumbuyf + pjf*supp_phiu(kg,gg1,gg2);
                sumbvxf = sumbvxf + pif*supp_phiv(kg,gg1,gg2);
                sumbvyf = sumbvyf + pjf*supp_phiv(kg,gg1,gg2);
                
                sumbxb = sumbxb + pib*supp_phi(kg,gg1,gg2);
                sumbyb = sumbyb + pjb*supp_phi(kg,gg1,gg2);
                sumbuxb = sumbuxb + pib*supp_phiu(kg,gg1,gg2);
                sumbuyb = sumbuyb + pjb*supp_phiu(kg,gg1,gg2);
                sumbvxb = sumbvxb + pib*supp_phiv(kg,gg1,gg2);
                sumbvyb = sumbvyb + pjb*supp_phiv(kg,gg1,gg2);
                
            end
            
            SBXf(gg1,gg2) = sumbxf;
            SBYf(gg1,gg2) = sumbyf;
            SBUXf(gg1,gg2) = sumbuxf;
            SBUYf(gg1,gg2) = sumbuyf;
            SBVXf(gg1,gg2) = sumbvxf;
            SBVYf(gg1,gg2) = sumbvyf;
            
            SBXb(gg1,gg2) = sumbxb;
            SBYb(gg1,gg2) = sumbyb;
            SBUXb(gg1,gg2) = sumbuxb;
            SBUYb(gg1,gg2) = sumbuyb;
            SBVXb(gg1,gg2) = sumbvxb;
            SBVYb(gg1,gg2) = sumbvyb;
        end
    end
    
    BIGXXf(1+(i-1)*xlen:i*xlen,1:6) = SBXf;
    BIGYYf(1+(i-1)*xlen:i*xlen,1:6) = SBYf;
    BIGMUXf(1+(i-1)*xlen:i*xlen,1:6) = SBUXf;
    BIGMUYf(1+(i-1)*xlen:i*xlen,1:6) = SBUYf;
    BIGMVXf(1+(i-1)*xlen:i*xlen,1:6) = SBVXf;
    BIGMVYf(1+(i-1)*xlen:i*xlen,1:6) = SBVYf;
    
    BIGXXb(1+(i-1)*xlen:i*xlen,1:6) = SBXb;
    BIGYYb(1+(i-1)*xlen:i*xlen,1:6) = SBYb;
    BIGMUXb(1+(i-1)*xlen:i*xlen,1:6) = SBUXb;
    BIGMUYb(1+(i-1)*xlen:i*xlen,1:6) = SBUYb;
    BIGMVXb(1+(i-1)*xlen:i*xlen,1:6) = SBVXb;
    BIGMVYb(1+(i-1)*xlen:i*xlen,1:6) = SBVYb;
end

BIGGXf(:,:,1) = BIGXXf;
BIGGXf(:,:,2) = BIGMUXf;
BIGGXf(:,:,3) = BIGMUYf;

BIGGYf(:,:,1) = BIGYYf;
BIGGYf(:,:,2) = BIGMVXf;
BIGGYf(:,:,3) = BIGMVYf;

BIGGXb(:,:,1) = BIGXXb;
BIGGXb(:,:,2) = BIGMUXb;
BIGGXb(:,:,3) = BIGMUYb;

BIGGYb(:,:,1) = BIGYYb;
BIGGYb(:,:,2) = BIGMVXb;
BIGGYb(:,:,3) = BIGMVYb;

BIGXf = BIGGXf;
BIGYf = BIGGYf;
BIGXb = BIGGXb;
BIGYb = BIGGYb;

end