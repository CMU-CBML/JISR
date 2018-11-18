function [RHSf, RHSb] = computeIntegrationFidelity(Btermf, Btermb, Bseg, Jm, PHI, Dm, H,Bvectf,Bvectb,par1,par2)

[~,w1] = ggquad(6);
xlen =6;
w2 = w1;
ac_ct = size(Jm,1);

for i = 1:ac_ct,
    
    term1f = Btermf(1+(i-1)*xlen:i*xlen,1:6,1);
    term2f = Btermf(1+(i-1)*xlen:i*xlen,1:6,2);
    
    term1b = Btermb(1+(i-1)*xlen:i*xlen,1:6,1);
    term2b = Btermb(1+(i-1)*xlen:i*xlen,1:6,2);
    
    term7 = Bseg(1+(i-1)*xlen:i*xlen,1:6,1);
    term8 = Bseg(1+(i-1)*xlen:i*xlen,1:6,2);
    term9 = Bseg(1+(i-1)*xlen:i*xlen,1:6,3);
    
    SB = Jm{i,1};
    supp_phi = PHI{i,1};
    supp_size = size(SB,1);
    
    
    for bg = 1:supp_size,
        %BBMvector = Bvect{SB(bg,2),1};
        valm1f = zeros(xlen,xlen);
        valm2f = zeros(xlen,xlen);
        
        valm1b = zeros(xlen,xlen);
        valm2b = zeros(xlen,xlen);
        
        for gg1 = 1:xlen,
            for gg2 = 1:xlen,
                
                phi_i = supp_phi(bg,gg1,gg2);
                
                valm1f(gg1,gg2) = par1*phi_i*term7(gg1,gg2)*term8(gg1,gg2) + par2*(phi_i)*(term1f(gg1,gg2));
                valm2f(gg1,gg2) = par1*phi_i*term7(gg1,gg2)*term9(gg1,gg2) + par2*(phi_i)*(term2f(gg1,gg2));
                valm1b(gg1,gg2) = par2*(phi_i)*(term1b(gg1,gg2));
                valm2b(gg1,gg2) = par2*(phi_i)*(term2b(gg1,gg2));
            end
        end
        
        h1 = H(i,1);
        h2 = H(i,2);
        
        val1f = w1' * valm1f * w2 * h1 * h2;
        val2f = w1' * valm2f * w2 * h1 * h2;
        
        val1b = w1' * valm1b * w2 * h1 * h2;
        val2b = w1' * valm2b * w2 * h1 * h2;
        
        btempp = Dm{SB(bg,2),1};
        bact_i = btempp{SB(bg,1),10};
        Bvectf(bact_i,1) = Bvectf(bact_i,1) + val1f;
        Bvectf(bact_i,2) = Bvectf(bact_i,2) + val2f;
        
        Bvectb(bact_i,1) = Bvectb(bact_i,1) + val1b;
        Bvectb(bact_i,2) = Bvectb(bact_i,2) + val2b;
    end
end

RHSf = Bvectf;
RHSb = Bvectb;
end