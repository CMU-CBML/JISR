function RHS = compute_Integ_Domainf_fid(Jm,Bseg,Bsegx,Bsegy,Bsegz,Bterm1f,Bterm2f,Bterm3f,RHSf,PHI1,par1,par2,w1,w2,w3,H)

ac_ct = size(Jm,1);
bf_ct = size(RHSf,1);
xlen = 4;

parfor i = 1:ac_ct
    
    RHSf1 = zeros(bf_ct,4);
    
    SB = Jm(i).nzsplines;
    supp_phi = PHI1(i).mat;
    supp_size = size(SB,1);

    hu = H(i,1);
    hv = H(i,2);
    hw = H(i,3);
    
    term1f = Bterm1f(1+(i-1)*xlen:i*xlen,:,:);
    term2f = Bterm2f(1+(i-1)*xlen:i*xlen,1:xlen,1:xlen);
    term3f = Bterm3f(1+(i-1)*xlen:i*xlen,1:xlen,1:xlen);

    term7 = Bseg(1+(i-1)*xlen:i*xlen,1:xlen,1:xlen);
    term8 = Bsegx(1+(i-1)*xlen:i*xlen,1:xlen,1:xlen);
    term9 = Bsegy(1+(i-1)*xlen:i*xlen,1:xlen,1:xlen);
    term10 = Bsegz(1+(i-1)*xlen:i*xlen,1:xlen,1:xlen);

    val1f = zeros(supp_size,1);
    val2f = zeros(supp_size,1);
    val3f = zeros(supp_size,1);
    
    valm1f = zeros(supp_size,xlen,xlen,xlen);
    valm2f = zeros(supp_size,xlen,xlen,xlen);
    valm3f = zeros(supp_size,xlen,xlen,xlen);
    
    for gg1 = 1:xlen
        for gg2 = 1:xlen
            for gg3 = 1:xlen
                phi_i  = supp_phi(:,gg1,gg2,gg3);

                valm1f(:,gg1,gg2,gg3) = par1*term7(gg1,gg2,gg3)*term8(gg1,gg2,gg3)*phi_i + par2*(phi_i)*(term1f(gg1,gg2,gg3));
                valm2f(:,gg1,gg2,gg3) = par1*term7(gg1,gg2,gg3)*term9(gg1,gg2,gg3)*phi_i + par2*(phi_i)*(term2f(gg1,gg2,gg3));
                valm3f(:,gg1,gg2,gg3) = par1*term7(gg1,gg2,gg3)*term10(gg1,gg2,gg3)*phi_i + par2*(phi_i)*(term3f(gg1,gg2,gg3));
                
                val1f(:,1) = val1f(:,1) + w1(gg1,1).*w2(gg2,1).*w3(gg3,1).*valm1f(:,gg1,gg2,gg3).*hu.*hv.*hw;
                val2f(:,1) = val2f(:,1) + w1(gg1,1).*w2(gg2,1).*w3(gg3,1).*valm2f(:,gg1,gg2,gg3).*hu.*hv.*hw;
                val3f(:,1) = val3f(:,1) + w1(gg1,1).*w2(gg2,1).*w3(gg3,1).*valm3f(:,gg1,gg2,gg3).*hu.*hv.*hw;
                
            end
        end
    end
    
    RHSf1(SB,1) = val1f;
    RHSf1(SB,2) = val2f;
    RHSf1(SB,3) = val3f;
    
    RHSf = RHSf  + RHSf1;
    
end

RHS = RHSf;
end

