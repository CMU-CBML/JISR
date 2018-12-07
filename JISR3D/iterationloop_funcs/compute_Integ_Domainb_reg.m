function RHS = compute_Integ_Domainb_reg(Jm,BIGX, BIGY, BIGZ,FBX, FBY,FBZ, FBXDX,FBXDY,FBXDZ, FBYDX, FBYDY, FBYDZ, FBZDX, FBZDY,FBZDZ, BFX,BFY,BFZ,BIGMUXf,BIGMUYf,BIGMUZf,BIGMVXf,BIGMVYf,BIGMVZf,BIGMWXf,BIGMWYf,BIGMWZf,RHSf,PHI1, PHIU1, PHIV1, PHIW1, PHIF1,lambda_2,lambda_3,mu,w1,w2,w3,H)

ac_ct = size(Jm,1);
bf_ct = size(RHSf,1);
xlen = size(BIGX,2);

parfor i = 1:ac_ct,
    
    RHSf1 = zeros(bf_ct,4);
    
    SB = Jm(i).nzsplines;
    supp_phi = PHI1(i).mat;
    supp_phiu = PHIU1(i).mat;
    supp_phiv = PHIV1(i).mat;
    supp_phiw = PHIW1(i).mat;
    supp_size = size(SB,1);
    
    supp_phif = PHIF1(i).mat;

    hu = H(i,1);
    hv = H(i,2);
    hw = H(i,3);

    term4f = BIGMUXf(1+(i-1)*xlen:i*xlen,1:xlen,1:xlen);
    term5f = BIGMUYf(1+(i-1)*xlen:i*xlen,1:xlen,1:xlen);
    term6f = BIGMUZf(1+(i-1)*xlen:i*xlen,1:xlen,1:xlen);
    
    term7f = BIGMVXf(1+(i-1)*xlen:i*xlen,1:xlen,1:xlen);
    term8f = BIGMVYf(1+(i-1)*xlen:i*xlen,1:xlen,1:xlen);
    term9f = BIGMVZf(1+(i-1)*xlen:i*xlen,1:xlen,1:xlen);
    
    term10f = BIGMWXf(1+(i-1)*xlen:i*xlen,1:xlen,1:xlen);
    term11f = BIGMWYf(1+(i-1)*xlen:i*xlen,1:xlen,1:xlen);
    term12f = BIGMWZf(1+(i-1)*xlen:i*xlen,1:xlen,1:xlen);

    gaussptx = BIGX(1+(i-1)*xlen:i*xlen,1:xlen,1:xlen);
    gausspty = BIGY(1+(i-1)*xlen:i*xlen,1:xlen,1:xlen);
    gaussptz = BIGZ(1+(i-1)*xlen:i*xlen,1:xlen,1:xlen);
    
    fbxx =  FBX(1+(i-1)*xlen:i*xlen,1:xlen,1:xlen);
    fbyy =  FBY(1+(i-1)*xlen:i*xlen,1:xlen,1:xlen);
    fbzz =  FBZ(1+(i-1)*xlen:i*xlen,1:xlen,1:xlen);

    bfxx =  BFX(1+(i-1)*xlen:i*xlen,1:xlen,1:xlen);
    bfyy =  BFY(1+(i-1)*xlen:i*xlen,1:xlen,1:xlen);
    bfzz =  BFZ(1+(i-1)*xlen:i*xlen,1:xlen,1:xlen);
    
    fbxdx = FBXDX(1+(i-1)*xlen:i*xlen,1:xlen,1:xlen);
    fbxdy = FBXDY(1+(i-1)*xlen:i*xlen,1:xlen,1:xlen);
    fbxdz = FBXDZ(1+(i-1)*xlen:i*xlen,1:xlen,1:xlen);
    
    fbydx = FBYDX(1+(i-1)*xlen:i*xlen,1:xlen,1:xlen);
    fbydy = FBYDY(1+(i-1)*xlen:i*xlen,1:xlen,1:xlen);
    fbydz = FBYDZ(1+(i-1)*xlen:i*xlen,1:xlen,1:xlen);
    
    fbzdx = FBZDX(1+(i-1)*xlen:i*xlen,1:xlen,1:xlen);
    fbzdy = FBZDY(1+(i-1)*xlen:i*xlen,1:xlen,1:xlen);
    fbzdz = FBZDZ(1+(i-1)*xlen:i*xlen,1:xlen,1:xlen);
        
    val1f = zeros(supp_size,1);
    val2f = zeros(supp_size,1);
    val3f = zeros(supp_size,1);
    
    valm1f = zeros(supp_size,xlen,xlen,xlen);
    valm2f = zeros(supp_size,xlen,xlen,xlen);
    valm3f = zeros(supp_size,xlen,xlen,xlen);
    
    for gg1 = 1:xlen,
        for gg2 = 1:xlen,
            for gg3 = 1:xlen
                phi_i  = supp_phi(:,gg1,gg2,gg3);
                phi_ui = supp_phiu(:,gg1,gg2,gg3);
                phi_vi = supp_phiv(:,gg1,gg2,gg3);
                phi_wi = supp_phiw(:,gg1,gg2,gg3);
                
                phi_if = supp_phif(:,gg1,gg2,gg3);

                xfb = [gaussptx(gg1,gg2,gg3)-fbxx(gg1,gg2,gg3), gausspty(gg1,gg2,gg3)-fbyy(gg1,gg2,gg3),gaussptz(gg1,gg2,gg3)-fbzz(gg1,gg2,gg3)];
                xbf = [-gaussptx(gg1,gg2,gg3)+bfxx(gg1,gg2,gg3), -gausspty(gg1,gg2,gg3)+bfyy(gg1,gg2,gg3),-gaussptz(gg1,gg2,gg3)+bfzz(gg1,gg2,gg3)];
                
                Jfb = [fbxdx(gg1,gg2,gg3), fbxdy(gg1,gg2,gg3), fbxdz(gg1,gg2,gg3); fbydx(gg1,gg2,gg3), fbydy(gg1,gg2,gg3), fbydz(gg1,gg2,gg3);fbzdx(gg1,gg2,gg3), fbzdy(gg1,gg2,gg3), fbzdz(gg1,gg2,gg3)];
                
                regfb = 2*mu*xfb*Jfb;
                
                reg1fb = regfb(1,1).*phi_i + 2*mu*xbf(1,1).*phi_if;
                reg2fb = regfb(1,2).*phi_i + 2*mu*xbf(1,2).*phi_if;
                reg3fb = regfb(1,3).*phi_i + 2*mu*xbf(1,3).*phi_if;
                
                valm1f(:,gg1,gg2,gg3) = 2*lambda_2*(term4f(gg1,gg2,gg3).*phi_ui+term5f(gg1,gg2,gg3).*phi_vi+term6f(gg1,gg2,gg3).*phi_wi) +  lambda_3.*reg1fb;
                valm2f(:,gg1,gg2,gg3) = 2*lambda_2*(term7f(gg1,gg2,gg3).*phi_ui+term8f(gg1,gg2,gg3).*phi_vi+term9f(gg1,gg2,gg3).*phi_wi) +  lambda_3.*reg2fb;
                valm3f(:,gg1,gg2,gg3) = 2*lambda_2*(term10f(gg1,gg2,gg3).*phi_ui+term11f(gg1,gg2,gg3).*phi_vi+term12f(gg1,gg2,gg3).*phi_wi) + lambda_3.*reg3fb;
                
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

