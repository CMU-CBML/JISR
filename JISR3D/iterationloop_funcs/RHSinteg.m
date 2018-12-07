VV1 = cell(ac_ct,1);
VV2 = cell(ac_ct,1);
VV3 = cell(ac_ct,3);

for i = 1:ac_ct,
    
    RHSf1 = zeros(bf_ct,4);
    RHSb1 = zeros(bf_ct,4);
    
    term1f = Bterm1f(1+(i-1)*xlen:i*xlen,1:xlen,1:xlen);
    term2f = Bterm2f(1+(i-1)*xlen:i*xlen,1:xlen,1:xlen);
    term3f = Bterm3f(1+(i-1)*xlen:i*xlen,1:xlen,1:xlen);
    
    term4f = BIGMUXf(1+(i-1)*xlen:i*xlen,1:xlen,1:xlen);
    term5f = BIGMUYf(1+(i-1)*xlen:i*xlen,1:xlen,1:xlen);
    term6f = BIGMUZf(1+(i-1)*xlen:i*xlen,1:xlen,1:xlen);
    
    term7f = BIGMVXf(1+(i-1)*xlen:i*xlen,1:xlen,1:xlen);
    term8f = BIGMVYf(1+(i-1)*xlen:i*xlen,1:xlen,1:xlen);
    term9f = BIGMVZf(1+(i-1)*xlen:i*xlen,1:xlen,1:xlen);
    
    term10f = BIGMWXf(1+(i-1)*xlen:i*xlen,1:xlen,1:xlen);
    term11f = BIGMWYf(1+(i-1)*xlen:i*xlen,1:xlen,1:xlen);
    term12f = BIGMWZf(1+(i-1)*xlen:i*xlen,1:xlen,1:xlen);
    
    term1b = Bterm1b(1+(i-1)*xlen:i*xlen,1:xlen,1:xlen);
    term2b = Bterm2b(1+(i-1)*xlen:i*xlen,1:xlen,1:xlen);
    term3b = Bterm3b(1+(i-1)*xlen:i*xlen,1:xlen,1:xlen);
    
    term4b = BIGMUXb(1+(i-1)*xlen:i*xlen,1:xlen,1:xlen);
    term5b = BIGMUYb(1+(i-1)*xlen:i*xlen,1:xlen,1:xlen);
    term6b = BIGMUZb(1+(i-1)*xlen:i*xlen,1:xlen,1:xlen);
    
    term7b = BIGMVXb(1+(i-1)*xlen:i*xlen,1:xlen,1:xlen);
    term8b = BIGMVYb(1+(i-1)*xlen:i*xlen,1:xlen,1:xlen);
    term9b = BIGMVZb(1+(i-1)*xlen:i*xlen,1:xlen,1:xlen);
    
    term10b = BIGMWXb(1+(i-1)*xlen:i*xlen,1:xlen,1:xlen);
    term11b = BIGMWYb(1+(i-1)*xlen:i*xlen,1:xlen,1:xlen);
    term12b = BIGMWZb(1+(i-1)*xlen:i*xlen,1:xlen,1:xlen);
    
    term7 = Bseg(1+(i-1)*xlen:i*xlen,1:xlen,1:xlen);
    term8 = Bsegx(1+(i-1)*xlen:i*xlen,1:xlen,1:xlen);
    term9 = Bsegy(1+(i-1)*xlen:i*xlen,1:xlen,1:xlen);
    term10 = Bsegz(1+(i-1)*xlen:i*xlen,1:xlen,1:xlen);
    
    SB = Jm(i).nzsplines;
    supp_phi = PHI{i,1};
    supp_phiu = PHIU{i,1};
    supp_phiv = PHIV{i,1};
    supp_phiw = PHIW{i,1};
    
    supp_size = size(SB,1);
    
    gaussptx = BIGX(1+(i-1)*xlen:i*xlen,1:xlen,1:xlen);
    gausspty = BIGY(1+(i-1)*xlen:i*xlen,1:xlen,1:xlen);
    gaussptz = BIGZ(1+(i-1)*xlen:i*xlen,1:xlen,1:xlen);
    
    fbxx =  FBX(1+(i-1)*xlen:i*xlen,1:xlen,1:xlen);
    fbyy =  FBY(1+(i-1)*xlen:i*xlen,1:xlen,1:xlen);
    fbzz =  FBZ(1+(i-1)*xlen:i*xlen,1:xlen,1:xlen);
    
    bfxx = BFX(1+(i-1)*xlen:i*xlen,1:xlen,1:xlen);
    bfyy = BFY(1+(i-1)*xlen:i*xlen,1:xlen,1:xlen);
    bfzz = BFZ(1+(i-1)*xlen:i*xlen,1:xlen,1:xlen);
    
    fbxdx = FBXDX(1+(i-1)*xlen:i*xlen,1:xlen,1:xlen);
    fbxdy = FBXDY(1+(i-1)*xlen:i*xlen,1:xlen,1:xlen);
    fbxdz = FBXDZ(1+(i-1)*xlen:i*xlen,1:xlen,1:xlen);
    
    fbydx = FBYDX(1+(i-1)*xlen:i*xlen,1:xlen,1:xlen);
    fbydy = FBYDY(1+(i-1)*xlen:i*xlen,1:xlen,1:xlen);
    fbydz = FBYDZ(1+(i-1)*xlen:i*xlen,1:xlen,1:xlen);
    
    fbzdx = FBZDX(1+(i-1)*xlen:i*xlen,1:xlen,1:xlen);
    fbzdy = FBZDY(1+(i-1)*xlen:i*xlen,1:xlen,1:xlen);
    fbzdz = FBZDZ(1+(i-1)*xlen:i*xlen,1:xlen,1:xlen);
    
    bfxdx = BFXDX(1+(i-1)*xlen:i*xlen,1:xlen,1:xlen);
    bfxdy = BFXDY(1+(i-1)*xlen:i*xlen,1:xlen,1:xlen);
    bfxdz = BFXDZ(1+(i-1)*xlen:i*xlen,1:xlen,1:xlen);
    
    bfydx = BFYDX(1+(i-1)*xlen:i*xlen,1:xlen,1:xlen);
    bfydy = BFYDY(1+(i-1)*xlen:i*xlen,1:xlen,1:xlen);
    bfydz = BFYDZ(1+(i-1)*xlen:i*xlen,1:xlen,1:xlen);
    
    bfzdx = BFZDX(1+(i-1)*xlen:i*xlen,1:xlen,1:xlen);
    bfzdy = BFZDY(1+(i-1)*xlen:i*xlen,1:xlen,1:xlen);
    bfzdz = BFZDZ(1+(i-1)*xlen:i*xlen,1:xlen,1:xlen);
    
    supp_phif = PHIF{i,1};
    supp_phiuf = PHIUF{i,1};
    supp_phivf = PHIVF{i,1};
    supp_phiwf = PHIWF{i,1};
    
    supp_phib = PHIB{i,1};
    supp_phiub = PHIUB{i,1};
    supp_phivb = PHIVB{i,1};
    supp_phiwb = PHIWB{i,1};
    
    vaaal1 = zeros(supp_size,1);
    vaaal2 = zeros(supp_size,1);
    vaaal3 = zeros(supp_size,1);
    
    for bg = 1:supp_size
        %BBMvector = Bvect{SB(bg,2),1};
        valm1f = zeros(xlen,xlen,xlen);
        valm2f = zeros(xlen,xlen,xlen);
        valm3f = zeros(xlen,xlen,xlen);
        
        valm1b = zeros(xlen,xlen,xlen);
        valm2b = zeros(xlen,xlen,xlen);
        valm3b = zeros(xlen,xlen,xlen);
        
        for gg1 = 1:xlen,
            for gg2 = 1:xlen,
                for gg3 = 1:xlen
           
                    phi_i = supp_phi(bg,gg1,gg2,gg3);
                    phi_ui = supp_phiu(bg,gg1,gg2,gg3);
                    phi_vi = supp_phiv(bg,gg1,gg2,gg3);
                    phi_wi = supp_phiw(bg,gg1,gg2,gg3);
           
                    phi_if = supp_phif(bg,gg1,gg2,gg3);
                    phi_uif = supp_phiuf(bg,gg1,gg2,gg3);
                    phi_vif = supp_phivf(bg,gg1,gg2,gg3);
                    phi_wif = supp_phiwf(bg,gg1,gg2,gg3);
                    
                    phi_ib = supp_phib(bg,gg1,gg2,gg3);
                    phi_uib = supp_phiub(bg,gg1,gg2,gg3);
                    phi_vib = supp_phivb(bg,gg1,gg2,gg3);
                    phi_wib = supp_phiwb(bg,gg1,gg2,gg3);
                    
                    arxf = 0;%term4f(gg1,gg2)*(term5f(gg1,gg2).^2 + term6f(gg1,gg2).^2) - (term3f(gg1,gg2)*term5f(gg1,gg2) + term4f(gg1,gg2)*term6f(gg1,gg2))*term5f(gg1,gg2);
                    brxf = 0;%term7f(gg1,gg2)*(term3f(gg1,gg2).^2 + term4f(gg1,gg2).^2) - (term3f(gg1,gg2)*term5f(gg1,gg2) + term4f(gg1,gg2)*term6f(gg1,gg2))*term3f(gg1,gg2);
                    crxf = 0;%term10f(gg1,gg2)*(term3f(gg1,gg2).^2 + term4f(gg1,gg2).^2) - (term3f(gg1,gg2)*term5f(gg1,gg2) + term4f(gg1,gg2)*term6f(gg1,gg2))*term3f(gg1,gg2);
                    
                    aryf = 0;%term4f(gg1,gg2)*(term5f(gg1,gg2).^2 + term6f(gg1,gg2).^2) - (term3f(gg1,gg2)*term5f(gg1,gg2) + term4f(gg1,gg2)*term6f(gg1,gg2))*term6f(gg1,gg2);
                    bryf = 0;%term6f(gg1,gg2)*(term3f(gg1,gg2).^2 + term4f(gg1,gg2).^2) - (term3f(gg1,gg2)*term5f(gg1,gg2) + term4f(gg1,gg2)*term6f(gg1,gg2))*term4f(gg1,gg2);
                    brxf = 0;%term5f(gg1,gg2)*(term3f(gg1,gg2).^2 + term4f(gg1,gg2).^2) - (term3f(gg1,gg2)*term5f(gg1,gg2) + term4f(gg1,gg2)*term6f(gg1,gg2))*term3f(gg1,gg2);
                    
                    arxb = 0;%term3b(gg1,gg2)*(term5b(gg1,gg2).^2 + term6b(gg1,gg2).^2) - (term3b(gg1,gg2)*term5b(gg1,gg2) + term4b(gg1,gg2)*term6b(gg1,gg2))*term5b(gg1,gg2);
                    brxb = 0;%term5b(gg1,gg2)*(term3b(gg1,gg2).^2 + term4b(gg1,gg2).^2) - (term3b(gg1,gg2)*term5b(gg1,gg2) + term4b(gg1,gg2)*term6b(gg1,gg2))*term3b(gg1,gg2);
                    
                    aryb = 0;%term4b(gg1,gg2)*(term5b(gg1,gg2).^2 + term6b(gg1,gg2).^2) - (term3b(gg1,gg2)*term5b(gg1,gg2) + term4b(gg1,gg2)*term6b(gg1,gg2))*term6b(gg1,gg2);
                    bryb = 0;%term6b(gg1,gg2)*(term3b(gg1,gg2).^2 + term4b(gg1,gg2).^2) - (term3b(gg1,gg2)*term5b(gg1,gg2) + term4b(gg1,gg2)*term6b(gg1,gg2))*term4b(gg1,gg2);
                    
                    xfb = [gaussptx(gg1,gg2,gg3)-fbxx(gg1,gg2,gg3), gausspty(gg1,gg2,gg3)-fbyy(gg1,gg2,gg3),gaussptz(gg1,gg2,gg3)-fbzz(gg1,gg2,gg3)];
                    xbf = [gaussptx(gg1,gg2,gg3)-bfxx(gg1,gg2,gg3), gausspty(gg1,gg2,gg3)-bfyy(gg1,gg2,gg3),gaussptz(gg1,gg2,gg3)-bfzz(gg1,gg2,gg3)];
                    
                    Jfb = [fbxdx(gg1,gg2,gg3), fbxdy(gg1,gg2,gg3), fbxdz(gg1,gg2,gg3); fbydx(gg1,gg2,gg3), fbydy(gg1,gg2,gg3), fbydz(gg1,gg2,gg3);fbzdx(gg1,gg2,gg3), fbzdy(gg1,gg2,gg3), fbzdz(gg1,gg2,gg3)];
                    Jbf = [bfxdx(gg1,gg2,gg3), bfxdy(gg1,gg2,gg3), bfxdz(gg1,gg2,gg3); bfydx(gg1,gg2,gg3), bfydy(gg1,gg2,gg3), bfydz(gg1,gg2,gg3);bfzdx(gg1,gg2,gg3), bfzdy(gg1,gg2,gg3), bfzdz(gg1,gg2,gg3)];
                    
                    regfb = 2*param.mu*xfb*Jfb;
                    regbf = 2*param.mu*xbf*Jbf;
                    
                    reg1fb = regfb(1,1)*phi_i + 2*param.mu*xfb(1,1)*phi_if;
                    reg2fb = regfb(1,2)*phi_i + 2*param.mu*xfb(1,2)*phi_if;
                    reg3fb = regfb(1,3)*phi_i + 2*param.mu*xfb(1,3)*phi_if;
                    
                    reg1bf = 2*param.mu*xbf(1,1)*phi_ib + regbf(1,1)*phi_i;
                    reg2bf = 2*param.mu*xbf(1,2)*phi_ib + regbf(1,2)*phi_i;
                    reg3bf = 2*param.mu*xbf(1,3)*phi_ib + regbf(1,3)*phi_i;
                    
                    valm1f(gg1,gg2,gg3) = param.par1*phi_i*term7(gg1,gg2,gg3)*term8(gg1,gg2,gg3) + param.par2*(phi_i)*(term1f(gg1,gg2,gg3))+ param.lambda_1*(term4f(gg1,gg2,gg3)*phi_ui+term5f(gg1,gg2,gg3)*phi_vi+term6f(gg1,gg2,gg3)*phi_wi)+param.lambda_2*(arxf*phi_ui+brxf*phi_vi)+ param.lambda_3*reg1fb;
                    valm2f(gg1,gg2,gg3) = param.par1*phi_i*term7(gg1,gg2,gg3)*term9(gg1,gg2,gg3) + param.par2*(phi_i)*(term2f(gg1,gg2,gg3))+ param.lambda_1*(term7f(gg1,gg2,gg3)*phi_ui+term8f(gg1,gg2,gg3)*phi_vi+term9f(gg1,gg2,gg3)*phi_wi)+param.lambda_2*(aryf*phi_ui+bryf*phi_vi)+ param.lambda_3*reg2fb;
                    valm3f(gg1,gg2,gg3) = param.par1*phi_i*term7(gg1,gg2,gg3)*term10(gg1,gg2,gg3) + param.par2*(phi_i)*(term3f(gg1,gg2,gg3))+ param.lambda_1*(term10f(gg1,gg2,gg3)*phi_ui+term11f(gg1,gg2,gg3)*phi_vi+term12f(gg1,gg2,gg3)*phi_wi)+param.lambda_2*(aryf*phi_ui+bryf*phi_vi)+ param.lambda_3*reg2fb;
                    
                    valm1b(gg1,gg2,gg3) = param.par2*phi_i*term1b(gg1,gg2,gg3)+ param.lambda_1*(term4b(gg1,gg2,gg3)*phi_ui+term5b(gg1,gg2,gg3)*phi_vi+term6b(gg1,gg2,gg3)*phi_wi)+param.lambda_2*(arxb*phi_ui+brxb*phi_vi)+ param.lambda_3*reg1bf;
                    valm2b(gg1,gg2,gg3) = param.par2*phi_i*term2b(gg1,gg2,gg3)+ param.lambda_1*(term7b(gg1,gg2,gg3)*phi_ui+term8b(gg1,gg2,gg3)*phi_vi+term9b(gg1,gg2,gg3)*phi_wi)+param.lambda_2*(aryb*phi_ui+bryb*phi_vi)+ param.lambda_3*reg2bf;
                    valm3b(gg1,gg2,gg3) = param.par2*phi_i*term3b(gg1,gg2,gg3)+ param.lambda_1*(term10b(gg1,gg2,gg3)*phi_ui+term11b(gg1,gg2,gg3)*phi_vi+term12b(gg1,gg2,gg3)*phi_wi)+param.lambda_2*(aryb*phi_ui+bryb*phi_vi)+ param.lambda_3*reg2bf;
                end
            end
        end
            
        val1f = 0;
        val2f = 0;
        val3f = 0;
        
        val1b = 0;
        val2b = 0;
        val3b = 0;
        
        for aa=1:xlen
            for bb = 1:xlen
                for cc = 1:xlen
                    val1f = val1f + w1(aa) * valm1f(aa,bb,cc) * w2(bb) * w3(cc) * H(i,1) * H(i,2) * H(i,3);
                    val2f = val2f + w1(aa) * valm2f(aa,bb,cc) * w2(bb) * w3(cc) * H(i,1) * H(i,2) * H(i,3);
                    val3f = val3f + w1(aa) * valm3f(aa,bb,cc) * w2(bb) * w3(cc) * H(i,1) * H(i,2) * H(i,3);
                    
                    val1b = val1b + w1(aa) * valm1b(aa,bb,cc) * w2(bb) * w3(cc) * H(i,1) * H(i,2) * H(i,3);
                    val2b = val2b + w1(aa) * valm2b(aa,bb,cc) * w2(bb) * w3(cc) * H(i,1) * H(i,2) * H(i,3);
                    val3b = val3b + w1(aa) * valm3b(aa,bb,cc) * w2(bb) * w3(cc) * H(i,1) * H(i,2) * H(i,3);
                end
            end
        end
        
        vaaal1(bg,1) = val1f;
        vaaal2(bg,1) = val2f;
        vaaal3(bg,1) = val3f;
    end
    

    
    RHSf1(SB,1) = vaaal1;
    RHSf1(SB,2) = vaaal2;
    RHSf1(SB,3) = vaaal3;
    
    RHSb1(SB,1) = val1b;
    RHSb1(SB,2) = val2b;
    RHSb1(SB,3) = val3b;
    
    RHSf = RHSf  + RHSf1;
    RHSb = RHSb  + RHSb1;
    
    VV1{i,1} = SB;
    VV2{i,1} = RHSf1;
    VV3{i,1} = vaaal1;
    VV3{i,2} = vaaal2;
    VV3{i,3} = vaaal3;
end