parfor i = 1:ac_ct,
    
    RHSf1 = zeros(bf_ct,4);
    RHSb1 = zeros(bf_ct,4);
    
    SB = Jm(i).nzsplines;
    supp_phi = PHI1(i).mat;
    supp_phiu = PHIU1(i).mat;
    supp_phiv = PHIV1(i).mat;
    supp_phiw = PHIW1(i).mat;
    supp_size = size(SB,1);
    
    supp_phif = PHIF1(i).mat;
    supp_phiuf = PHIFU1(i).mat;
    supp_phivf = PHIFV1(i).mat;
    supp_phiwf = PHIFW1(i).mat;
    
    supp_phib = PHIB1(i).mat;
    supp_phiub = PHIBU1(i).mat;
    supp_phivb = PHIBV1(i).mat;
    supp_phiwb = PHIBW1(i).mat;
    
    hu = H(i,1);
    hv = H(i,2);
    hw = H(i,3);
    
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
    
    fu_normf = term4f.^2 + term5f.^2 + term6f.^2;
    fv_normf = term7f.^2 + term8f.^2 + term9f.^2;
    fw_normf = term10f.^2 + term11f.^2 + term12f.^2;
    
    fu_fvf = term4f.*term7f+term5f.*term8f+term6f.*term9f;
    fv_fwf = term7f.*term10f+term8f.*term11f+term9f.*term12f;
    fw_fuf = term10f.*term4f+term11f.*term5f+term12f.*term6f;
    
    fu_normb = term4b.^2 + term5b.^2 + term6b.^2;
    fv_normb = term7b.^2 + term8b.^2 + term9b.^2;
    fw_normb = term10b.^2 + term11b.^2 + term12b.^2;
    
    fu_fvb = term4b.*term7b+term5b.*term8b+term6b.*term9b;
    fv_fwb = term7b.*term10b+term8b.*term11b+term9b.*term12b;
    fw_fub = term10b.*term4b+term11b.*term5b+term12b.*term6b;
    
    arxf = 2.*(term4f.*(fv_normf) - term7f.*fu_fvf + term4f.*fw_normf - term10f.*fw_fuf);
    brxf = 2.*(term7f.*(fu_normf) - term4f.*fu_fvf + term7f.*fw_normf - term10f.*fv_fwf);
    crxf = 2.*(term10f.*(fv_normf) - term7f.*fv_fwf + term10f.*fu_normf - term4f.*fw_fuf);
    
    aryf = 2.*(term5f.*(fv_normf) - term8f.*fu_fvf + term5f.*fw_normf - term11f.*fw_fuf);
    bryf = 2.*(term8f.*(fu_normf) - term5f.*fu_fvf + term8f.*fw_normf - term11f.*fv_fwf);
    cryf = 2.*(term11f.*(fv_normf) - term8f.*fv_fwf + term11f.*fu_normf - term5f.*fw_fuf);
    
    arzf = 2.*(term6f.*(fv_normf) - term9f.*fu_fvf + term6f.*fw_normf - term12f.*fw_fuf);
    brzf = 2.*(term9f.*(fu_normf) - term6f.*fu_fvf + term9f.*fw_normf - term12f.*fv_fwf);
    crzf = 2.*(term12f.*(fv_normf) - term9f.*fv_fwf + term12f.*fu_normf - term6f.*fw_fuf);
    
    arxb = 2.*(term4b.*(fv_normb) - term7b.*fu_fvb + term4b.*fw_normb - term10b.*fw_fub);
    brxb = 2.*(term7b.*(fu_normb) - term4b.*fu_fvb + term7b.*fw_normb - term10b.*fv_fwb);
    crxb = 2.*(term10b.*(fv_normb) - term7b.*fv_fwb + term10b.*fu_normb - term4b.*fw_fub);
    
    aryb = 2.*(term5b.*(fv_normb) - term8b.*fu_fvb + term5b.*fw_normb - term11b.*fw_fub);
    bryb = 2.*(term8b.*(fu_normb) - term5b.*fu_fvb + term8b.*fw_normb - term11b.*fv_fwb);
    cryb = 2.*(term11b.*(fv_normb) - term8b.*fv_fwb + term11b.*fu_normb - term5b.*fw_fub);
    
    arzb = 2.*(term6b.*(fv_normb) - term9b.*fu_fvb + term6b.*fw_normb - term12b.*fw_fub);
    brzb = 2.*(term9b.*(fu_normb) - term6b.*fu_fvb + term9b.*fw_normb - term12b.*fv_fwb);
    crzb = 2.*(term12b.*(fv_normb) - term9b.*fv_fwb + term12b.*fu_normb - term6b.*fw_fub);
    
    val1f = zeros(supp_size,1);
    val2f = zeros(supp_size,1);
    val3f = zeros(supp_size,1);
    
    val1b = zeros(supp_size,1);
    val2b = zeros(supp_size,1);
    val3b = zeros(supp_size,1);
    
    valm1f = zeros(supp_size,xlen,xlen,xlen);
    valm2f = zeros(supp_size,xlen,xlen,xlen);
    valm3f = zeros(supp_size,xlen,xlen,xlen);
    
    valm1b = zeros(supp_size,xlen,xlen,xlen);
    valm2b = zeros(supp_size,xlen,xlen,xlen);
    valm3b = zeros(supp_size,xlen,xlen,xlen);
    
    for gg1 = 1:xlen,
        for gg2 = 1:xlen,
            for gg3 = 1:xlen
                phi_i  = supp_phi(:,gg1,gg2,gg3);
                phi_ui = supp_phiu(:,gg1,gg2,gg3);
                phi_vi = supp_phiv(:,gg1,gg2,gg3);
                phi_wi = supp_phiw(:,gg1,gg2,gg3);
                
                phi_if = supp_phif(:,gg1,gg2,gg3);
                phi_uif = supp_phiuf(:,gg1,gg2,gg3);
                phi_vif = supp_phivf(:,gg1,gg2,gg3);
                phi_wif = supp_phiwf(:,gg1,gg2,gg3);
                
                phi_ib = supp_phib(:,gg1,gg2,gg3);
                phi_uib = supp_phiub(:,gg1,gg2,gg3);
                phi_vib = supp_phivb(:,gg1,gg2,gg3);
                phi_wib = supp_phiwb(:,gg1,gg2,gg3);
                
                xfb = [gaussptx(gg1,gg2,gg3)-fbxx(gg1,gg2,gg3), gausspty(gg1,gg2,gg3)-fbyy(gg1,gg2,gg3),gaussptz(gg1,gg2,gg3)-fbzz(gg1,gg2,gg3)];
                xbf = [gaussptx(gg1,gg2,gg3)-bfxx(gg1,gg2,gg3), gausspty(gg1,gg2,gg3)-bfyy(gg1,gg2,gg3),gaussptz(gg1,gg2,gg3)-bfzz(gg1,gg2,gg3)];
                
                Jfb = [fbxdx(gg1,gg2,gg3), fbxdy(gg1,gg2,gg3), fbxdz(gg1,gg2,gg3); fbydx(gg1,gg2,gg3), fbydy(gg1,gg2,gg3), fbydz(gg1,gg2,gg3);fbzdx(gg1,gg2,gg3), fbzdy(gg1,gg2,gg3), fbzdz(gg1,gg2,gg3)];
                Jbf = [bfxdx(gg1,gg2,gg3), bfxdy(gg1,gg2,gg3), bfxdz(gg1,gg2,gg3); bfydx(gg1,gg2,gg3), bfydy(gg1,gg2,gg3), bfydz(gg1,gg2,gg3);bfzdx(gg1,gg2,gg3), bfzdy(gg1,gg2,gg3), bfzdz(gg1,gg2,gg3)];
                
                
                regfb = 2*param.mu*xfb*Jfb;
                regbf = 2*param.mu*xbf*Jbf;
                
                reg1fb = regfb(1,1).*phi_i + 2*param.mu*xfb(1,1).*phi_if;
                reg2fb = regfb(1,2).*phi_i + 2*param.mu*xfb(1,2).*phi_if;
                reg3fb = regfb(1,3).*phi_i + 2*param.mu*xfb(1,3).*phi_if;
                
                reg1bf = 2*param.mu*xbf(1,1).*phi_ib + regbf(1,1).*phi_i;
                reg2bf = 2*param.mu*xbf(1,2).*phi_ib + regbf(1,2).*phi_i;
                reg3bf = 2*param.mu*xbf(1,3).*phi_ib + regbf(1,3).*phi_i;
                
                valm1f(:,gg1,gg2,gg3) = param.par1*phi_i*term7(gg1,gg2,gg3)*term8(gg1,gg2,gg3) + param.par2*(phi_i)*(term1f(gg1,gg2,gg3))+ param.lambda_1*(term4f(gg1,gg2,gg3).*phi_ui+term5f(gg1,gg2,gg3).*phi_vi+term6f(gg1,gg2,gg3).*phi_wi) + param.lambda_2*(arxf(gg1,gg2,gg3).*phi_ui+brxf(gg1,gg2,gg3).*phi_vi+crxf(gg1,gg2,gg3).*phi_wi)+ param.lambda_3.*reg1fb;
                valm2f(:,gg1,gg2,gg3) = param.par1*phi_i*term7(gg1,gg2,gg3)*term9(gg1,gg2,gg3) + param.par2*(phi_i)*(term2f(gg1,gg2,gg3))+ param.lambda_1*(term7f(gg1,gg2,gg3).*phi_ui+term8f(gg1,gg2,gg3).*phi_vi+term9f(gg1,gg2,gg3).*phi_wi) + param.lambda_2*(aryf(gg1,gg2,gg3).*phi_ui+bryf(gg1,gg2,gg3).*phi_vi+cryf(gg1,gg2,gg3).*phi_wi)+ param.lambda_3.*reg2fb;
                valm3f(:,gg1,gg2,gg3) = param.par1*phi_i*term7(gg1,gg2,gg3)*term10(gg1,gg2,gg3) + param.par2*(phi_i)*(term3f(gg1,gg2,gg3))+ param.lambda_1*(term10f(gg1,gg2,gg3).*phi_ui+term11f(gg1,gg2,gg3).*phi_vi+term12f(gg1,gg2,gg3).*phi_wi) + param.lambda_2*(arzf(gg1,gg2,gg3).*phi_ui+brzf(gg1,gg2,gg3).*phi_vi+crzf(gg1,gg2,gg3).*phi_wi)+ param.lambda_3.*reg3fb;
                
                valm1b(:,gg1,gg2,gg3) = param.par2*(phi_i)*(term1b(gg1,gg2,gg3))+ param.lambda_1*(term4b(gg1,gg2,gg3).*phi_ui+term5b(gg1,gg2,gg3).*phi_vi+term6b(gg1,gg2,gg3).*phi_wi) + param.lambda_2*(arxb(gg1,gg2,gg3).*phi_ui+brxb(gg1,gg2,gg3).*phi_vi+crxb(gg1,gg2,gg3).*phi_wi)+ param.lambda_3.*reg1fb;
                valm2b(:,gg1,gg2,gg3) = param.par2*(phi_i)*(term2b(gg1,gg2,gg3))+ param.lambda_1*(term7b(gg1,gg2,gg3).*phi_ui+term8b(gg1,gg2,gg3).*phi_vi+term9b(gg1,gg2,gg3).*phi_wi) + param.lambda_2*(aryb(gg1,gg2,gg3).*phi_ui+bryb(gg1,gg2,gg3).*phi_vi+cryb(gg1,gg2,gg3).*phi_wi)+ param.lambda_3.*reg2fb;
                valm3b(:,gg1,gg2,gg3) = param.par2*(phi_i)*(term3b(gg1,gg2,gg3))+ param.lambda_1*(term10b(gg1,gg2,gg3).*phi_ui+term11b(gg1,gg2,gg3).*phi_vi+term12b(gg1,gg2,gg3).*phi_wi) + param.lambda_2*(arzb(gg1,gg2,gg3).*phi_ui+brzb(gg1,gg2,gg3).*phi_vi+crzb(gg1,gg2,gg3).*phi_wi)+ param.lambda_3.*reg3fb;
                
                val1f(:,1) = val1f(:,1) + w1(gg1,1).*w2(gg2,1).*w3(gg3,1).*valm1f(:,gg1,gg2,gg3).*hu.*hv.*hw;
                val2f(:,1) = val2f(:,1) + w1(gg1,1).*w2(gg2,1).*w3(gg3,1).*valm2f(:,gg1,gg2,gg3).*hu.*hv.*hw;
                val3f(:,1) = val3f(:,1) + w1(gg1,1).*w2(gg2,1).*w3(gg3,1).*valm3f(:,gg1,gg2,gg3).*hu.*hv.*hw;
                
                val1b(:,1) = val1b(:,1) + w1(gg1,1).*w2(gg2,1).*w3(gg3,1).*valm1b(:,gg1,gg2,gg3).*hu.*hv.*hw;
                val2b(:,1) = val2b(:,1) + w1(gg1,1).*w2(gg2,1).*w3(gg3,1).*valm2b(:,gg1,gg2,gg3).*hu.*hv.*hw;
                val3b(:,1) = val3b(:,1) + w1(gg1,1).*w2(gg2,1).*w3(gg3,1).*valm3b(:,gg1,gg2,gg3).*hu.*hv.*hw;
            end
        end
    end
    
    RHSf1(SB,1) = val1f;
    RHSf1(SB,2) = val2f;
    RHSf1(SB,3) = val3f;
    
    RHSb1(SB,1) = val1b;
    RHSb1(SB,2) = val2b;
    RHSb1(SB,3) = val3b;
    
    RHSf = RHSf  + RHSf1;
    RHSb = RHSb  + RHSb1;
    
end


