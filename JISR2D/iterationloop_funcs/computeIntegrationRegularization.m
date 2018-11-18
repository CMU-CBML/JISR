function [RHSff,RHSbb,PHIF] = computeIntegrationRegularization(parameters,H,BIGXf,BIGYf,BIGXb,BIGYb,Jm,Dm, PHII,BIGX,BIGY, FBX,BFX, PHIFD,PHIBD,RHSf,RHSb)

ac_ct = size(Jm,1);
xlen = 6;
[s1,w1] = ggquad(xlen);
w2 = w1;
BIGMUXf = BIGXf(:,:,2);
BIGMUYf = BIGXf(:,:,3);
BIGMVXf = BIGYf(:,:,2);
BIGMVYf = BIGYf(:,:,3);

BIGMUXb = BIGXb(:,:,2);
BIGMUYb = BIGXb(:,:,3);
BIGMVXb = BIGYb(:,:,2);
BIGMVYb = BIGYb(:,:,3);

PHIF = PHIFD{1,1};
PHIB=PHIBD{1,1};
PHI = PHII{1,1};
PHUI = PHII{1,2};
PHVI = PHII{1,3};
for i = 1:ac_ct,
    
    term3f = BIGMUXf(1+(i-1)*xlen:i*xlen,1:6);
    term4f = BIGMUYf(1+(i-1)*xlen:i*xlen,1:6);
    term5f = BIGMVXf(1+(i-1)*xlen:i*xlen,1:6);
    term6f = BIGMVYf(1+(i-1)*xlen:i*xlen,1:6);
    
    
    term3b = BIGMUXb(1+(i-1)*xlen:i*xlen,1:6);
    term4b = BIGMUYb(1+(i-1)*xlen:i*xlen,1:6);
    term5b = BIGMVXb(1+(i-1)*xlen:i*xlen,1:6);
    term6b = BIGMVYb(1+(i-1)*xlen:i*xlen,1:6);
    
    
    SB = Jm{i,1};
    supp_phi = PHI{i,1};
    supp_phiu = PHUI{i,1};
    supp_phiv = PHVI{i,1};
    supp_size = size(SB,1);
    
    gaussptx = BIGX(1+(i-1)*xlen:i*xlen,1:6);
    gausspty = BIGY(1+(i-1)*xlen:i*xlen,1:6);
    
    fbxx =  FBX(1+(i-1)*xlen:i*xlen,1:6,2);
    fbyy =  FBX(1+(i-1)*xlen:i*xlen,1:6,1);
    
    bfxx = BFX(1+(i-1)*xlen:i*xlen,1:6,2);
    bfyy = BFX(1+(i-1)*xlen:i*xlen,1:6,1);
    
    fbxdx = FBX(1+(i-1)*xlen:i*xlen,1:6,3);
    fbxdy = FBX(1+(i-1)*xlen:i*xlen,1:6,4);
    
    fbydx = FBX(1+(i-1)*xlen:i*xlen,1:6,5);
    fbydy = FBX(1+(i-1)*xlen:i*xlen,1:6,6);
    
    bfxdx = BFX(1+(i-1)*xlen:i*xlen,1:6,3);
    bfxdy = BFX(1+(i-1)*xlen:i*xlen,1:6,4);
    
    bfydx = BFX(1+(i-1)*xlen:i*xlen,1:6,5);
    bfydy = BFX(1+(i-1)*xlen:i*xlen,1:6,6);
    
    supp_phif = PHIF{i,1};
    supp_phib = PHIB{i,1};

    for bg = 1:supp_size,
        %BBMvector = Bvect{SB(bg,2),1};
        valm1f = zeros(xlen,xlen);
        valm2f = zeros(xlen,xlen);
        
        valm1b = zeros(xlen,xlen);
        valm2b = zeros(xlen,xlen);
        
        for gg1 = 1:xlen,
            for gg2 = 1:xlen,
                
                phi_i = supp_phi(bg,gg1,gg2);
                phi_ui = supp_phiu(bg,gg1,gg2);
                phi_vi = supp_phiv(bg,gg1,gg2);
                
                phi_if = supp_phif(bg,gg1,gg2);
                phi_ib = supp_phib(bg,gg1,gg2);

                arxf = term3f(gg1,gg2)*(term5f(gg1,gg2).^2 + term6f(gg1,gg2).^2) - (term3f(gg1,gg2)*term5f(gg1,gg2) + term4f(gg1,gg2)*term6f(gg1,gg2))*term5f(gg1,gg2);
                brxf = term5f(gg1,gg2)*(term3f(gg1,gg2).^2 + term4f(gg1,gg2).^2) - (term3f(gg1,gg2)*term5f(gg1,gg2) + term4f(gg1,gg2)*term6f(gg1,gg2))*term3f(gg1,gg2);
                aryf = term4f(gg1,gg2)*(term5f(gg1,gg2).^2 + term6f(gg1,gg2).^2) - (term3f(gg1,gg2)*term5f(gg1,gg2) + term4f(gg1,gg2)*term6f(gg1,gg2))*term6f(gg1,gg2);
                bryf = term6f(gg1,gg2)*(term3f(gg1,gg2).^2 + term4f(gg1,gg2).^2) - (term3f(gg1,gg2)*term5f(gg1,gg2) + term4f(gg1,gg2)*term6f(gg1,gg2))*term4f(gg1,gg2);
                
                arxb = term3b(gg1,gg2)*(term5b(gg1,gg2).^2 + term6b(gg1,gg2).^2) - (term3b(gg1,gg2)*term5b(gg1,gg2) + term4b(gg1,gg2)*term6b(gg1,gg2))*term5b(gg1,gg2);
                brxb = term5b(gg1,gg2)*(term3b(gg1,gg2).^2 + term4b(gg1,gg2).^2) - (term3b(gg1,gg2)*term5b(gg1,gg2) + term4b(gg1,gg2)*term6b(gg1,gg2))*term3b(gg1,gg2);
                aryb = term4b(gg1,gg2)*(term5b(gg1,gg2).^2 + term6b(gg1,gg2).^2) - (term3b(gg1,gg2)*term5b(gg1,gg2) + term4b(gg1,gg2)*term6b(gg1,gg2))*term6b(gg1,gg2);
                bryb = term6b(gg1,gg2)*(term3b(gg1,gg2).^2 + term4b(gg1,gg2).^2) - (term3b(gg1,gg2)*term5b(gg1,gg2) + term4b(gg1,gg2)*term6b(gg1,gg2))*term4b(gg1,gg2);
                
                xfb = [gaussptx(gg1,gg2)-fbxx(gg1,gg2), gausspty(gg1,gg2)-fbyy(gg1,gg2)];
                xbf = [gaussptx(gg1,gg2)-bfxx(gg1,gg2), gausspty(gg1,gg2)-bfyy(gg1,gg2)];
                
                Jfb = [fbxdx(gg1,gg2), fbxdy(gg1,gg2); fbydx(gg1,gg2), fbydy(gg1,gg2)];
                Jbf = [bfxdx(gg1,gg2), bfxdy(gg1,gg2); bfydx(gg1,gg2), bfydy(gg1,gg2)];
                
                regfb = 2*xfb*Jfb;
                regbf = 2*xbf*Jbf;
                
                reg1fb = regfb(1,1)*phi_i + 2*xfb(1,1)*phi_if;
                reg2fb = regfb(1,2)*phi_i + 2*xfb(1,2)*phi_if;
                
                reg1bf = 2*xbf(1,1)*phi_ib + regbf(1,1)*phi_i;
                reg2bf = 2*xbf(1,2)*phi_ib + regbf(1,2)*phi_i;
                
                valm1f(gg1,gg2) = parameters.lambda1*(term3f(gg1,gg2)*phi_ui+term5f(gg1,gg2)*phi_vi)+parameters.lambda2*(arxf*phi_ui+brxf*phi_vi)+ parameters.lambda3*reg1fb;
                valm2f(gg1,gg2) = parameters.lambda1*(term4f(gg1,gg2)*phi_ui+term6f(gg1,gg2)*phi_vi)+parameters.lambda2*(aryf*phi_ui+bryf*phi_vi)+ parameters.lambda3*reg2fb;
                
                valm1b(gg1,gg2) = parameters.lambda1*(term3b(gg1,gg2)*phi_ui+term5b(gg1,gg2)*phi_vi)+parameters.lambda2*(arxb*phi_ui+brxb*phi_vi)+ parameters.lambda3*reg1bf;
                valm2b(gg1,gg2) = parameters.lambda1*(term4b(gg1,gg2)*phi_ui+term6b(gg1,gg2)*phi_vi)+parameters.lambda2*(aryb*phi_ui+bryb*phi_vi)+ parameters.lambda3*reg2bf;
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
        RHSf(bact_i,1) = RHSf(bact_i,1) + val1f;
        RHSf(bact_i,2) = RHSf(bact_i,2) + val2f;
        
        RHSb(bact_i,1) = RHSb(bact_i,1) + val1b;
        RHSb(bact_i,2) = RHSb(bact_i,2) + val2b;
    end
end

RHSff = RHSf;
RHSbb = RHSb;
end