while(abs(rs_final-rs_initial)>tol)
    
    tic
    rs_initial = rs_final;
    
    Pmf = Pmold;
    Pmb = Pmold;
 
    ACPf = ACPold;
    ACPb = ACPold;
    ACCf = ACPold;
    ACCb = ACPold;
    
    RHSf = RHS_init;
    RHSb = RHS_init;
    Bvectf = RHS_init;
    Bvectb = RHS_init;
    
    H_phi = regular_Heiviside_fun3D(phi);
    
    c1 = (sum(sum(sum(F_orig.*H_phi))) * volume)./(sum(sum(sum(H_phi))) * volume);
    c2 = (sum(sum(sum(F_orig.*(1-H_phi)))) * volume)./(sum(sum(sum(1-H_phi))) * volume);
    
    [D1E] = Delta_fun3D(phi);
    vecD = D1E(:)';
    coef_D = img2coef3D(vecD,nx,ny,nz);
    CID = coef_D(:)';
    coef_matD = reshape(CID, nx+img_pU,ny+img_pV, nz+img_pW);
    
    iterct = iterct +1;
    iteration = iteration+1;
    disp('iterct:');
    disp(iterct);
    toc
    tic
    
    [cI1f, cDI1Yf, cDI1Xf, cDI1Zf] = BsplineComposeImage3D_single(BIGY, BIGX, BIGZ, coef_matM, size(BIGX,1), size(BIGX,2), size(BIGX,3));
    [cI1b, cDI1Yb, cDI1Xb, cDI1Zb] = BsplineComposeImage3D_single(BIGY, BIGX, BIGZ, coef_matF, size(BIGX,1), size(BIGX,2), size(BIGX,3));
    disp('cI1');
    toc
    
    tic
    [ctempp,cDphiY, cDphiX,cDphiZ] = BsplineComposeImage3D_single(BIGY, BIGX, BIGZ, coef_matP, size(BIGX,1), size(BIGX,2), size(BIGX,3));
    [dDelta, dt, dt1, dt2] = BsplineComposeImage3D_single(BIGY, BIGX, BIGZ, coef_matD, size(BIGX,1), size(BIGX,2), size(BIGX,3));
    disp('cphi');
    toc
    
    tic
    denominatef = sqrt(cDI1Xf.^2 + cDI1Yf.^2 + cDI1Zf.^2 + smallNumber);
    denominateb = sqrt(cDI1Xb.^2 + cDI1Yb.^2 + cDI1Zb.^2 + smallNumber);
    
    Bterm1f = (cI1f - cI0f).*2.*cDI1Yf./denominatef;
    Bterm2f = (cI1f - cI0f).*2.*cDI1Xf./denominatef;
    Bterm3f = (cI1f - cI0f).*2.*cDI1Zf./denominatef;
    
    Bterm1b = (cI1b - cI0b).*2.*cDI1Yb./denominateb;
    Bterm2b = (cI1b - cI0b).*2.*cDI1Xb./denominateb;
    Bterm3b = (cI1b - cI0b).*2.*cDI1Zb./denominateb;
    
    %denominate_phi = sqrt((cDphiX.^2) + (cDphiY.^2) + smallNumber);
    
    Bseg = ((-c1 + cI0f).^2 - (-c2 + cI0f).^2).*dDelta;
    Bsegx = cDphiY;
    Bsegy = cDphiX;
    Bsegz = cDphiZ;
    
    disp('Bterm1');
    toc
    tic
    
    RHSf = compute_Integ_Domainf_fid(Jm,Bseg, Bsegx,Bsegy,Bsegz,Bterm1f,Bterm2f,Bterm3f,RHSf,PHI1,param.par1,param.par2,Wu,Wv,Ww,H);
    RHSb = compute_Integ_Domainb_fid(Jm,Bterm1b,Bterm2b,Bterm3b, RHSb,PHI1,param.par2,Wu,Wv,Ww,H);
    toc
    tic
    RHSf = bcondition3D(RHSf);
    RHSb = bcondition3D(RHSb);
    disp('integration done')
    toc
    tic
    ACCf(:,1:3) = ACCf(:,1:3) - timestep.*RHSf(:,1:3);
    ACCb(:,1:3) = ACCb(:,1:3) - timestep.*RHSb(:,1:3);
    tic
    [BIGXXf, BIGYYf, BIGZZf, BIGMUXf, BIGMUYf, BIGMUZf, BIGMVXf, BIGMVYf, BIGMVZf, BIGMWXf, BIGMWYf, BIGMWZf] = computenewPoints(Jm, ACCf, PHI1, PHIU1, PHIV1, PHIW1, orderGauss);
    [BIGXXb, BIGYYb, BIGZZb, BIGMUXb, BIGMUYb, BIGMUZb, BIGMVXb, BIGMVYb, BIGMVZb, BIGMWXb, BIGMWYb, BIGMWZb] = computenewPoints(Jm, ACCb, PHI1, PHIU1, PHIV1, PHIW1, orderGauss);
    disp('BIGXX computed');
    toc
    
    [s1,w1] = ggquad(orderGauss);
    [s2,w2] = ggquad(orderGauss);
    [s3,w3] = ggquad(orderGauss);
    
    [PHIF1,PHIFU1,PHIFV1,PHIFW1] = GaussPhi_fb(ac,BIGXXf, BIGYYf, BIGZZf, Jm,knotvectorU,knotvectorV,knotvectorW, Coeff,param);
    PHIF1 = cell2struct(PHIF1,'mat',2);
    PHIFU1 = cell2struct(PHIFU1,'mat',2);
    PHIFV1 = cell2struct(PHIFV1,'mat',2);
    PHIFW1 = cell2struct(PHIFW1,'mat',2);
    
    [PHIB1,PHIBU1,PHIBV1,PHIBW1] = GaussPhi_fb(ac,BIGXXb, BIGYYb, BIGZZb, Jm,knotvectorU,knotvectorV,knotvectorW, Coeff,param);
    PHIB1 = cell2struct(PHIB1,'mat',2);
    PHIBU1 = cell2struct(PHIBU1,'mat',2);
    PHIBV1 = cell2struct(PHIBV1,'mat',2);
    PHIBW1 = cell2struct(PHIBW1,'mat',2);
    
    disp('PHIFB calculated');
    toc
    tic
    [FBX, FBY, FBZ, FBXDX, FBXDY, FBXDZ, FBYDX, FBYDY, FBYDZ, FBZDX, FBZDY, FBZDZ] = computenewPoints(Jm, ACCf, PHIF1, PHIFU1, PHIFV1, PHIFW1, orderGauss);
    [BFX, BFY, BFZ, BFXDX, BFXDY, BFXDZ, BFYDX, BFYDY, BFYDZ, BFZDX, BFZDY, BFZDZ] = computenewPoints(Jm, ACCb, PHIB1, PHIBU1, PHIBV1, PHIBW1, orderGauss);
    disp('FBX, BFX calculated');
    tic
    
    
    Bvectf = compute_Integ_Domainf_reg(Jm,BIGX, BIGY, BIGZ,FBX, FBY,FBZ, FBXDX,FBXDY,FBXDZ, FBYDX, FBYDY, FBYDZ, FBZDX, FBZDY,FBZDZ, BFX, BFY, BFZ,BIGMUXf,BIGMUYf,BIGMUZf,BIGMVXf,BIGMVYf,BIGMVZf,BIGMWXf,BIGMWYf,BIGMWZf,Bvectf,PHI1, PHIU1, PHIV1, PHIW1, PHIF1,param.lambda_1,param.lambda_3,param.mu,Wu,Wv,Ww,H);
    Bvectb = compute_Integ_Domainb_reg(Jm,BIGX, BIGY, BIGZ,FBX, FBY,FBZ, FBXDX,FBXDY,FBXDZ, FBYDX, FBYDY, FBYDZ, FBZDX, FBZDY,FBZDZ, BFX,BFY,BFZ,BIGMUXf,BIGMUYf,BIGMUZf,BIGMVXf,BIGMVYf,BIGMVZf,BIGMWXf,BIGMWYf,BIGMWZf,Bvectb,PHI1, PHIU1, PHIV1, PHIW1, PHIF1,param.lambda_2,param.lambda_3,param.mu,Wu,Wv,Ww,H);
    clear('BIGXXf','BIGYYf','BIGZZf','BIGMUXf','BIGMUYf','BIGMUZf','BIGMVXf','BIGMVYf','BIGMVZf','BIGMWXf','BIGMWYf','BIGMWZf');
    clear('BIGXXb','BIGYYb','BIGZZb','BIGMUXb','BIGMUYb','BIGMUZb','BIGMVXb','BIGMVYb','BIGMVZb','BIGMWXb','BIGMWYb','BIGMWZb');
    toc
    tic
    Bvectf = bcondition3D(Bvectf);
    Bvectb = bcondition3D(Bvectb);
    disp('integration done')
    toc
    tic
    ACPf(:,1:3) = ACCf(:,1:3) - timestep.*Bvectf(:,1:3);
    ACPb(:,1:3) = ACCb(:,1:3) - timestep.*Bvectb(:,1:3);
    tic
    [pxxf,pyyf,pzzf] = tripleIterLoop(sizeImage, Pixel, Jm, ACPf);
    [pxxb,pyyb,pzzb] = tripleIterLoop(sizeImage, Pixel, Jm, ACPb);
    
    Vxf_temp = pxxf-Vx_ident;
    Vyf_temp = pyyf-Vy_ident;
    Vzf_temp = pzzf-Vz_ident;
    
    Vxb_temp = pxxb-Vx_ident;
    Vyb_temp = pyyb-Vy_ident;
    Vzb_temp = pzzb-Vz_ident;
    
    Vxf_old = Vxf_temp + Vx_ident;
    Vyf_old = Vyf_temp + Vy_ident;
    Vzf_old = Vzf_temp + Vz_ident;
    
    Vxb_old = Vxb_temp + Vx_ident;
    Vyb_old = Vyb_temp + Vy_ident;
    Vzb_old = Vzb_temp + Vz_ident;
    
    Vxf_old = reshape(Vxf_old,1,nx*ny*nz);
    Vyf_old = reshape(Vyf_old,1,nx*ny*nz);
    Vzf_old = reshape(Vzf_old,1,nx*ny*nz);
    
    Vxb_old = reshape(Vxb_old,1,nx*ny*nz);
    Vyb_old = reshape(Vyb_old,1,nx*ny*nz);
    Vzb_old = reshape(Vzb_old,1,nx*ny*nz);
    
    coef_xf_fun = reshape(coef_xf, 1, (nx+img_pU)*(ny+img_pV)*(nz+img_pW));
    coef_yf_fun = reshape(coef_yf, 1, (nx+img_pU)*(ny+img_pV)*(nz+img_pW));
    coef_zf_fun = reshape(coef_zf, 1, (nx+img_pU)*(ny+img_pV)*(nz+img_pW));
    
    coef_xb_fun = reshape(coef_xb, 1, (nx+img_pU)*(ny+img_pV)*(nz+img_pW));
    coef_yb_fun = reshape(coef_yb, 1, (nx+img_pU)*(ny+img_pV)*(nz+img_pW));
    coef_zb_fun = reshape(coef_zb, 1, (nx+img_pU)*(ny+img_pV)*(nz+img_pW));
    
    Vxf_old_inv = -Vxf_temp + Vx_ident;
    Vyf_old_inv = -Vyf_temp + Vy_ident;
    Vzf_old_inv = -Vzf_temp + Vz_ident;
    
    Vxb_old_inv = -Vxb_temp + Vx_ident;
    Vyb_old_inv = -Vyb_temp + Vy_ident;
    Vzb_old_inv = -Vzb_temp + Vz_ident;
    
    Vxf_old_inv = reshape(Vxf_old_inv,1,nx*ny*nz);
    Vyf_old_inv = reshape(Vyf_old_inv,1,nx*ny*nz);
    Vzf_old_inv = reshape(Vzf_old_inv,1,nx*ny*nz);
    
    Vxb_old_inv = reshape(Vxb_old_inv,1,nx*ny*nz);
    Vyb_old_inv = reshape(Vyb_old_inv,1,nx*ny*nz);
    Vzb_old_inv = reshape(Vzb_old_inv,1,nx*ny*nz);
    
    coef_xf_inv_fun = reshape(coef_xf_inv, 1, (nx+img_pU)*(ny+img_pV)*(nz+img_pW));
    coef_yf_inv_fun = reshape(coef_yf_inv, 1, (nx+img_pU)*(ny+img_pV)*(nz+img_pW));
    coef_zf_inv_fun = reshape(coef_zf_inv, 1, (nx+img_pU)*(ny+img_pV)*(nz+img_pW));
    
    coef_xb_inv_fun = reshape(coef_xb_inv, 1, (nx+img_pU)*(ny+img_pV)*(nz+img_pW));
    coef_yb_inv_fun = reshape(coef_yb_inv, 1, (nx+img_pU)*(ny+img_pV)*(nz+img_pW));
    coef_zb_inv_fun = reshape(coef_zb_inv, 1, (nx+img_pU)*(ny+img_pV)*(nz+img_pW));
    
    [VXf_new, VYf_new, VZf_new] = BsplineCompose3D(nx, ny, nz, Vxf_old, Vyf_old, Vzf_old, coef_xf_fun, coef_yf_fun, coef_zf_fun);
    [VXf_new_inv, VYf_new_inv, VZf_new_inv] = BsplineCompose3D( nx, ny, nz, Vxf_old_inv, Vyf_old_inv, Vzf_old_inv, coef_xf_inv_fun, coef_yf_inv_fun, coef_zf_inv_fun);
    
    [VXb_new, VYb_new, VZb_new] = BsplineCompose3D(nx, ny, nz, Vxb_old, Vyb_old, Vzb_old, coef_xb_fun, coef_yb_fun, coef_zb_fun);
    [VXb_new_inv, VYb_new_inv, VZb_new_inv] = BsplineCompose3D( nx, ny, nz, Vxb_old_inv, Vyb_old_inv, Vzb_old_inv, coef_xb_inv_fun, coef_yb_inv_fun, coef_zb_inv_fun);
    
    temp_coef_xf = imfilter(VXf_new, Bsplinekernel); %Vx  %doubt
    temp_coef_yf = imfilter(VYf_new, Bsplinekernel); %Vy
    temp_coef_zf = imfilter(VZf_new, Bsplinekernel);
    
    temp_coef_xb = imfilter(VXb_new, Bsplinekernel); %Vx  %doubt
    temp_coef_yb = imfilter(VYb_new, Bsplinekernel); %Vy
    temp_coef_zb = imfilter(VZb_new, Bsplinekernel);
    
    coef_xf(4:ny,4:nx, 4:nz) = temp_coef_xf(2:end-2,2:end-2,2:end-2);
    coef_yf(4:ny,4:nx, 4:nz) = temp_coef_yf(2:end-2,2:end-2,2:end-2);
    coef_zf(4:ny,4:nx, 4:nz) = temp_coef_zf(2:end-2,2:end-2,2:end-2);
    
    coef_xb(4:ny,4:nx, 4:nz) = temp_coef_xb(2:end-2,2:end-2,2:end-2);
    coef_yb(4:ny,4:nx, 4:nz) = temp_coef_yb(2:end-2,2:end-2,2:end-2);
    coef_zb(4:ny,4:nx, 4:nz) = temp_coef_zb(2:end-2,2:end-2,2:end-2);
    
    temp_coef_xf_inv = imfilter(VXf_new_inv, Bsplinekernel); %Vx  %doubt
    temp_coef_yf_inv = imfilter(VYf_new_inv, Bsplinekernel); %Vy
    temp_coef_zf_inv = imfilter(VZf_new_inv, Bsplinekernel);
    
    temp_coef_xb_inv = imfilter(VXb_new_inv, Bsplinekernel); %Vx  %doubt
    temp_coef_yb_inv = imfilter(VYb_new_inv, Bsplinekernel); %Vy
    temp_coef_zb_inv = imfilter(VZb_new_inv, Bsplinekernel);
    
    coef_xf_inv(4:ny,4:nx, 4:nz) = temp_coef_xf_inv(2:end-2,2:end-2,2:end-2);
    coef_yf_inv(4:ny,4:nx, 4:nz) = temp_coef_yf_inv(2:end-2,2:end-2,2:end-2);
    coef_zf_inv(4:ny,4:nx, 4:nz) = temp_coef_zf_inv(2:end-2,2:end-2,2:end-2);
    
    coef_xb_inv(4:ny,4:nx, 4:nz) = temp_coef_xb_inv(2:end-2,2:end-2,2:end-2);
    coef_yb_inv(4:ny,4:nx, 4:nz) = temp_coef_yb_inv(2:end-2,2:end-2,2:end-2);
    coef_zb_inv(4:ny,4:nx, 4:nz) = temp_coef_zb_inv(2:end-2,2:end-2,2:end-2);
    
    Vxf = VXf_new - Vx_ident;
    Vyf = VYf_new - Vy_ident;
    Vzf = VZf_new - Vz_ident;
    
    Vxb = VXb_new - Vx_ident;
    Vyb = VYb_new - Vy_ident;
    Vzb = VZb_new - Vz_ident;
    
    Vxf_inv = VXf_new_inv - Vx_ident;
    Vyf_inv = VYf_new_inv - Vy_ident;
    Vzf_inv = VZf_new_inv - Vz_ident;
    
    Vxb_inv = VXb_new_inv - Vx_ident;
    Vyb_inv = VYb_new_inv - Vy_ident;
    Vzb_inv = VZb_new_inv - Vz_ident;
    
    Vyf=-Vyf;
    Vyf_inv = -Vyf_inv;
    
    Vyb=-Vyb;
    Vyb_inv = -Vyb_inv;
    
    disp('update transformations')
    toc
    tic
    CIMM = coef_Morig(:)';
    CIPP = coef_Porig(:)';
    VXFF_new = VXf_new(:)';
    VYFF_new = VYf_new(:)';
    VZFF_new = VZf_new(:)';
    
    %M = BsplineComposeImage3D(nx,ny,nz,VXFF_new,VYFF_new,VZFF_new,CIMM);
    M = interp3(Vx_ident,Vy_ident,Vz_ident,M_orig,VXf_new,VYf_new, VZf_new);
    %M = interp3(Vx_ident,Vy_ident,Vz_ident,M,pyyf,pxxf, pzzf,'*linear');
    M(isnan(M)) = 0;
    %M = round(M);
    vecM = M(:)';
    coef_M= img2coef3D(vecM,nx,ny,nz);
    CIM = coef_M(:)';
    coef_matM = reshape(CIM, ny+img_pU, nx+img_pV, nz+img_pW);
    
    F = interp3(Vx_ident,Vy_ident,Vz_ident,F_orig,VXb_new,VYb_new, VZb_new);
    F(isnan(F)) = 0;
    vecF = F(:)';
    coef_F= img2coef3D(vecF,nx,ny,nz);
    CIF = coef_F(:)';
    coef_matF = reshape(CIF, ny+img_pU, nx+img_pV, nz+img_pW);
    
    %phi = BsplineComposeImage3D(nx,ny,nz,VXFF_new,VYFF_new,VZFF_new,CIPP);
    phi = interp3(Vx_ident,Vy_ident,Vz_ident,phi_orig,VXf_new,VYf_new, VZf_new);
    phi(isnan(phi)) = 0;
    vecP = phi(:)';
    coef_P = img2coef3D(vecP,nx,ny,nz);
    CIP = coef_P(:)';
    coef_matP = reshape(CIP, ny+img_pU, nx+img_pV, nz+img_pW);
    
    disp('new images')
    toc
    tic
    
    %Iplot = M;%interp3(X_vtk,Y_vtk,Z_vtk,M,pyyf,pxxf, pzzf);
    Iplot = resize3Dmatrix(size(M00,1), size(M00,2), size(M00,3),M/scale);
    
    if(level==1)
        Iplot = M;
    end
    
    save('Mbunny.mat','M');
    save('Iplotbunny.mat','Iplot');
    save('phibunny.mat','phi');
    Mlabel = ConvertIntensity3D_binary(Iplot);
    Flabel = ConvertIntensity3D_binary(F00);
    
    sum11 = (Iplot-F00).^2;
    sum1 = sum(sum11(:));
    msd = sum1/(200*200*200)
    [ Dice ] = DiceSimilarity3D(size(M00,1), size(M00,2), size(M00,3), Mlabel,Flabel ,2)
    PlotImage(iterct,M00,F00,Iplot);
    if(mod(iterct,1)==0)
        filename_vtk = sprintf('evolve%d.vtk',iterct);
        vtkwrite(filename_vtk, 'structured_grid',X,Y,Z,'scalars','Intensity',Iplot);
    end
    
    [MDX,MDY,MDZ] = gradient(M);
    [FDX,FDY,FDZ] = gradient(F);
    [DphiX,DphiY,DphiZ] = gradient(phi);
    
    initial_error = norm3d(M_orig-F_orig);
    final_error = norm3d(F_orig-M);
    disp('final error')
    disp(final_error)
    
    relative_residual = final_error/initial_error;
    disp('relative_residual')
    disp(relative_residual)
    
    rs_final = relative_residual;
    iteration_vec = [iteration_vec;iteration];
    rs_vec = [rs_vec;relative_residual];
    
    if(rs_final>rs_initial)
        break;
    end
    
    if(iteration==maxiteration)
        break;
    end
    
    disp('display images')
    toc
end