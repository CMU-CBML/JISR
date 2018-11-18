par1 = parameters.par1;
par2 = parameters.par2;

Btermf = zeros(ac_ct*orderGauss,orderGauss,2);
Btermb = zeros(ac_ct*orderGauss,orderGauss,2);
Bsegmat = zeros(ac_ct*orderGauss,orderGauss,3);

iteration=0;
while(abs(rs_final-rs_initial)>tol)
    
    Pmf = Pmold;
    Pmb = Pmold;
    Cmf = Pmold;
    Cmb = Pmold;
    
    rs_initial = rs_final;
    H_phi = regular_Heiviside_fun(phi);
    
    c1 = (sum(sum(F_orig.*H_phi)) * area)/(sum(sum(H_phi)) * area);
    c2 = (sum(sum(F_orig.*(1-H_phi))) * area)/(sum(sum(1-H_phi)) * area );
    
    [D1E] = Delta_fun(phi);
    vecD = D1E(:)';
    coef_D = img2coef2D(nx,ny, vecD);
    CID = coef_D(:)';
    coef_matD = reshape(CID, ny+3, nx+3);
    
    iterct = iterct +1;
    iteration=iteration+1;
    iterct
    
    Bvectf = zeros(bf_ct,2);
    Bvectb = zeros(bf_ct,2);
    
    [cI1f,cDI1Xf, cDI1Yf] = BsplineComposeImage2D_single(BIGY, BIGX, coef_matM, size(BIGX,1), size(BIGX,2));
    [cI1b,cDI1Xb, cDI1Yb] = BsplineComposeImage2D_single(BIGY, BIGX, coef_matF, size(BIGX,1), size(BIGX,2));
    [ctempp,cDphiX, cDphiY] = BsplineComposeImage2D_single(BIGY, BIGX, coef_matP, size(BIGX,1), size(BIGX,2));
    [dDelta, dt, dt1] = BsplineComposeImage2D_single(BIGY, BIGX, coef_matD, size(BIGX,1), size(BIGX,2));
    
    denominatef = sqrt(cDI1Xf.^2+cDI1Yf.^2 + smallNumber);
    denominateb = sqrt(cDI1Xb.^2+cDI1Yb.^2 + smallNumber);
    
    Bterm1f = (cI1f - cI0f).*2.*cDI1Xf./denominatef;
    Bterm2f = (cI1f - cI0f).*2.*cDI1Yf./denominatef;
    Btermf(:,:,1) = Bterm1f;
    Btermf(:,:,2) = Bterm2f;
    
    Bterm1b = (cI1b - cI0b).*2.*cDI1Xb./denominateb;
    Bterm2b = (cI1b - cI0b).*2.*cDI1Yb./denominateb;
    Btermb(:,:,1) = Bterm1b;
    Btermb(:,:,2) = Bterm2b;
    
    Bseg = ((-c1 + cI0f).^2 - (-c2 + cI0f).^2).*dDelta;
    Bsegx = cDphiX;
    Bsegy = cDphiY;
    Bsegmat(:,:,1) = Bseg;
    Bsegmat(:,:,2) = Bsegx;
    Bsegmat(:,:,3) = Bsegy;
    
    [Bvectf, Bvectb] = computeIntegrationFidelity(Btermf, Btermb, Bsegmat, Jm, PHI, Dm, HH,Bvectf,Bvectb,parameters.par1,parameters.par2);

    Bvectf = bcondition1(Bvectf,Dm,bf,nobU);
    Bvectb = bcondition1(Bvectb,Dm,bf,nobU);
    
    %Update the control points
    C1f = zeros(bf_ct,2);
    C1b = zeros(bf_ct,2);
    for i=1:bf_ct,
        bbc = bf(i,1:2);
        bf_lev = bf(i,3);
        bi = nobU(bf_lev,1)*(bbc(1,2)-1)+bbc(1,1);
        Cif = Cmf{bf_lev,1};
        Cib = Cmb{bf_lev,1};
        
        C1f(i,1) = Cif(bi,1);
        C1f(i,2) = Cif(bi,2);
        
        C1b(i,1) = Cib(bi,1);
        C1b(i,2) = Cib(bi,2);
    end
    timestep
    C1f = C1f - timestep*Bvectf; %update set of first control points
    C1b = C1b - timestep*Bvectb;
    
    for i = 1:bf_ct,
        bbc = bf(i,1:2);
        bf_lev = bf(i,3);
        bi = nobU(bf_lev,1)*(bbc(1,2)-1)+bbc(1,1);
        
        Cif = Cmf{bf_lev,1};
        Cif(bi,1) = C1f(i,1);
        Cif(bi,2) = C1f(i,2);
        Cmf{bf_lev,1} = Cif;
        
        Cib = Cmb{bf_lev,1};
        Cib(bi,1) = C1b(i,1);
        Cib(bi,2) = C1b(i,2);
        Cmb{bf_lev,1} = Cib;
    end
    
    %% now do the regularization
    RHSf = zeros(bf_ct,2);
    RHSb = zeros(bf_ct,2);
    
    [BIGXf,BIGYf,BIGXb,BIGYb] = computedeformation(Jm,PHI,PHIU,PHIV,Cmf,Cmb);
    
    PHII = cell(1,3);
    PHII{1,1} = PHI;
    PHII{1,2} = PHIU;
    PHII{1,3} = PHIV;
    
    BIGXXf = BIGXf(:,:,1);
    BIGMUXf = BIGXf(:,:,2);
    BIGMUYf = BIGXf(:,:,3);
    
    BIGYYf = BIGYf(:,:,1);
    BIGMVXf = BIGYf(:,:,2);
    BIGMVYf = BIGYf(:,:,3);
    
    BIGXXb = BIGXb(:,:,1);
    BIGMUXb = BIGXb(:,:,2);
    BIGMUYb = BIGXb(:,:,3);
    
    BIGYYb = BIGYb(:,:,1);
    BIGMVXb = BIGYb(:,:,2);
    BIGMVYb = BIGYb(:,:,3);
    
    BIGXXfb = BIGYYf;
    BIGYYfb = BIGXXf;
    BIGXXbf = BIGYYb;
    BIGYYbf = BIGXXb;
    
    [FBXX,BFXX, PHIF,PHIB,PHIUF,PHIVF,PHIUB,PHIVB]  = bidirectionupdate(ac, parameters, Jm, Coeff, Em, knotvectorU, knotvectorV, nobU, nobV, Cmf,Cmb,BIGXXfb,BIGYYfb,BIGXXbf,BIGYYbf);
    
    FBX = FBXX(:,:,1);
    FBY = FBXX(:,:,2);
    FBXDX = FBXX(:,:,3);
    FBXDY = FBXX(:,:,4);
    FBYDX = FBXX(:,:,5);
    FBYDY = FBXX(:,:,6);
    
    BFX = BFXX(:,:,1);
    BFY = BFXX(:,:,2);
    BFXDX = BFXX(:,:,3);
    BFXDY = BFXX(:,:,4);
    BFYDX = BFXX(:,:,5);
    BFYDY = BFXX(:,:,6);
    
    PHIFF = cell(1,3);
    PHIBB = cell(1,3);
    PHIFF{1,1} = PHIF;
    PHIFF{1,2} = PHIUF;
    PHIFF{1,3} = PHIVF;

    PHIBB{1,1} = PHIB;
    PHIBB{1,2} = PHIUB;
    PHIBB{1,3} = PHIVB;
    
    [RHSff,RHSbb,PHIF11] = computeIntegrationRegularization(parameters,HH,BIGXf,BIGYf, BIGXb,BIGYb, Jm,Dm, PHII,BIGX,BIGY, FBXX,BFXX, PHIFF,PHIBB,RHSf,RHSb);
                             
    RHSf = bcondition1(RHSf,Dm,bf,nobU);
    RHSb = bcondition1(RHSb,Dm,bf,nobU);
    
    %Update the control points
    P1f = zeros(bf_ct,2);
    P1b = zeros(bf_ct,2);
    
    for i=1:bf_ct,
        bbc = bf(i,1:2);
        bf_lev = bf(i,3);
        bi = nobU(bf_lev,1)*(bbc(1,2)-1)+bbc(1,1);
        Pif = Cmf{bf_lev,1};
        Pib = Cmb{bf_lev,1};
        
        P1f(i,1) = Pif(bi,1);
        P1f(i,2) = Pif(bi,2);
        
        P1b(i,1) = Pib(bi,1);
        P1b(i,2) = Pib(bi,2);
    end
    
    P1fd = P1f - timestep*RHSf;
    P1bd = P1b - timestep*RHSb;
    
    for i = 1:bf_ct,
        bbc = bf(i,1:2);
        bf_lev = bf(i,3);
        bi = nobU(bf_lev,1)*(bbc(1,2)-1)+bbc(1,1);
        
        Pif1 = Pmf{bf_lev,1};
        Pif1(bi,1) = P1fd(i,1);
        Pif1(bi,2) = P1fd(i,2);
        Pmf{bf_lev,1} = Pif1;
        
        Pib1 = Pmb{bf_lev,1};
        Pib1(bi,1) = P1bd(i,1);
        Pib1(bi,2) = P1bd(i,2);
        Pmb{bf_lev,1} = Pib1;
    end

    [PXF,PXB] = computenewpoints(M,Pixel,Jm,Pmf,Pmb);
    pxxf = PXF(:,:,1);
    pyyf = PXF(:,:,2);
    pxxb = PXB(:,:,1);
    pyyb = PXB(:,:,2);

    Vxf_temp = pxxf-Vx_ident;
    Vyf_temp = pyyf-Vy_ident;
    
    Vxb_temp = pxxb-Vx_ident;
    Vyb_temp = pyyb-Vy_ident;
    
    Vxf_old = Vxf_temp + Vx_ident;
    Vyf_old = Vyf_temp + Vy_ident;
    
    Vxb_old = Vxb_temp + Vx_ident;
    Vyb_old = Vyb_temp + Vy_ident;
    
    Vxf_old = reshape(Vxf_old,1,nx*ny);
    Vyf_old = reshape(Vyf_old,1,nx*ny);
    
    Vxb_old = reshape(Vxb_old,1,nx*ny);
    Vyb_old = reshape(Vyb_old,1,nx*ny);
    
    coef_xf_fun = reshape(coef_xf, 1, (nx+pU)*(ny+pV));
    coef_yf_fun = reshape(coef_yf, 1, (nx+pU)*(ny+pV));
    
    coef_xb_fun = reshape(coef_xb, 1, (nx+pU)*(ny+pV));
    coef_yb_fun = reshape(coef_yb, 1, (nx+pU)*(ny+pV));
    
    Vxf_old_inv = -Vxf_temp + Vx_ident;
    Vyf_old_inv = -Vyf_temp + Vy_ident;
    
    Vxb_old_inv = -Vxb_temp + Vx_ident;
    Vyb_old_inv = -Vyb_temp + Vy_ident;
    
    Vxf_old_inv = reshape(Vxf_old_inv,1,nx*ny);
    Vyf_old_inv = reshape(Vyf_old_inv,1,nx*ny);
    
    Vxb_old_inv = reshape(Vxb_old_inv,1,nx*ny);
    Vyb_old_inv = reshape(Vyb_old_inv,1,nx*ny);
    
    coef_xf_inv_fun = reshape(coef_xf_inv, 1, (nx+pU)*(ny+pV));
    coef_yf_inv_fun = reshape(coef_yf_inv, 1, (nx+pU)*(ny+pV));
    
    coef_xb_inv_fun = reshape(coef_xb_inv, 1, (nx+pU)*(ny+pV));
    coef_yb_inv_fun = reshape(coef_yb_inv, 1, (nx+pU)*(ny+pV));
    
    [VXf_new, VYf_new] = BsplineCompose2D( Vxf_old, Vyf_old, coef_xf_fun, coef_yf_fun, nx, ny);
    [VXf_new_inv, VYf_new_inv] = BsplineCompose2D( Vxf_old_inv, Vyf_old_inv, coef_xf_inv_fun, coef_yf_inv_fun, nx, ny);
    
    [VXb_new, VYb_new] = BsplineCompose2D( Vxb_old, Vyb_old, coef_xb_fun, coef_yb_fun, nx, ny);
    [VXb_new_inv, VYb_new_inv] = BsplineCompose2D( Vxb_old_inv, Vyb_old_inv, coef_xb_inv_fun, coef_yb_inv_fun, nx, ny);
    
    temp_coef_xf = imfilter(VXf_new, M_filt); %Vx  %doubt
    temp_coef_yf = imfilter(VYf_new, M_filt); %Vy
    
    temp_coef_xb = imfilter(VXb_new, M_filt); %Vx  %doubt
    temp_coef_yb = imfilter(VYb_new, M_filt); %Vy
    
    coef_xf(4:ny,4:nx) = temp_coef_xf(2:end-2,2:end-2);
    coef_yf(4:ny,4:nx) = temp_coef_yf(2:end-2,2:end-2);
    
    coef_xb(4:ny,4:nx) = temp_coef_xb(2:end-2,2:end-2);
    coef_yb(4:ny,4:nx) = temp_coef_yb(2:end-2,2:end-2);
    
    temp_coef_xf_inv = imfilter(VXf_new_inv, M_filt); %Vx  %doubt
    temp_coef_yf_inv = imfilter(VYf_new_inv, M_filt); %Vy
    
    temp_coef_xb_inv = imfilter(VXb_new_inv, M_filt); %Vx  %doubt
    temp_coef_yb_inv = imfilter(VYb_new_inv, M_filt); %Vy
    
    coef_xf_inv(4:ny,4:nx) = temp_coef_xf_inv(2:end-2,2:end-2);
    coef_yf_inv(4:ny,4:nx) = temp_coef_yf_inv(2:end-2,2:end-2);
    
    coef_xb_inv(4:ny,4:nx) = temp_coef_xb_inv(2:end-2,2:end-2);
    coef_yb_inv(4:ny,4:nx) = temp_coef_yb_inv(2:end-2,2:end-2);
    
    Vxf = VXf_new - Vx_ident;
    Vyf = VYf_new - Vy_ident;
    
    Vxb = VXb_new - Vx_ident;
    Vyb = VYb_new - Vy_ident;
    
    Vxf_inv = VXf_new_inv - Vx_ident;
    Vyf_inv = VYf_new_inv - Vy_ident;
    
    Vxb_inv = VXb_new_inv - Vx_ident;
    Vyb_inv = VYb_new_inv - Vy_ident;
    
    Vyf=-Vyf;
    Vyf_inv = -Vyf_inv;
    
    Vyb=-Vyb;
    Vyb_inv = -Vyb_inv;
    
    %M = BsplineComposeImage2D(VXf_new(:)',VYf_new(:)', CIMorig, nx, ny);
    M = interp2(Vx_ident,Vy_ident,M_orig,VXf_new,VYf_new);
    M(isnan(M)) = 0;
    vecM = M(:)';
    coef_M= img2coef2D(nx,ny, vecM);
    CIM = coef_M(:)';
    coef_matM = reshape(CIM, ny+3, nx+3);
    
    %F = BsplineComposeImage2D(VXb_new(:)',VYb_new(:)', CIForig, nx, ny);
    F = interp2(Vx_ident,Vy_ident,F_orig,VXb_new,VYb_new);
    F(isnan(F)) = 0;
    vecF = F(:)';
    coef_F= img2coef2D(nx,ny, vecF);
    CIF = coef_F(:)';
    coef_matF = reshape(CIF, ny+3, nx+3);
    
    %phi = BsplineComposeImage2D(VXf_new(:)',VYf_new(:)', CIPorig, nx, ny);
    phi = interp2(Vx_ident,Vy_ident,phi_orig,VXf_new,VYf_new);
    phi(isnan(phi)) = 0;
    vecP = phi(:)';
    coef_P = img2coef2D(nx,ny, vecP);
    CIP = coef_P(:)';
    coef_matP = reshape(CIP, ny+3, nx+3);
    
    [MDX,MDY] = gradient(M);
    [FDX,FDY] = gradient(F);
    [DphiX,DphiY] = gradient(phi);
    
    In_gridX=zeros(size(FgridX));
    In_gridY=zeros(size(FgridY));
    
    for ii=1:ngridX
        for jj=1:ngridY
            
            j = gridY0(jj,ii);
            i = gridX0(jj,ii);
            
            VGx(jj,ii) = Vxf_inv(j,i);
            VGy(jj,ii) = Vyf_inv(j,i);
            
            In_gridX(jj,ii) = round(i+VGx(jj,ii));
            In_gridY(jj,ii) = round(j-VGy(jj,ii));
            
        end
    end
    
    for ii=1:ngridX
        for jj=1:ngridY
            
            j = gridY0(jj,ii);
            i = gridX0(jj,ii);
            
            VGx(jj,ii) = Vxb_inv(j,i);
            VGy(jj,ii) = Vyb_inv(j,i);
            
            FgridX(jj,ii) = round(i+VGx(jj,ii));
            FgridY(jj,ii) = round(j-VGy(jj,ii));
            
        end
    end
    
    initial_error = norm(M_orig-F_orig);
    final_error = norm(F_orig-M);
    disp('final error')
    disp(final_error)
    
    relative_residual = final_error/initial_error;
    disp('relative_residual')
    disp(relative_residual)
    
    rs_final = relative_residual;
    if(rs_final>rs_initial)
        M = M_1;
        break;
    end
    M_1 = M;
    
    clf;
    colormap gray;
    subplot(3,3,1);
    imagesc(M);
    title('Source Image');
    axis([1 nx 1 ny])
    %set(gca,'position',[0 0 1 1],'units','normalized')
    
    colormap gray;
    subplot(3,3,2);
    imagesc(F_orig);
    title('Target Image')
    axis([1 nx 1 ny])
    %set(gca,'position',[0 0 1 1],'units','normalized')
    
    colormap gray;
    subplot(3,3,3);
    imagesc(M-F_orig);
    title('Image difference');
    axis([1 nx 1 ny])
    %set(gca,'position',[0 0 1 1],'units','normalized')
    
    colormap gray;
    subplot(3,3,4);
    imagesc(phi);
    phi_img = 255.*ones(nx,ny);
    phi_img(phi<=mean(phi(:))) = -255;
    hold on
    contour(phi_img,[0,0],'g','Linewidth',1.24);
    hold off
    title('Segmented Image');
    axis([1 nx 1 ny])
    %set(gca,'position',[0 0 1 1],'units','normalized')
    
    th2=subplot(3,3,5);
    cla(th2);
    plotGrid(In_gridX, In_gridY)
    axis image
    title(['Inverse deformed grid at step: ',num2str(k) ])
    colormap('gray')
    set(gca,'YDir','reverse');
    
    
    th1=subplot(3,3,6);
    cla(th1);
    plotGrid(FgridX, FgridY)
    axis image
    title(['Deformed grid at step: ',num2str(k) ])
    colormap('gray')
    set(gca,'YDir','reverse');
    
    quiver_matrix_x = Vxf(plot_pix_index_row, plot_pix_index_column);
    quiver_matrix_y = Vyf(plot_pix_index_row, plot_pix_index_column);
    
    th5=subplot(3,3,7);
    cla(th5);
    quiver(-quiver_matrix_x(end:-1:1,:),-quiver_matrix_y(end:-1:1,:),2)
    title(['Transformation field at step ', num2str(k)])
    axis tight
    drawnow
    
    Jdet = computeJacobian(parameters,ac,Coeff,Em,Jm,knotvectorU,knotvectorV,Pmf);
    
    min(Jdet)
    minDet(iterct,1) = min(Jdet);
    
    if(minDet(iterct,1)<=0.1)
        break
    end
    
    if(iteration==parameters.maxiteration)
        break
    end
end