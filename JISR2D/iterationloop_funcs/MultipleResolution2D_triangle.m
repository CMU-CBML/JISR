function [M,F,phi,VF,VB,In_gridX, In_gridY,FgridX, FgridY] = MultipleResolution2D_triangle(M,F,phi,parameters)

M_orig = M;
F_orig = F;
phi_orig = phi;

phi00 = phi;
F00 = F;
M00 = M;

nx = size(F,1);
ny = size(F,2);

[X,Y] = meshgrid(0.5:size(M,2)-0.5,0.5:size(M,1)-0.5);

initial_error = norm(F-M);
disp('initial error')
disp(initial_error);

[MDX,MDY] = gradient(M);
[FDX,FDY] = gradient(F);
[DphiX,DphiY] = gradient(phi);

iterct = 0;
smallNumber =  1e-12;

%% Basis value computed by substituting corresponding midpoint
Mb = [1/24, 11/24, 11/24, 1/24];

%% Bspline kernel
M_filt = Mb'*Mb;

%% Multilevel framework
VXLF=X(:)';
VYLF=Y(:)';

VXLB=X(:)';
VYLB=Y(:)';

coef_xf = img2coef2D(nx,ny, VXLF);   %coefficients for position x
coef_yf = img2coef2D(nx,ny, VYLF);   %coefficients for position y

coef_xf_fun = reshape(coef_xf, 1, (nx+parameters.pU)*(ny+parameters.pV));
coef_yf_fun = reshape(coef_yf, 1, (nx+parameters.pU)*(ny+parameters.pV));

coef_xb = img2coef2D(nx,ny, VXLB);   %coefficients for position x
coef_yb = img2coef2D(nx,ny, VYLB);   %coefficients for position y

coef_xb_fun = reshape(coef_xb, 1, (nx+parameters.pU)*(ny+parameters.pV));
coef_yb_fun = reshape(coef_yb, 1, (nx+parameters.pU)*(ny+parameters.pV));

%% Start multilevel
VXf = X;
VYf = Y;
VXb = X;
VYb = Y;

for level = parameters.maxlevel:-1:1
    
    disp(['Registration level:' num2str(parameters.maxlevel-level+1)]);
    
    %% downsample image
    scale = 2^-(level-1);  % image scale
    F = imresize(F,scale);
    npx = size(F,1); %size of scaled image
    npy = size(F,2);
    vecF = F(:)';
    coef_F= img2coef2D(npx,npy, vecF);
    CIF = coef_F(:)';
    coef_matF = reshape(CIF, npy+parameters.pU, npx+parameters.pV); % this is an unnecessary step, check later!!!!
    
    M = imresize(M,scale);
    vecM = M(:)';
    coef_M= img2coef2D(npx,npy, vecM);
    CIM = coef_M(:)';
    coef_matM = reshape(CIM, npy+parameters.pU, npx+parameters.pV);
    
    phi = imresize(phi,scale); % does not matter if phi00 or phi, because this only counts
    % for the first level which is anyway phi00.
    vecP = phi(:)';
    coef_P= img2coef2D(npx,npy, vecP);
    CIP = coef_P(:)';
    coef_matP = reshape(CIP, npy+parameters.pU, npx+parameters.pV);
    
    M_orig = imresize(M00,scale);
    vecMorig = M_orig(:)';
    coef_Morig= img2coef2D(npx,npy, vecMorig);
    CIMorig = coef_Morig(:)';
    coef_matMorig = reshape(CIMorig, npy+parameters.pU, npx+parameters.pV);
    
    F_orig = imresize(F00,scale);
    vecForig = F_orig(:)';
    coef_Forig= img2coef2D(npx,npy, vecForig);
    CIForig = coef_Forig(:)';
    coef_matForig = reshape(CIForig, npy+parameters.pU, npx+parameters.pV);
    
    phi_orig = imresize(phi00,scale);
    vecphiorig = phi_orig(:)';
    coef_phiorig= img2coef2D(npx,npy, vecphiorig);
    CIPorig = coef_phiorig(:)';
    coef_matPorig = reshape(CIPorig, npy+parameters.pU, npx+parameters.pV);
    
    [MDX,MDY] = gradient(M);
    [FDX,FDY] = gradient(F);
    [DphiX,DphiY] = gradient(phi);
    
    if(level==parameters.maxlevel)
        
        nx = size(M,2);
        ny = size(M,1);
        
        [Vx_ident, Vy_ident] = meshgrid((0.5:nx-0.5),(0.5:ny-0.5));
        
        Vxlf = Vx_ident;
        Vylf = Vy_ident;
        
        Vxlb = Vx_ident;
        Vylb = Vy_ident;
        
        Vxlf_inv = Vx_ident;
        Vylf_inv = Vy_ident;
        
        Vxlb_inv = Vx_ident;
        Vylb_inv = Vy_ident;
        
        area = (nx-1)*(ny-1);
        
    else
        
        nx = size(M,2);
        ny = size(M,1);
        
        [Vx_ident, Vy_ident] = meshgrid((0.5:nx-0.5),(0.5:ny-0.5));
        
        Vxlf = imresize(Vxf*scale,scale);
        Vylf = imresize(Vyf*scale,scale);
        
        Vxlb = imresize(Vxb*scale,scale);
        Vylb = imresize(Vyb*scale,scale);
        
        Vxlf_inv = imresize(Vxf_inv*scale,scale);
        Vylf_inv = imresize(Vyf_inv*scale,scale);
        
        Vxlb_inv = imresize(Vxb_inv*scale,scale);
        Vylb_inv = imresize(Vyb_inv*scale,scale);
        
        Vxlf = Vxlf+Vx_ident;
        Vylf = -Vylf+Vy_ident;
        
        Vxlb = Vxlb+Vx_ident;
        Vylb = -Vylb+Vy_ident;
        
        Vxlf_inv = Vxlf_inv+Vx_ident;
        Vylf_inv = -Vylf_inv+Vy_ident;
        
        Vxlb_inv = Vxlb_inv+Vx_ident;
        Vylb_inv = -Vylb_inv+Vy_ident;
        
        II0M_ori = M_orig(:)';
        II0F_ori = F_orig(:)';
        II0phi_ori = phi_orig(:)';
        
        coef_IM = img2coef2D(nx,ny, II0M_ori);
        coef_IF = img2coef2D(nx,ny, II0F_ori);
        coef_Iphi = img2coef2D(nx,ny, II0phi_ori);
        
        CI0M = coef_IM(:)';
        CI0F = coef_IF(:)';
        CI0phi = coef_Iphi(:)';
        
        VXLf=Vxlf(:)';
        VYLf=Vylf(:)';
        
        VXLb=Vxlb(:)';
        VYLb=Vylb(:)';
        
        M = BsplineComposeImage2D(VXLf, VYLf, CI0M, nx, ny);
        vecM = M(:)';
        coef_M= img2coef2D(nx,ny, vecM);
        CIM = coef_M(:)';
        coef_matM = reshape(CIM, ny+parameters.pU, nx+parameters.pV);
        
        F = BsplineComposeImage2D(VXLb, VYLb, CI0F, nx, ny);
        vecF = F(:)';
        coef_F = img2coef2D(nx,ny, vecF);
        CIF = coef_F(:)';
        coef_matF = reshape(CIF, ny+parameters.pU, nx+parameters.pV);
        
        phi = BsplineComposeImage2D(VXLf, VYLf, CI0phi, nx, ny);
        vecP = phi(:)';
        coef_P = img2coef2D(nx,ny, vecP);
        CIP = coef_P(:)';
        coef_matP = reshape(CIP, ny+parameters.pU, nx+parameters.pV);
        
        [MDX,MDY] = gradient(M);  % might get obsolete
        [FDX,FDY] = gradient(F);
        [DphiX,DphiY] = gradient(phi);
        
        area = (nx-1)*(ny-1);
        
    end
    
    ngridX = 40;
    ngridY = 40;
    
    VGx = zeros(ngridY,ngridX);
    VGy = zeros(ngridY,ngridX);
    
    [FgridX, FgridY, plot_pix_index_column, plot_pix_index_row] = makeGrid(ngridX,ngridY,nx,ny);
    
    gridX0 = FgridX;
    gridY0 = FgridY;
    
    VXLf = Vxlf(:)';
    VYLf = Vylf(:)';
    
    VXLb = Vxlb(:)';
    VYLb = Vylb(:)';
    
    coef_xf = img2coef2D(nx,ny,VXLf);
    coef_yf = img2coef2D(nx,ny,VYLf);
    
    coef_xb = img2coef2D(nx,ny,VXLb);
    coef_yb = img2coef2D(nx,ny,VYLb);
    
    VXLf_inv = Vxlf_inv(:)';
    VYLf_inv = Vylf_inv(:)';
    
    VXLb_inv = Vxlb_inv(:)';
    VYLb_inv = Vylb_inv(:)';
    
    coef_xf_inv = img2coef2D(nx,ny,VXLf_inv);
    coef_yf_inv = img2coef2D(nx,ny,VYLf_inv);
    
    coef_xb_inv = img2coef2D(nx,ny,VXLb_inv);
    coef_yb_inv = img2coef2D(nx,ny,VYLb_inv);
    
    %% Construct B-spline grid
    maxlev = parameters.maxlevel-level+1;
    Bvectf = cell(maxlev,1);
    Bvectb = cell(maxlev,1);
    
    [Dm,Pm,Em,Bvect,knotvectorU,knotvectorV,nobU,nobV,nelemU] = setBsplineGrid(maxlev,parameters,F);
    Pmold = Pm;
    for multilev = 0:1:maxlev-1
        if(multilev>0)
            for j =1:bf_ct
                bbc = bf(j,1:2);
                bf_lev = bf(j,3);
                EEM = Em{bf_lev,1};
                BEM = Dm{bf_lev,1};
                bind = nobU(bf_lev,1)*(bbc(1,2)-1)+bbc(1,1);
                supp_cells = BEM{bind,6};
                grad = 0;
                supp_ct = 0;
                for i =1:size(supp_cells,2)
                    if(supp_cells(1,i)~=0)
                        supp_ct = supp_ct + 1;
                        ac_ind = EEM{supp_cells(1,i),11};
                        grad  = grad + Cell_grad(ac_ind,1);
                    end
                end
                
                grad = grad/supp_ct;
                
                %Refinement to create next level
                rho = parameters.rho(multilev+1);
                if(grad>=(rho*meangrad))
                    [Dmold,Emold,Pmold] =  Refine2D(bbc(1,1),bbc(1,2),bf_lev,Dm,Em,Pmold,knotvectorU,knotvectorV,pU,pV);
                    [Dm,Em,Pm] =  Refine2D(bbc(1,1),bbc(1,2),bf_lev,Dm,Em,Pm,knotvectorU,knotvectorV,pU,pV);
                end
            end
        end
        
        Pmold = Pm;
        
        ac_ct = 0;
        bf_ct = 0;
        ac = zeros(1,2);
        bf = zeros(1,3);
        
        for lev = 1:(multilev+1)
            
            EE = Em{lev,1};
            BE = Dm{lev,1};
            sizee = size(EE,1);
            sizeb = size(BE,1);
            
            for i = 1:sizee
                if(EE{i,4}==1)
                    ac_ct = ac_ct+1;
                    ac(ac_ct,1) = EE{i,1};
                    ac(ac_ct,2) = lev;
                    EE{i,11} = ac_ct;
                end
            end
            
            for j = 1:sizeb,
                if(BE{j,3}==1),
                    bf_ct = bf_ct + 1;
                    bf(bf_ct,1:2) = BE{j,1};
                    bf(bf_ct,3) = lev;
                    BE{j,10} = bf_ct;
                end
            end
            
            Em{lev,1} = EE;
            Dm{lev,1} = BE;
            
            
        end
        
        [Jm,Coeff,Pixel,HH,PHI,PHIU,PHIV,BIGX,BIGY] = constructAdaptiveGrid(ac,parameters,Dm,Em,M,knotvectorU,knotvectorV,multilev,nobU,nobV,nelemU);
        
        cell_co = zeros(ac_ct,2);
        for i = 1:ac_ct
            cell_id = ac(i,1);
            cell_le = ac(i,2);
            EEM = Em{cell_le,1};
            cell_co(i,1) = EEM{cell_id,8};
            cell_co(i,2) = EEM{cell_id,9};
        end
        
        M_255 = M.*255;
        F_255 = F_orig.*255;
        [MDX,MDY] = gradient(M_255-F_255);
        Idiff = sqrt((MDX).^2 + (MDY).^2);
        Cell_grad = interp2(Vx_ident,Vy_ident,Idiff,cell_co(:,2),cell_co(:,1));
        meangrad = mean2(Idiff);
        Pmf = Pm;
        Pmb = Pm;

        %displayAdaptiveGrid(ac,Coeff,Em,knotvectorU,knotvectorV,Jm,Pmold,parameters,nx,ny);
    end
    
    timestep = parameters.timestep(parameters.maxlevel-level+1);
    
    [cI0f,tempu1,tempu2] = BsplineComposeImage2D_single(BIGY, BIGX, coef_matForig, size(BIGX,1), size(BIGX,2));
    [cI0b,tempu3,tempu4] = BsplineComposeImage2D_single(BIGY, BIGX, coef_matMorig, size(BIGX,1), size(BIGX,2));
    
    figure
    tol = 10^-4;
    rs_initial = 1;
    rs_final = 2;
    
    iterationloop_triangle
    
    figure
    imagesc(phi);
    phi_img = 255.*ones(nx,ny);
    phi_img(phi<=mean(phi(:))) = -255;
    hold on
    contour(phi_img,[0,0],'g','Linewidth',3.24);
    hold off
    set(gca,'position',[0 0 1 1],'units','normalized')
    colormap gray
    
    figure
    imagesc(M);
    colormap gray
    set(gca,'position',[0 0 1 1],'units','normalized')
    
    figure
    plotGrid(In_gridX, In_gridY)
    axis ([1 nx 1 ny])
    colormap('gray')
    set(gca,'position',[0 0 1 1],'units','normalized')
    set(gca,'YDir','reverse');
    
    figure
    plotGrid(FgridX, FgridY)
    axis ([1 nx 1 ny])
    colormap('gray')
    set(gca,'position',[0 0 1 1],'units','normalized')
    set(gca,'YDir','reverse');
    
    figure
    imagesc(M-F_orig);
    set(gca,'position',[0 0 1 1],'units','normalized')
    colormap gray
    
    Vxf = imresize(Vxf/scale,size(M00));
    Vyf = imresize(Vyf/scale,size(M00));
    
    Vxb = imresize(Vxb/scale,size(M00));
    Vyb = imresize(Vyb/scale,size(M00));
    
    Vxf_inv = imresize(Vxf_inv/scale,size(M00));
    Vyf_inv = imresize(Vyf_inv/scale,size(M00));
    
    Vxb_inv = imresize(Vxb_inv/scale,size(M00));
    Vyb_inv = imresize(Vyb_inv/scale,size(M00));
    
    M = imresize(M/scale,size(M00));
    F = imresize(F/scale,size(F00));
    phi = imresize(phi/scale,size(phi00));
    
    
end
VF = struct('VXf_new',VXf_new,'VXf_new_inv',VXf_new_inv,'VYf_new',VYf_new,'VYf_new_inv',VYf_new_inv);
VB = struct('VXb_new',VXb_new,'VXb_new_inv',VXb_new_inv,'VYb_new',VYb_new,'VYb_new_inv',VYb_new_inv);
end
