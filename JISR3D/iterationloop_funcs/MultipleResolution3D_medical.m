function [M,F,phi,VF,VB] = MultipleResolution3D_medical(M,F,phi,param)

[X,Y,Z] = meshgrid(0.5:size(F,1)-0.5,0.5:size(F,2)-0.5,0.5:size(F,3)-0.5);

initial_error = norm3d(F-M);
Idiff = F-M;
Idiff2 = (Idiff).^2;
initial_error = sqrt(sum(sum(sum(Idiff2,3),2),1));
disp('initial error')
disp(initial_error)

img_pU = 3;
img_pV = 3;
img_pW = 3;

toc
tic
%%
[MDX,MDY,MDZ] = gradient(M);
[FDX,FDY,FDZ] = gradient(F);
[DphiX,DphiY,DphiZ] = gradient(phi);

M_orig = M;
F_orig = F;
phi_orig = phi;
phi00 = phi;
F00 = F;
M00 = M;

%Basis value computed by substituting corresponding midpoint
%Bspline kernel
[Bsplinekernel] = BsplineKernel3D;

% Display intial config
orderGauss = 4;
iterct = 0;
xlen = 4;
smallNumber =  1e-12;
maxiteration = 20;

nx = size(F,1);
ny = size(F,2);
nz = size(F,3);

%% Multilevel framework
VXLF=X(:)';
VYLF=Y(:)';
VZLF=Z(:)';

VXLB=X(:)';
VYLB=Y(:)';
VZLB=Z(:)';

coef_xf = img2coef3D(VXLF,nx,ny,nz);   %coefficients for position x
coef_yf = img2coef3D(VYLF,nx,ny,nz);   %coefficients for position y
coef_zf = img2coef3D(VZLF,nx,ny,nz);   %coefficients for position z

coef_xf_fun = reshape(coef_xf, 1, (nx+img_pU)*(ny+img_pV)*(nz+img_pW));
coef_yf_fun = reshape(coef_yf, 1, (nx+img_pU)*(ny+img_pV)*(nz+img_pW));
coef_zf_fun = reshape(coef_zf, 1, (nx+img_pU)*(ny+img_pV)*(nz+img_pW));

coef_xb = img2coef3D(VXLB,nx,ny,nz);   %coefficients for position x
coef_yb = img2coef3D(VYLB,nx,ny,nz);   %coefficients for position y
coef_zb = img2coef3D(VZLB,nx,ny,nz);   %coefficients for position z

coef_xb_fun = reshape(coef_xb, 1,  (nx+img_pU)*(ny+img_pV)*(nz+img_pW));
coef_yb_fun = reshape(coef_yb, 1,  (nx+img_pU)*(ny+img_pV)*(nz+img_pW));
coef_zb_fun = reshape(coef_zb, 1,  (nx+img_pU)*(ny+img_pV)*(nz+img_pW));

% Start multilevel
VXf = X;
VYf = Y;
VZf = Z;

VXb = X;
VYb = Y;
VZb = Z;

ori_size_x=size(F,2);
ori_size_y=size(F,1);
ori_size_z=size(F,3);

disp('coefficients of initial vectors');
toc

for level = param.maxlevel:-1:1
    tic
    disp(['Registration level:' num2str(param.maxlevel-level+1)]);
    
    %% downsample image
    scale = 2^-(level-1);  % image scale
    
    new_size_x=round(ori_size_x*scale);
    new_size_y=round(ori_size_y*scale);
    new_size_z=round(ori_size_z*scale);
    npx = new_size_x; %size of scaled image
    npy = new_size_y;
    npz = new_size_z;
    
    F = resize3Dmatrix(new_size_x,new_size_y,new_size_z,F);
    vecF = F(:)';
    coef_F= img2coef3D(vecF,npx,npy,npz);
    CIF = coef_F(:)';
    coef_matF = reshape(CIF, npy+img_pU, npx+img_pV, npz+img_pW); % this is an unnecessary step, check later!!!!
    
    M =  resize3Dmatrix(new_size_x,new_size_y,new_size_z,M);
    vecM = M(:)';
    coef_M = img2coef3D(vecM,npx,npy,npz);
    CIM = coef_M(:)';
    coef_matM = reshape(CIM,npy+img_pU, npx+img_pV, npz+img_pW);
    
    phi = resize3Dmatrix(new_size_x,new_size_y,new_size_z,phi); % does not matter if phi00 or phi, because this only counts
    % for the first level which is anyway phi00.
    vecP = phi(:)';
    coef_P = img2coef3D(vecP,npx,npy,npz);
    CIP = coef_P(:)';
    coef_matP = reshape(CIP,npy+img_pU, npx+img_pV, npz+img_pW);
    
    M_orig = resize3Dmatrix(new_size_x,new_size_y,new_size_z,M00);
    vecMorig = M_orig(:)';
    coef_Morig = img2coef3D(vecMorig,npx,npy,npz);
    CIMorig = coef_Morig(:)';
    coef_matMorig = reshape(CIMorig,npy+img_pU, npx+img_pV, npz+img_pW);
    
    F_orig = resize3Dmatrix(new_size_x,new_size_y,new_size_z,F00);
    vecForig = F_orig(:)';
    coef_Forig= img2coef3D(vecForig,npx,npy,npz);
    CIForig = coef_Forig(:)';
    coef_matForig = reshape(CIForig,npy+img_pU, npx+img_pV, npz+img_pW);
    
    phi_orig = resize3Dmatrix(new_size_x,new_size_y,new_size_z,phi00);
    
    [MDX,MDY,MDZ] = gradient(M);
    [FDX,FDY,FDZ] = gradient(F);
    [DphiX,DphiY,DphiZ] = gradient(phi);
    
    disp('coefficient of scaled images');
    toc
    tic
    
    if(level==1)
        maxiteration = 100;
    end
    
    if(level==param.maxlevel)
        
        nx = size(M,2);
        ny = size(M,1);
        nz = size(M,3);
        
        [Vx_ident,Vy_ident,Vz_ident] = meshgrid((0.5:nx-0.5),(0.5:ny-0.5),(0.5:nz-0.5));
        [X_vtk,Y_vtk,Z_vtk] = meshgrid(1:nx,1:ny,1:nz);
        
        Vxlf = Vx_ident;
        Vylf = Vy_ident;
        Vzlf = Vz_ident;
        
        Vxlb = Vx_ident;
        Vylb = Vy_ident;
        Vzlb = Vz_ident;
        
        Vxlf_inv = Vx_ident;
        Vylf_inv = Vy_ident;
        Vzlf_inv = Vz_ident;
        
        Vxlb_inv = Vx_ident;
        Vylb_inv = Vy_ident;
        Vzlb_inv = Vz_ident;
        
        volume = (nx-1)*(ny-1)*(nz-1);
        
    else
        
        nx = size(M,2);
        ny = size(M,1);
        nz = size(M,3);
        
        [Vx_ident, Vy_ident,Vz_ident] = meshgrid((0.5:nx-0.5),(0.5:ny-0.5),(0.5:nz-0.5));
        [X_vtk,Y_vtk,Z_vtk] = meshgrid(1:nx,1:ny,1:nz);
        
        Vxlf = resize3Dmatrix(new_size_x,new_size_y,new_size_z,Vxf*scale);
        Vylf = resize3Dmatrix(new_size_x,new_size_y,new_size_z,Vyf*scale);
        Vzlf = resize3Dmatrix(new_size_x,new_size_y,new_size_z,Vzf*scale);
        
        Vxlb = resize3Dmatrix(new_size_x,new_size_y,new_size_z,Vxb*scale);
        Vylb = resize3Dmatrix(new_size_x,new_size_y,new_size_z,Vyb*scale);
        Vzlb = resize3Dmatrix(new_size_x,new_size_y,new_size_z,Vzb*scale);
        
        Vxlf_inv = resize3Dmatrix(new_size_x,new_size_y,new_size_z,Vxf_inv*scale);
        Vylf_inv = resize3Dmatrix(new_size_x,new_size_y,new_size_z,Vyf_inv*scale);
        Vzlf_inv = resize3Dmatrix(new_size_x,new_size_y,new_size_z,Vzf_inv*scale);
        
        Vxlb_inv = resize3Dmatrix(new_size_x,new_size_y,new_size_z,Vxb_inv*scale);
        Vylb_inv = resize3Dmatrix(new_size_x,new_size_y,new_size_z,Vyb_inv*scale);
        Vzlb_inv = resize3Dmatrix(new_size_x,new_size_y,new_size_z,Vzb_inv*scale);
        
        Vxlf = Vxlf+Vx_ident;
        Vylf = -Vylf+Vy_ident;
        Vzlf = Vzlf+Vz_ident;
        
        Vxlb = Vxlb+Vx_ident;
        Vylb = -Vylb+Vy_ident;
        Vzlb = Vzlb+Vz_ident;
        
        Vxlf_inv = Vxlf_inv+Vx_ident;
        Vylf_inv = -Vylf_inv+Vy_ident;
        Vzlf_inv = Vzlf_inv + Vz_ident;
        
        Vxlb_inv = Vxlb_inv+Vx_ident;
        Vylb_inv = -Vylb_inv+Vy_ident;
        Vzlb_inv = Vzlb_inv+Vz_ident;
        
        II0M_ori = M_orig(:)';
        II0F_ori = F_orig(:)';
        II0phi_ori = phi_orig(:)';
        
        %% coef is in the vector form, see cdmffd code
        coef_IM = img2coef3D(II0M_ori, nx,ny, nz);
        coef_IF = img2coef3D(II0F_ori, nx,ny, nz);
        coef_Iphi = img2coef3D(II0phi_ori,nx,ny,nz);
        
        CI0M = coef_IM(:)';
        CI0F = coef_IF(:)';
        CI0phi = coef_Iphi(:)';
        
        VXLf = Vxlf(:)';
        VYLf = Vylf(:)';
        VZLf = Vzlf(:)';
        
        VXLb = Vxlb(:)';
        VYLb = Vylb(:)';
        VZLb = Vzlb(:)';
        
        M = BsplineComposeImage3D(nx,ny,nz,VXLf, VYLf,VZLf, CI0M);
        vecM = M(:)';
        coef_M= img2coef3D(vecM,nx,ny,nz);
        CIM = coef_M(:)';
        coef_matM = reshape(CIM, npy+img_pU, npx+img_pV, npz+img_pW);
        
        F = BsplineComposeImage3D(nx,ny,nz,VXLb,VYLb,VZLb, CI0F);
        vecF = F(:)';
        coef_F = img2coef3D(vecF,nx,ny,nz);
        CIF = coef_F(:)';
        coef_matF = reshape(CIF,npy+img_pU, npx+img_pV, npz+img_pW);
        
        phi = BsplineComposeImage3D(nx,ny,nz,VXLf, VYLf, VZLf, CI0phi);
        vecP = phi(:)';
        coef_P = img2coef3D(vecP, nx, ny,nz);
        CIP = coef_P(:)';
        coef_matP = reshape(CIP,npy+img_pU, npx+img_pV, npz+img_pW);
        
        [MDX,MDY,MDZ] = gradient(M);  % might get obsolete
        [FDX,FDY,FDZ] = gradient(F);
        [DphiX,DphiY,DphiZ] = gradient(phi);
        
        volume = (nx-1)*(ny-1)*(nz-1);
    end
    
    %removed the grid plotting functions
    VXLf = Vxlf(:)';
    VYLf = Vylf(:)';
    VZLf = Vzlf(:)';
    
    VXLb = Vxlb(:)';
    VYLb = Vylb(:)';
    VZLb = Vzlb(:)';
    
    coef_xf = img2coef3D(VXLf,nx,ny,nz);
    coef_yf = img2coef3D(VYLf,nx,ny,nz);
    coef_zf = img2coef3D(VZLf,nx,ny,nz);
    
    coef_xb = img2coef3D(VXLb,nx,ny,nz);
    coef_yb = img2coef3D(VYLb,nx,ny,nz);
    coef_zb = img2coef3D(VZLb,nx,ny,nz);
    
    VXLf_inv = Vxlf_inv(:)';
    VYLf_inv = Vylf_inv(:)';
    VZLf_inv = Vzlf_inv(:)';
    
    VXLb_inv = Vxlb_inv(:)';
    VYLb_inv = Vylb_inv(:)';
    VZLb_inv = Vzlb_inv(:)';
    
    coef_xf_inv = img2coef3D(VXLf_inv,nx,ny,nz);
    coef_yf_inv = img2coef3D(VYLf_inv,nx,ny,nz);
    coef_zf_inv = img2coef3D(VZLf_inv,nx,ny,nz);
    
    coef_xb_inv = img2coef3D(VXLb_inv,nx,ny,nz);
    coef_yb_inv = img2coef3D(VYLb_inv,nx,ny,nz);
    coef_zb_inv = img2coef3D(VZLb_inv,nx,ny,nz);
    
    disp('set initial vectors for the level');
    toc
    tic
    %% Construct B-spline grid
    maxlev = param.maxlevel-level+1;
    
    setBsplineGrid3D
    
    ActiveNodes = [];
    Node = [];
    
    for levels = 1:maxlev,
        [nx1,ny1,nz1] = meshgrid(uknotvectorV{levels,1},uknotvectorU{levels,1},uknotvectorW{levels,1});
        Node = [Node;ny1(:),nx1(:),nz1(:)];
    end
    
    disp('set bspline grid');
    toc
    tic
    %Loop over each refinement level
    disp('Loop over each refinement level...');
    for multil = 0:1:maxlev-1
        
        iteration = 0;
        
        fprintf('Refinement at level %i...\n',multil+1);
        if(multil>0)
            [Em,Dm,Pm,ActiveNodes] = THB_Refinement(Em,Dm,Pm,knotvectorU, knotvectorV,knotvectorW,bf,CellGrad,meanGrad,param,multil,ActiveNodes);
        end
        
        disp('Collecting active elements, control points and basis functions...');
        [ac, bf, ACP, RHS,Em,Dm,ActiveNodes] = storeActiveElem(Em,Dm,Pm,multil,ActiveNodes);
        
        ac_ct = size(ac,1);
        bf_ct = size(bf,1);
        
        ActiveNodes = unique(ActiveNodes);
        ACPf = ACP;
        ACPb = ACP;
        ACPold = ACP;
        
        cell_co = zeros(ac_ct,3);
        for i = 1:ac_ct
            cell_id = ac(i,1);
            cell_le = ac(i,2);
            cell_co(i,:) = Em(cell_le).cell_centre(cell_id,:);
        end
        
        Idiff = sqrt((MDX).^2 + (MDY).^2);
        CellGrad = interp3(Vx_ident,Vy_ident,Vz_ident,Idiff,cell_co(:,2),cell_co(:,1),cell_co(:,3));
        meanGrad = mean2(Idiff);
    end
    
    Pmf = Pm;
    Pmb = Pm;
    Pmold = Pm;
    
    [Jm, Coeff] = computeNonZeroSplines(ac, param, Em, Dm,multil);
    
    disp('Computing the basis functions at pixel coordinates...');
    numPixels = int64(prod(sizeImage));
    pix =  [Vy_ident(:),Vx_ident(:),Vz_ident(:)];
    [Pixel, Pix2] = storePixelPhi(numPixels, multil,pix, knotvectorU, knotvectorV, knotvectorW, Em, Coeff, param);
    for p_ind = 1:numPixels
        Pixel(p_ind).phi = Pix2{p_ind};
    end
    clear Pix2
    
    disp('Computing the basis functions at gaussian points...');
    %compute the gaussian points and weights of the given gauss order
    [Gu,Wu] = ggquad(param.orderGauss);
    [Gv,Wv] = ggquad(param.orderGauss);
    [Gw,Ww] = ggquad(param.orderGauss);
    
    [PHI,PHIU,PHIV,PHIW,BIGX,BIGY,BIGZ,H] = GaussPhi(ac,Em,knotvectorU,knotvectorV,knotvectorW,Coeff,param,maxlev);
    % interpolate the intensity values of the target image at the gauss
    % points stored in BIGX, BIGY, BIGZ
    BIGXF = BIGX;
    BIGYF = BIGY;
    BIGZF = BIGZ;
    
    %cI0fd = interp3(X_vtk,Y_vtk,Z_vtk,F_orig,BIGY,BIGX,BIGZ,'*linear',min(M(:)));
    [cI0f, tempu1, tempu2, tempu3] = BsplineComposeImage3D_single(BIGY, BIGX, BIGZ, coef_matForig, size(BIGX,1), size(BIGX,2), size(BIGX,3));
    
    [cI0b, tempu4, tempu5, tempu6] = BsplineComposeImage3D_single(BIGY, BIGX, BIGZ, coef_matMorig, size(BIGX,1), size(BIGX,2), size(BIGX,3));
    
    Node(:,4) = 0;
    for i=1:size(ActiveNodes,1)
        Node(ActiveNodes(i,1),4) = i;
    end
    
    PlotGrid
    
    disp('refinement done')
    toc
    tic
    %figure
    %scatter3(ACP(:,1),ACP(:,2),ACP(:,3),'filled');
    
    %THBcheck
    
    % start the iteration loop
    %% Update the iteration loop here
    timestep = param.timestep(level,1);
    
    RHS_init = RHS;
    
    figure
    tol  = 10^-6;
    rs_initial = 1;
    rs_final = 2;
    
    PHI1 = cell2struct(PHI,'mat',2);
    PHIU1 = cell2struct(PHIU,'mat',2);
    PHIV1 = cell2struct(PHIV,'mat',2);
    PHIW1 = cell2struct(PHIW,'mat',2);
    
    disp('set iteration parameters');
    toc;
    iteration_vec = [];
    rs_vec = [];
    
    VXf_new = Vx_ident;
    VYf_new = Vy_ident;
    VZf_new = Vz_ident;
    ACCf = ACPold;
    ACCb = ACPold;
    
    iterationloop_medical
    
    tic
    Vxf = resize3Dmatrix(size(M00,1), size(M00,2), size(M00,3),Vxf/scale);
    Vyf = resize3Dmatrix(size(M00,1), size(M00,2), size(M00,3),Vyf/scale);
    Vzf = resize3Dmatrix(size(M00,1), size(M00,2), size(M00,3),Vzf/scale);
    
    Vxb = resize3Dmatrix(size(M00,1), size(M00,2), size(M00,3),Vxb/scale);
    Vyb = resize3Dmatrix(size(M00,1), size(M00,2), size(M00,3),Vyb/scale);
    Vzb = resize3Dmatrix(size(M00,1), size(M00,2), size(M00,3),Vzb/scale);
    
    Vxf_inv = resize3Dmatrix(size(M00,1), size(M00,2), size(M00,3),Vxf_inv/scale);
    Vyf_inv = resize3Dmatrix(size(M00,1), size(M00,2), size(M00,3),Vyf_inv/scale);
    Vzf_inv = resize3Dmatrix(size(M00,1), size(M00,2), size(M00,3),Vzf_inv/scale);
    
    Vxb_inv = resize3Dmatrix(size(M00,1), size(M00,2), size(M00,3),Vxb_inv/scale);
    Vyb_inv = resize3Dmatrix(size(M00,1), size(M00,2), size(M00,3),Vyb_inv/scale);
    Vzb_inv = resize3Dmatrix(size(M00,1), size(M00,2), size(M00,3),Vzb_inv/scale);
    
    M = resize3Dmatrix(size(M00,1), size(M00,2), size(M00,3),M/scale);
    F = resize3Dmatrix(size(M00,1), size(M00,2), size(M00,3),F/scale);
    
    phi = resize3Dmatrix(size(M00,1), size(M00,2), size(M00,3),phi/scale);
    
    disp('scale up images')
    toc
    
end

end