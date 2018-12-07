%% Create configuration object of class 'coder.MexCodeConfig'.
cfg = coder.config('mex');
cfg.GenerateReport = true;
%cfg.ReportPotentialDifferences = false;
cfg.GlobalDataSyncMethod = 'NoSync';

%% Define argument types for entry-point 'img2coef3D'.
ARGS1 = cell(1,1);
ARGS1{1} = cell(1,4);
ARGS1{1}{1} = coder.typeof(0,[1 Inf],[0 1]);
ARGS1{1}{2} = coder.typeof(0);
ARGS1{1}{3} = coder.typeof(0);
ARGS1{1}{4} = coder.typeof(0);

%% Invoke MATLAB Coder.
codegen -config cfg img2coef3D -args ARGS1{1} -o img2coef3D

%% Define argument types for entry-point 'BsplineComposeImage3D'.
ARGS2 = cell(1,1);
ARGS2{1} = cell(7,1);
ARGS2{1}{1} = coder.typeof(0);
ARGS2{1}{2} = coder.typeof(0);
ARGS2{1}{3} = coder.typeof(0);
ARGS2{1}{4} = coder.typeof(0,[1 Inf],[0 1]);
ARGS2{1}{5} = coder.typeof(0,[1 Inf],[0 1]);
ARGS2{1}{6} = coder.typeof(0,[1 Inf],[0 1]);
ARGS2{1}{7} = coder.typeof(0,[1 Inf],[0 1]);

%% Invoke MATLAB Coder.
codegen -config cfg BsplineComposeImage3D -args ARGS2{1} -o BsplineComposeImage3D

%% Define argument types for entry-point 'BsplineComposeImage3D_single'.
ARGS3 = cell(1,1);
ARGS3{1} = cell(7,1);
ARGS3{1}{1} = coder.typeof(single(0), [Inf,4,4],[0,1]);
ARGS3{1}{2} = coder.typeof(single(0), [Inf,4,4],[0,1]);
ARGS3{1}{3} = coder.typeof(single(0), [Inf,4,4],[0,1]);
ARGS3{1}{4} = coder.typeof(0, [Inf,Inf,Inf],[0,1]);
ARGS3{1}{5} = coder.typeof(0);
ARGS3{1}{6} = coder.typeof(0);
ARGS3{1}{7} = coder.typeof(0);

%% Invoke MATLAB Coder.
codegen -config cfg BsplineComposeImage3D_single -args ARGS3{1} -o BsplineComposeImage3D_single

%% Define argument types for entry-point 'computenewPoints'.
ARGS4 = cell(1,1);
ARGS4{1} = cell(7,1);

ARGS4{1}{1} = struct;
ARGS4{1}{1}.nzsplines = coder.typeof(int64(0), [Inf,1],[1,0]);
ARGS4{1}{1} = coder.typeof(ARGS4{1}{1}, [Inf,1],[1,0]);

ARGS4{1}{2} = coder.typeof(double(0), [Inf,4],[1,0]);

ARGS4{1}{3} = struct;
ARGS4{1}{3}.mat = coder.typeof(single(0), [Inf,4,4,4],[1,0,0,0]);
ARGS4{1}{3} = coder.typeof(ARGS4{1}{3}, [Inf,1],[1,0]);

ARGS4{1}{4} = struct;
ARGS4{1}{4}.mat = coder.typeof(single(0), [Inf,4,4,4],[1,0,0,0]);
ARGS4{1}{4} = coder.typeof(ARGS4{1}{4}, [Inf,1],[1,0]);

ARGS4{1}{5} = struct;
ARGS4{1}{5}.mat = coder.typeof(single(0), [Inf,4,4,4],[1,0,0,0]);
ARGS4{1}{5} = coder.typeof(ARGS4{1}{5}, [Inf,1],[1,0]);

ARGS4{1}{6} = struct;
ARGS4{1}{6}.mat = coder.typeof(single(0), [Inf,4,4,4],[1,0,0,0]);
ARGS4{1}{6} = coder.typeof(ARGS4{1}{6}, [Inf,1],[1,0]);

ARGS4{1}{7} = coder.typeof(0);

%% Invoke MATLAB Coder.
codegen -config cfg computenewPoints -args ARGS4{1} -o computenewPoints

%% Define argument types for entry-point 'tripleIterLoop'.
ARGS5 = cell(1,1);
ARGS5{1} = cell(4,1);
ARGS5{1}{1} = coder.typeof(0,[1 3]);

ARGS5{1}{2} = struct;
ARGS5{1}{2}.active_cell = coder.typeof(0);
ARGS5{1}{2}.phi = coder.typeof(single(0),[Inf  1],[1 0]);
ARGS5{1}{2} = coder.typeof(ARGS5{1}{2},[Inf  1],[1 0]);

ARGS5{1}{3} = struct;
ARGS5{1}{3}.nzsplines = coder.typeof(int64(0),[Inf  1],[1 0]);
ARGS5{1}{3} = coder.typeof(ARGS5{1}{3},[Inf  1],[1 0]);
ARGS5{1}{4} = coder.typeof(0,[Inf  4],[1 0]);

%% Invoke MATLAB Coder.
codegen -config cfg tripleIterLoop -args ARGS5{1}

%% Define argument types for entry-point 'compute_Integ_Domainf'.
ARGS6 = cell(1,1);
ARGS6{1} = cell(50,1);
ARGS6{1}{1} = struct;
ARGS6{1}{1}.nzsplines = coder.typeof(int64(0),[Inf  1],[1 0]);
ARGS6{1}{1} = coder.typeof(ARGS6{1}{1},[Inf  1],[1 0]);

ARGS6{1}{2} = coder.typeof(0,[Inf  4  4],[1 0 0]);
ARGS6{1}{3} = coder.typeof(0,[Inf  4  4],[1 0 0]);
ARGS6{1}{4} = coder.typeof(0,[Inf  4  4],[1 0 0]);
ARGS6{1}{5} = coder.typeof(0,[Inf  4  4],[1 0 0]);

ARGS6{1}{6} = coder.typeof(single(0),[Inf  4  4],[1 0 0]);
ARGS6{1}{7} = coder.typeof(single(0),[Inf  4  4],[1 0 0]);
ARGS6{1}{8} = coder.typeof(single(0),[Inf  4  4],[1 0 0]);

ARGS6{1}{9} = coder.typeof(single(0),[Inf  4  4],[1 0 0]);
ARGS6{1}{10} = coder.typeof(single(0),[Inf  4  4],[1 0 0]);
ARGS6{1}{11} = coder.typeof(single(0),[Inf  4  4],[1 0 0]);

ARGS6{1}{12} = coder.typeof(single(0),[Inf  4  4],[1 0 0]);
ARGS6{1}{13} = coder.typeof(single(0),[Inf  4  4],[1 0 0]);
ARGS6{1}{14} = coder.typeof(single(0),[Inf  4  4],[1 0 0]);

ARGS6{1}{15} = coder.typeof(single(0),[Inf  4  4],[1 0 0]);
ARGS6{1}{16} = coder.typeof(single(0),[Inf  4  4],[1 0 0]);
ARGS6{1}{17} = coder.typeof(single(0),[Inf  4  4],[1 0 0]);

ARGS6{1}{18} = coder.typeof(single(0),[Inf  4  4],[1 0 0]);
ARGS6{1}{19} = coder.typeof(single(0),[Inf  4  4],[1 0 0]);
ARGS6{1}{20} = coder.typeof(single(0),[Inf  4  4],[1 0 0]);

ARGS6{1}{21} = coder.typeof(single(0),[Inf  4  4],[1 0 0]);
ARGS6{1}{22} = coder.typeof(single(0),[Inf  4  4],[1 0 0]);
ARGS6{1}{23} = coder.typeof(single(0),[Inf  4  4],[1 0 0]);

ARGS6{1}{24} = coder.typeof(single(0),[Inf  4  4],[1 0 0]);
ARGS6{1}{25} = coder.typeof(single(0),[Inf  4  4],[1 0 0]);
ARGS6{1}{26} = coder.typeof(single(0),[Inf  4  4],[1 0 0]);

ARGS6{1}{27} = coder.typeof(single(0),[Inf  4  4],[1 0 0]);
ARGS6{1}{28} = coder.typeof(single(0),[Inf  4  4],[1 0 0]);
ARGS6{1}{29} = coder.typeof(single(0),[Inf  4  4],[1 0 0]);

ARGS6{1}{30} = coder.typeof(single(0),[Inf  4  4],[1 0 0]);
ARGS6{1}{31} = coder.typeof(single(0),[Inf  4  4],[1 0 0]);
ARGS6{1}{32} = coder.typeof(single(0),[Inf  4  4],[1 0 0]);

ARGS6{1}{33} = coder.typeof(single(0),[Inf  4  4],[1 0 0]);
ARGS6{1}{34} = coder.typeof(single(0),[Inf  4  4],[1 0 0]);
ARGS6{1}{35} = coder.typeof(single(0),[Inf  4  4],[1 0 0]);

ARGS6{1}{36} = coder.typeof(0,[Inf  4],[1 0]);

ARGS6{1}{37} = struct;
ARGS6{1}{37}.mat = coder.typeof(single(0),[Inf  4  4  4],[1 0 0 0]);
ARGS6{1}{37} = coder.typeof(ARGS6{1}{37},[Inf  1],[1 0]);

ARGS6{1}{38} = struct;
ARGS6{1}{38}.mat = coder.typeof(single(0),[Inf  4  4  4],[1 0 0 0]);
ARGS6{1}{38} = coder.typeof(ARGS6{1}{38},[Inf  1],[1 0]);

ARGS6{1}{39} = struct;
ARGS6{1}{39}.mat = coder.typeof(single(0),[Inf  4  4  4],[1 0 0 0]);
ARGS6{1}{39} = coder.typeof(ARGS6{1}{39},[Inf  1],[1 0]);

ARGS6{1}{40} = struct;
ARGS6{1}{40}.mat = coder.typeof(single(0),[Inf  4  4  4],[1 0 0 0]);
ARGS6{1}{40} = coder.typeof(ARGS6{1}{40},[Inf  1],[1 0]);

ARGS6{1}{41} = struct;
ARGS6{1}{41}.mat = coder.typeof(single(0),[Inf  4  4  4],[1 0 0 0]);
ARGS6{1}{41} = coder.typeof(ARGS6{1}{41},[Inf  1],[1 0]);

ARGS6{1}{42} = coder.typeof(0);
ARGS6{1}{43} = coder.typeof(0);
ARGS6{1}{44} = coder.typeof(0);

ARGS6{1}{45} = coder.typeof(0);
ARGS6{1}{46} = coder.typeof(0);

ARGS6{1}{47} = coder.typeof(0,[Inf  1],[1 0]);
ARGS6{1}{48} = coder.typeof(0,[Inf  1],[1 0]);
ARGS6{1}{49} = coder.typeof(0,[Inf  1],[1 0]);

ARGS6{1}{50} = coder.typeof(single(0),[Inf  3],[1 0]);
%% Invoke MATLAB Coder.
codegen -config cfg compute_Integ_Domainf_mod -args ARGS6{1}

%% Define argument types for entry-point 'compute_Integ_Domainb_mod'.
ARGS7 = cell(1,1);
ARGS7{1} = cell(45,1);
ARGS7{1}{1} = struct;
ARGS7{1}{1}.nzsplines = coder.typeof(int64(0),[Inf  1],[1 0]);
ARGS7{1}{1} = coder.typeof(ARGS7{1}{1},[Inf  1],[1 0]);

ARGS7{1}{2} = coder.typeof(single(0),[Inf  4  4],[1 0 0]);
ARGS7{1}{3} = coder.typeof(single(0),[Inf  4  4],[1 0 0]);
ARGS7{1}{4} = coder.typeof(single(0),[Inf  4  4],[1 0 0]);

ARGS7{1}{5} = coder.typeof(single(0),[Inf  4  4],[1 0 0]);
ARGS7{1}{6} = coder.typeof(single(0),[Inf  4  4],[1 0 0]);
ARGS7{1}{7} = coder.typeof(single(0),[Inf  4  4],[1 0 0]);

ARGS7{1}{8} = coder.typeof(single(0),[Inf  4  4],[1 0 0]);
ARGS7{1}{9} = coder.typeof(single(0),[Inf  4  4],[1 0 0]);
ARGS7{1}{10} = coder.typeof(single(0),[Inf  4  4],[1 0 0]);

ARGS7{1}{11} = coder.typeof(single(0),[Inf  4  4],[1 0 0]);
ARGS7{1}{12} = coder.typeof(single(0),[Inf  4  4],[1 0 0]);
ARGS7{1}{13} = coder.typeof(single(0),[Inf  4  4],[1 0 0]);

ARGS7{1}{14} = coder.typeof(single(0),[Inf  4  4],[1 0 0]);
ARGS7{1}{15} = coder.typeof(single(0),[Inf  4  4],[1 0 0]);
ARGS7{1}{16} = coder.typeof(single(0),[Inf  4  4],[1 0 0]);

ARGS7{1}{17} = coder.typeof(single(0),[Inf  4  4],[1 0 0]);
ARGS7{1}{18} = coder.typeof(single(0),[Inf  4  4],[1 0 0]);
ARGS7{1}{19} = coder.typeof(single(0),[Inf  4  4],[1 0 0]);

ARGS7{1}{20} = coder.typeof(single(0),[Inf  4  4],[1 0 0]);
ARGS7{1}{21} = coder.typeof(single(0),[Inf  4  4],[1 0 0]);
ARGS7{1}{22} = coder.typeof(single(0),[Inf  4  4],[1 0 0]);

ARGS7{1}{23} = coder.typeof(single(0),[Inf  4  4],[1 0 0]);
ARGS7{1}{24} = coder.typeof(single(0),[Inf  4  4],[1 0 0]);
ARGS7{1}{25} = coder.typeof(single(0),[Inf  4  4],[1 0 0]);

ARGS7{1}{26} = coder.typeof(single(0),[Inf  4  4],[1 0 0]);
ARGS7{1}{27} = coder.typeof(single(0),[Inf  4  4],[1 0 0]);
ARGS7{1}{28} = coder.typeof(single(0),[Inf  4  4],[1 0 0]);

ARGS7{1}{29} = coder.typeof(single(0),[Inf  4  4],[1 0 0]);
ARGS7{1}{30} = coder.typeof(single(0),[Inf  4  4],[1 0 0]);
ARGS7{1}{31} = coder.typeof(single(0),[Inf  4  4],[1 0 0]);

ARGS7{1}{32} = coder.typeof(0,[Inf  4],[1 0]);

ARGS7{1}{33} = struct;
ARGS7{1}{33}.mat = coder.typeof(single(0),[Inf  4  4  4],[1 0 0 0]);
ARGS7{1}{33} = coder.typeof(ARGS7{1}{33},[Inf  1],[1 0]);

ARGS7{1}{34} = struct;
ARGS7{1}{34}.mat = coder.typeof(single(0),[Inf  4  4  4],[1 0 0 0]);
ARGS7{1}{34} = coder.typeof(ARGS7{1}{34},[Inf  1],[1 0]);

ARGS7{1}{35} = struct;
ARGS7{1}{35}.mat = coder.typeof(single(0),[Inf  4  4  4],[1 0 0 0]);
ARGS7{1}{35} = coder.typeof(ARGS7{1}{35},[Inf  1],[1 0]);

ARGS7{1}{36} = struct;
ARGS7{1}{36}.mat = coder.typeof(single(0),[Inf  4  4  4],[1 0 0 0]);
ARGS7{1}{36} = coder.typeof(ARGS7{1}{36},[Inf  1],[1 0]);

ARGS7{1}{37} = struct;
ARGS7{1}{37}.mat = coder.typeof(single(0),[Inf  4  4  4],[1 0 0 0]);
ARGS7{1}{37} = coder.typeof(ARGS7{1}{37},[Inf  1],[1 0]);

ARGS7{1}{38} = coder.typeof(0);
ARGS7{1}{39} = coder.typeof(0);

ARGS7{1}{40} = coder.typeof(0);

ARGS7{1}{41} = coder.typeof(0);

ARGS7{1}{42} = coder.typeof(0,[Inf  1],[1 0]);
ARGS7{1}{43} = coder.typeof(0,[Inf  1],[1 0]);
ARGS7{1}{44} = coder.typeof(0,[Inf  1],[1 0]);

ARGS7{1}{45} = coder.typeof(single(0),[Inf  3],[1 0]);
%% Invoke MATLAB Coder.
codegen -config cfg compute_Integ_Domainb_mod -args ARGS7{1}

%% Define argument types for entry-point 'BsplineCompose3D'.
ARGS8 = cell(1,1);
ARGS8{1} = cell(9,1);

ARGS8{1}{1} = coder.typeof(0);
ARGS8{1}{2} = coder.typeof(0);
ARGS8{1}{3} = coder.typeof(0);

ARGS8{1}{4} = coder.typeof(0,[1 Inf],[0 1]);
ARGS8{1}{5} = coder.typeof(0,[1 Inf],[0 1]);
ARGS8{1}{6} = coder.typeof(0,[1 Inf],[0 1]);

ARGS8{1}{7} = coder.typeof(0,[1 Inf],[0 1]);
ARGS8{1}{8} = coder.typeof(0,[1 Inf],[0 1]);
ARGS8{1}{9} = coder.typeof(0,[1 Inf],[0 1]);

%% Invoke MATLAB Coder.
codegen -config cfg BsplineCompose3D -args ARGS8{1}

%% Define argument types for entry-point 'storePixelPhi'.
ARGS = cell(1,1);
ARGS{1} = cell(9,1);
ARGS{1}{1} = coder.typeof(int64(0));
ARGS{1}{2} = coder.typeof(0);
ARGS{1}{3} = coder.typeof(0,[Inf  3],[1 0]);
ARG = coder.typeof(0,[1 Inf],[0 1]);
ARGS{1}{4} = coder.typeof({ARG}, [Inf  1],[1 0]);
ARG = coder.typeof(0,[1 Inf],[0 1]);
ARGS{1}{5} = coder.typeof({ARG}, [Inf  1],[1 0]);
ARG = coder.typeof(0,[1 Inf],[0 1]);
ARGS{1}{6} = coder.typeof({ARG}, [Inf  1],[1 0]);
ARGS{1}{7} = struct;
ARGS{1}{7}.knot_ind = coder.typeof(0,[Inf  3  2],[1 0 0]);
ARGS{1}{7}.flag_active = coder.typeof(0,[Inf  1],[1 0]);
ARGS{1}{7}.IEN = coder.typeof(0,[Inf  27],[1 0]);
ARGS{1}{7}.chdElem = coder.typeof(0,[Inf  8],[1 0]);
ARGS{1}{7}.cell_centre = coder.typeof(0,[Inf  3],[1 0]);
ARGS{1}{7}.node = coder.typeof(0);
ARGS{1}{7}.parElem = coder.typeof(0,[Inf  1],[1 0]);
ARGS{1}{7}.actE = coder.typeof(0,[Inf  1],[1 0]);
ARGS{1}{7} = coder.typeof(ARGS{1}{7},[Inf  1],[1 0]);
ARGS{1}{8} = struct;
ARGS{1}{8}.mat = coder.typeof(single(0),[Inf  27],[1 0]);
ARGS{1}{8} = coder.typeof(ARGS{1}{8},[Inf  1],[1 0]);
ARGS{1}{9} = struct;
ARGS{1}{9}.pU = coder.typeof(0);
ARGS{1}{9}.pV = coder.typeof(0);
ARGS{1}{9}.pW = coder.typeof(0);
ARGS{1}{9}.maxlevel = coder.typeof(0);
ARGS{1}{9}.nelemU = coder.typeof(0,[Inf  1],[1 0]);
ARGS{1}{9}.nelemV = coder.typeof(0,[Inf  1],[1 0]);
ARGS{1}{9}.nelemW = coder.typeof(0,[Inf  1],[1 0]);
ARGS{1}{9}.orderGauss = coder.typeof(0);
ARGS{1}{9}.kU = coder.typeof(0,[Inf  1],[1 0]);
ARGS{1}{9}.kV = coder.typeof(0,[Inf  1],[1 0]);
ARGS{1}{9}.kW = coder.typeof(0,[Inf  1],[1 0]);
ARGS{1}{9}.nobU = coder.typeof(0,[Inf  1],[1 0]);
ARGS{1}{9}.nobV = coder.typeof(0,[Inf  1],[1 0]);
ARGS{1}{9}.nobW = coder.typeof(0,[Inf  1],[1 0]);
ARGS{1}{9}.rho = coder.typeof(0,[3 1]);
ARGS{1}{9}.timestep = coder.typeof(0,[Inf  1],[1 0]);
ARGS{1}{9}.smallNumber = coder.typeof(0);
ARGS{1}{9}.lambda_1 = coder.typeof(0);
ARGS{1}{9}.lambda_2 = coder.typeof(0);
ARGS{1}{9}.lambda_3 = coder.typeof(0);
ARGS{1}{9}.mu= coder.typeof(0);
ARGS{1}{9}.par1= coder.typeof(0);
ARGS{1}{9}.par2= coder.typeof(0);
ARGS{1}{9} = coder.typeof(ARGS{1}{9});

%% Invoke MATLAB Coder.
codegen -config cfg storePixelPhi -args ARGS{1}

%% Define argument types for entry-point 'GaussPhi'.
ARGS = cell(1,1);
ARGS{1} = cell(8,1);
ARGS{1}{1} = coder.typeof(0,[Inf  2],[1 0]);
ARGS{1}{2} = struct;
ARGS{1}{2}.knot_ind = coder.typeof(0,[Inf  3  2],[1 0 0]);
ARGS{1}{2}.flag_active = coder.typeof(0,[Inf  1],[1 0]);
ARGS{1}{2}.IEN = coder.typeof(0,[Inf  27],[1 0]);
ARGS{1}{2}.chdElem = coder.typeof(0,[Inf  8],[1 0]);
ARGS{1}{2}.cell_centre = coder.typeof(0,[Inf  3],[1 0]);
ARGS{1}{2}.node = coder.typeof(0,[Inf  1],[1 0]);
ARGS{1}{2}.parElem = coder.typeof(0,[Inf  1],[1 0]);
ARGS{1}{2}.actE = coder.typeof(0,[Inf  1],[1 0]);
ARGS{1}{2} = coder.typeof(ARGS{1}{2},[Inf  1],[1 0]);
ARG = coder.typeof(0,[1 Inf],[0 1]);
ARGS{1}{3} = coder.typeof({ARG}, [Inf  1],[1 0]);
ARG = coder.typeof(0,[1 Inf],[0 1]);
ARGS{1}{4} = coder.typeof({ARG}, [Inf  1],[1 0]);
ARG = coder.typeof(0,[1 Inf],[0 1]);
ARGS{1}{5} = coder.typeof({ARG}, [Inf  1],[1 0]);
ARGS{1}{6} = struct;
ARGS{1}{6}.mat = coder.typeof(single(0),[Inf  27],[1 0]);
ARGS{1}{6} = coder.typeof(ARGS{1}{6},[Inf  1],[1 0]);
ARGS{1}{7} = struct;
ARGS{1}{7}.pU = coder.typeof(0);
ARGS{1}{7}.pV = coder.typeof(0);
ARGS{1}{7}.pW = coder.typeof(0);
ARGS{1}{7}.maxlevel = coder.typeof(0);
ARGS{1}{7}.nelemU = coder.typeof(0,[Inf  1],[1 0]);
ARGS{1}{7}.nelemV = coder.typeof(0,[Inf  1],[1 0]);
ARGS{1}{7}.nelemW = coder.typeof(0,[Inf  1],[1 0]);
ARGS{1}{7}.orderGauss = coder.typeof(0);
ARGS{1}{7}.kU = coder.typeof(0,[Inf  1],[1 0]);
ARGS{1}{7}.kV = coder.typeof(0,[Inf  1],[1 0]);
ARGS{1}{7}.kW = coder.typeof(0,[Inf  1],[1 0]);
ARGS{1}{7}.nobU = coder.typeof(0,[Inf  1],[1 0]);
ARGS{1}{7}.nobV = coder.typeof(0,[Inf  1],[1 0]);
ARGS{1}{7}.nobW = coder.typeof(0,[Inf  1],[1 0]);
ARGS{1}{7}.rho = coder.typeof(0,[Inf  1],[1 0]);
ARGS{1}{7}.timestep = coder.typeof(0,[Inf  1],[1 0]);
ARGS{1}{7}.smallNumber = coder.typeof(0);
ARGS{1}{7}.lambda_1 = coder.typeof(0);
ARGS{1}{7}.lambda_2 = coder.typeof(0);
ARGS{1}{7}.lambda_3 = coder.typeof(0);
ARGS{1}{7}.mu= coder.typeof(0);
ARGS{1}{7}.par1= coder.typeof(0);
ARGS{1}{7}.par2= coder.typeof(0);
ARGS{1}{7} = coder.typeof(ARGS{1}{7});
ARGS{1}{8} = coder.typeof(0);
%% Invoke MATLAB Coder.
codegen -config cfg GaussPhi -args ARGS{1}

%% Define argument types for entry-point 'GaussPhi_fb'.
ARGS = cell(1,1);
ARGS{1} = cell(8,1);
ARGS{1}{1} = coder.typeof(0,[Inf  2],[1 0]);
ARGS{1}{2} = coder.typeof(single(0),[Inf  4  4],[1 0]);
ARGS{1}{3} = coder.typeof(single(0),[Inf  4  4],[1 0]);
ARGS{1}{4} = coder.typeof(single(0),[Inf  4  4],[1 0]);
ARGS{1}{5} = struct;
ARGS{1}{5}.nzsplines = coder.typeof(int64(0), [Inf,1],[1,0]);
ARGS{1}{5} = coder.typeof(ARGS4{1}{1}, [Inf,1],[1,0]);
ARGS{1}{6} = coder.typeof({ARG}, [Inf  1],[1 0]);
ARG = coder.typeof(0,[1 Inf],[0 1]);
ARGS{1}{7} = coder.typeof({ARG}, [Inf  1],[1 0]);
ARG = coder.typeof(0,[1 Inf],[0 1]);
ARGS{1}{8} = coder.typeof({ARG}, [Inf  1],[1 0]);
ARG = coder.typeof(0,[1 Inf],[0 1]);
ARGS{1}{9} = struct;
ARGS{1}{9}.mat = coder.typeof(single(0),[Inf  27],[1 0]);
ARGS{1}{9} = coder.typeof(ARGS{1}{9},[Inf  1],[1 0]);
ARGS{1}{10} = struct;
ARGS{1}{10}.pU = coder.typeof(0);
ARGS{1}{10}.pV = coder.typeof(0);
ARGS{1}{10}.pW = coder.typeof(0);
ARGS{1}{10}.maxlevel = coder.typeof(0);
ARGS{1}{10}.nelemU = coder.typeof(0,[Inf  1],[1 0]);
ARGS{1}{10}.nelemV = coder.typeof(0,[Inf  1],[1 0]);
ARGS{1}{10}.nelemW = coder.typeof(0,[Inf  1],[1 0]);
ARGS{1}{10}.orderGauss = coder.typeof(0);
ARGS{1}{10}.kU = coder.typeof(0,[Inf  1],[1 0]);
ARGS{1}{10}.kV = coder.typeof(0,[Inf  1],[1 0]);
ARGS{1}{10}.kW = coder.typeof(0,[Inf  1],[1 0]);
ARGS{1}{10}.nobU = coder.typeof(0,[Inf  1],[1 0]);
ARGS{1}{10}.nobV = coder.typeof(0,[Inf  1],[1 0]);
ARGS{1}{10}.nobW = coder.typeof(0,[Inf  1],[1 0]);
ARGS{1}{10}.rho = coder.typeof(0,[Inf  1],[1 0]);
ARGS{1}{10}.timestep = coder.typeof(0,[Inf  1],[1 0]);
ARGS{1}{10}.smallNumber = coder.typeof(0);
ARGS{1}{10}.lambda_1 = coder.typeof(0);
ARGS{1}{10}.lambda_2 = coder.typeof(0);
ARGS{1}{10}.lambda_3 = coder.typeof(0);
ARGS{1}{10}.mu= coder.typeof(0);
ARGS{1}{10}.par1= coder.typeof(0);
ARGS{1}{10}.par2= coder.typeof(0);
ARGS{1}{10} = coder.typeof(ARGS{1}{10});
%% Invoke MATLAB Coder.
codegen -config cfg GaussPhi_fb -args ARGS{1}