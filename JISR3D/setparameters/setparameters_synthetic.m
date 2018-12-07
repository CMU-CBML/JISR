function param = setparameters_synthetic()
%In this function, set the required parameters for the input to the
%registration code

%degree of B-splines
pU = 2;
pV = 2;
pW = 2;

% maximum number of levels
maxlevel = 3;
maxlevelg = 3;

% order of gauss quadrature
orderGauss = 4;

%initial grid: number of elements in x,y,z direction for level 1
m_var = 5;
n_var = 5;
o_var = 5;

%regularization parameters 
lambda_1 = 0.1;
lambda_2 = 0;
lambda_3 = 10;
mu = 0.001;

par1 = 50;
par2 = 1;

%gamma in g(x)term
smallNumber = 10^-12;

%number of elements in each direction at each level
nelemU = zeros(maxlevelg,1);
nelemV = zeros(maxlevelg,1);
nelemW = zeros(maxlevelg,1);

%number of splines in each direction at each level
nobU = zeros(maxlevelg,1);
nobV = zeros(maxlevelg,1);
nobW = zeros(maxlevelg,1);

%number of knot vectors in each direction at each level
kU = zeros(maxlevelg,1); 
kV = zeros(maxlevelg,1);
kW = zeros(maxlevelg,1);

rho = zeros(maxlevelg,1); %refinement parameter
timestep = zeros(maxlevelg,1); %timestep for each refinement level

%number of elements in each direction
for level = 1:maxlevelg
nelemU(level,1) = m_var*2^(level-1);
nelemV(level,1) = n_var*2^(level-1);
nelemW(level,1) = o_var*2^(level-1);

kU(level,1) = nelemU(level,1)+2*pU+1;
kV(level,1) = nelemV(level,1)+2*pU+1;
kW(level,1) = nelemW(level,1)+2*pU+1;

nobU(level,1) = kU(level,1) - pU - 1;
nobV(level,1) = kV(level,1) - pV - 1;
nobW(level,1) = kW(level,1) - pW - 1;
end

rho(1,1) = 0.5; %level 2 refinement
rho(2,1) = 1; %level 3 refinement
rho(3,1) = 1;

timestep(1,1) = 0.08;
timestep(2,1) = 0.15;
timestep(3,1) = 0.25;

%make a struct variable param, with all the parameters
param = struct('pU',pU,'pV',pV,'pW',pW,'maxlevel',maxlevel,'nelemU',nelemU,'nelemV',nelemV,'nelemW',nelemW,...
    'orderGauss',orderGauss,'kU',kU,'kV',kV,'kW',kW,'nobU',nobU,...
    'nobV',nobV,'nobW',nobW,'rho',rho,'timestep',timestep,'smallNumber',smallNumber,'lambda_1',lambda_1,'lambda_2',lambda_2,'lambda_3',lambda_3,'mu',mu,'par1',par1,'par2',par2);
end