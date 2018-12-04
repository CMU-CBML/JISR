function parameters = setparameters_triangle()

rho(1,1) = 0.001;
rho(2,1) = 0.05;
rho(3,1) = 0.10;

orderGauss = 6;

maxlevel = 3;

lambda1 = 0.1;
lambda2 = 0.1;
lambda3 = 0.01;

maxiteration = 100;

timestep(1,1) = 0.35;
timestep(2,1) = 0.35;
timestep(3,1) = 0.35;

nelemx = 25;
nelemy = 25;

pU = 3;
pV = 3;

par1 = 100;
par2 = 1;

parameters = struct('pU',pU,'pV',pV,'maxlevel',maxlevel,'nelemx',nelemx,'nelemy',nelemy,...
    'orderGauss',orderGauss,'rho',rho,'timestep',timestep,'lambda1',lambda1,'lambda2',lambda2,'lambda3',lambda3,'par1',par1,'par2',par2,'maxiteration',maxiteration);
end