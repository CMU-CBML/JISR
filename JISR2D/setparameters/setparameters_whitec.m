function parameters = setparameters_whitec()

rho(1,1) = 0.25;
rho(2,1) = 0.5;
rho(3,1) = 1;

orderGauss = 6;

maxlevel = 3;

lambda1 = 0.01;
lambda2 = 0.01;
lambda3 = 0.01;

maxiteration = 100;

timestep(1,1) = 0.40;
timestep(2,1) = 0.40;
timestep(3,1) = 0.40;

nelemx = 30;
nelemy = 30;

pU = 3;
pV = 3;

par1 = 100;
par2 = 1;

parameters = struct('pU',pU,'pV',pV,'maxlevel',maxlevel,'nelemx',nelemx,'nelemy',nelemy,...
    'orderGauss',orderGauss,'rho',rho,'timestep',timestep,'lambda1',lambda1,'lambda2',lambda2,'lambda3',lambda3,'par1',par1,'par2',par2,'maxiteration',maxiteration);
end