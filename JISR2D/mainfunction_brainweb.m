clc
clear all
close all

addpath('../JISR2Dcopy/readimg');
addpath('../JISR2Dcopy/images2D');
addpath('../JISR2Dcopy/bsplinelevelset');
addpath('../JISR2Dcopy/composition');
addpath('../JISR2Dcopy/setparameters');
addpath('../JISR2Dcopy/thbspline');
addpath('../JISR2Dcopy/iterationloop_funcs');
addpath('../JISR2Dcopy/levelsetseg_funcs');
addpath('../JISR2Dcopy/postprocessing');

parameters = setparameters_brainweb();

%% Read Image, Initialize
filenameI0='phantom_1.0mm_msles3_crisp.rawb';
filenamephi1 = 'phantom_1.0mm_msles3_wht.rawb';
[I1,I2] = read3DImage_brainWeb(filenameI0,filenameI0);
[phi1,phi1] = read3DImage_brainWeb_phi(filenamephi1,filenamephi1);


F = I1(:,:,113);
M = I1(:,:,111);
phi = phi1(:,:,111);
phi = (phi-min(phi(:)))./(max(phi(:))-min(phi(:)));
M = (M-min(M(:)))./(max(M(:))-min(M(:)));
F = (F-min(F(:)))./(max(F(:))-min(F(:)));


F00 = F;
M00 = M;
phi00 = phi;

nx = size(F,1);
ny = size(F,2);

%% Display intial config
figure
axis equal
imagesc(M)
colormap gray;
set(gca,'position',[0 0 1 1],'units','normalized')

figure
axis equal
imagesc(F)
colormap gray;
set(gca,'position',[0 0 1 1],'units','normalized')

figure
axis equal
imagesc(F-M)
colormap gray;
set(gca,'position',[0 0 1 1],'units','normalized')

figure
axis equal
imagesc(phi)
colormap gray;
set(gca,'position',[0 0 1 1],'units','normalized')
phi_img = 255.*ones(nx,ny);
phi_img(phi<mean(phi(:))) = -255;
hold on
contour(phi_img,[0,0],'g','Linewidth',3.24);
hold off

tic
[M,F,phi,VF,VB,In_gridX, In_gridY,FgridX, FgridY] = MultipleResolution2D_brainweb(M,F,phi,parameters);
toc

phi2 = phi1(:,:,113);
phi2 =(phi2-min(phi2(:)))/(max(phi2(:))-min(phi2(:)));
mlabel = ConvertIntensity2D_binary(phi);
flabel = ConvertIntensity2D_binary(phi2);
dice = DiceSimilarity2D(nx,ny,mlabel,flabel,2);

sum11 = (M-F00).^2;
sum1 = sum(sum11(:));
msd = sum1./(255*255);