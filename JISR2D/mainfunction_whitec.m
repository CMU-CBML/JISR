% Main function
clc
clear all
close all

addpath('../JISR2D/readimg');
addpath('../JISR2D/images2D');
addpath('../JISR2D/bsplinelevelset');
addpath('../JISR2D/composition');
addpath('../JISR2D/setparameters');
addpath('../JISR2D/thbspline');
addpath('../JISR2D/iterationloop_funcs');
addpath('../JISR2D/levelsetseg_funcs');
addpath('../JISR2D/postprocessing');

parameters = setparameters_whitec();

%% Read Image, Initialize
filename1 = 'whitecircle.png';
filename2 = 'whitec.png';
[Img1,Img2] = readImg(filename1,filename2);

F = Img2;
M = Img1;
phi = M;

F00 = Img2;
M00 = Img1;
phi00 = M;

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
[M,F,phi,VF,VB,In_gridX, In_gridY,FgridX, FgridY] = MultipleResolution2D_whitec(M,F,phi,parameters);
toc

mlabel = ConvertIntensity2D_binary(M);
flabel = ConvertIntensity2D_binary(F00);
dice = DiceSimilarity2D(nx,ny,mlabel,flabel,2);

sum11 = (M-F00).^2;
sum1 = sum(sum11(:));
msd = sum1./(255*255);