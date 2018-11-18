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

parameters = setparameters_lung();

%% Read Image, Initialize
fid1 = fopen('case1_T00_s.img');
fid2 = fopen('case1_T50_s.img');

Img1 = fread(fid1,'*int16');
Img1 = double(Img1);
Img1=(Img1-min(Img1(:)))/(max(Img1(:))-min(Img1(:)));
Img1 = reshape(Img1,[256,256,94]);
Img1 = resize3Dmatrix(256,256,256,Img1);
Img1 = squeeze(Img1(:,128,:));
Img1 = Img1';
fclose(fid1);

Img2 = fread(fid2,'*int16');
Img2 = double(Img2);
Img2=(Img2-min(Img2(:)))/(max(Img2(:))-min(Img2(:)));
Img2 = reshape(Img2,[256,256,94]);
Img2 = resize3Dmatrix(256,256,256,Img2);
Img2 = squeeze(Img2(:,128,:));
Img2 = Img2';
fclose(fid2);

[X,Y] = meshgrid(0.5:size(Img1,2)-0.5,0.5:size(Img1,1)-0.5);
load 'phi_source_lung.mat';

F = Img2;
M = Img1;
phi = phi_img;

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
[M,F,phi,VF,VB,In_gridX, In_gridY,FgridX, FgridY] = MultipleResolution2D(M,F,phi,parameters);
toc

mlabel = ConvertIntensity2D_binary(M);
flabel = ConvertIntensity2D_binary(F00);
dice = DiceSimilarity2D(nx,ny,mlabel,flabel,2);