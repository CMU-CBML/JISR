% Main function
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

parameters = setparameters_triangle();

%% Read Image, Initialize
load 'triangle.mat';
load 'blackcircle.mat';

F = triangle;
M = blackcircle;
phi = M;

nx = size(F,1);
ny = size(F,2);

phi = (phi-min(phi(:)))./(max(phi(:))-min(phi(:)));
M = (M-min(M(:)))./(max(M(:))-min(M(:)));
F = (F-min(F(:)))./(max(F(:))-min(F(:)));

F00 = F;
M00 = M;
phi00 = M;
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
[M,F,phi,VF,VB,In_gridX, In_gridY,FgridX, FgridY] = MultipleResolution2D_triangle(M,F,phi,parameters);
toc

load 'triangle.mat';
T = triangle;
T = (T-min(T(:)))./(max(T(:))-min(T(:)));

mlabel = ConvertIntensity2D_binary(M);
flabel = ConvertIntensity2D_binary(T);
dice = DiceSimilarity2D(nx,ny,mlabel,flabel,2);

sum11 = (M-T).^2;
sum1 = sum(sum11(:));
msd = sum1./(255*255);

figure
imagesc(F00)
colormap gray
hold on
phi_img = 255.*ones(255,255);
phi_img(phi<=mean(phi(:))) = -255;
contour(phi_img,[0,0],'g','Linewidth',3.24);
set(gca,'position',[0 0 1 1],'units','normalized')