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

parameters = setparameters_brainweb();

%% Read Image, Initialize
load 'Img11.mat';
load 'Img22.mat';
load 'phi11.mat';

M = Img11(:,:,90);
F = Img11(:,:,92);
phiq = zeros(180,180,180);
phiq(phi11~=128) = 1;
phi2 = phiq(:,:,92);
phi = phiq(:,:,90);

%Img=(Img-min(Img(:)))/(max(Img(:))-min(Img(:)));
%Img_seg=(Img_seg-min(Img_seg(:)))/(max(Img_seg(:))-min(Img_seg(:)));

%F = Img(:,:,103);
%M = Img(:,:,100);
%phi = Img_seg(:,:,100);
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
[M,F,phi,VF,VB,In_gridX, In_gridY,FgridX, FgridY] = MultipleResolution2D(M,F,phi,parameters);
toc

mlabel = ConvertIntensity2D_binary(M);
flabel = ConvertIntensity2D_binary(F00);
dice = DiceSimilarity2D(nx,ny,mlabel,flabel,2);