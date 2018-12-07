% Main function
clc
clear all
close all

%% Read Image, Initialize
tic
% if computing Dice Similarity for medical images from Brainweb, setflagDS
% = 1;
setflagDS = 0;
saveVTK = 0;
plotImage = 1; %plot image for display and save the image in .png format
saveImage = 1; %save image result in .mat format


disp('Added paths....');
addpath('../JISR3D/thbspline');
addpath('../JISR3D/setparameters');
addpath('../JISR3D/readimg');
addpath('../JISR3D/postprocessing');
addpath('../JISR3D/levelsetseg');
addpath('../JISR3D/iterationloop_funcs');
addpath('../JISR3D/images3D');
addpath('../JISR3D/bsplinelevelset');
%add the paths of subfolders of the software


%set the parameters for running registration process
disp('Setting parameters...');
param = setparameters_synthetic();

%load images
load Sphere_200.mat;
load bunny_img_200.mat;

Img1 = Sphere_200;
Img2 = bunny_img_200;

Img1 = resize3Dmatrix(50,50,50,Img1);
Img2 = resize3Dmatrix(50,50,50,Img2);

F = Img2;
M = Img1;
phi = M;

%Normalize image intensities
phi = (phi-min(phi(:)))./(max(phi(:))-min(phi(:)));
M = (M-min(M(:)))./(max(M(:))-min(M(:)));
F = (F-min(F(:)))./(max(F(:))-min(F(:)));

nx = size(M,1);
ny = size(M,2);
nz = size(M,3);

figure
axis equal
imagesc(M(:,:,floor(nz/2)))
colormap gray;
set(gca,'position',[0 0 1 1],'units','normalized')

figure
axis equal
imagesc(F(:,:,floor(nz/2)))
colormap gray;
set(gca,'position',[0 0 1 1],'units','normalized')

figure
axis equal
imagesc(F(:,:,floor(nz/2))-M(:,:,floor(nz/2)))
colormap gray;
set(gca,'position',[0 0 1 1],'units','normalized')

figure
axis equal
imagesc(phi(:,:,floor(nz/2)))
colormap gray;
set(gca,'position',[0 0 1 1],'units','normalized')
phi_img = 255.*ones(nx,ny,nz);
phi_img(phi<mean(phi(:))) = -255;
hold on
contour(phi_img(:,:,floor(nz/2)),[0,0],'g','Linewidth',1.24);
hold off

disp('display images');
toc
tic

%Start the multiresolution framework
[M,F,phi,VF,VB] = MultipleResolution3D_synthetic(M,F,phi,param);