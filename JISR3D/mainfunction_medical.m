% Main function
clc
clear all
close all

%% Read Image, Initialize
tic
addpath('../3D_adaptive_image_registration73/images3D');
% if computing Dice Similarity for medical images from Brainweb, setflagDS
% = 1;
setflagDS = 0;
saveVTK = 0;
plotImage = 1; %plot image for display and save the image in .png format
saveImage = 1; %save image result in .mat format

output_file='Dice_similarity_sub20_38_2step.txt';
fid_out = fopen(output_file,'w');
fprintf(fid_out, 'Level Iteration Background  CSF  GrayMatter  WhiteMatter  Fat  Muscle Muscle/Skin  Skull  vessels  AroundFat  DuraMatter  BoneMarrow\r\n');
disp('Added paths....');
%add the paths of subfolders of the software
%addpaths();
%set the parameters for running registration process
disp('Setting parameters...');
param = setparameters_brainweb_2038();

filenameI0='subject20_crisp_v.rawb';
filenameI1='subject38_crisp_v.rawb';
filenamephi1 = 'subject20_gm_v.rawb';
filenamephi2 = 'subject38_gm_v.rawb';
[I1,I2] = read3DImage_brainWeb(filenameI0,filenameI1);
[phi1,phi2] = read3DImage_brainWeb_phi(filenamephi1,filenamephi2);

I1_in = I1; %store initial moving image
I2_in = I2; %store initial target image

F = permute(I2,[1,3,2]);
M = permute(I1,[1,3,2]);
phi = permute(phi1,[1,3,2]);

phi = (phi-min(phi(:)))./(max(phi(:))-min(phi(:)));
M = (M-min(M(:)))./(max(M(:))-min(M(:)));
F = (F-min(F(:)))./(max(F(:))-min(F(:)));

nx = size(F,1);
ny = size(F,2);
nz = size(F,3);

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
phi_img = 255.*ones(ny,nx,nz);
phi_img(phi<mean(phi(:))) = -255;
hold on
contour(phi_img(:,:,floor(nz/2)),[0,0],'g','Linewidth',1.24);
set(gca,'position',[0 0 1 1],'units','normalized')
hold off

disp('display images');
toc

[M,F,phi,VF,VB] = MultipleResolution3D_medical(M,F,phi,param);




