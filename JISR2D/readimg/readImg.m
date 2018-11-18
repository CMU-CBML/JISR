function [Img1, Img2] = readImg(fn1,fn2)

Img11 = imread(fn1);
Img22 = imread(fn2);

Img1 = im2double(Img11(:,:,1));
Img2 = im2double(Img22(:,:,1));

% Img1 = double(Img11(:,:,1));
% Img2 = double(Img22(:,:,1));
end