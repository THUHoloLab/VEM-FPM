clc
clear



% load('results\data_saved_Red.mat');
% img_rgb_final(:,:,1) = abs(wavefront1).^2;
% load('results\data_saved_Green.mat');
% img_rgb_final(:,:,2) = abs(wavefront1).^2;
% load('results\data_saved_blue.mat');
% img_rgb_final(:,:,3) = abs(wavefront1).^2;

img_rgb_final(:,:,1) = double(imread('R\0\3\1.png'));
img_rgb_final(:,:,2) = double(imread('G\0\2\1.png'));
img_rgb_final(:,:,3) = double(imread('B\0\2\1.png'));

img_rgb_final = img_rgb_final / max(img_rgb_final(:));
%%%% 对齐和光强矫正
f = img_rgb_final;
% f_r = f(:,:,1);
% f_b = f(:,:,3);
% usfac = 1500;

f(:,:,1) = 0.7*f(:,:,1);
f(:,:,2) = 0.9*f(:,:,2);
f(:,:,3) = 0.7*f(:,:,3);

% f = img_rgb_raw;
f = gather(f);
title('请选择背景区域')
[temp,rect] = imcrop(f);
if rem(size(temp,1),2) == 1
    rect(4) = rect(4) - 1;
end
if rem(size(temp,2),2) == 1
    rect(3) = rect(3) - 1;
end
pix = fix((rect(4) + rect(3))/2);
pix = pix + mod(pix,2);
rect = fix(rect);
% close all

area = f(rect(2):rect(2)+pix-1,rect(1):rect(1)+pix-1,:);

r_sum = mean(mean(area(:,:,1)));
g_sum = mean(mean(area(:,:,2)));
b_sum = mean(mean(area(:,:,3)));

white = [230,230,230]/255;
f_ic = f;

f_ic(:,:,1) = f(:,:,1) .* white(1) / r_sum;
% f_ic(:,:,2) = imfilter(f(:,:,1),fspecial('gaussian',5,1.5),'symmetric') .* white(2) / g_sum;

% f_ic(:,:,2) = f(:,:,2) .* white(2) / g_sum;
f_ic(:,:,2) = imfilter(f(:,:,2),fspecial('gaussian',5,1.5),'symmetric') .* white(2) / g_sum;

% f_ic(:,:,3) = f(:,:,3) .* white(3) / b_sum;
f_ic(:,:,3) = imfilter(f(:,:,3),fspecial('gaussian',15,3),'symmetric') .* white(3) / b_sum;

figure();imshow(f_ic,[])

% save rgb_img f_ic
imwrite(f_ic,'data_ware_20231204_raw.png');



function out = it(in,ref)
    max1 = max(ref(:));
    min1 = min(ref(:));

    out = (max1 - min1) * mat2gray(in) + min1;
end