%{
    VEM-FPM for fpm reconstruction
    Using median filter for image denoising.
%}

clear
clc

addpath(genpath('func_ddfpm'));
path = {'datas/red/','datas/green/','datas/blue/'};

led_num = [1,8,12,16,24,32];
led_total = sum(led_num(:));
rot_ang = 0 / 180 * pi;

pix = 512;
rect = [1024-256,1024-256];


img_rgb_raw = zeros(pix,pix,3);


learning_rate = [0.01,0.01,0.05];  % learning rate for r,g,b
denoise_str = [0.01,0.01,0.01];    % denoising strength for r,g,b


for color_index = 1:3

% load loc_pos.mat
imRaw_new = zeros(pix,pix,led_total);

% Load data, crop image
for num_of_image = 1:led_total
    clc
    disp(num_of_image);
    img = single(imread([path{color_index},num2str(num_of_image),'.tif'], ...
                          'PixelRegion',{[rect(2),rect(2)+pix-1],...
                                         [rect(1),rect(1)+pix-1]}));

    imRaw_new(:,:,num_of_image) = mean(img,3);
end
% normalization data
imRaw_new = imRaw_new - min(imRaw_new(:));
imRaw_new = imRaw_new / max(imRaw_new(:));
imRaw_new = gpuArray(single(sqrt(imRaw_new)));

img_rgb_raw(:,:,color_index) = imRaw_new(:,:,1);

[f_pos_set_true,pratio,Pupil0] = misc.init_recon(color_index,...
                                                    pix, ...
                                                    led_num, ...
                                                    rot_ang);



f_pos_set_true();

fpm_cube = combine(arrayDatastore(f_pos_set_true, 'IterationDimension',1),...
                   arrayDatastore(imRaw_new, 'IterationDimension',3));

% (0 ~ 93) set mini-batch size a total of 93 images for FPM recon
batchSize = 24; 

fpm_cube = minibatchqueue(fpm_cube,...
            'MiniBatchSize',     batchSize,...
            'MiniBatchFormat',   ["",""],...
            'OutputEnvironment', {'gpu'},...
            'OutputAsDlarray',   false,...
            'OutputCast',        'single');

numEpochs = 50;
numIterationsPerEpoch  = size(imRaw_new,3) / batchSize;
numIterations = numEpochs * numIterationsPerEpoch;

epoch = 0;



%% The iterative recovery process for FP
disp('initializing parameters')


oI = (imresize(mean(imRaw_new(:,:,1),3),pratio)); 

dtd = abs(psf2otf([-1,1],[size(oI,1),size(oI,2)])).^2;
dtd = dtd + abs(psf2otf([-1;1],[size(oI,1),size(oI,2)])).^2;
dtd = gpuArray(single(fftshift(dtd)));

wavefront1 = gpuArray(fftshift(fft2(oI)));
wavefront2 = gpuArray(Pupil0);   

deconv_data.hto = wavefront1 * 0;  
deconv_data.hth = deconv_data.hto;
deconv_data.oth = 0;
deconv_data.oto = 0;

oI = gpuArray(imresize(mean(imRaw_new(:,:,1),3),pratio)); 

disp('begin solving-----')

v = 0;
u = 0;

error_bef = inf;

while epoch < numEpochs
    epoch = epoch + 1;

    fpm_cube.reset();

    tic
    deconv_data.hth = deconv_data.hth .* 0;
    deconv_data.hto = deconv_data.hto .* 0;
    deconv_data.oth = deconv_data.oth .* 0;
    deconv_data.oto = deconv_data.oto .* 0;
    
    iteration = 0;
    
    %% E-step: get latent image
    while fpm_cube.hasdata()
        iteration = iteration + 1;
        this_ratio = round(iteration/numIterationsPerEpoch * 1000)/1000;

        disp(['at ',num2str(epoch),'-epoch, loading for ---- ',...
                              num2str(this_ratio * 100),'%'])

        [leds,dY_obs] = fpm_cube.next();

        [loss,deconv_data] = helpers.E_step(wavefront1, ...
                                            wavefront2 , ...
                                            deconv_data, ...
                                            leds, ...
                                            dY_obs, ...
                                            pratio, ...
                                            learning_rate(color_index));
    end
    clc

    %% M-step: deconvolution learning the parameters
    % [wavefront1,oI,v,u,error_now] = helpers.M_step(deconv_data,...
    %                                                v,...
    %                                                u, ...
    %                                                dtd, ...
    %                                                denoise_str, ...
    %                                                oR, ...
    %                                                'retinex');

    [wavefront1,oI,error_now] = helpers.M_step_HQS(deconv_data,...
                                                   dtd, ...
                                                   denoise_str(color_index), ...
                                                   oI, ...
                                                   'retinex');


    wavefront2 = deconv_data.oth./(deconv_data.oto + 1e-4) .* Pupil0;
    wavefront2 = min(max(abs(wavefront2),0.8),1.2) .* ...
                                               sign(wavefront2) .* Pupil0;
    
    if abs(error_bef - error_now)/error_bef < 0.05
        learning_rate(color_index) = learning_rate(color_index) * 0.75;
        denoise_str(color_index) = denoise_str(color_index) * 2;
        if (learning_rate(color_index) < 1e-6) || (denoise_str(color_index) > 1e4)
            break;
        end
    end
    error_bef = error_now;

    toc

    if mod(epoch,1) == 0
        figure(7);
        subplot(1,2,1)
        img_spe = log(abs(wavefront1) + 1);
        mm = max(max(log(abs(wavefront1)+1)))/2;
        img_spe(img_spe>mm) = mm;
        img_spe(img_spe<0) = 0;
        img_spe = mat2gray(img_spe);
        imshow(img_spe,[])
        title('Fourier spectrum');
        drawnow;

        % Show the reconstructed amplitude
        % subplot(1,3,2)
        % imshow(((angle(oI))),[]);colorbar;
        % title(['Iteration No. = ',int2str(epoch), '  \alpha = ',num2str(tv_max)])

        subplot(1,2,2)
        imshow(((abs(oI))),[]);
        title(['Iteration No. = ',int2str(epoch), '  \alpha = ',num2str(learning_rate(color_index))])
        drawnow;
    end

end
img_rgb_final(:,:,color_index) = abs(oI).^2; 

end


img_rgb_raw = imresize(img_rgb_raw,pratio,'box');
f = img_rgb_final;

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
f_ic(:,:,2) = f(:,:,2) .* white(2) / g_sum;
f_ic(:,:,3) = f(:,:,3) .* white(3) / b_sum;

figure();imshow(f_ic,[])

% save rgb_img f_ic
imwrite(f_ic,'rgb_VEM_FPM.png');
imwrite(img_rgb_raw,'rgb_raw.png');