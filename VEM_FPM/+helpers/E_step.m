function [loss,deconv_data] = E_step(wavefront1, ...
                                     wavefront2, ...
                                     deconv_data,...
                                     b_ledpos, ...
                                     dY_obs, ...
                                     pratio, ...
                                     learning_rate)

loss = 0;
sub_wavefront = dY_obs * 0;

% forward inference
for data_con = 1:size(dY_obs,3)
    kt = b_ledpos(data_con,3);
    kb = b_ledpos(data_con,4);
    kl = b_ledpos(data_con,1);
    kr = b_ledpos(data_con,2);
    sub_wavefront(:,:,data_con) = wavefront1(kt:kb,kl:kr);
end
clear wavefront1;

%calculate latent variable
latent_z = (misc.ifft2_ware(bsxfun(@times,sub_wavefront,wavefront2),true) / pratio^2);
[loss,dm] = misc.ret_loss(abs(latent_z) - dY_obs,'isotropic',1);

latent_z = bsxfun(@times,...
                        abs(latent_z) - learning_rate .* dm,...
                        sign(latent_z));

clear dm;

latent_z = misc.fft2_ware(latent_z, true) .* pratio^2;


for data_con = 1:size(dY_obs,3)
    kt = b_ledpos(data_con,3);
    kb = b_ledpos(data_con,4);
    kl = b_ledpos(data_con,1);
    kr = b_ledpos(data_con,2);

    deconv_data.hto(kt:kb,kl:kr) = deconv_data.hto(kt:kb,kl:kr) + ...
                                latent_z(:,:,data_con) .* conj(wavefront2);

    deconv_data.hth(kt:kb,kl:kr) = deconv_data.hth(kt:kb,kl:kr) +...
                                                        abs(wavefront2).^2;
end

deconv_data.oth = deconv_data.oth + sum(latent_z .* conj(sub_wavefront),3);
deconv_data.oto = deconv_data.oto + sum(abs(sub_wavefront).^2,3);

end