function [out,o,v,u,err] = M_step(dec_o,v,u,dtd,denoise_str,oR,type)
    o_old = oR;
    v_old = v;
    u_old = u;

    % solving g-subproblem
    switch type
        case 'normal'
            fenzi = dec_o.hto .* (1) + fftshift(fft2(denoise_str * (v - u)));
            fenmu = dec_o.hth .* (1) + denoise_str + 0.00001;
        case 'retinex'
            fenzi = dec_o.hto .* dtd + fftshift(fft2(denoise_str * (v - u) + 0.01*oR));
            fenmu = dec_o.hth .* dtd + denoise_str + 0.01;
        case 'none'
            fenzi = dec_o.hto .* dtd + fftshift(fft2(0.01*oR));
            fenmu = dec_o.hth .* dtd + 0.01;
        otherwise
            error('type must be normal or retinex')
    end
    out = fenzi ./ fenmu;
    o = ifft2(ifftshift(out));

    % solving g-subproblem
    % v = complex_TV(o + u,0.001/tv_max,'isotropic');
    v = medfilt2(real(o + u),[3,3]) + 1i * medfilt2(imag(o + u),[3,3]);
    % solving admm variables
    u = u + (o - v);


    N = size(o,1) * size(o,2);
    err_o = (1/sqrt(N))*(sqrt(sum(sum((o - o_old).^2))));
    err_v = (1/sqrt(N))*(sqrt(sum(sum((v - v_old).^2))));
    err_u = (1/sqrt(N))*(sqrt(sum(sum((u - u_old).^2))));
    err = err_o + err_v + err_u;

end
