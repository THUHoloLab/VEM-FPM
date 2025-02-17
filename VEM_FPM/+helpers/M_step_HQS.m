function [out,o,err] = M_step_HQS(dec_o,dtd,de_str,oR,type)
    o_old = oR;

    v = denoiser(oR);

    % solving g-subproblem
    switch type
        case 'normal'
            fenzi = dec_o.hto .* (1) + fftshift(fft2(de_str * v));
            fenmu = dec_o.hth .* (1) + de_str + 0.00001;
        case 'retinex'
            fenzi = dec_o.hto .* dtd + fftshift(fft2(de_str * v));
            fenmu = dec_o.hth .* dtd + de_str;
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
    
    % solving admm variables

    N = size(o,1) * size(o,2);
    err = (1/sqrt(N))*(sqrt(sum(sum((o - o_old).^2))));
end

function o = denoiser(o)
    o = medfilt2(real(o)) + 1i * medfilt2(imag(o));
end