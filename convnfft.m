%Convolution of two vectors in the Fourier domain
function R = convnfft(u,v)
    nsize = length(u) + length(v) - 1;
    vf = fft(v, nsize);
    uf = fft(u, nsize);
    R  = ifft(vf .* uf);
end
