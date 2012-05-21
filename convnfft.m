%Convolution of two vectors
%uses the same algo as Bruno Luong's convnfft
function R = convnfft(u,v)
    nsize = length(u) + length(v) - 1;
    vf = fft(v, nsize);
    uf = fft(u, nsize);
    R  = ifft(vf .* uf);
end
