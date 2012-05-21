%Convolution of a vector and a matrix
%uses the same algo as Bruno Luong's convnfft
function R = myconvnfft(M,v)
    nsize = size(M,2)+length(v)-1;
    vf = fft(v,nsize);
    Mf = fft(M,nsize,2);
    R  = ifft(bsxfun(@times,vf,Mf),[],2);
end