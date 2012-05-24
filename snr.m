function snr = snr(sig, noise)
%compute snr estimate of signal and noise vectors

    L = length(sig);
    
    if length(noise) ~= L
        error('sig and noise should have the same length');
    end
    
    snr = 10 * log10(sum(sig.^2) / sum(noise.^2));