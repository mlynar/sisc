function snr = snr(sig, noise)
%compute snr estimate of signal and noise vectors

    L = length(sig);
    
    if length(noise) ~= L
        error('sig and noise should have the same length');
    end
    
    snr = 20 * log10(sqrt(sum(sig.^2) / L) / sqrt(sum(noise.^2) / L));