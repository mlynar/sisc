function [f g gc] = phi_objective(phi, x, s)

L = length(x);
B = size(s,1);
N = length(phi) / B;

phi = reshape(phi, N, B);

hL = floor(N / 2);

EI = zeros(L,1);
for i = 1:size(phi,2)
    EI = EI + conv(s(i,:), phi(:,i), 'same')';
end

E = x' - EI;

f = 0.5 * sum(E.^2);

dphi = zeros(size(phi));
for t = 1:L
    srt = t - hL; 
    fin = t + hL; 
    
    tind = srt:fin; 
    ind = 1:2*hL+1; ind = ind(tind > 0 & tind <= L);
    tind = tind(tind > 0 & tind <= L);
    
    dphi(ind,:) = dphi(ind,:) - s(:,t) * E(tind);
end

dphic = zeros(size(phi));
for t = -hL:hL
    str = max(1, t);
    fin = min(L, L+t);
    
    tind = str:fin;
    
    sind = tind - t;
    sind = sind(sind >= 1 & sind <= L);
    
    size(s(:,sind))
    size(E(tind))
    
    dphic(t+hL+1,:) = dphic(t+hL+1,:) - s(:,sind) * E(tind);
end

g = dphi(:);
gc = dphic(:);