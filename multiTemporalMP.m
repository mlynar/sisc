function [w r] = multiTemporalMP(y, B, noneg, maxiter, mindelta, snrTrs, spkTrs)

deadzone = 0;

nK = size(B, 3);
nTp = size(B, 2);
nCh = size(B, 1);

yLen = size(y, 2);

if nCh ~= size(y, 1)
    error('Dictionary and signal have different number of channels');
end

if ~mod(nTp, 2)
    error('Kernel lengths should be odd');
end

%if even number of channels add phony zeros
modFlag = 0;
if ~mod(nCh, 2)
    B(nCh+1, :, :) = zeros(1, nTp, nK);
    y(nCh+1,:) = zeros(1, yLen);
    nCh = nCh + 1;
    modFlag = 1;
    
    fprintf('Changing dimensionality from %d to %d\n', nCh-1, nCh);
end

%reshape dictionary
rKl = nTp * nCh;
rB = zeros(rKl, nK);

for k = 1:nK
    rB(:,k) = reshape(B(:,:,k), rKl, 1);
end

%reshape signal
rY = reshape(y, 1, yLen * nCh);

%run 1d MP skipping each nCh time points
[wr rr] = mdimTemporalMP(rY, rB, noneg, maxiter, mindelta, deadzone, ...
    snrTrs, spkTrs, nCh);

%reshape stuff back to original dimensions
w = zeros(length(y), nK);
[xp yp] = find(wr);
for i = 1:length(xp)
    xnr = xp(i) - ceil(nCh / 2) + 1;
    w((xnr - 1) / nCh + 1, yp(i)) = wr(xp(i), yp(i)); 
end

r = reshape(rr, size(y,1), size(y,2));
if modFlag
    r = r(1:end-1,:);
end
