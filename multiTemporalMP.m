function [w r] = multiTemporalMP(y, B, noneg, maxiter, mindelta, snrTrs, spkTrs)

deadzone = 0;

nK = size(B, 3);
nTp = size(B, 2);
nCh = size(B, 1);

if nCh ~= size(y, 1)
    error('Dictionary and signal have different number of channels');
end

rKl = nTp * nCh;
rB = zeros(rKl, nK);

for k = 1:nK
    rB(:,k) = reshape(B(:,:,k), rKl, 1);
end


yLen = size(y, 2);
rY = reshape(y, 1, yLen * nCh);

[wr rr] = mdimTemporalMP(rY, rB, noneg, maxiter, mindelta, deadzone, ...
    snrTrs, spkTrs, nCh);


%=====================TO BE DELETED=============================
figure;
subplot(4,1,1);
wrl = zeros(size(wr)); wrl(wr~=0) = 1;
imagesc(wrl'); title('wr from multitemporalmp');
subplot(4,1,2);
imagesc(rr'); title('rr from multitemporalmp');
subplot(4,1,3)
imagesc(rY); title('original ry')
subplot(4,1,4)
wl=wr';
rc = zeros(size(rr));
for i = 1:nK
    rc = rc + conv(wl(i,:), rB(:,i), 'same')';
end
imagesc(rc'); title('reconstruction from wr');
%==============================================================

w = zeros(length(y), nK);
[xp yp] = find(wr);
for i = 1:length(xp)
    xnr = xp(i) - ceil(nCh / 2) + 1;
    w((xnr - 1) / nCh + 1, yp(i)) = wr(xp(i), yp(i)); 
end
r = reshape(rr, size(y));