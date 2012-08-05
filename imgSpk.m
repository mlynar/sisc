function [lim] = imgSpk(w, lm, fqs, cmap)
%if highest frequency component is the first column of phi,
%then it's plotted as the topmost row of spikegram

if nargin < 3
    fqs = 1:size(w,2);
end
if nargin < 4
    cmap = french(512);
end

[fqs ind] = sort(fqs, 'ascend');
w = w(:,ind);

l = ceil(0.1 * length(w))
rg = (-l:l)';
v = 100;

img = conv2(w,exp(-rg.^2/2/v^2),'same')';

if nargin < 2
    lim = 0.8 * max(abs(img(:)));
else
    lim = lm;
end


imagesc(img, [-lim lim]); colormap(cmap);
set(gca, 'YTick', 1:floor(0.2*length(fqs)):length(fqs));
set(gca, 'YTickLabel', fqs(end:-floor(0.2*length(fqs)):1));

ylabel('kernel freq [Hz]');
xlabel('time [sample]');

