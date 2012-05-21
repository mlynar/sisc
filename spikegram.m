function spikegram(W, phi, scm, wMax)
%W - n * t (n - basis functions, t - time points)
%phi - basis functions in columns
%scm - int flag - if 0 - plot spikes at position of maximum or if 1 begining
%                of kernel otherwise or if 2 middle of the kernel
%TODO: plot center frequencies, add a flag to plot kernel itself - line 66
%TODO: plot kernels ordered by frequency

if nargin < 3
    scm = 0;
    wMax = max(abs(W(:)));
elseif nargin == 3
    wMax = max(abs(W(:)));
end

wMax

if scm > 2 || scm < 0
    error('scm can be only {0,1,2}. Type help spikegram');
end

nB = size(phi,2);
bL = size(phi,1);
hL = floor(bL/2); 
if ~mod(bL, 2)
    hL = hL+1;
end

offset = zeros(1, nB);
for i = 1:nB
    if scm == 0
        offset(i) = hL - com(phi(:,i)); %TODO: inverted or not
    elseif scm == 1
        offset(i) = hL;
    else
        offset(i) = 0;
    end
end

nRows = size(W,1); 
if nRows ~= nB
    error('w and phi should share a dimension');
end
nCols = size(W,2);

circleSize = 0.01;

%figure;
plot(0,0,'w.');
hold on

dX = 1 / nCols;
dY = 1 / nRows;

usdKernels = [];

yLevel = 0;
for i = 1:nRows
    yLevel = yLevel + dY;
    for j = 1:nCols;
        if W(i,j) ~= 0
            sc = abs(W(i,j)) / wMax;
            
            if W(i,j) > 0
                clr = sc * [1 0 0];
            else
                clr = sc * [0 0 1];
            end
            
            xP = (j - offset(i)) * dX;
            
            %TODO: add a flag to plot kernel
            if false
                plot(dX * (j-hL:j+hL), yLevel + 0.2 * phi(:,i));
            end
            
            circle(xP, yLevel, sc * circleSize, clr);
            usdKernels = [usdKernels i];
        end
    end
end
usdKernels = unique(usdKernels);

set(gca, 'XTick', 0:0.2:1); set(gca, 'XTickLabel', (0:0.2:1) / dX);
set(gca, 'YTick', usdKernels * dY); set(gca, 'YTickLabel', usdKernels);
tStr = sprintf('Spikegram - spike no: %d', nnz(W(:)));
title(tStr);
ylabel('kernel number');
xlabel('time');

xlim([0 1+2*circleSize])
ylim([0 1+2*circleSize])