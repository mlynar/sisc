smp_hz = 6000;
fcoefs = MakeERBFilters(smp_hz, 32, 200);
gammaTones = ERBFilterBank([1 zeros(1,256)], fcoefs)';
gammaTones = normalize_matrix(gammaTones);
freqs = ERBSpace(200, smp_hz/2, 32);
%%
clc;
x = zeros(1, 2000);
x(100:100+256) = gammaTones(:,30) * 3.2;

%%
x = zeros(1,2000);
for i = 1:50
    t = randi(1700);
    x(t:t+256) = x(t:t+256)' + gammaTones(:, 3 + 29) * rand() * 3;
end
plot(x)
%%
lambda = 0.01;
x0 = zeros(13, length(x));
[xmin fx] = minimize(x0(:), 's_objective', 10, x, gammaTones(:,20:32), lambda);
xmin = reshape(xmin, 13, length(x));

%%
clc
close all;
niter = 20;
fhist = zeros(1, niter);
s = zeros(13, length(x));
lr = 0.005;
for i = 1:niter
    [f g] = s_objective(s, x, gammaTones(:,20:32), lambda);
    fprintf('%d - value %.4f\n', i, f);
    fhist(i) = f;
    
    
    g = reshape(g, 13, length(x));
    s = s - lr * g;
end

figure;
plot(fhist, 'LineWidth', 2); 
hold on; plot(fx, 'r', 'LineWidth',2);
grid()

%%
q = 0.9;
xt = xmin;

for i=1:size(xt,1)
    xt(i, abs(xt(i,:)) < quantile(abs(xt(i,:)), q)) = 0;
end

s = reconstructSignal(xmin, gammaTones(:,20:32));
sr = reconstructSignal(xt, gammaTones(:,20:32));

fprintf('SNR: %.3f\n', snr(x, x-s));
fprintf('SNR PRUNNED: %.3f\n', snr(x, x-sr));

figure;
subplot(2,1,1)
imagesc(xt); colorbar()
subplot(2,1,2)
plot(x); hold on;
plot(s, 'r'); plot(sr, 'k');
legend('signal', 'estimate', 'est prunned');

