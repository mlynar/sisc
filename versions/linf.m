phi = rand(101,1);
ntrials = 1000;

B = size(phi,2);
N = size(phi,1);

lambda = 0.01;

lr = 0.0001;

%x = sampleToyData(gammaTones(:,30),  5000, 2, 0.2);
%plot(x)

L = length(x);

for n = 15:ntrials
    
    %inference
    s0 = zeros(B, L);
    [s1 fs] = minimize(s0(:), 's_objective', 10, x, phi, lambda);
    s1 = reshape(s1, B, L);
    
    %reconstruction
    EI = reconstructSignal(s1, phi);
    
    %learning
%     [obj0,g] = phi_objective(phi(:), x, s1);
%     dphi = reshape(g, N, B);
% 
%     phi1 = phi - lr * dphi;
% 
%     [obj1,g] = phi_objective(phi1(:), x, s1);
% 
%     if obj1 > obj0
%         fprintf('objective function increased\n');
%         lr = lr * 0.9;
%     else
%         phi = phi1;
%     end

    [phi fp fpc] = minimize(phi(:), 'phi_objective', 10, x, s1);

    phi = normalize_matrix(phi);
    [obj2 g] = phi_objective(phi(:), x, s1);
    
    fprintf('Iteration %d error: %.6f lr: %.5f \n',n, obj2, lr);
    
    if 1
        clf
        figure(4);
        subplot(2,1,1)
        imagesc(s1); colorbar();
        subplot(2,1,2)
        plot(phi);
        drawnow;
    end
end
    

    
    