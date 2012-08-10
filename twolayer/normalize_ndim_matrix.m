function A = normalize_ndim_matrix(W)
    A = W;
    nx = size(A,1);
    ny = size(A,2);
    
    for i = 1:size(W, 3)
        if abs(sum(A(:,:,i))) > 0
            tmp = reshape(A(:,:,i), 1, nx * ny);
            tmp = tmp ./ norm(tmp, 2);
            A(:,:,i) = reshape(tmp, nx, ny);
        end
    end