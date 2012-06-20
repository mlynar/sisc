function A = normalize_matrix(W)
    A = W;
    for i = 1:size(A, 2)
        if abs(sum(A(:,i))) > 0
            A(:,i) = A(:,i) ./ norm(A(:,i), 2);
        end
    end
    