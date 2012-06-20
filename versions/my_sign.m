function sgn = my_sign(x, eps)

if nargin < 2
    eps = 1e-6;
end

sgn = zeros(size(x));
ind = abs(x) > eps;
sgn(ind) = sign(x(ind));
