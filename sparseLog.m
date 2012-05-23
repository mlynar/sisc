function s = sparseLog(x)
%sparse penalty d/dx log(1 + x^2)
    s = 2 * x ./ (1 + x.^2);