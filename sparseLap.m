function s = sparseLap(x)
%sparse penalty d/dx exp(-|x|/b)
%b = 1
    s = sign(x);