function x = rndlap(mu, b, rng)
%draw a sample from a laplace distribution over rng interval
%with mu mean and b scaling

if nargin < 1
    mu = 0;
    b = 1;
    rng = -5:0.1:5;
end

if nargin < 2
    b = 1;
    rng = -5:0.1:5;
end

if nargin < 3
    rng = -5:0.1:5;
end

l = rng(end) - rng(1);

x = rng(1) + l * rand();
y = rand();
while (exp(-abs(x - mu) / b) / 2 * b) < y
    x = rng(1) + l * rand();
    y = rand();
end
