function circle(x, y, r, c)
%fill circle
%	x,y - center coordinates
%	r - radius
%	c - color, black by default

if nargin < 4
    c = 'k';
end

edge = 0:0.01:2*pi; 
xp = r * cos(edge);
yp = r * sin(edge);

fill(x + xp, y + yp, c, 'EdgeColor', 'none');


