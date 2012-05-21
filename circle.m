function circle(x,y,r,c)
%x and y are the coordinates of the center of the circle
%r is the radius of the circle
%0.01 is the angle step, bigger values will draw the circle faster but
%you might notice imperfections (not very smooth)
if nargin < 4
    c = 'k';
end

ang=0:0.01:2*pi; 
xp=r*cos(ang);
yp=r*sin(ang);
fill(x+xp,y+yp,c, 'EdgeColor', 'none');
end