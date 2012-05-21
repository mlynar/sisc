%the same as myconvnfft but in temporal domain
function convM = myconv(m, v)

nm = size(m,1);
convM = zeros(nm, size(v,2));
for i = 1:nm
    convM(i,:) = conv(v, m(i,:), 'same');
end
