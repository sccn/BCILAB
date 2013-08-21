function y = MRImask(n,beams)
%
% produces the fan MRI mask, of size n*n,
% beams is the number of angles 
%
m = ceil(sqrt(2)*n);
aux = zeros(m,m); ima = aux;
aux(round(m/2+1),:) = 1;
angle = 180/beams;
angles = (0:angle:180-angle);
for a = 1:length(angles)
    ang = angles(a);
    temp = imrotate(aux,ang,'crop');
    ima = ima + temp;
end
ima = ima(round(m/2+1) - n/2:round(m/2+1) + n/2-1,...
          round(m/2+1) - n/2:round(m/2+1) + n/2-1);

y = (ima > 0);
