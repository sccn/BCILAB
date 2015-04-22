function C = geodesic_mean(COV,metric)

I = size(COV,3);
C = COV(:,:,1);
for i=2:I
    C = geodesic(C,COV(:,:,i),1/i,metric);
end
