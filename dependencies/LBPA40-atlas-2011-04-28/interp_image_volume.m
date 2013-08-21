function signal = interp_image_volume(V,points,d)
% Uses the spm_bsplinc and spm_bsplins SPM routines to interp a set of
% points in a image.

% V = spm_vol('/home/alejandro/Documents/MATLAB/testing/lpba40.spm5.avg152T1.maxprob.nii');

x = points(:,1);
y = points(:,2);
z = points(:,3);

c = spm_bsplinc(V,d);
signal = spm_bsplins(c,x,y,z,d);