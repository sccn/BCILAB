function [patches, mean_patch, V] = preprocess(patches)

[patches, mean_patch] = removeDC(patches,2);
[patches,V,E,D] = zca2(patches);
