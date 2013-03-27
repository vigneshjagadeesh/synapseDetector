%% Auxillary
function [threshImg retainSizes] = cullThresh(threshImg, numRetain)
CC = bwconncomp(threshImg);
numPixels = cellfun(@numel,CC.PixelIdxList);
if( numRetain > CC.NumObjects ), retainSizes = zeros(numRetain,1); return; end
[biggest,idx] = sort(numPixels);
for iter = 1:numel(idx)-numRetain
    threshImg(CC.PixelIdxList{idx(iter)}) = 0;
end
retainSizes = biggest(numel(idx)-numRetain+1:end);