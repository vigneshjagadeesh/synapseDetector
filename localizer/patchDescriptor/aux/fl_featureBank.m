function [featVal featName] = f_featureBank(I,F)

addpath(genpath('/home/vignesh/Dropbox/scratchSpace/synapseProject/interfaces/featureComputation/imrender/'));

if( nargin < 1 )
    display('I need an input argument to begin ... processing a dummy image');
    clear; close all; clc;
    I = imread('cameraman.tif');
    ctr =1;
    for filtSup = [49]
        F{ctr} = makeSubFilters(filtSup);
        ctr = ctr + 1;
    end
end

[noRows noCols noSlices] = size(I);

% Threshold Features
display('Computing Threshold Features');
thFeat = cell(1);
fiFeat = cell(1);
suFeat = cell(1);
ctr = 1;
for thIter = .2:.2:.8
    Ith = im2bw( I, thIter );
    CC = bwconncomp(Ith);
    numObjects = CC.NumObjects/500;
    numPixels = cellfun(@numel,CC.PixelIdxList);
    [desc, idx] = sort(numPixels, 'descend');
    topSize = 0;
    for iter = 1:min(numel(numPixels), 3)
        topSize = topSize + numel( CC.PixelIdxList{idx(iter)} );
    end
    topSize = topSize / (3* noRows * noCols);
    thFeat{1} = [thFeat{1}; topSize; numObjects];
    thName{ctr} = ['topTh' num2str(thIter)]; ctr = ctr + 1;
    thName{ctr} = ['topTh' num2str(thIter)]; ctr = ctr + 1;
end

% Filter Features
display('Computing Filter Features');
ctr = 1;
for filtSup = 1:numel(F)
    for filtIter = 1:size(F{filtSup},3)
        FR = imfilter( double(I), F{filtSup}(:,:,filtIter), 'same', 'conv' );
        fiFeat{1} = [fiFeat{1}; mean(FR(:)); std(FR(:))];
        fiName{ctr} = ['meanFilt ' num2str(filtSup) ' ' num2str(filtIter)]; ctr = ctr + 1;
        fiName{ctr} = ['stdFilt ' num2str(filtSup) ' ' num2str(filtIter)]; ctr = ctr + 1;
    end
end

% Superpixel Features
display('Computing Superpixel Features');
msFeat = [5 3 8 5 13 11];
ctr = 1;
%for iter = 1:3
%    S = vgg_segment_ms(uint8(I), msFeat(2*iter-1), msFeat(2*iter), 300);
%    sFeat = numel( unique(S) )/1000;
%    suFeat{1} = [suFeat{1}; sFeat];
%    suName{ctr} = ['suMS ' num2str(msFeat(2*iter-1)) ' ' num2str(msFeat(2*iter))]; ctr = ctr + 1;
%end
gbFeat = [.8 200 1.1 400 1.1 600];
for iter = 1:3
    S = vgg_segment_gb(uint8(I), gbFeat(2*iter-1), gbFeat(2*iter), 500);
    sFeat = numel( unique(S) )/1000;
    suFeat{1} = [suFeat{1}; sFeat];
    suName{ctr} = ['suGB ' num2str(gbFeat(2*iter-1)) ' ' num2str(gbFeat(2*iter))]; ctr = ctr + 1;
end

featName = [thName fiName suName];
featVal = [thFeat{1}; fiFeat{1}; suFeat{1}];