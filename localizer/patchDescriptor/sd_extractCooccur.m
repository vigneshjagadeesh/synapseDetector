function [cleftFeat1 ribbonFeat1 logFeat1] = sd_extractCooccur(I,F,fil,VIS, bboxes)

if( nargin < 1 )
clear; close all; clc;
%% Test Co-Occur
addpath( genpath( '../../../../interfaces/' ) );
[F FFreq] = makeRFSfilters(); F = F(:,:,31:36);
fil = fspecial('log', 49, 7);
dataDir = '/media/OS/researchCode/dataBase/synapseDataset/synapseDataset/ribbon/ribTrain/';
allFiles = dir( [dataDir '/*.png'] );

for fileIter = 1:numel(allFiles)
I = imread( [dataDir allFiles(fileIter).name] );
I = double( I(:,:,1) );
[cleftFeat1, ribbonFeat1, logFeat1] = extractCooccur(I,F,fil);
pause(.01);
end
end
%%% Keep Updating this function in the interfaces folder, so that the
%%% latest version is used for benchmarking ...

c1 = I; 
[border1 mask1 ] = genBorderMask(c1 );
masker = zeros( size(I,1), size(I,2) );
masker(40:end-40, 40:end-40 ) = 1;

%% Cleft Features ... threshold=7 with equalization and 3 without
% S1: Adaptive Histogram Equalization
% S2: Smoothen Image 
% S3: Filter with Steered Filter Bank
% S4: Threshold with 7
% S5: Retain only top 20 Connected Components
% S6: Compute cleft features for entire image
% S7: Compute cleft features for pyramidal patches
origImage = I; 
I = double( adapthisteq( uint8(I) ) );                                      %S1
I = imfilter(I, fspecial('gauss', [49 49], 7) );                            %S2
for iter =1:size(F,3),  
    FR(:,:,iter) = imfilter(double(I), F(:,:,iter), 'same', 'conv') .* double( masker ); 
end
FR1 = reshape( max( reshape( FR(:,:,1:6), size(I,1)*size(I,2), [] ), [], 2 ), [size(I,1) size(I,2)] );
threshImg = FR1 > 7;                                                        %S4
threshImg = cullThresh(threshImg, 20) .* FR1;                               %S5
%cleftFeat1 = sum( threshImg(:) )/numel(threshImg);                         %S6
cleftFeat1 = sum( threshImg(:) );
for scaleIter = 1:4                                                         %S7
    cleftFeat1( scaleIter + 1 ) = sum( sum( mask1(:,:,scaleIter) .* threshImg ) )/sum( sum( mask1(:,:,scaleIter) ) );
    %cleftFeat1( scaleIter + 1 ) = sum( sum( mask1(:,:,scaleIter) .* threshImg ) );
end
cleftBinary = (1-imdilate(threshImg>0, strel('disk',3)));

%% Vesicle Features Heuristic threshold for current image size set to 1.5 top 500 retentions
% S1: Perform Median Filtering - Kernel 5 * 5
% S2: Adaptive Histogram Equalization
% S3: Perform Laplacian Gaussian Filtering
% S4: Threshold Weak Responses (set to 1.5)
% S5: Detect Peaks in 2D Map
% S6: Remove Border Artifacts
% S7: Retain top 500 responses by ranking them
% S8: Compute full Image Feature and Dispersion
% S9: Compute Pyramid Features of Dispersion

c1 = c1;%medfilt2( double(c1), [5 5] );                                         %S1
cEq = double( adapthisteq( uint8(c1) ) );                                   %S2
log1 = imfilter( cEq, fil, 'same', 'conv');                                 %S3
log1( log1 < 1.5) = 0;   %was 1.5                                                   %S4
logMax1 = log1 > imdilate(log1, [1 1 1; 1 0 1; 1 1 1]);                     %S5
logMax1 = logMax1 .* border1 .* cleftBinary;                                               %S6
logMax1 = retainTop( log1, logMax1, 500);                                   %S7
fl_detectBox(origImage, log1.*logMax1, 500, 3, VIS);
[yF xF] = find( logMax1~=0 );                                               %S8
logFeat1(1) = sum( logMax1(:) .* log1(:) .* cleftBinary );
logDisp1(1) = median( pdist([xF yF],'euclidean') );
for scaleIter = 1:4                                                         %S9
    currMask = logMax1 .* mask1(:,:,scaleIter);
    [y x] = find( currMask~= 0 );
    if( ~isempty(x) )
        logFeat1( scaleIter + 1 ) = sum( sum( currMask .* log1 ) );
        logDisp1( scaleIter + 1 ) = median( pdist([x y],'euclidean') );
    else
        logFeat1( scaleIter + 1 ) = 0;
        logDisp1( scaleIter + 1 ) = 0;
    end
end
logFeat1 = [logFeat1(:); logDisp1(:)];

%% Ribbon Features
% S1: Adaptive Histogram Equalization
% S2: Thrshold with high threshold factor
% S3: Open the Image with Disk of Size 3
% S4: Remove all small connected components

Ieq = adapthisteq( uint8(c1) );                                             %S1
threshImg = im2bw(255-Ieq, 0.9);                                            %S2
threshImg = imopen(threshImg, strel('disk',3));                             %S3
[threshImg ribbonFeat1] = cullThresh(threshImg, 3);                         %S4
ribbonFeat1 = ribbonFeat1 / numel(threshImg); 
ribbonBinary = threshImg;

if(VIS), 
    %figure(1); imshow(uint8(c1)); 
    hold on; contour( cleftBinary, [0 0], 'c' ); contour( ribbonBinary, 'g'); 
    plot(xF,yF,'m*'); hold off;
    title(['LoG Strength' num2str(logFeat1(1)) 'LoG Dispersion' num2str(logDisp1(1)) ]);
end
cleftFeat1 = cleftFeat1(:);
ribbonFeat1 = ribbonFeat1(:);
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

function binMaskTop = retainTop( response, binMask, NUM)
binMaskTop = zeros( size(binMask,1), size(binMask,2) );
[y x]= find(binMask ~= 0 );
allVals = response( binMask ~= 0 );
[sortedRes sortedIds] = sort( allVals, 'descend' );
NUM = min( NUM, numel(y) );
for iter = 1:NUM
    binMaskTop( y( sortedIds(iter) ), x( sortedIds(iter) ) ) = 1;
end

function [border mask ] = genBorderMask(c )

[noRows noCols noSlices] = size(c);

% Creating Borders Alone from the Given Input Image
border = zeros( noRows, noCols ); 
rowSplitter = ceil( noRows/2);
colSplitter = ceil( noCols/2);
border(3:end-3, 3:end-3) = 1; 

mask = zeros( noRows, noCols);
mask1 = mask; mask1(1:rowSplitter, 1:colSplitter) = 1;
mask2 = mask; mask2(rowSplitter+1:end, 1:colSplitter) = 1;
mask3 = mask; mask3(1:rowSplitter, colSplitter+1:end) = 1;
mask4 = mask; mask4(rowSplitter+1:end, colSplitter+1:end) = 1;

mask = cat(3, mask1, mask2, mask3, mask4);

