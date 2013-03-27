clear; close all; clc;
% Detect Ribbons ...
if( ismac )
addpath(genpath('/Users/vigneshjagadeesh/Dropbox/synapseProject/interfaces/'));
dataDir = '/research/datasets/synapseDataset/synapseDataset/ribbon/ribSynS1/';
else
    addpath(genpath('C:\Users\vignesh\Dropbox\synapseProject\interfaces\'));
    dataDir = 'C:\researchCode\dataBase\synapseDataset\synapseDataset\ribbon\ribSynS1\';
end
ribFiles = dir( dataDir );

%I = imread('sampleRibbon.png');
for fileIter = 3:10
    I = imread( [dataDir ribFiles(fileIter).name] );
    I = I(:,:,1); Iorig = I;
    I = imfilter(I, fspecial('gauss', [49 49], 7) );    
    [F FFreq] = makeRFSfilters();
    F = cat(3, F(:,:,13:18), F(:,:,31:36) );
    % FFreq = FFreq(:,:,19:24);
    % FR = vrl_imfilter(I, FFreq);
    masker = zeros( size(I,1), size(I,2) ); masker(40:end-40, 40:end-40 ) = 1;
    for iter =1:size(F,3), 
        FR(:,:,iter) = imfilter(double(I), F(:,:,iter), 'same', 'conv') .* double( masker );
    end    
    % This is the initial formulation and needs to be changed ... 
    FR1 = reshape( max( reshape( FR(:,:,1:6), size(I,1)*size(I,2), [] ), [], 2 ), [size(I,1) size(I,2)] );
    FR2 = reshape( max( reshape( FR(:,:,7:12), size(I,1)*size(I,2), [] ), [], 2 ), [size(I,1) size(I,2)] );    
    figure(1); imshow(Iorig); hold on;
    % contour( FR1 > 10, [0 0], 'r' );
    threshImg = FR2 > 3;
    CC = bwconncomp(threshImg);
    numPixels = cellfun(@numel,CC.PixelIdxList);
    [biggest,idx] = sort(numPixels);
    for iter = 1:numel(idx)-10
        threshImg(CC.PixelIdxList{idx(iter)}) = 0;
    end    
    contour( threshImg, [0 0], 'g' );
    hold off;
end

% Alg 2
% function cellWallDetect()
clear; close all; clc;
%load emcAnno126;
if( ispc )
    emcPath = 'C:\researchCode\dataBase\emChallenge\514_ipl_00_219\';
    addpath( genpath( 'C:\Users\vignesh\Dropbox\synapseDetect\interfaces\' ) );
else
    emcPath = '/Users/vigneshjagadeesh/Dropbox/synapseDetect/';
    addpath( genpath( '/Users/vigneshjagadeesh/Dropbox/synapseDetect/interfaces/' ) );
end

load emcSynAnno126;
load emcNonAnno126;
FILTERINFO.mapping=getmapping(8,'u2');

I = imread([emcPath '514_ipl_01_126.tif']);
I = medfilt2( I, [11 11] );
% figure; imshow(I);

for synIter = 1:14
    I = synPatch{synIter};
    
    Ismooth = imfilter( double(I), fspecial('gauss', [15 15], 3) , 'conv', 'same'  );
    [Ix Iy] = gradient(Ismooth);
    [Ixx Ixy] = gradient(Ix);
    [Ixy Iyy] = gradient(Iy);
    gradMag = sqrt(Ix.^2 + Iy.^2);
    sdMag   = sqrt(Ixx.^2 + Iyy.^2);
    edgeInd = 1./(1+exp(-gradMag/3));
    sdInd   = 1./(1+exp(-sdMag/10) );
    figure(1); imagesc( edgeInd );
    figure(2); imagesc(sdInd); 
    
    
    BW = edgeInd > .8;
    BW(1,:) = 0; BW(end, :) = 0; BW(:,1) = 0; BW(:, end) = 0;
    BW = imopen(BW, strel('disk', 5) );
    figure; imshow( BW );
    CC = bwconncomp(BW);
    numPixels = cellfun(@numel,CC.PixelIdxList);
    [biggest,idx] = sort(numPixels, 'descend');
    
    SNODE = 200;
    if( numel(idx) > SNODE )
        for iter = SNODE:numel(idx)
            BW(CC.PixelIdxList{idx(iter)}) = 0;
        end
    end
    figure; imshow(BW);
    figure; imshow( uint8(imdilate(BW, strel('disk', 3))) .* I );
end
