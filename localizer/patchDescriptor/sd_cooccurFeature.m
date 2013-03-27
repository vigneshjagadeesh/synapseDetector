function [fullFeat] = sd_cooccurFeature(I, F, VIS, bboxes)

addpath('patchDescriptor');
if( nargin < 1 )
    clear; close all; clc;
    %% Test Co-Occur
    [F FFreq] = makeRFSfilters(); F1 = F(:,:,31:36);
    fil = fspecial('log', 49, 7);
    dataDir = '/cluster/home/retinamap/synapseProject/synapseDataset/ribbon/ribSynS1/';
    allFiles = dir( [dataDir '/*.png'] );
    F = cat(3, F1,fil);
    bboxes = [1 1 1024 1024];
    VIS = 1;
    
    I = imread( [dataDir allFiles(10).name] );
    I = double( I(:,:,1) );
    
end

if( nargin == 3 )
    bboxes = [1 1 size(I,2) size(I,1)];
end

% Set default parameters for thresholds
cleftChannel = 1:6;
vesicleChannel = 7;
numFilters = size( F, 3 );

% Tunable Parameters
cleftTHRESHOLD = 0.006; %7
cleftRETAINCC  = 20;
vesicleTHRESHOLD = 0.008;%1.5
vesicleRETAINMAX = 500;


% Preprocess and Filter the images

Iorig = uint8(I);
I = double( adapthisteq( uint8(I) ) );
Ismooth = imfilter(I, fspecial('gauss', [49 49], 7) );
% filterCornerMask = zeros( size(I,1), size(I,2) );
% filterCornerMask(40:end-40, 40:end-40) = 1;
% filterCornerMask = repmat( filterCornerMask, [1 1 7] );
% filtResponseCleft = vrl_imfilter(double(Ismooth), F(:,:,cleftChannel));
% filtResponseVesicle = vrl_imfilter(double(I), F(:,:,vesicleChannel));

FR = aniFiltering(Ismooth, I);

%% Get Cleft Response and Cleft Mask
cleftResponse = reshape( max( reshape( FR(:,:,cleftChannel), size(I,1)*size(I,2), [] ), [], 2 ), [size(I,1) size(I,2)] );
cleftMaxima   = fl_cullThresh( cleftResponse > cleftTHRESHOLD, cleftRETAINCC );
cleftBinary = (1-imdilate(cleftMaxima>0, strel('disk',3)));

%% Get Vescile Response and Vesicle Mask
vesicleResponse = FR(:,:,vesicleChannel);
vesicleResponse( vesicleResponse < vesicleTHRESHOLD ) = 0;
vesicleMaxima = vesicleResponse > imdilate( vesicleResponse, [1 1 1; 1 0 1; 1 1 1] );
vesicleMaxima = vesicleMaxima .* cleftBinary;
vesicleMaxima = fl_retainTop( vesicleResponse, vesicleMaxima, vesicleRETAINMAX);
[yVes xVes] = find(vesicleMaxima > 0);

if(VIS),
    
    %% Get Ribbon Response
    ribbonResponse = im2bw(255-uint8(I), 0.9);
    ribbonResponse = imopen(ribbonResponse, strel('disk',3));
    [ribbonResponse ribbonFeat] = fl_cullThresh(ribbonResponse, 3);
    ribbonFeat = ribbonFeat / numel(ribbonResponse);
    ribbonBinary = ribbonResponse;
    
    fl_detectBox(I, vesicleResponse.*vesicleMaxima, ceil(size(I,1)/10), 3, VIS);
    hold on; contour( cleftBinary, [0 0], 'c' ); contour( ribbonBinary, 'g');
    plot(xVes,yVes,'m*'); hold off;
end
%% Extract Features only on interest points
fullFeat = zeros( size(bboxes,1), 18 );
for bIter = 1:size(bboxes,1)
    % Cleft and Vesicle Features
    x1 = bboxes( bIter, 1 );
    y1 = bboxes( bIter, 2 );
    x2 = bboxes( bIter, 3 );
    y2 = bboxes( bIter, 4 );
    croppedImage = I( y1:y2, x1:x2, : );
    currCleftMaxima = cleftMaxima( y1:y2, x1:x2 );
    currVesicleMaxima = vesicleMaxima( y1:y2, x1:x2 );
    currCleftResponse = cleftResponse( y1:y2, x1:x2 );
    currVesicleResponse = vesicleResponse( y1:y2, x1:x2 );
    [foo spatialMasks ] = fl_genBorderMask( croppedImage );
    
    % Ribbon Feature
    ribbonResponse = im2bw(255-uint8(croppedImage), 0.9);
    ribbonResponse = imopen(ribbonResponse, strel('disk',3));
    [ribbonResponse ribbonFeat] = fl_cullThresh(ribbonResponse, 3);
    ribbonFeat = ribbonFeat / numel(ribbonResponse);
    
    % Cleft and Vesicle Features
    for maskIter = 1:size( spatialMasks, 3)
        currCleftMask = currCleftMaxima .* currCleftResponse .* spatialMasks(:,:,maskIter);
        cleftFeat(maskIter) = sum( currCleftMask(:) ); % changed this from sum
        
        currVesicleMask = currVesicleMaxima .* currVesicleResponse .* spatialMasks(:,:,maskIter);
        [y x] = find( currVesicleMask ~= 0 );
        logFeat(maskIter) = sum( currVesicleMask(:) );
        if( numel(x) > 1 )
            logFeat(maskIter + size(spatialMasks,3) ) = median( pdist([x y],'euclidean') )/size(I,1);
        else
            logFeat(maskIter + size(spatialMasks,3) ) = 0;
        end
    end
    fullFeat( bIter, : ) = [ribbonFeat(:)' logFeat(:)' cleftFeat(:)'];
end