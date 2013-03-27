function fl_detectBox(img, I, detSize, numDet, VIS)


if( nargin < 1 )
    clear; close all; clc;
    I = zeros(10,10);
    img = I;
    I(2:4,2:4) = 1;
    I(6,6) = 1;
    I(7,1) = 1;
    detSize = 3;
    numDet = 3;
    VIS = true;
end

if(VIS)
    figure(1);
    imshow( uint8(img) ); %subplot(211); imagesc(intImage); subplot(212); imagesc(I);
    hold on;
end

resBoxes = [];
for detIter = 1:numDet
    intImage = integralImage(I);
    [noRows noCols noSlices] = size(intImage);
    detMask = zeros( size(intImage,1), size(intImage,2) );
    detMask(1:noRows-detSize, 1:noCols-detSize) = 1;
    [startI startJ] = find(detMask == 1 );
    detCandidates = [startI startJ startI+detSize-1 startJ+detSize-1];
    detScores = intImage( sub2ind([noRows noCols], detCandidates(:,3), detCandidates(:,4)) ) + ...
        intImage( sub2ind([noRows noCols], detCandidates(:,1), detCandidates(:,2)) ) - ...
        intImage( sub2ind([noRows noCols], detCandidates(:,1), detCandidates(:,4)) ) - ...
        intImage( sub2ind([noRows noCols], detCandidates(:,3), detCandidates(:,2)) );
    detImage = reshape(detScores, [noRows-detSize noCols-detSize]);
    
    [maxVal idx] = max(detScores);
    resBoxes = [resBoxes; detCandidates( idx, :)];
    I( resBoxes(end,1):resBoxes(end,3), resBoxes(end,2):resBoxes(end,4) ) = 0;
    
    iStart = resBoxes(end,1); iEnd = resBoxes(end,3);
    jStart = resBoxes(end,2); jEnd = resBoxes(end,4);
    if( VIS )
        rectangle('Position', [jStart iStart jEnd-jStart iEnd-iStart], 'LineWidth', 5, 'EdgeColor', [1 0 0] );
        text(jStart, iStart, ['Detection Score ' num2str(detScores(idx))], 'HorizontalAlignment','left', 'BackgroundColor',[.7 .9 .7]);
    end
end
return;
%nmMask = double( detImage > imdilate(detImage, [1 1 1; 1 0 1; 1 1 1]) );
%detImage = nmMask .* detImage; detScores = detImage(:);
%[maxCand idx] = sort(detScores, 'descend');

function intImage = integralImage(I)
%integralImage Compute integral image.
%   J = integralImage(I) computes integral image of an intensity image I.
%   The output integral image, J, is zero padded on top and left, resulting
%   in size(J) = size(I) + 1. This facilitates easy computation of pixel
%   sums along all image boundaries. Integral image, J, is essentially
%   a padded version of CUMSUM(CUMSUM(I)).
%
%   Class Support
%   -------------
%   Intensity image I can be any numeric class. The class of output
%   integral image, J, is double.
%
%   Example
%   -------
%   % Compute the integral image and use it to compute sum of pixels
%   % over a rectangular region in I.
%   I = magic(5)
%
%   % define rectangular region as
%   % [startingRow, startingColumn, endingRow, endingColumn]
%   [sR sC eR eC] = deal(1, 3, 2, 4);
%
%   % compute the sum over the region using the integral image
%   J = integralImage(I);
%   regionSum = J(eR+1,eC+1) - J(eR+1,sC) - J(sR,eC+1) + J(sR,sC)
%
%   See also integralFilter, CUMSUM

%   Copyright 2010 The MathWorks, Inc.
%   $Revision: 1.1.6.1 $  $Date: 2011/10/31 06:52:39 $

%   References:
%      P.A. Viola and M.J. Jones. Rapid object detection using boosted
%      cascade of simple features. In CVPR (1), pages 511-518, 2001.
%
%#codegen
%#ok<*EMCLS>
%#ok<*EMCA>

validateattributes(I, {'numeric','logical'}, {'2d', 'nonsparse', 'real'},...
    'integralImage', 'I');

if ~isempty(I)
    outputSize = size(I) + 1;
    
    intImage = zeros(outputSize);
    intImage(2:end, 2:end) = cumsum(cumsum(double(I),1),2);
else
    intImage = [];
end

