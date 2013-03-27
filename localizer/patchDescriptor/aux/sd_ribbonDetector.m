clear; close all; clc;
% Detect Ribbons ...
dataDir = '/research/datasets/synapseDataset/synapseDataset/ribbon/ribSynS1/';
ribFiles = dir( dataDir );

% I = imread('sampleRibbon.png');
% figure(1); imshow(I);

for fileIter = 3:numel(ribFiles)
    I = imread( [dataDir ribFiles(fileIter).name] );
    I = adapthisteq(I(:,:,1));
    figure(1); imshow(I);
    
    threshImg = im2bw(255-I, 0.9);
    threshImg = imopen(threshImg, strel('disk',3));
    
    CC = bwconncomp(threshImg);
    numPixels = cellfun(@numel,CC.PixelIdxList);
    [biggest,idx] = sort(numPixels);
    for iter = 1:numel(idx)-3
        threshImg(CC.PixelIdxList{idx(iter)}) = 0;
    end
    for ribIter = 1:3
    [yy xx] = ind2sub([size(I,1), size(I,2)], CC.PixelIdxList{idx(iter+ribIter)});
    yy = yy ./ size( I,1 );
    xx = xx ./ size( I,2 );
    XX1 = [xx(:)'; yy(:)'];
    
    m1 = mean(XX1, 2);
    XX1 = XX1 - repmat( m1, [ 1 size(XX1, 2) ] );
    
    c1 = cov(XX1');
    [u1, s1, v1] = svd(c1);
    ribbonFeat( ribIter, : ) = [1000*numel(yy)/(size(I,1)*size(I,2) ) ]; 
    end
    load ribbonGMM;
    softProb = pdf(obj,ribbonFeat);
    [vals picker] = sort(softProb, 'descend');
    masker = zeros( size(I,1), size(I,2) );
    if(vals(1) > 0.09 )
        masker( CC.PixelIdxList{idx(iter+picker(1))} ) = 1;
        hold on; contour( masker, [0 0], 'r' ); hold off; title( [num2str(fileIter-2)] );
        pause;
    else
        hold on; contour( threshImg, [0 0], 'r' ); hold off; title( [num2str(fileIter-2)] );
        pause;
        title('Ribbon Not Found');
    end
end




% Alg 2
function validateRibbon()

%%% Validate Ribbon Synapse

ribPath = 'C:\Users\vignesh\Desktop\ribbonTester\posExamples\';
ribFiles = dir(ribPath);
load (['C:\Users\vignesh\Desktop\ribbonTester\' 'ribGt.mat']);
for ribIter = 3:numel( ribFiles )
    I =  ( rgb2gray ( imread([ribPath ribFiles( ribIter ).name]) ) );
    % figure(1); imshow(I);
    % try
    %    hold on; contour( ribbonGt(:,:,ribIter-2 ), 'r' ); hold off; pause(0.5);
    % catch me
    %    ribbonGt(:,:,ribIter-2) = roipoly(I);
    % end
    ribSize(ribIter-2) = sum( sum( ribbonGt(:,:,ribIter-2 ) ) );
    svRat( ribIter - 2)  = computeMoments(ribbonGt(:,:,ribIter-2 ));
    
end
% save([ribPath 'ribGt.mat'], 'ribbonGt');
figure(1); subplot(211); plot( 1:numel(ribFiles)-2, ribSize ); subplot(212); plot( 1:numel(ribFiles)-2, svRat )

for ribIter = 3:numel( ribFiles )
    I =  ( rgb2gray ( imread([ribPath ribFiles( ribIter ).name]) ) );
    ctr = 1;
    for thCtr = .2:.2:.8
        bw(:,:,ctr) = 1 - im2bw(I, thCtr);
        bw(:,:,ctr) = imerode( bw(:,:,ctr), strel('disk', 3) ); % imclose ...
        ctr = ctr + 1;
    end
    stableRegions = prod( bw , 3);
    figure(1);
    imshow(I);
    % hold on;     contour( stableRegions, [0 0], 'r');     hold off;
    CC = bwconncomp(stableRegions);
    numPixels = cellfun(@numel,CC.PixelIdxList);
    [tempo idx] = sort(numPixels, 'descend');
    for idxIter = 4:numel(numPixels)
        stableRegions(CC.PixelIdxList{ idx(idxIter) }) = 0;
    end
    CC = bwconncomp(stableRegions);
    hold on;     contour( stableRegions, [0 0], 'b');  hold off;
    pause(1);
end

function svRat = computeMoments(im1)

[y x] = find(im1);

yM = mean(y);
xM = mean(x);

Cxx = sum( ( x - xM ).^2 );
Cyy = sum( ( y - yM ).^2 );
Cxy = sum( (x-xM) .* (y-yM) );

C = [Cxx Cxy; Cxy Cyy];
[u s v] = svd(C);
svRat = s(1)/s(4);
