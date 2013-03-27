clear; close all; clc;

% Load Ribbon Model .. featurs are [ numPixels eig1 eig2 ] ... ribbon model
% has been learnt
load savedRibbons;

for iter = 1:size(currRibbon,3)
    ribbonImg = double( currRibbon(:,:,iter) );
    CC = bwconncomp( ribbonImg );
    [yy xx] = find( ribbonImg );
    yy = yy ./ size( currRibbon,1 );
    xx = xx ./ size( currRibbon,2 );
    XX1 = [xx(:)'; yy(:)'];
    
    m1 = mean(XX1, 2);
    XX1 = XX1 - repmat( m1, [ 1 size(XX1, 2) ] );
    
    c1 = cov(XX1');
    [u1, s1, v1] = svd(c1);
    ribbonFeat( iter, : ) = [1000*numel(yy)/numel(ribbonImg)];    
end

figure;
%scatter(ribbonFeat(:,1),ribbonFeat(:,2),10,'*')
hold on;
options = statset('Display','final');
obj = gmdistribution.fit(ribbonFeat,5,'Options',options);
save('ribbonGMM', 'obj'); %return;
% h = ezcontour(@(x,y)pdf(obj,[x y]),[0 max(ribbonFeat(:,1) )],[0 max(ribbonFeat(:,2) ) ]);
hold off;

softProb = pdf(obj,ribbonFeat);

for iter = 1:size(currRibbon,3)
    ribbonImg = double( currRibbon(:,:,iter) );    
    imshow(ribbonImg); hold on; contour( currRibbon(:,:,iter), 'r' ); hold off;
    %title( [ 'Size: ' num2str(ribbonFeat(iter,1)) 'Shape: '  num2str(ribbonFeat(iter,2)) 'PDF: ' num2str(softProb(iter)) ] ); pause;
    title( [ 'Size: ' num2str(ribbonFeat(iter,1)) 'PDF: ' num2str(softProb(iter)) ] ); pause;
end
