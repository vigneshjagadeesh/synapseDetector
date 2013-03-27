function [lbpFeat, grayFeat, hogFeat] = compFeatures( I, splitter )

% Compute LBP

% [noRows noCols noSlices] = size(I);
% midX = ceil(noCols/2); midY = ceil(noRows/2);
% if( noRows > 512 & noCols > 512 )
%     I = I( midY - 255: midY + 256, midX - 255: midX + 256, 1 );
% end

mapping=getmapping(8,'riu2');
lbpFeat = [];
grayFeat = [];
hogFeat = [];
[Ix Iy] = gradient( double(I) );
if( splitter == 2 )
    [Is  Isg Isa] = splitI(I, Ix, Iy);
else
    Is{1} = I;
    Ise{1} = edge( Is{1}, 'canny' );
    Isg{1} = sqrt( Ix.^2 + Iy.^2 );
    Isa{1} = atan2( Iy, Ix ) * 180 /pi;
end

for i = 1:numel(Is)
    currLbp=lbp(Is{i},1,8,mapping,'h');
    currFeat = currLbp ./ sum( currLbp);
    lbpFeat = [lbpFeat currLbp(:)'];
    
    currGray = histc(Is{i}(:), linspace(0,255,17) );
    currGray = currGray(1:end-1);
    if( sum(currGray > 0  ) )
        currGray = currGray / sum(currGray);
    else
        max(currGray(:))
    end
    grayFeat = [grayFeat currGray(:)' mean(currGray(:)) entropy(currGray(:))];
    
    hogFeat = [hogFeat hog(Isg{i}, Isa{i}, Ise{i})];
    
end

function [Is Isg Isa Ise] = splitI(I, Ix, Iy)

[noRows noCols noSlices] = size(I);
SX = ceil( size(I,2) /2 );
SY = ceil( size(I,1) / 2 );
Is{1} = I(1:SY         , 1:SX         );
Is{2} = I((SY+1):noRows, 1:SX         );
Is{3} = I(1:SY         , (SX+1):noCols);
Is{4} = I((SY+1):noRows, (SX+1):noCols);

Ise{1} = edge( Is{1}, 'canny' );
Ise{2} = edge( Is{2}, 'canny' );
Ise{3} = edge( Is{3}, 'canny' );
Ise{4} = edge( Is{4}, 'canny' );

Ig = sqrt( Ix.^2 + Iy.^2 );
Isg{1} = Ig(1:SY         , 1:SX         );
Isg{2} = Ig((SY+1):noRows, 1:SX         );
Isg{3} = Ig(1:SY         , (SX+1):noCols);
Isg{4} = Ig((SY+1):noRows, (SX+1):noCols);

Isa{1} = atan2( Iy(1:SY         , 1:SX         ), Ix(1:SY         , 1:SX         ) ) * 180 /pi;
Isa{2} = atan2( Iy((SY+1):noRows, 1:SX         ), Ix((SY+1):noRows, 1:SX         ) ) * 180 /pi;
Isa{3} = atan2( Iy(1:SY         , (SX+1):noCols), Ix(1:SY         , (SX+1):noCols) ) * 180 /pi;
Isa{4} = atan2( Iy((SY+1):noRows, (SX+1):noCols), Ix((SY+1):noRows, (SX+1):noCols) ) * 180 /pi;

function feat = hog(gradMag, ang, edgeMap)
noBins = 16;
feat = zeros( noBins, 1 );
gradMag( edgeMap == 0 ) = [];
ang( edgeMap == 0 ) = [];
ang( ang < 0 ) = ang( ang < 0 )  + 360;
binSpace = linspace(0, 360, noBins+1);
for iter = 1:numel(binSpace)-1
    tempInd = ang > binSpace(iter) & ang < binSpace(iter+1);
    feat(iter) = sum( gradMag( tempInd(:) ) );
end
feat = feat(:)';