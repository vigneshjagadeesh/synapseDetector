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

% I = imread([emcPath '514_ipl_01_126.tif']);
% I = medfilt2( I, [11 11] );
sizeVals = cell(14,1);
eVal     = cell(14,1);
for synIter = 1:14
    if( synIter < 8 )
        I = double( ( synPatch{synIter} ) );
    else
        I = double( ( nonPatch{synIter-7} ) );
    end
    I = 255 .* ( ( I - min(I(:)) ) ./ max( I(:) ) );
    
    %[outIm,whatScale,Direction] = FrangiFilter2D(im2double(I)); figure(1); imagesc(outIm); pause; continue;
    BW = I < 40;
    BW(1,:) = 0; BW(end, :) = 0; BW(:,1) = 0; BW(:, end) = 0;
    % BW = imopen(BW, strel('disk', 5) );
    CC = bwconncomp(BW);
    numPixels = cellfun(@numel,CC.PixelIdxList);
    [biggest,idx] = sort(numPixels, 'descend');
    
    SNODE = 10;
    SIZEP = 500;
    ELLPT  = 10;
    if( numel(idx) > SNODE )
        for iter = SNODE:numel(idx)
            BW(CC.PixelIdxList{idx(iter)}) = 0;
        end
    end
    ctr = 0;
    for iter = 1:SNODE
        if( numel( CC.PixelIdxList{idx(iter)} ) < SIZEP )
            BW(CC.PixelIdxList{idx(iter)}) = 0;
        else
            [y x] = ind2sub([size(I,1) size(I,2)], CC.PixelIdxList{idx(iter)} );
            svRat = computeSVMoments(y, x);            
            if( svRat > ELLPT )
                ctr = ctr + 1;
            else
                BW(CC.PixelIdxList{idx(iter)}) = 0;
            end
        end
    end
    
    BW = imclose( BW, strel('disk', 5) );
    CC = bwconncomp(BW);
    for iter = 1:CC.NumObjects
                [y x] = ind2sub([size(I,1) size(I,2)], CC.PixelIdxList{iter} );
                svRat = computeSVMoments(y, x);  
                if(  svRat > ELLPT  )
                    eVal{synIter} = [eVal{synIter}; svRat];
                    sizeVals{synIter} = [sizeVals{synIter}; numel( CC.PixelIdxList{iter} )];
                else
                        BW(CC.PixelIdxList{iter}) = 0;
                end
    end
    figure(1); imshow( uint8(I) ); hold on; contour( BW, [0 0], 'r'); hold off;
    
    if( synIter < 8 )
        figure(2); hold on; plot(sizeVals{synIter}, 'r'); hold off;
        figure(3); hold on; plot(sort( eVal{synIter}, 'descend' ), 'r'); hold off;
    else
        figure(2); hold on; plot(sizeVals{synIter}, 'b'); hold off;
        figure(3); hold on; plot(sort( eVal{synIter}, 'descend'), 'b'); hold off;
    end
end