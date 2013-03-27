function connectomeCaller(xStart, xEnd, sliceNumber)
if( nargin < 2 )
clear; close all; clc;
xStart = 1;
xEnd = 30;
sliceNumber  = 145;
else
    xStart = str2double(xStart);
    xEnd = str2double(xEnd);
    sliceNumber = str2double( sliceNumber );
end
% This is a script to perform synapse detection on an entire slice of RC1 ... and can operate in three modes,
% random where it picks random locations for testing
% single where it processes
% sequential where is scans brute force sequentially through the data
model.minSize = 2667;
model.maxSize = 15862;
model.meanInt = 107.3966;
model.stdInt  = 13.7951;
%addpath('/cluster/home/retinamap/synapseProject/source/synapseDetector/benchmark');
%addpath( genpath( ['/cluster/home/retinamap/synapseProject/interfaces'] ) );

% Initialize Parameters
% load synapseGroundTruths; % loadSynGT;
opMode = 'sequential'; % 'sequential' and 'single'

downSampleLevel = 1;
Xval = []; Yval = [];
rc1Path = sprintf(['/cluster/retinamap/RABBIT/%.4d/mosaic/%.3d/'], sliceNumber, downSampleLevel);
chAPath = sprintf(['/cluster/retinamap/RABBIT/%.4d/ChannelA/%.3d/'], sliceNumber, downSampleLevel);
xmlFile = sprintf([rc1Path '%.4d.xml'], sliceNumber);
try,
sliceDetails = xml2struct('/cluster/retinamap/RABBIT/0143/mosaic/001/0143.xml');
sliceWidth = str2num(sliceDetails.Attributes(4).Value);
sliceHeight = str2num(sliceDetails.Attributes(5).Value);
catch me,
    sliceWidth = 458;
    sliceHeight = 460;
end
Xval = 262;
Yval = 206;
if( isempty(Xval) ),
    Xval = ceil(rand*sliceWidth);
end
if( isempty(Yval) ),
    Yval = ceil(rand*sliceHeight);
end
GridRange = 5;
if( strcmp(opMode, 'random') || strcmp(opMode, 'single') )
    xrange = Xval;
    yrange = Yval;
    VIS = true;
elseif( strcmp(opMode, 'sequential') )
    offset = GridRange+1;
    xrange = (offset):(2*GridRange+1):(sliceWidth-offset);
    yrange = (offset):(2*GridRange+1):(sliceHeight-offset);
    xrange( xrange < xStart ) = [];
    xrange( xrange > xEnd ) = [];
    VIS = false;
end

% Create Dirs
opDir = ['/cluster/home/retinamap/datasetCreator/utahDataBase/AA_largeScaleDetector/' num2str(sliceNumber) '_Output'];
opFileName = [opDir 'synMarginal%.3d_%.3d.mat'];
if( ~isdir(opDir) )
mkdir(opDir);
end

imgName = [rc1Path '%.4d_X%.3d_Y%.3d.png'];
heatImgName = [chAPath '%.4d_X%.3d_Y%.3d.png'];
imgSize = [256*(2*GridRange+1) 256*(2*GridRange+1)];

%indStoreName = ['/cluster/home/retinamap/synZ_' num2str(sliceNumber) '_X_' num2str(xrange(1)) '.mat'];
%h = waitbar(0,'Please wait...');
fullCtr = 0;
miniCtr = 0;
totImg = numel(xrange) * numel(yrange);
%fileStore = cell(totImg*(2*GridRange+1)^2,1);
%fileLabel = zeros(totImg*(2*GridRange+1)^2,1);
 for Xsing = xrange
    for Ysing = yrange
        fullCtr = fullCtr + 1;
        %waitbar(fullCtr/totImg,h);
        display(['Processing ratio ' num2str(fullCtr/totImg) ]);
        opFile = sprintf(opFileName,Xsing, Ysing);
        currImg = uint8(zeros(imgSize));
        xCtr = 0;
        for xIter = (Xsing - GridRange) : (Xsing + GridRange)
            yCtr = 0;
            for yIter = (Ysing - GridRange) : (Ysing + GridRange)
            miniCtr = miniCtr + 1;
		try,
                    I = imread(sprintf(imgName, sliceNumber, xIter, yIter));
                catch me,
                    I = zeros(256, 256);
                end
                currImg((yCtr*256+1):((yCtr+1)*256),(xCtr*256+1):((xCtr+1)*256)) = I;
                yCtr = yCtr + 1;
            end
            xCtr = xCtr + 1;
        end
        tic, [detMask finalMask256] = simpleDetect(currImg, model, 0, VIS); toc
        xCtr = 0;
        for xIter = (Xsing - GridRange) : (Xsing + GridRange)
            yCtr = 0;
            for yIter = (Ysing - GridRange) : (Ysing + GridRange)
                currImName = sprintf(heatImgName, sliceNumber, xIter, yIter);
                heatImg = finalMask256(( (256*yCtr) + 1 ):( 256*(yCtr+1) ), ( (256*xCtr) + 1 ):( 256*(xCtr+1) ) );
                heatImg = uint8( 255.*heatImg);
                imwrite(heatImg, currImName);
		%fileStore{fullCtr} = currImName;
		%fileLabel(fullCtr) = (heatImg(128,128)>0);
                yCtr = yCtr + 1;
            end
            xCtr = xCtr + 1;
        end
        %save(opFile, 'synMarginal');
    end
end
%save(indStoreName, 'fileStore', 'fileLabel');
exit(0);
