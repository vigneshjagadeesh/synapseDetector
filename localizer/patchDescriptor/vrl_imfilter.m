function [filtResponse feat]= vrl_imfilter(I, F)

%% Function to perform filtering in the frequency domain ...
% PLEASE INPUT IMAGE AND FILTER IN THE SPATIAL DOMAIN
% Input Parameters:
% I - Input Image
% F - Filter Response
% Output Parameters:
% filtResponse - Filtered Image
% 
% tic,
% SFR = zeros(noRows, noCols, noFilters);
% for iter = 1:size(F,3)
%     SFR(:,:,iter) = imfilter( I, F(:,:,iter), 'same', 'conv');
% end
% toc

% Preprocess for variables ... 
I = double(I);
[noRows noCols noSlices] = size(I);
[filtRows filtCols noFilters] = size(F);
Fpad = zeros(noRows, noCols);
filtResponse = zeros( noRows, noCols, noFilters);
feat = zeros(noFilters, 1);
imgFft = fft2(I);

% Perform the filtering
begRow = floor(noRows / 2) - floor(filtRows / 2);
begCol = floor(noCols / 2) - floor(filtCols / 2);
for filtIter = 1:noFilters
    Fpad(begRow:begRow+filtRows-1, begCol:begCol+filtCols-1) = F(:,:,filtIter);
    % Take FFT
    filFft = fft2(Fpad);
    % Take inverse to get the filter response
    D = real ( fftshift( ifft2( imgFft .* filFft ) ) );
    feat(filtIter) = mean(D(:));
    feat(noFilters + filtIter) = sqrt(mean(mean(( D - feat(filtIter) * ones(noRows,noRows) ).^2)));
    filtResponse(:,:,filtIter) = D;
end

feat(1:end/2) = feat(1:end/2) - mean(feat(1:end/2));
feat(1:end/2) = feat(1:end/2)./std(feat(1:end/2));
feat(end/2+1:end) = feat(end/2+1:end) - mean(feat(end/2+1:end));
feat(end/2+1:end) = feat(end/2+1:end)./std(feat(end/2+1:end));