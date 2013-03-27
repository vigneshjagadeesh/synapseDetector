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

mask = cat(3, 1-mask, mask1, mask2, mask3, mask4);


