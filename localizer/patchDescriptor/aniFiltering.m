function [filtRes] = aniFiltering(im, im1)
if( nargin <  1 ),
    im = imread('cameraman.tif');
end
%computes MR8 filterbank using recursive Gaussian filters
%input: intensity image
%output: MR8 feature vector
dropper = zeros( size(im,1), size(im,2) );
dropper(40:end-40, 40:end-40 ) = 1;

in = double(im) - mean(mean(im));
in = in ./ sqrt(mean(mean(in .^ 2)));

inVes = double(im1) - mean(mean(im1));
inVes = inVes ./ sqrt(mean(mean(inVes .^ 2)));

i=1;

sfac = 1.0;
mulfac = 2.0;
filtRes = zeros( size(im,1), size(im,2), 7 );
%s1 = 3*sfac; s2 = 1*sfac;
s1 = 24; s2 = 8;
for j=0:0,
    for k=0:5,
        phi = (k/6.0)*180.0;
        
        im2 = anigauss(in, s1, s2, phi, 0, 2);

        % take max of abs response for first order derivative
        % Varma&Zisserman also take abs max of second order...
        %im2 = abs(im2);
        filtRes(:,:,k+1) = im2 .* dropper;
        if (k==0)
            maxim2 = im2;
        else
            maxim2 = max(maxim2, im2);
        end
    end
    maxCleftRes = maxim2;     
end

sigma = 10.0*sfac;
im1 = anigauss(inVes, sigma, sigma, 0.0, 2, 0);
im2 = anigauss(inVes, sigma, sigma, 0.0, 0, 2);
filtRes(:,:,7) = (im1+im2) .* dropper;
%filtRes(:,:,7) = imfilter( inVes, fspecial('log', [49 49], 7), 'conv', 'same' );

return;