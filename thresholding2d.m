function out = thresholding2d(in, threshold)
% Description
%       Set all signal components under threshold to zero
%
% Synopsis
%		out = thresholding2d(in, threshold)
%
% Inputs ([]s are optional)
%		(vector)    in      input signal of arbitrary length,
%                           power of 2 prefered
%		(scalar)    threshold
%                           ignore all signal components under this value
%
% Outputs ([]s are optional)
%		(vector)    out     signal after thresholding
%
% Examples
%		img = imread('zelda.bmp');
%		in = img(:,36);
%		out = thresholding(in, 128);
%
out = in;
size_in = size(in);
for i = 1:size_in(1)
    for j = 1:size_in(2)
        if (abs(out(i, j)) <= threshold)
            out(i, j) = 0;
        end
    end
end