function out = thresholding(in, threshold)
% Description
%       Set all signal components under threshold to zero
%
% Synopsis
%		out = thresholding(in, threshold)
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
i = 1;
out = in;
while (i <= length(out))
    if (abs(out(i)) <= threshold)
        out(i) = 0;
    end
    i = i + 1;
end