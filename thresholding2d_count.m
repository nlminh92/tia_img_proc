function [out, retained_point_count] = thresholding2d_count(in, threshold)
%       Set all signal components under threshold to zero
%
% Inputs ([]s are optional)
%		(vector)    in      input signal of arbitrary length,
%                           power of 2 prefered
%		(scalar)    threshold
%                           ignore all signal components under this value
%
% Outputs ([]s are optional)
%		(vector)    out     signal after thresholding
%       (scalar)    retained_point_count
%
out = in;
size_in = size(in);
retained_point_count = numel(in);
for i = 1:size_in(1)
    for j = 1:size_in(2)
        if (abs(out(i, j)) <= threshold)
            retained_point_count = retained_point_count - 1;
            out(i, j) = 0;
        end
    end
end