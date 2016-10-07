function [meanError, meanSquareError] = compute_errors(before, after)
% Description
%       Compute mean errors
%
% Inputs
%		(matrix)    before      original signal
%		(matrix)    after       reconstructed signal
%
% Outputs
%		(array 1x2) out         mean absolute error and mean square error
%
inSize = size(before);
before = double(before);
after = double(after);
%meanError = norm(before - after, 1) / (inSize(1)*inSize(2));
%meanSquareError = norm(before - after, 2)^2 / (inSize(1)*inSize(2));
meanError = mean2(abs(before - after));
meanSquareError = mean2((before - after).^2);