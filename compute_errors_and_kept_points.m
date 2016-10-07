function [me, mse, kept_count] = compute_errors_and_kept_points(in, thres)
% Compute mean error, mean square error, and retained point count after
% thresholding
out = quadratic_2d(in);
[outFiltered, kept_count] = thresholding2d_count(out, thres);
inFiltered = quadratic_2d_inverse(outFiltered);
[me, mse] = compute_errors(in, inFiltered);