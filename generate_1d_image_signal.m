% Generate 1D input signal from image
function out = generate_1d_image_signal()
img = double(imread('peppers.pgm'));
out = img(:, 128);