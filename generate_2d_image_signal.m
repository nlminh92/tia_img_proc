function [image, map] = generate_2d_image_signal()
[image, map] = imread('peppers.pgm');
image = double(image);