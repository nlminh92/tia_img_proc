% Generate 1D incremental input signal
function out = generate_1d_increment_signal()
img = double(imread('peppers.pgm'));
img_size = size(img);
signal_length = img_size(1);
out = zeros(signal_length, 1);
for i = 1:signal_length
    out(i) = i;
end