% Generate 1D sinusoidal input signal
function out = generate_1d_trigo_signal()
img = double(imread('peppers.pgm'));
img_size = size(img);
signal_length = img_size(1);
out = zeros(signal_length, 1);
for i = 1:signal_length
    if i <= signal_length/2
        out(i) = sin(2*pi*i);
    else
        out(i) = 1/2 + sin(2*pi*i);
    end
end