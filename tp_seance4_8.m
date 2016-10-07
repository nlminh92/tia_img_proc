function tp_seance4_8()
disp('Exercice 7: Erreur en fonction du seuil.');
in = generate_2d_image_signal();
n = 20;
me = zeros(n, 3);
mse = zeros(n, 3);
T = zeros(n, 1);
for i = 1:n
    T(i) = 0.01*i;
    [me(i, 1), mse(i, 1), ~] = compute_transformation_errors(in, T(i), 'c');
    [me(i, 2), mse(i, 2), ~] = compute_transformation_errors(in, T(i), 'd');
    [me(i, 3), mse(i, 3), ~] = compute_transformation_errors(in, T(i), 'g');
end
figure
plot(T, me(:, 1), '--rs', 'LineWidth', 2,...
                'MarkerEdgeColor', 'k',...
                'MarkerFaceColor', 'g',...
                'MarkerSize', 5);
hold on;
plot(T, mse(:, 1), '--ro', 'LineWidth', 2,...
                'MarkerEdgeColor', 'k',...
                'MarkerFaceColor', 'g',...
                'MarkerSize', 5);
plot(T, me(:, 2), '--gs', 'LineWidth', 2,...
                'MarkerEdgeColor', 'k',...
                'MarkerFaceColor', 'g',...
                'MarkerSize', 5);
plot(T, mse(:, 2), '--go', 'LineWidth', 2,...
                'MarkerEdgeColor', 'k',...
                'MarkerFaceColor', 'g',...
                'MarkerSize', 5);
plot(T, me(:, 3), '--bs', 'LineWidth', 2,...
                'MarkerEdgeColor', 'k',...
                'MarkerFaceColor', 'g',...
                'MarkerSize', 5);
plot(T, mse(:, 3), '--bo', 'LineWidth', 2,...
                'MarkerEdgeColor', 'k',...
                'MarkerFaceColor', 'g',...
                'MarkerSize', 5);
legend('Transform c_{min}, ME', 'Transform c_{min}, MSE',...
    'Transform d_{min}, ME', 'Transform d_{min}, MSE',...
    'Transform g_{min}, ME', 'Transform g_{min}, MSE');
xlabel('Threshold')
ylabel('Errors')

function [me, mse, kept_count] = compute_transformation_errors(in, thres, transform_id)
% Compute mean error, mean square error, and retained point count after
% thresholding
% Input
%       (matrix)    in          input signal of length power-of-2
%       (scalar)    thres       threshold
%       (string)    transform_id ('c' or 'd' or 'g')
%                               id of 1D nonlinear quadratic transform
% Output
%       (matrix)    out         transformed signal of identical length
%
if ~exist('transform_id', 'var') || isempty(transform_id)
    transform_id = 'c';
end
if transform_id == 'c'
    disp(['Transform with c_min, threshold = ', num2str(thres)]);
elseif transform_id == 'd'
    disp(['Transform with d_min, threshold = ', num2str(thres)]);
else
    disp(['Transform with g_min, threshold = ', num2str(thres)]);
end
out = nonlinear_quadratic_2d(in, transform_id);
[outFiltered, kept_count] = thresholding2d_count(out, thres);
inFiltered = nonlinear_quadratic_2d_inverse(outFiltered, transform_id);
[me, mse] = compute_errors(in, inFiltered);

function out = nonlinear_quadratic_2d(in, transform_id)
% Implementation of 2D nonlinear quadratic transform
% Input
%       (matrix)    in          input signal of length power-of-2
% Output
%       (matrix)    out         transformed signal of identical length
%
if ~exist('transform_id', 'var') || isempty(transform_id)
    transform_id = 'c';
end
in_size = size(in);
temp = zeros(in_size);
out = zeros(in_size);
% apply 1D nonlinear quadratic transform to rows
for i = 1:in_size(1)
    temp(i,:) = nonlinear_quadratic_1d(in(i,:), transform_id);
end
% apply 1D nonlinear quadratic transform to columns
for j = 1:in_size(2)
    out(:,j) = nonlinear_quadratic_1d(temp(:,j), transform_id);
end

function out = nonlinear_quadratic_2d_inverse(in, transform_id)
% Implementation of 2D inverse nonlinear quadratic transform
% Input([] is optional)
%       (matrix)    in          input signal of length power-of-2
%       (string)    transform_id ('c' or 'd' or 'g')
%                               id of 1D nonlinear quadratic transform
% Output
%       (matrix)    out         transformed signal of identical length
%
if ~exist('transform_id', 'var') || isempty(transform_id)
    transform_id = 'c';
end
in_size = size(in);
temp = zeros(in_size);
out = zeros(in_size);
% apply 1D inverse nonlinear quadratic transform to columns
for j = 1:in_size(2)
    temp(:,j) = nonlinear_quadratic_1d_inverse(in(:,j), transform_id);
end
% apply 1D inverse nonlinear quadratic transform to rows
for i = 1:in_size(1)
    out(i,:) = nonlinear_quadratic_1d_inverse(temp(i,:), transform_id);
end

function out = nonlinear_quadratic_1d(in, transform_id)
% nonlinear_quadratic_1d computes the forward quadratic transform
% of a vector.
% Input([] is optional)
%       (vector)    in      input signal of length power-of-2
%       (string)    transform_id ('c' or 'd' or 'g')
%                           id of 1D nonlinear quadratic transform
% Output
%       (vector)    output signal of length power-of-2
%
if ~exist('transform_id', 'var') || isempty(transform_id)
    transform_id = 'c';
end
signal_length = length(in);
signal_length_l2 = log2(signal_length);
out = in;
temp = in;
m = signal_length;
for k = signal_length_l2:-1:3
    m = m/2;
    for i = 1:m
        temp(i) = 0.5*(out(2*i - 1) + out(2*i));
    end
    % Border value, i = 1
    i = 1;
    c = abs(temp(m) - temp(i)) + abs(temp(i) - temp(i + 1));
    g = abs(temp(m) - temp(m - 1)) + abs(temp(i) - temp(m));
    d = abs(temp(i + 1) - temp(i)) + abs(temp(i + 2) - temp(i + 1));
    if transform_id == 'c'
        temp(m + i) = out(2*i - 1) - (temp(i) + (temp(m) - temp(i + 1))/8);
    elseif transform_id == 'd'
        temp(m + i) = out(2*i - 1) - (5*temp(i)/8 + temp(i + 1)/2 - temp(i + 2)/8);
    else
        temp(m + i) = out(2*i - 1) - (11*temp(i)/8 - temp(m)/2 + temp(m - 1)/8);
    end
    
    % Border value, i = 2
    i = 2;
    c = abs(temp(i - 1) - temp(i)) + abs(temp(i) - temp(i + 1));
    g = abs(temp(i - 1) - temp(m)) + abs(temp(i) - temp(i - 1));
    d = abs(temp(i + 1) - temp(i)) + abs(temp(i + 2) - temp(i + 1));
    if transform_id == 'c'
        temp(m + i) = out(2*i - 1) - (temp(i) + (temp(i - 1) - temp(i + 1))/8);
    elseif transform_id == 'd'
        temp(m + i) = out(2*i - 1) - (5*temp(i)/8 + temp(i + 1)/2 - temp(i + 2)/8);
    else
        temp(m + i) = out(2*i - 1) - (11*temp(i)/8 - temp(i - 1)/2 + temp(m)/8);
    end
    
    % Border value, i = m - 1
    i = m - 1;
    c = abs(temp(i - 1) - temp(i)) + abs(temp(i) - temp(i + 1));
    g = abs(temp(i - 1) - temp(i - 2)) + abs(temp(i) - temp(i - 1));
    d = abs(temp(i + 1) - temp(i)) + abs(temp(1) - temp(i + 1));
    if transform_id == 'c'
        temp(m + i) = out(2*i - 1) - (temp(i) + (temp(i - 1) - temp(i + 1))/8);
    elseif transform_id == 'd'
        temp(m + i) = out(2*i - 1) - (5*temp(i)/8 + temp(i + 1)/2 - temp(1)/8);
    else
        temp(m + i) = out(2*i - 1) - (11*temp(i)/8 - temp(i - 1)/2 + temp(i - 2)/8);
    end
    
    % Border value, i = m
    i = m;
    c = abs(temp(i - 1) - temp(i)) + abs(temp(i) - temp(1));
    g = abs(temp(i - 1) - temp(i - 2)) + abs(temp(i) - temp(i - 1));
    d = abs(temp(i + 1) - temp(i)) + abs(temp(2) - temp(1));
    if transform_id == 'c'
        temp(m + i) = out(2*i - 1) - (temp(i) + (temp(i - 1) - temp(1))/8);
    elseif transform_id == 'd'
        temp(m + i) = out(2*i - 1) - (5*temp(i)/8 + temp(1)/2 - temp(2)/8);
    else
        temp(m + i) = out(2*i - 1) - (11*temp(i)/8 - temp(i - 1)/2 + temp(i - 2)/8);
    end
    
    for i = 3:(m - 2)
        c = abs(temp(i - 1) - temp(i)) + abs(temp(i) - temp(i + 1));
        g = abs(temp(i - 1) - temp(i - 2)) + abs(temp(i) - temp(i - 1));
        d = abs(temp(i + 1) - temp(i)) + abs(temp(i + 2) - temp(i + 1));
        if transform_id == 'c'
            temp(m + i) = out(2*i - 1) - (temp(i) + (temp(i - 1) - temp(i + 1))/8);
        elseif transform_id == 'd'
            temp(m + i) = out(2*i - 1) - (5*temp(i)/8 + temp(i + 1)/2 - temp(i + 2)/8);
        else
            temp(m + i) = out(2*i - 1) - (11*temp(i)/8 - temp(i - 1)/2 + temp(i - 2)/8);
        end
    end
    out = temp;
end

function out = nonlinear_quadratic_1d_inverse(in, transform_id)
% nonlinear_quadratic_1d_inverse computes
% the inverse nonlinear quadratic transform of a vector.
% Input([] is optional)
%       (vector)    in      input signal of length power-of-2
%       (string)    transform_id ('c' or 'd' or 'g')
%                           id of 1D nonlinear quadratic transform
% Output
%       (vector)    output signal of length power-of-2
%
if ~exist('transform_id', 'var') || isempty(transform_id)
    transform_id = 'c';
end
signal_length = length(in);
signal_length_l2 = log2(signal_length);
out = in;
temp = in;
m = 4;
for k = 3:1:signal_length_l2
        
    % Border value, i = 1
    i = 1;
    c = abs(out(m) - out(i)) + abs(out(i) - out(i + 1));
    g = abs(out(m) - out(m - 1)) + abs(out(i) - out(m));
    d = abs(out(i + 1) - out(i)) + abs(out(i + 2) - out(i + 1));
    if transform_id == 'c'
        temp(2*i - 1) = out(m + i) + (out(i) + (out(m) - out(i + 1))/8);
    elseif transform_id == 'd'
        temp(2*i - 1) = out(m + i) + (5*out(i)/8 + out(i + 1)/2 - out(i + 2)/8);
    else
        temp(2*i - 1) = out(m + i) + (11*out(i)/8 - out(m)/2 + out(m - 1)/8);
    end
    
    % Border value, i = 2
    i = 2;
    c = abs(out(i - 1) - out(i)) + abs(out(i) - out(i + 1));
    g = abs(out(i - 1) - out(m)) + abs(out(i) - out(i - 1));
    d = abs(out(i + 1) - out(i)) + abs(out(i + 2) - out(i + 1));
    if transform_id == 'c'
        temp(2*i - 1) = out(m + i) + (out(i) + (out(i - 1) - out(i + 1))/8);
    elseif transform_id == 'd'
        temp(2*i - 1) = out(m + i) + (5*out(i)/8 + out(i + 1)/2 - out(i + 2)/8);
    else
        temp(2*i - 1) = out(m + i) + (11*out(i)/8 - out(i - 1)/2 + out(m)/8);
    end

    % Border value, i = m - 1
    i = m - 1;
    c = abs(out(i - 1) - out(i)) + abs(out(i) - out(i + 1));
    g = abs(out(i - 1) - out(i - 2)) + abs(out(i) - out(i - 1));
    d = abs(out(i + 1) - out(i)) + abs(out(1) - out(i + 1));
    if transform_id == 'c'
        temp(2*i - 1) = out(m + i) + (out(i) + (out(i - 1) - out(i + 1))/8);
    elseif transform_id == 'd'
        temp(2*i - 1) = out(m + i) + (5*out(i)/8 + out(i + 1)/2 - out(1)/8);
    else
        temp(2*i - 1) = out(m + i) + (11*out(i)/8 - out(i - 1)/2 + out(i - 2)/8);
    end
    
    % Border value, i = m
    i = m;
    c = abs(out(i - 1) - out(i)) + abs(out(i) - out(1));
    g = abs(out(i - 1) - out(i - 2)) + abs(out(i) - out(i - 1));
    d = abs(out(i + 1) - out(i)) + abs(out(2) - out(1));
    if transform_id == 'c'
        temp(2*i - 1) = out(m + i) + (out(i) + (out(i - 1) - out(1))/8);
    elseif transform_id == 'd'
        temp(2*i - 1) = out(m + i) + (5*out(i)/8 + out(1)/2 - out(2)/8);
    else
        temp(2*i - 1) = out(m + i) + (11*out(i)/8 - out(i - 1)/2 + out(i - 2)/8);
    end
    
    for i = 3:(m - 2)
        c = abs(out(i - 1) - out(i)) + abs(out(i) - out(i + 1));
        g = abs(out(i - 1) - out(i - 2)) + abs(out(i) - out(i - 1));
        d = abs(out(i + 1) - out(i)) + abs(out(i + 2) - out(i + 1));
        if transform_id == 'c'
            temp(2*i - 1) = out(m + i) + (out(i) + (out(i - 1) - out(i + 1))/8);
        elseif transform_id == 'd'
            temp(2*i - 1) = out(m + i) + (5*out(i)/8 + out(i + 1)/2 - out(i + 2)/8);
        else
            temp(2*i - 1) = out(m + i) + (11*out(i)/8 - out(i - 1)/2 + out(i - 2)/8);
        end
    end
    for i = 1:m
       temp(2*i) =  2*out(i) - temp(2*i - 1);
    end
    out = temp;
    m = m * 2;
end