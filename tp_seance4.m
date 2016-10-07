function tp_seance4()
clc;
clear all;
close all;
exercice = input('Choose an exercice to run\n 1 / 2 / 5 / 6 / 7 / 8\n', 's');
if (exercice == '1')
    exercice12();
elseif (exercice == '2')
    exercice12();
elseif (exercice == '5')
    exercice5();
elseif (exercice == '6')
    exercice6();
elseif (exercice == '7')
    exercice7();
elseif (exercice == '8')
    tp_seance4_8();
else
    return;
end

function exercice12()
disp('Exercice 1 et 2: appliquer la transformation quadratique non-lineaire 1D');
in1 = generate_1d_increment_signal();
in2 = generate_1d_image_signal();
in3 = generate_1d_trigo_signal();
in4 = 100*randn(1, 16)

out1 = nonlinear_quadratic_1d(in1);
out2 = nonlinear_quadratic_1d(in2);
out3 = nonlinear_quadratic_1d(in3);
out4 = nonlinear_quadratic_1d(in4)

in1T = nonlinear_quadratic_1d_inverse(out1);
in2T = nonlinear_quadratic_1d_inverse(out2);
in3T = nonlinear_quadratic_1d_inverse(out3);
in4T = nonlinear_quadratic_1d_inverse(out4)

disp('Erreurs de signal cumulatif avant et apres transformation:');
[ME1 MSE1] = compute_errors(in1, in1T)
disp('Erreurs de signal d''image avant et apres transformation:');
[ME2 MSE2] = compute_errors(in2, in2T)
disp('Erreurs de signal trigonométrique avant et apres transformation:');
[ME3 MSE3] = compute_errors(in3, in3T)
disp('Erreurs de signal trigonométrique avant et apres transformation:');
[ME4 MSE4] = compute_errors(in4, in4T)

function exercice5()
disp('Exercice 5: appliquer la transformation quadratique non-lineaire 2D');
[in, map] = generate_2d_image_signal();
out = nonlinear_quadratic_2d(in);
inT = nonlinear_quadratic_2d_inverse(out);
disp('Erreurs d''images avant et apres transformation quadratique non-lineaire:');
[ME MSE] = compute_errors(in, inT)
figure;
subplot(1,2,1);imshow(in, map);title('image original avant trans. quadratique non-lineaire');
subplot(1,2,2);imshow(inT, map);title('image apres reconstruction');

function exercice6()
disp('Exercice 6: appliquer la transformation quadratique non-lineaire 2D avec seuillage pour debruitage');
T = 0;
[in, map] = generate_2d_image_signal();
in = double(in);
inSize = size(in);
sigma = 0.5 * (max(in(:)) - min(in(:))); % noise level
inNoisy = in + sigma * randn(inSize(1), inSize(2)); % noisy image
outNoisy = nonlinear_quadratic_2d(inNoisy);
outFiltered = thresholding2d(outNoisy, T);
inFiltered = nonlinear_quadratic_2d_inverse(outFiltered);
disp(['Erreurs d''images avant et apres transformation quadratique + seuillage, seuil = ', num2str(T)]);
[ME MSE] = compute_errors(in, inFiltered)
figure;
subplot(1,3,1); imshow(in, map); title('original image');
subplot(1,3,2); imshow(inNoisy, map); title(['noisy image, deviation = ', num2str(sigma)]);
subplot(1,3,3); imshow(inFiltered, map); title(['reconstructed image, threshold = ', num2str(T)]);

function [me, mse, kept_count] = compute_transformation_errors(in, thres)
% Compute mean error, mean square error, and retained point count after
% thresholding
out = nonlinear_quadratic_2d(in);
[outFiltered, kept_count] = thresholding2d_count(out, thres);
inFiltered = nonlinear_quadratic_2d_inverse(outFiltered);
[me, mse] = compute_errors(in, inFiltered);

function exercice7()
disp('Exercice 7: Erreur en fonction du seuil.');
in = generate_2d_image_signal();
n = 20;
me = zeros(n, 1);
mse = zeros(n, 1);
T = zeros(n, 1);
for i = 1:n
    T(i) = 0.01*i;
    [me(i), mse(i), ~] = compute_transformation_errors(in, T(i));
end
figure
plot(T, me, '--bs', 'LineWidth', 2,...
                'MarkerEdgeColor', 'k',...
                'MarkerFaceColor', 'g',...
                'MarkerSize', 5);
hold on;
plot(T, mse, '--bo', 'LineWidth', 2,...
                'MarkerEdgeColor', 'k',...
                'MarkerFaceColor', 'g',...
                'MarkerSize', 5);
legend('ME', 'MSE');
xlabel('Threshold')
ylabel('Errors')

function out = nonlinear_quadratic_2d(in)
% Implementation of 2D nonlinear quadratic transform
% Input
%       (matrix)    in          input signal of length power-of-2
% Output
%       (matrix)    out         transformed signal of identical length
%
in_size = size(in);
temp = zeros(in_size);
out = zeros(in_size);
% apply 1D nonlinear quadratic transform to rows
for i = 1:in_size(1)
    temp(i,:) = nonlinear_quadratic_1d(in(i,:));
end
% apply 1D nonlinear quadratic transform to columns
for j = 1:in_size(2)
    out(:,j) = nonlinear_quadratic_1d(temp(:,j));
end

function out = nonlinear_quadratic_2d_inverse(in)
% Implementation of 2D inverse nonlinear quadratic transform
% Input
%       (matrix)    in          input signal of length power-of-2
% Output
%       (matrix)    out         transformed signal of identical length
%
in_size = size(in);
temp = zeros(in_size);
out = zeros(in_size);
% apply 1D inverse nonlinear quadratic transform to columns
for j = 1:in_size(2)
    temp(:,j) = nonlinear_quadratic_1d_inverse(in(:,j));
end
% apply 1D inverse nonlinear quadratic transform to rows
for i = 1:in_size(1)
    out(i,:) = nonlinear_quadratic_1d_inverse(temp(i,:));
end

function out = nonlinear_quadratic_1d(in)
% nonlinear_quadratic_1d computes the forward quadratic transform
% of a vector.
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
    if c == min([c g d])
        temp(m + i) = out(2*i - 1) - (temp(i) + (temp(m) - temp(i + 1))/8);
    elseif g == min([c g d])
        temp(m + i) = out(2*i - 1) - (5*temp(i)/8 + temp(i + 1)/2 - temp(i + 2)/8);
    else
        temp(m + i) = out(2*i - 1) - (11*temp(i)/8 - temp(m)/2 + temp(m - 1)/8);
    end
    
    % Border value, i = 2
    i = 2;
    c = abs(temp(i - 1) - temp(i)) + abs(temp(i) - temp(i + 1));
    g = abs(temp(i - 1) - temp(m)) + abs(temp(i) - temp(i - 1));
    d = abs(temp(i + 1) - temp(i)) + abs(temp(i + 2) - temp(i + 1));
    if c == min([c g d])
        temp(m + i) = out(2*i - 1) - (temp(i) + (temp(i - 1) - temp(i + 1))/8);
    elseif g == min([c g d])
        temp(m + i) = out(2*i - 1) - (5*temp(i)/8 + temp(i + 1)/2 - temp(i + 2)/8);
    else
        temp(m + i) = out(2*i - 1) - (11*temp(i)/8 - temp(i - 1)/2 + temp(m)/8);
    end
    
    % Border value, i = m - 1
    i = m - 1;
    c = abs(temp(i - 1) - temp(i)) + abs(temp(i) - temp(i + 1));
    g = abs(temp(i - 1) - temp(i - 2)) + abs(temp(i) - temp(i - 1));
    d = abs(temp(i + 1) - temp(i)) + abs(temp(1) - temp(i + 1));
    if c == min([c g d])
        temp(m + i) = out(2*i - 1) - (temp(i) + (temp(i - 1) - temp(i + 1))/8);
    elseif g == min([c g d])
        temp(m + i) = out(2*i - 1) - (5*temp(i)/8 + temp(i + 1)/2 - temp(1)/8);
    else
        temp(m + i) = out(2*i - 1) - (11*temp(i)/8 - temp(i - 1)/2 + temp(i - 2)/8);
    end
    
    % Border value, i = m
    i = m;
    c = abs(temp(i - 1) - temp(i)) + abs(temp(i) - temp(1));
    g = abs(temp(i - 1) - temp(i - 2)) + abs(temp(i) - temp(i - 1));
    d = abs(temp(i + 1) - temp(i)) + abs(temp(2) - temp(1));
    if c == min([c g d])
        temp(m + i) = out(2*i - 1) - (temp(i) + (temp(i - 1) - temp(1))/8);
    elseif g == min([c g d])
        temp(m + i) = out(2*i - 1) - (5*temp(i)/8 + temp(1)/2 - temp(2)/8);
    else
        temp(m + i) = out(2*i - 1) - (11*temp(i)/8 - temp(i - 1)/2 + temp(i - 2)/8);
    end
    
    for i = 3:(m - 2)
        c = abs(temp(i - 1) - temp(i)) + abs(temp(i) - temp(i + 1));
        g = abs(temp(i - 1) - temp(i - 2)) + abs(temp(i) - temp(i - 1));
        d = abs(temp(i + 1) - temp(i)) + abs(temp(i + 2) - temp(i + 1));
        if c == min([c g d])
            temp(m + i) = out(2*i - 1) - (temp(i) + (temp(i - 1) - temp(i + 1))/8);
        elseif g == min([c g d])
            temp(m + i) = out(2*i - 1) - (5*temp(i)/8 + temp(i + 1)/2 - temp(i + 2)/8);
        else
            temp(m + i) = out(2*i - 1) - (11*temp(i)/8 - temp(i - 1)/2 + temp(i - 2)/8);
        end
    end
    out = temp;
end

function out = nonlinear_quadratic_1d_inverse(in)
% nonlinear_quadratic_1d_inverse computes
% the inverse nonlinear quadratic transform of a vector.
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
    if c == min([c g d])
        temp(2*i - 1) = out(m + i) + (out(i) + (out(m) - out(i + 1))/8);
    elseif g == min([c g d])
        temp(2*i - 1) = out(m + i) + (5*out(i)/8 + out(i + 1)/2 - out(i + 2)/8);
    else
        temp(2*i - 1) = out(m + i) + (11*out(i)/8 - out(m)/2 + out(m - 1)/8);
    end
    
    % Border value, i = 2
    i = 2;
    c = abs(out(i - 1) - out(i)) + abs(out(i) - out(i + 1));
    g = abs(out(i - 1) - out(m)) + abs(out(i) - out(i - 1));
    d = abs(out(i + 1) - out(i)) + abs(out(i + 2) - out(i + 1));
    if c == min([c g d])
        temp(2*i - 1) = out(m + i) + (out(i) + (out(i - 1) - out(i + 1))/8);
    elseif g == min([c g d])
        temp(2*i - 1) = out(m + i) + (5*out(i)/8 + out(i + 1)/2 - out(i + 2)/8);
    else
        temp(2*i - 1) = out(m + i) + (11*out(i)/8 - out(i - 1)/2 + out(m)/8);
    end

    % Border value, i = m - 1
    i = m - 1;
    c = abs(out(i - 1) - out(i)) + abs(out(i) - out(i + 1));
    g = abs(out(i - 1) - out(i - 2)) + abs(out(i) - out(i - 1));
    d = abs(out(i + 1) - out(i)) + abs(out(1) - out(i + 1));
    if c == min([c g d])
        temp(2*i - 1) = out(m + i) + (out(i) + (out(i - 1) - out(i + 1))/8);
    elseif g == min([c g d])
        temp(2*i - 1) = out(m + i) + (5*out(i)/8 + out(i + 1)/2 - out(1)/8);
    else
        temp(2*i - 1) = out(m + i) + (11*out(i)/8 - out(i - 1)/2 + out(i - 2)/8);
    end
    
    % Border value, i = m
    i = m;
    c = abs(out(i - 1) - out(i)) + abs(out(i) - out(1));
    g = abs(out(i - 1) - out(i - 2)) + abs(out(i) - out(i - 1));
    d = abs(out(i + 1) - out(i)) + abs(out(2) - out(1));
    if c == min([c g d])
        temp(2*i - 1) = out(m + i) + (out(i) + (out(i - 1) - out(1))/8);
    elseif g == min([c g d])
        temp(2*i - 1) = out(m + i) + (5*out(i)/8 + out(1)/2 - out(2)/8);
    else
        temp(2*i - 1) = out(m + i) + (11*out(i)/8 - out(i - 1)/2 + out(i - 2)/8);
    end
    
    for i = 3:(m - 2)
        c = abs(out(i - 1) - out(i)) + abs(out(i) - out(i + 1));
        g = abs(out(i - 1) - out(i - 2)) + abs(out(i) - out(i - 1));
        d = abs(out(i + 1) - out(i)) + abs(out(i + 2) - out(i + 1));
        if c == min([c g d])
            temp(2*i - 1) = out(m + i) + (out(i) + (out(i - 1) - out(i + 1))/8);
        elseif g == min([c g d])
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

