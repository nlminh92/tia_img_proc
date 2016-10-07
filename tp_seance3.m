function tp_seance3()
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
    exercice8();
else
    return;
end

function exercice12()
disp('Exercice 1 et 2: appliquer la transformation quadratique multi-echelle 1D');
in1 = generate_1d_increment_signal();
in2 = generate_1d_image_signal();
in3 = generate_1d_trigo_signal();

out1 = quadratic_1d(in1);
out2 = quadratic_1d(in2);
out3 = quadratic_1d(in3);

in1T = quadratic_1d_inverse(out1);
in2T = quadratic_1d_inverse(out2);
in3T = quadratic_1d_inverse(out3);

disp('Erreurs de signal cumulatif avant et apres transformation:');
[ME1 MSE1] = compute_errors(in1, in1T)
disp('Erreurs de signal d''image avant et apres transformation:');
[ME2 MSE2] = compute_errors(in2, in2T)
disp('Erreurs de signal trigonométrique avant et apres transformation:');
[ME3 MSE3] = compute_errors(in3, in3T)

function exercice5()
disp('Exercice 5: appliquer la transformation quadratique multi-echelle 2D');
[in, map] = generate_2d_image_signal();
out = quadratic_2d(in);
inT = quadratic_2d_inverse(out);
disp('Erreurs d''images avant et apres transformation quadratique:');
[ME MSE] = compute_errors(in, inT)
figure;
subplot(1,2,1);imshow(in, map);title('image original avant transformation quadratique');
subplot(1,2,2);imshow(inT, map);title('image apres reconstruction');

function exercice6()
disp('Exercice 6: appliquer la transformation quadratique multi-echelle 2D avec seuillage pour debruitage');
T = 10;
[in, map] = generate_2d_image_signal();
in = double(in);
inSize = size(in);
sigma = 0.1 * (max(in(:)) - min(in(:))); % noise level
inNoisy = in + sigma * randn(inSize(1), inSize(2)); % noisy image
outNoisy = quadratic_2d(inNoisy);
outFiltered = thresholding2d(outNoisy, T);
inFiltered = quadratic_2d_inverse(outFiltered);
disp(['Erreurs d''images avant et apres transformation quadratique + seuillage, seuil = ', num2str(T)]);
[ME MSE] = compute_errors(in, inFiltered)
figure;
subplot(1,3,1); imshow(in, map); title('original image');
subplot(1,3,2); imshow(inNoisy, map); title('noisy image');
subplot(1,3,3); imshow(inFiltered, map); title('reconstructed image');

function [me, mse, kept_count] = compute_transformation_errors(in, thres)
% Compute mean error, mean square error, and retained point count after
% thresholding
out = quadratic_2d(in);
[outFiltered, kept_count] = thresholding2d_count(out, thres);
inFiltered = quadratic_2d_inverse(outFiltered);
[me, mse] = compute_errors(in, inFiltered);

function exercice7()
disp('Exercice 7: Erreur en fonction du seuil.');
in = generate_2d_image_signal();
n = 20;
me = zeros(n, 1);
mse = zeros(n, 1);
T = generate_thresholds(n); % size 1 x n
for i = 1:n
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

function exercice8()
n = 20;
T = generate_thresholds(n);
nPoints = zeros(1, n);
me = zeros(n, 1);
mse = zeros(n, 1);
in = generate_2d_image_signal();
for i = 1:n
    [me(i), mse(i), nPoints(i)] = compute_transformation_errors(in, T(i));
end
figure
xlabel('Retained points')
ylabel('Errors')
plot(nPoints, me, '--bs', 'LineWidth', 2,...
                'MarkerEdgeColor', 'k',...
                'MarkerFaceColor', 'g',...
                'MarkerSize', 5);
hold on;
plot(nPoints, mse, '--bo', 'LineWidth', 2,...
                'MarkerEdgeColor', 'k',...
                'MarkerFaceColor', 'g',...
                'MarkerSize', 5);
legend('ME', 'MSE');
xlabel('Retained points');
ylabel('Errors');

function out = quadratic_2d(in)
% Implementation of 2D multiscale quadratic transform
% Input
%       (matrix)    in          input signal of length power-of-2
% Output
%       (matrix)    out         transformed signal of identical length
%
in_size = size(in);
temp = zeros(in_size);
out = zeros(in_size);
% apply quadratic transform to rows
for i = 1:in_size(1)
    temp(i,:) = quadratic_1d(in(i,:));
end
% apply quadratic transform to columns
for j = 1:in_size(2)
    out(:,j) = quadratic_1d(temp(:,j));
end

function out = quadratic_2d_inverse(in)
% Implementation of 2D inverse multiscale quadratic transform
% Input
%       (matrix)    in          input signal of length power-of-2
% Output
%       (matrix)    out         transformed signal of identical length
%
in_size = size(in);
temp = zeros(in_size);
out = zeros(in_size);
% apply inverse quadratic transform to rows
for i = 1:in_size(1)
    temp(i,:) = quadratic_1d_inverse(in(i,:));
end
% apply inverse quadratic transform to columns
for j = 1:in_size(2)
    out(:,j) = quadratic_1d_inverse(temp(:,j));
end

function out = quadratic_1d(in)
% Implementation of 1D multiscale quadratic transform
% Input
%       (vector)    in          input signal of length power-of-2
% Output
%       (vector)    out         transformed signal of identical length
%
in_length = length(in);
in_length_2 = log2(in_length);
out = in;
temp = in;
m = in_length;
for k = in_length_2:-1:2
    m = m/2;
    for i = 1:m
        temp(i) = 0.5*(out(2*i - 1) + out(2*i));
    end
    % Border values, i = 1, i = m
    i = 1;
    temp(m + i) = out(2*i - 1) - (temp(i) + (temp(m) - temp(i + 1))/8);
    i = m;
    temp(m + i) = out(2*i - 1) - (temp(i) + (temp(i - 1) - temp(1))/8);
    for i = 2:(m - 1)
        temp(m + i) = out(2*i - 1) - (temp(i) + (temp(i - 1) - temp(i + 1))/8);
    end
    out = temp;
end

function out = quadratic_1d_inverse(in)
% Implementation of 1D inverse multiscale quadratic transform
% Input
%       (vector)    in          input signal of length power-of-2
% Output
%       (vector)    out         transformed signal of identical length
%
in_length = length(in);
in_length_2 = log2(in_length);
out = in;
temp = in;
m = 2;
for k = 2:1:in_length_2
    % Border values, i = 1, i = m
    i = 1; 
    temp(2*i - 1) = out(i) + (out(m) - out(i + 1))/8 + out(m + i);
    i = m;
    temp(2*i - 1) = out(i) + (out(i - 1) - out(1))/8 + out(m + i);
    for i=2:(m - 1)
        temp(2*i - 1) = out(i) + (out(i - 1) - out(i + 1))/8 + out(m + i);
    end
    for i = 1:m
       temp(2*i) =  2*out(i) - temp(2*i - 1);
    end
    out = temp;
    m = m * 2;
end
