function tp_seance2()
clc;
clear all;
close all;
exercice = input('Choose an exercice to run\n 1 / 2 / 6 / 7 / 8\n', 's');
if (exercice == '1')
    exercice1();
elseif (exercice == '2')
    exercice2();
elseif (exercice == '6')
    exercice6();
elseif (exercice == '7')
    exercice7();
elseif (exercice == '8')
    exercice8();
else
    return;
end

function exercice1()
disp('Exercice 1: reconstruction de signaux en utilisant la transformation de Haar 1d avec seuillage:');
T = 5;
str = sprintf('Avec seuillage, valeur de seuil = %d', T);
disp(str);
in1 = generate_1d_increment_signal();
in2 = generate_1d_image_signal();
in3 = generate_1d_trigo_signal();
disp('Erreurs de signal cumulatif avant et apres transformation: ');
[ME MSE] = reconstruct_haar_1d(in1, T)
disp('Erreurs de signal d''image avant et apres transformation: ');
[ME MSE] = reconstruct_haar_1d(in2, T)
disp('Erreurs de signal trigonométrique avant et apres transformation: ');
[ME MSE] = reconstruct_haar_1d(in3, T)

function exercice2()
disp('Exercice 2: la qualite de reconstruction par rapport aux seuils utilisees.');
in1 = generate_1d_increment_signal();
in2 = generate_1d_image_signal();
in3 = generate_1d_trigo_signal();

n = 20; % nombre de seuils
errors1 = zeros(n, 2);
errors2 = zeros(n, 2);
errors3 = zeros(n, 2);
T = zeros(1, n);
for i = 1:n
    T(i) = 0.2 * i;
    [errors1(i, 1), errors1(i, 2)] = reconstruct_haar_1d(in1, T(i));
    [errors2(i, 1), errors2(i, 2)] = reconstruct_haar_1d(in2, T(i));
    [errors3(i, 1), errors3(i, 2)] = reconstruct_haar_1d(in3, T(i));
end

figure;
plot(T, errors1(:,1), '--rs', 'LineWidth', 2,...
                'MarkerEdgeColor', 'k',...
                'MarkerFaceColor', 'g',...
                'MarkerSize', 5);
hold on;
plot(T, errors1(:,2), '--ro', 'LineWidth', 2,...
                'MarkerEdgeColor', 'k',...
                'MarkerFaceColor', 'g',...
                'MarkerSize', 5);
plot(T, errors2(:,1), '--bs', 'LineWidth', 2,...
                'MarkerEdgeColor', 'k',...
                'MarkerFaceColor', 'g',...
                'MarkerSize', 5);
plot(T, errors2(:,2), '--bo', 'LineWidth', 2,...
                'MarkerEdgeColor', 'k',...
                'MarkerFaceColor', 'g',...
                'MarkerSize', 5);
plot(T, errors3(:,1), '--gs', 'LineWidth', 2,...
                'MarkerEdgeColor', 'k',...
                'MarkerFaceColor', 'g',...
                'MarkerSize', 5);
plot(T, errors3(:,2), '--go', 'LineWidth', 2,...
                'MarkerEdgeColor', 'k',...
                'MarkerFaceColor', 'g',...
                'MarkerSize', 5);
legend('ME, incremental signal', 'MSE, incremental signal',...
        'ME, image vector', 'MSE, image vector',...
        'ME, sinusoidal signal', 'MSE, sinusoidal signal')
xlabel('Threshold');
ylabel('Errors');

function [ME MSE] = reconstruct_haar_1d(in, threshold)
out = haar_1d(in);
outT = thresholding(out, threshold);
inT = inverse_haar_1d(outT);
[ME MSE] = compute_errors(in, inT);

function exercice6()
disp('Exercice 6: errors d''signaux des images avant et apres la transformation 2D de Haar');
[in, map] = generate_2d_image_signal();
out = haar_2d(in);
inT = haar_2d_inverse(out);
[ME, MSE] = compute_errors(inT, in)
figure;
subplot(1,2,1);imshow(in,map);title('image original');
subplot(1,2,2);imshow(inT,map);title('image apres reconstruction');

function [me, mse, inT] = reconstruct_haar_2d(in, threshold) %norm_exercice7
% Generate reconstructed image and errors
% Output:
%       (scalar)    me              mean absolute error
%       (scalar)    mse             mean square error
%       (matrix)    inT             reconstructed image
% Input:
%       (matrix)    in              input image
%       (scalar)    threshold       extreme value for thresholding
%
out = haar_2d(in);
outT = thresholding2d(out, threshold);
inT = haar_2d_inverse(outT);
[me, mse] = compute_errors(inT, in);

function exercice7()
T = 3;
disp('Exercice 7: Processus d''image avec transformationd 2d de Haar');
str = sprintf('Et seuillage, avec seuil = %d', T);
disp(str);
[in, map] = generate_2d_image_signal();
[ME, MSE, inT] = reconstruct_haar_2d(in, T);
ME
MSE
figure;
subplot(1,2,1);imshow(in, map);title('image original');
subplot(1,2,2);imshow(inT, map);title(['image apres reconstruction, threshold = ', num2str(T)]);

    
function exercice8()
disp('Exercice 8: Qualite de reconstruction par rapport aux seuils utilisees');
in = generate_2d_image_signal();
n = 20; % nombre de threshold
T = generate_thresholds(n); % size 1 x n
me = zeros(1, n);
mse = zeros(1, n);
for i = 1:n
    [me(i), mse(i), ~] = reconstruct_haar_2d(in, T(i));
end

figure
hold on;
plot(T, me, '--bs', 'LineWidth', 2,...
                'MarkerEdgeColor', 'k',...
                'MarkerFaceColor', 'g',...
                'MarkerSize', 5);
plot(T, mse, '--bo', 'LineWidth', 2,...
                'MarkerEdgeColor', 'k',...
                'MarkerFaceColor', 'g',...
                'MarkerSize', 5);
legend('ME, image signal', 'MSE, image signal');
xlabel('Threshold')
ylabel('Errors')

function out = haar_2d(in, use_while_loop)
if ~exist('use_while_loop', 'var') || isempty(use_while_loop)
    use_while_loop = true;
end
% size x y are identical
size_in = size(in);
out = zeros(size_in);
temp = zeros(size_in);
if (use_while_loop)
    k = 1;
    %Transformation on rows
    while (k <= size_in(1))
        temp(k,:) = haar_1d(in(k,:), true);
        k = k + 1;
    end
    k = 1;
    %Transformation on columns
    while (k <= size_in(2))
        out(:,k) = haar_1d(temp(:,k), true);
        k = k + 1;
    end
else
    %Transformation on rows
    for i = 1:size_in(1)
        temp(i,:) = haar_1d(in(i,:), false);
    end
    %Transformation on columns
    for j = 1:size_in(2)
        out(:,j) = haar_1d(temp(:,j), false);
        k = k + 1;
    end
end

function out = haar_2d_inverse(in, use_while_loop)
if ~exist('use_while_loop', 'var') || isempty(use_while_loop)
    use_while_loop = true;
end
% size x y are identical
size_in = size(in);
out = zeros(size_in);
temp = zeros(size_in);
if (use_while_loop)
    k = 1;
    %Transformation on rows
    while (k <= size_in(1))
        temp(k,:) = inverse_haar_1d(in(k,:), true);
        k = k + 1;
    end
    k = 1;
    %Transformation on columns
    while (k <= size_in(2))
        out(:,k) = inverse_haar_1d(temp(:,k), true);
        k = k + 1;
    end
else
    %Transformation on rows
    for i = 1:size_in(1)
        temp(i,:) = inverse_haar_1d(in(i,:), false);
    end
    %Transformation on columns
    for j = 1:size_in(2)
        out(:,j) = inverse_haar_1d(temp(:,j), false);
        k = k + 1;
    end
end

function out = haar_1d(in, use_while_loop)
% Description
%       1D Haar multilevel transform
%
% Synopsis
%		out = haar_1d(in, [use_while_loop])
%
% Inputs ([]s are optional)
%		(vector)    in      input signal of arbitrary length,
%                           power of 2 prefered
%		(boolean)   [use_while_loop]
%                           if true, calculate using while loops,
%                           otherwise, code using for loops are in action
%
% Outputs ([]s are optional)
%		(vector)    out     output signal of normalized length (power of 2).
%
% Examples
%		img = imread('zelda.bmp');
%		in = img(:,36);
%		out = haar_1d(in);
%       out = haar_1d(in, false);
%
% Requirements
%		haar_1d_one_level_while_loop
%       haar_1d_one_level_for_loop

if ~exist('use_while_loop', 'var') || isempty(use_while_loop)
    use_while_loop = true;
end
temp = in;

signal_length = length(temp);
scale = log2(signal_length);
additional_length = 0;

% Begin normalization
% Adding supplementary item to input signal so that signal length is a power of 2
while (scale - floor(scale) ~= 0)
    additional_length = additional_length + 1;
    scale = log2(signal_length + additional_length);
end

if (additional_length > 0)
    additional_length_int = floor(additional_length / signal_length);
    additional_length_remain = additional_length - additional_length_int;
    if (additional_length_int > 0)
        for i = 1:additional_length_int
            temp = [temp; in];
        end
    end
    for j = 1:additional_length_remain
        temp = [temp; in(j)];
    end
    signal_length = signal_length + additional_length_int + additional_length_remain;
end
% End normalization

m = signal_length;
out = temp;

if (use_while_loop)
    % execute the transformation using while loops
    while (m > 1)
        out(1:m) = haar_1d_one_level_while_loop(temp(1:m));
        temp(1:m) = out(1:m);
        m = m / 2;
    end
else
    % execute the transformation using for loops
    e = log2(m);
    for i = 1:e
        out(1:m) = haar_1d_one_level_for_loop(temp(1:m));
        temp(1:m) = out(1:m);
        m = m / 2;
    end
end

function out = haar_1d_one_level_while_loop(in)
% Description
%       1D Haar one-level transform using while loops
%
% Synopsis
%		out = haar_1d_one_level_while_loop(in)
%
% Inputs ([]s are optional)
%		(vector)    in      input signal of even length
%
% Outputs ([]s are optional)
%		(vector)    out     result of transformation
%
% Examples
%		img = imread('zelda.bmp');
%		in = img(:,36);
%		out = haar_1d_one_level_while_loop(in);
%
signal_length = length(in); % doit etre pair
out = in;
m = 1;
while (m <= signal_length / 2)
    out(m) = (in(2*m - 1) + in(2*m)) / 2;
    out(m + signal_length / 2) = (in(2*m - 1) - in(2*m)) / 2;
    m = m + 1;
end

% transformation de Haar dans un bas, avec la boucle for
function out = haar_1d_one_level_for_loop(in)
% Description
%       1D Haar one-level transform using for loops
%
% Synopsis
%		out = haar_1d_one_level_for_loop(in)
%
% Inputs ([]s are optional)
%		(vector)    in      input signal of even length
%
% Outputs ([]s are optional)
%		(vector)    out     result of transformation
%
% Examples
%		img = imread('zelda.bmp');
%		in = img(:,36);
%		out = haar_1d_one_level_for_loop(in);
%
signal_length = length(in); % doit etre pair
out = in;
for m = 1 : signal_length / 2
    out(m) = (in(2*m - 1) + in(2*m)) / 2;
    out(m + signal_length / 2) = (in(2*m - 1) - in(2*m)) / 2;
end

% Transformee inverse de la transformee de Haar multiechelles
% en utilisant le boucle while et for
function out = inverse_haar_1d(in, use_while_loop)
% Description
%       1D Haar inverse multilevel transform
%
% Synopsis
%		out = inverse_haar_1d(in, [use_while_loop])
%
% Inputs ([]s are optional)
%		(vector)    in      input signal of length power of 2
%		(boolean)   [use_while_loop]
%                           if true, calculate using while loops,
%                           otherwise, code using for loops are in action
%
% Outputs ([]s are optional)
%		(vector)    out     output signal of length power of 2
%
% Examples
%		img = imread('zelda.bmp');
%		in = img(:,36);
%		out = haar_1d(in);
%       inverse_out = inverse_haar_1d(out);
%       inverse_out = inverse_haar_1d(out, false);
%
% Requirements
%		inverse_haar_1d_one_level_for_loop
%       inverse_haar_1d_one_level_while_loop
if ~exist('use_while_loop', 'var') || isempty(use_while_loop)
    use_while_loop = true;
end
signal_length = length(in);
m = 1;
temp = in;
out = temp;

if (use_while_loop)
    % executer la transformation
    while (m < signal_length)
        m = m * 2;
        temp(1:m) = out(1:m);
        out(1:m) = inverse_haar_1d_one_level_while_loop(temp(1:m));
    end
else
    % executer la transformation, avec la boucle for
    e = log2(signal_length);
    for i = 1:e
        m = m * 2;
        temp(1:m) = out(1:m);
        out(1:m) = inverse_haar_1d_one_level_for_loop(temp(1:m));
    end
end

function out = inverse_haar_1d_one_level_while_loop(in)
% Description
%       1D Haar inverse one-level transform using while loops
%
% Synopsis
%		out = inverse_haar_1d_one_level_while_loop(in)
%
% Inputs ([]s are optional)
%		(vector)    in      input signal of even length
%
% Outputs ([]s are optional)
%		(vector)    out     result of transformation
%
% Examples
%		img = imread('zelda.bmp');
%		in = img(:,36);
%		out = inverse_haar_1d_one_level_while_loop(in);
%
signal_length = length(in);
if (mod(signal_length, 2) == 1)
    disp('Input signal should be of even length');
    return
end
out = in;
m = 1;
while (m <= signal_length / 2)
    out(2*m - 1) = in(m) + in(m + signal_length / 2);
    out(2*m) = in(m) - in(m + signal_length / 2);
    m = m + 1;
end

function out = inverse_haar_1d_one_level_for_loop(in)
% Description
%       1D Haar inverse one-level transform using for loops
%
% Synopsis
%		out = inverse_haar_1d_one_level_for_loop(in)
%
% Inputs ([]s are optional)
%		(vector)    in      input signal of even length
%
% Outputs ([]s are optional)
%		(vector)    out     result of transformation
%
% Examples
%		img = imread('zelda.bmp');
%		in = img(:,36);
%		out = inverse_haar_1d_one_level_for_loop(in);
%
signal_length = length(in);
if (mod(signal_length, 2) == 1)
    disp('Input signal should be of even length');
    return
end
out = in;
for m = 1:signal_length / 2
    out(2*m - 1) = in(m) + in(m + signal_length / 2);
    out(2*m) = in(m) - in(m + signal_length / 2);
end