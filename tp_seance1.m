%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%           EXERCICE 1 & 3 & 5             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function tp_seance1()
clc;
clear all;
close all;
exercice = input('Choose an exercice to run\n 1 / 2 / 3 / 4 / 5\n', 's');
if (exercice == '1')
    exercice12();
elseif (exercice == '2')
    exercice12();
elseif (exercice == '3')
    exercice3();
elseif (exercice == '4')
    exercice4();
elseif (exercice == '5')
    exercice5();
else
    return;
end

function exercice12()
disp('Exercice 1 et 2: processus des signaux donnees (en utilisant la boucle while):');
in1 = generate_1d_increment_signal();
in2 = generate_1d_image_signal();
in3 = generate_1d_trigo_signal();
out1 = haar_1d(in1, true);
out2 = haar_1d(in2, true);
out3 = haar_1d(in3, true);
in1T = inverse_haar_1d(out1, true);
in2T = inverse_haar_1d(out2, true);
in3T = inverse_haar_1d(out3, true);
disp('Erreurs de signal cumulatif avant et apres transformation: ');
[ME, MSE] = compute_errors(in1T, in1)
disp('Erreurs de signal d''image avant et apres transformation: ');
[ME, MSE] = compute_errors(in2T, in2)
disp('Erreurs de signal trigonométrique avant et apres transformation: ');
[ME, MSE] = compute_errors(in3T, in3)

function exercice3()
disp('Exercice 3: processus de signal x1 - x3, en utilisant la boucle while');
in1 = generate_1d_increment_signal();
in3 = generate_1d_trigo_signal();
in = in1 - in3;
out = haar_1d(in, true);
inT = inverse_haar_1d(out, true);
disp('Erreurs de signal x1 - x3 avant et apres transformation: ');
[ME, MSE] = compute_errors(in, inT)

function exercice4()
disp('Exercice 4: processus des signaux donnees (en utilisant la boucle for):');
in1 = generate_1d_increment_signal();
in2 = generate_1d_image_signal();
in3 = generate_1d_trigo_signal();
out1 = haar_1d(in1, false);
out2 = haar_1d(in2, false);
out3 = haar_1d(in3, false);
in1T = inverse_haar_1d(out1, false);
in2T = inverse_haar_1d(out2, false);
in3T = inverse_haar_1d(out3, false);
disp('Erreurs de signal cumulatif avant et apres transformation: ');
[ME, MSE] = compute_errors(in1T, in1)
disp('Erreurs de signal d''image avant et apres transformation: ');
[ME, MSE] = compute_errors(in2T, in2)
disp('Erreurs de signal trigonométrique avant et apres transformation: ');
[ME, MSE] = compute_errors(in3T, in3)

function exercice5()
disp('Exercice 5: processus de signal x1 - x3, en utilisant la boucle for');
in1 = generate_1d_increment_signal();
in3 = generate_1d_trigo_signal();
in = in1 - in3;
out = haar_1d(in, false);
inT = inverse_haar_1d(out, false);
disp('Erreurs de signal x1 - x3 avant et apres transformation:');
[ME, MSE] = compute_errors(in, inT)

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