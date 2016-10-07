function thress = generate_thresholds(n)
% Generate array of threshold values for testing
% output
%       (array)     thress  array of size n
% input
%       (scalar)    n       number of threshold values
%
T = zeros(1, n);
for i = 1:n
    if i == 1
        T(i) = 0.05;
    elseif i == 2
        T(i) = 0.1;
    else
        T(i) = T(i - 1) + T(i - 2)/2;
    end
end
thress = T;