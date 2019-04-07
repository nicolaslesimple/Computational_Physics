function result = myfft_slow(inputarray)
% This funciton called myfft computes discrete Fourier transform with a radix 2.
%
% Arguments :
%   - input array (1D array complex, size 2?N, N>0): data to transform ;
%
% Returns : 
%   - 1D array complex, transformed data of the same shape as an input array .

p=nextpow2(length(inputarray)); % the function nextpow2 returns the smallest power of two that is greater than or equal to the absolute value of A. (That is, p that satisfies 2^p >= abs(A)).
inputarray=[inputarray zeros(1,(2^p)-length(inputarray))]; % we add zeros to the inputarray if the size is not of lenght n, where n should be the power of 2
N=length(inputarray); 
size_decrease=N/2;
for step=1:log2(N)
    for j=0:(N/(2^(step-1))):(N-1)
        for n=0:(size_decrease-1)
            position=n+j+1;
            power_part=(2^(step-1))*n; % power part of the complex multiplier
            w = exp((-1i)*(2*pi)*power_part/N); % complex multiplier
            a = inputarray(position) + inputarray(position + size_decrease);
            b = (inputarray(position) - inputarray(position + size_decrease)).* w;
            inputarray(position)=a;   % saving computation of the 1-st part
            inputarray(position+size_decrease)=b;% saving computation of the 2-nd part
        end
    end
size_decrease=size_decrease/2; % computing the next "Half" value
end
result=bitrevorder(inputarray);   % performing bit-reverse operation and returning the result from function
return