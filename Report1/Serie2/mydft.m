
function result = mydft(input_array)
% This function called mydft computes the discrete Fourier transform of a 1D signal
% Arguments :
%   - input_array (1D array complex) : data to transform ;
% Returns : 
%   - 1D array complex, transformed data of the same shape as an input array .

% We check the format of the input to allow line or column vector
input_size = size(input_array);
if input_size(1)> 1
    input_array = transpose(input_array);
end

% We check if the input array is 1D. If not, an error is throw.
nD = size(input_size);
if (min(input_size)>1)~=0 || (nD(2)>2)~=0
    msg = 'The input array is in 2D or more.';
    error(msg);
end

% We first declare the variables needed to solve the problem
N = length(input_array);
result = zeros(1, N);
sum=0;
% We apply the calculation algorithm of the DFT
for k=1:N
    for l=1:N
        sum=sum+input_array(l)*exp(-2*pi*1i*(l-1)*(k-1)/N);
    end
    result(k)=sum;
sum=0;
end

% If the input was not in the optimal format, we change it during the
% process and thus there we give him back is original form.
if input_size(1)> 1
    result = transpose(result);
end
return


