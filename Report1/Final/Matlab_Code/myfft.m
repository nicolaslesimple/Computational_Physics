function result = myfft(input_array)
  % This funciton called myfft computes discrete Fourier transform with a radix 2.
  %
  % Arguments :
  %   - input array (1D array complex, size 2N, N>0): data to transform ;
  %
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
  % We check if N is a power of two
  if mod(log2(input_size(2)),2)~=0 && mod(log2(input_size(2)),2)~=1
      msg = 'The number of data is not a power of two.';
      error(msg);
  end
  
  n = length(input_array);
  N = pow2(ceil(log2(n)));
  input_array = [input_array, zeros(1, N - n)]; % Create the new arry to store data
  result = myfft_module(input_array); % Call the module declared below
  result = result(1:n); % Fill up the results
  
  % If the input was not in the optimal format, we change it during the
  % process and thus there we give him back is original form.
  if input_size(1)> 1
    result = transpose(result);
  end
  
return

function result_module = myfft_module(input_array_module)
  % This functions allows the recursive call of myfft function
  % It uses another function called W()
  
  n = length(input_array_module); 
  % When n is 1 we have nearly finish
  if (n == 1)
    result_module = input_array_module;
  else
    % Odd and even number are calculate separtly to optimize the
    % computation time 
    f_even = input_array_module(1:2:n);
    f_odd = input_array_module(2:2:n);
    % The recursive part of the algorithm is made thank to the call to the
    % function myfft again
    X1 = myfft(f_even);
    X2 = myfft(f_odd).*Wn(n);  
    % Calculate and store the results
    F1 = X1 + X2;
    F2 = X1 - X2;
    result_module = [F1 F2];
  end
return

function w = Wn(n)
  % This function allows the calculation of specific exponential terms
  % added to the result
  m = n/2;
  w = exp(-2*pi*1i.*(0:1:m-1)/n);
return

