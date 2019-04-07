clear all;
close all;
clc;

% This script tests myfft.m for correctness
tic

fprintf('Test 1: Gaussian ...');
sample = exp(-linspace(-4,4,16).^2);
assert(all(abs(myfft(sample) - fft(sample))<1e-10));
fprintf('\tpassed\n');

fprintf('Test 1.1: Gaussian with a different Matlab dimension ...');
sample = exp(-linspace(-4,4,16).^2)';
assert(all(abs(myfft(sample) - fft(sample))<1e-10));
fprintf('\tpassed\n');

fprintf('Test 1.2: Gaussian complex ...');
sample = 1i*exp(-linspace(-4,4,128).^2);
assert(all(abs(myfft(sample) - fft(sample))<1e-10));
fprintf('\tpassed\n');

fprintf('Test 2: sawtooth ...');
sample = linspace(-1,1,128);
assert(all(abs(myfft(sample) - fft(sample))<1e-10));
fprintf('\tpassed\n');

fprintf('Test 3: sin and sin2 ...');
sample = sin(linspace(-pi,pi,128));
assert(all(abs(myfft(sample) - fft(sample))<1e-10));
sample = sin(2*linspace(-pi,pi,128));
assert(all(abs(myfft(sample) - fft(sample))<1e-10));
fprintf('\tpassed\n');

fprintf('Test 4: 2D input\n');
sample = sin(linspace(-pi,pi,10));
sample = [sample, sample];
testCase = matlab.unittest.TestCase.forInteractiveUse;
verifyError(testCase,@() dft(sample),'MATLAB:UndefinedFunction');

fprintf('Test 5: Power of 2\n');
sample = sin(linspace(-pi,pi,100));
sample = [sample, sample];
testCase = matlab.unittest.TestCase.forInteractiveUse;
verifyError(testCase,@() dft(sample),'MATLAB:UndefinedFunction');

fprintf('All tests passed!\n');
timeElapsed = toc;
fprintf('Time to pass all the tests : ')
disp(timeElapsed)
