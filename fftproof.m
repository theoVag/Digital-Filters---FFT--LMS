% "Proof" by MATLAB
% A simple technique to develop and verify the steps of a proof 
% using random data input
%
% N P P
% Cornell U 
% Sept 1992
%

clear

n = 8; % any even

% input 
x = randn(n,1) + 1i*randn(n,1);
% correct answer
ys = fft(x);

% root of unity
w = @(n,e) exp(-2*pi*1i.*e/n);

k = (0:n-1)';

% DFT proof steps
y = zeros(n,1);

for j = 0:n-1
  y(j +1) = sum(w(n,j*k) .* x(k +1));
end

fprintf('DFT : %e\n', norm(y - ys))

% split output top bottom
y = zeros(n,1);

for j = 0:n/2-1
  y(j +1) = sum(w(n,j*k) .* x(k +1));
end
for j = n/2:n-1
  y(j +1) = sum(w(n,j*k) .* x(k +1));
end

fprintf('split output top bottom : %e\n', norm(y - ys))

% split input even odd
y = zeros(n,1);

k = (0:n/2-1)';
for j = 0:n/2-1
  y(j +1) = sum(w(n,j*2*k) .* x(2*k +1)) + sum(w(n,j*(2*k+1)) .* x(2*k+1 +1));
end
for j = n/2:n-1
  y(j +1) = sum(w(n,j*2*k) .* x(2*k +1)) + sum(w(n,j*(2*k+1)) .* x(2*k+1 +1));
end

fprintf('split input even odd : %e\n', norm(y - ys))

% apply w identities
% etc
% ...
% ...

% to complete the proof

fe = fft(x((0:2:n-1) +1));
fo = fft(x((1:2:n-1) +1));

wfo = w(n,(0:n/2-1)') .* fo; 

y = [fe + wfo; fe - wfo];

%%fft proof2
sfe=zeros(n/2,1);
sfo=zeros(n/2,1);
for k=0:n/2-1
    for t = 0:n/2-1
      wn=exp((-1i)*2*2*pi*k*t/n);
      sfe(k+1)=sfe(k+1) + x(2*t+2)*wn;
      sfo(k+1)=sfo(k+1) + x(2*t+1)*wn;
    end 
end
w2 = @(n,e)exp(-2*pi*1i.*e/n);
wfe2 = w2(n,(0:n/2-1)') .* sfe; 
y2 = [sfo + wfe2; sfo - wfe2];
fprintf('FFT algorithm : %e\n', norm(y2 - ys))