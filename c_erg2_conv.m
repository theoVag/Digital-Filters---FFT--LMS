
clear
n1 = 1024; 
n2 = 1024;
% input 
x = randn(n1,1) + 1i*randn(n1,1);
y = randn(n2,1) + 1i*randn(n2,1);

%a' convolution with conv
g1=conv(x,y);

%b' convolution with Yx
Y=toeplitz([y;zeros(n2-1,1)],[y(1);zeros(n2-1,1)]);
g2=Y*x;
fprintf('conv b error : %e\n', norm(g1 - g2))

%c' convolution with Cx
yr=[y(1) zeros(1,n2-1) y(end:-1:2).'];
%C=circulant(yr,1);
C=gallery('circul',yr);
g3=C*[x ; zeros(n1-1,1)];
fprintf('conv c error : %e\n', norm(g1 - g3))

%d' convolution with fourier
yp=[y; zeros(n2-1,1)];
xp=[x; zeros(n1-1,1)];
g4=ifft(fft(yp).*fft(xp));
fprintf('conv d error : %e\n', norm(g1 - g4))