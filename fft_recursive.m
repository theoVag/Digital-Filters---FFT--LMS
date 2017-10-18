function X = fft_recursive(x)

%only works if N = 2^k
%x1=even , x2=odd
N = length(x);
x1 = x(1:2:end); %x even (matlab starts counting from 1)
x2 = x(2:2:end); %x odd 
if N>1 
    X1 = fft_recursive(x1); 
    X2 = fft_recursive(x2); 
    X = zeros(N,1);
    Wn = exp(-1i*2*pi*((0:N/2-1)')/N);
    temp = Wn .* X2;
    X = [(X1 + temp);(X1 -temp)];

else
    if N==1 
        X=x;
    end
end
end
