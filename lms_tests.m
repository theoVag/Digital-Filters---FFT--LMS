
clear all
close all

n = 100000;
trials = 10;
sigma = 0.57;

v  = zeros(trials,  n);
for t = 1:trials
    v(t,:) = sqrt(sigma)*randn(n,1); v(t,:) = v(t,:) - mean(v(t,:));
end

uin = zeros(trials, n);
a = 0.34;
for t = 1:trials
   
    uin(t,1) = v(t, 1);
    for i=2:n
        uin(t,i) = -a * uin(t, i-1) + v(t, i);
    end
    
end

d = plant(uin);

%% adaptation
mu = 0.0005;
M = 2^10; % taps
mu_unconstrained=0.0004;
%%Block LMS algorithms
J = zeros(n, 1);
Jm = J; Js2=J;Js3=J;
for t=1:trials
  
    
  %% nested loops version
  y = zeros(n,1);
  e = zeros(n,1);
  u = uin(t,:)';
  
  w = zeros(M, 1);

  for i = 2*M:M:n
    g = zeros(M,1);
    for j=M-1:-1:0
      r = i-j;
      
      y(r) = w'*u((r-M+1):r);
      e(r) = d(t, r) - y(r);
      
      g = g + u((r-M+1):(r))*e(r);
    end
    w = w + mu*g;
    
  end
  J = J + e.^2;
  
  %% loop + matrix version
  ym = zeros(n,1);
  em = zeros(n,1);
  wm = zeros(M,1);
  
  for i = 2*M:M:n
    
    S = toeplitz(u((i-M+1):i), u((i-M+1):-1:(i-2*M+2)));
    ym(i-M+1:i) = S * wm;
    em(i-M+1:i) = d(t, i-M+1:i) - ym(i-M+1:i)';
    
    g = S' * em(i-M+1:i);
    wm = wm + mu * g;
  end
  Jm = Jm + em.^2;
  
%% Constrained Fast lms (Frequency Domain)
  ys = zeros(n,1);
  es = zeros(n,1);
  ws = zeros(M,1);
  W = zeros(2*M,1);

 for i = 2*M:M:n
    
    v = fft(u(i-2*M+1:i));
    C = ifft(W .* v);
    ys(i-M+1:i) = C(M+1:end);
    es(i-M+1:i) = d(t, i-M+1:i) - ys(i-M+1:i)';
    
    g = ifft(fft([zeros(M,1);es(i-M+1:i)]) .* conj(v));
    g=g(1:M);
    g=fft([g; zeros(M,1)]);
    W = W + mu * g;
  end

  Js2 = Js2 + es.^2;

  %% unconstrained fast lms (Frequency Domain)

  ys = zeros(n,1);
  es = zeros(n,1);
  ws = zeros(M,1);
  W = zeros(2*M,1);

 for i = 2*M:M:n
    
    v = fft(u(i-2*M+1:i));
    C = ifft(W.* v);
    ys(i-M+1:i) = C(M+1:end);
    es(i-M+1:i) = d(t, i-M+1:i) - ys(i-M+1:i)';
    
    g = fft([zeros(M,1);es(i-M+1:i)]) .* conj(v);

    W = W + mu_unconstrained *g ;

  end

  Js3 = Js3 + es.^2;
  
end

J = J / trials;
Jm = Jm / trials;
Js2 = Js2 / trials;
Js3 = Js3 / trials;

%%plots
figure(1)
plot([J - Jm])
xlabel('Time steps')
ylabel('error');
title('Block LMS and FAST LMS error difference')

figure(3)
plot([J - Js2])
xlabel('Time steps')
ylabel('error');
title('Block LMS and constrained FAST LMS(Frequency Domain) error difference')

figure(4)
plot([J - Js3])
xlabel('Time steps')
ylabel('error');
title('Block LMS and unconstrained FAST LMS(Frequency Domain) error difference')

figure(5)
plot([J Jm Js2 Js3])
xlabel('Time steps')
ylabel('error');
title('Comparison of LMS')

figure(9)
plot([J])
xlabel('Time steps')
ylabel('error');
title('2 nested Loops')

figure(6)
plot([Jm])
xlabel('Time steps')
ylabel('error');
title('Loop + matrix LMS')

figure(7)
plot([Js2])
xlabel('Time steps')
ylabel('error');
title('Constrained Fast LMS - Frequency Domain')

figure(8)
plot([Js3])
xlabel('Time steps')
ylabel('error');
title('Unconstrained Fast LMS - Frequency Domain')

fprintf('J-Jm error : %e\n', norm(J - Jm))
fprintf('J-JS2 error : %e\n', norm(J - Js2))
