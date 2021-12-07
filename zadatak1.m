%% Zadatak 1.
clear all;
close all;
clc;

%% Inicijalizacija i ucitavanje signala

I = mod((2017*(0+0+1+7)),4); 
J = mod(17,3); 
K = mod((17*(2+0+1+7)),5); 
dat = load(['data',num2str(I),'.mat']);
x = dat.x;
N = length(x);
a_an = dat.a;

%% Odredjivanje analiticke SGS procesa x[n] za date parametre AR modela

N_an = 256;
u_an = randn(1,N_an);
p_an  = length(a_an)-1;
x_an = filter(1, a_an, u_an);

[H_an, W_an] = freqz(1, a_an);
f_an = W_an/(2*pi); 

a_an_pom = zeros(1,length(f_an));

for k = 0:p_an
    a_an_pom = a_an_pom + a_an(k+1)*exp(-1i*2*k*pi*f_an);
end
Pxx_an = (abs(1./a_an_pom)).^2;
Pxx_an = Pxx_an/max(Pxx_an);

%% Izbor optimalnog reda modela
p = 1:10; 
FPE = zeros(1, length(p));
AIC = zeros(1,length(p));
sigma_sq = zeros(1, length(p));

for k = 1:length(p)
    [a, sigma_sq(k)] = m_covar(x, p(k));
    FPE(k) = ((N + p(k))/(N - p(k)))*sigma_sq(k);   
    AIC(k) = N*log(sigma_sq(k)) + 2*p(k); 
end
figure
    subplot(211)
    plot(p, FPE);
    xlabel('red modela p')
    ylabel('FPE(p)')
    grid on
    title('Final prediction error')
    
    subplot(212)
    plot(p, AIC)
    xlabel('red modela p')
    ylabel('AIC(p)')
    grid on
    title('Akaike information criterion')


i_FPE_opt = find(FPE == min(FPE)); 
disp(['Optimalni red modela dobijen FPE metodom : ',...
      num2str(p(i_FPE_opt))]);
i_AIC_opt = find(AIC == min(AIC)); 
disp(['Optimalni red modela dobijen AIC metodom : ',...
      num2str(p(i_AIC_opt))]);

%%
p_opt = 3; %prva vrednost na kojoj opadne primetno vrednost izraza koji se minimizuje
[a_1, sigma_sq_1] = m_covar(x, p_opt);

a_1_pom = zeros(1,length(f_an));

for k = 0:p_opt
    a_1_pom = a_1_pom + a_1(k+1)*exp(-1i*2*k*pi*f_an);
end
Pxx_1 = sigma_sq_1*(abs(1./a_1_pom)).^2;
Pxx_1 = Pxx_1/max(Pxx_1);

[Pxx_2, W] = pmcov(x, p_opt, 2*pi*f_an);
Pxx_2 = Pxx_2/max(Pxx_2);

figure
    hold all
    plot(f_an,10*log10(Pxx_1))
    plot(f_an,10*log10(Pxx_2))
    plot(f_an,10*log10(Pxx_an))
    hold off
    xlabel('f [1/odb]')
    ylabel('SGS [dB]')
    legend('modifikovana kovarijantna','ugradjena m. kovarijantna','analiticka');
    grid on
    
%% Autokorelaciona funkcija procene pobudnog procesa
u_est = zeros(1, N); 
for n = 1:N
    if n <= p_opt
        u_est(n) = sum(a_1(1:n).*x(n:(-1):n-mod(n-1, p_opt)));
    else
        u_est(n) = sum(a_1.*x(n:(-1):n-p_opt)); 
    end
end

ruu = akf(u_est);

figure
    stem(-N+1:N-1, ruu)
    xlim([-10 10]) 
    xlabel('k')
    ylabel('r_{uu}[k]')
    title('Nepomerena AKF procene pobudnog procesa'); 
    grid on
    
%% Polozaj nula estimiranog polinoma A(z)
% polozaj polova funkcije prenosa AR modelovanog filtra
z = tf('z');
A_1 = 0;
for k = 0:p_opt
    A_1 = A_1 + a_1(k+1)*z^(-k);
end

figure
    pzmap(A_1);
    axis equal

