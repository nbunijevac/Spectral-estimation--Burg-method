%% Zadatak 2.
clear all;
close all;
clc;
%% Ucitavanje signala i deljenje sekvence

[y, Fs]= audioread('nbA.wav');
x = y(4001:12000);
x = x';
t = 1; % [s]

tx = 0.02; % [s] 
x_20 = zeros(t/tx,tx*Fs);

for i = 1:50
    x_20(i,:) = x(((i-1)*160+1):i*160);
end

%% biranje reda modela

N = 160;
p = 1:10; 
FPE = zeros(1, length(p));
AIC = zeros(1,length(p));
sigma_sq = zeros(1, length(p));

figure(500)
hold all
for i = 1:50
for k = 1:length(p)
    [a, sigma_sq(k)] = m_covar(x_20(i,:), p(k));   
    AIC(k) = N*log(sigma_sq(k)) + 2*p(k); 
end
    plot(p, AIC)
    xlabel('red modela p')
    ylabel('AIC(p)')
    grid on
    title('AIC')
end
hold off

%% Estimacija SGS Burgovom metodom

br_sekv = 5;
sekv = sort(randi(t/tx,br_sekv,1));
maxp = 10; 
f = linspace(0,0.5,1000);
p = [2 3 4 5 6 7 8 9 10];

for j = 1:length(p)
figure(j)
    hold all
    for i = 1:length(sekv)
        [A, var] = burgova(x_20(sekv(i),:),maxp);
        a_1 = A(p(j)+1, 1:p(j)+1);

        a_1_pom = zeros(1,length(f));
        for k = 0:p(j)
            a_1_pom = a_1_pom + a_1(k+1)*exp(-1i*2*k*pi*f);
        end
        Pxx_1 = var(p(j)+1)*(abs(1./a_1_pom)).^2;
        Pxx_1 = Pxx_1/max(Pxx_1);

        [Pxx_2, W] = pburg(x_20(sekv(i),:),p(j),2*pi*f);
        Pxx_2 = Pxx_2/max(Pxx_2);

        plot(f,10*log10(Pxx_1))
        plot(f,10*log10(Pxx_2),'--','LineWidth',2)
    end
    hold off
    xlabel('f [1/odb]')
    ylabel('SGS [dB]')
    title(['SGS procesa x[n] dobijena Burgovom metodom za red modela ',...
          'p=',num2str(p(j))])
    legend([num2str(sekv(1)), '. sekvenca'],...
           [num2str(sekv(1)), '. sekvenca (pburg)'],...
           [num2str(sekv(2)), '. sekvenca'],...
           [num2str(sekv(2)), '. sekvenca (pburg)'],...
           [num2str(sekv(3)), '. sekvenca'],...
           [num2str(sekv(3)), '. sekvenca (pburg)'],...
           [num2str(sekv(4)), '. sekvenca'],...
           [num2str(sekv(4)), '. sekvenca (pburg)'],...
           [num2str(sekv(5)), '. sekvenca'],...
           [num2str(sekv(5)), '. sekvenca (pburg)'],...
           'Location','SouthWest')
    grid on
end    

%% Promena estimiranih parametara modela u vremenu
p_opt = 2; %za njega sam primetila najveci pad kad sam proveravala
a_2 = zeros(t/tx,p_opt+1);
t_osa = tx:tx:t;
for i = 1:t/tx
    [A, var] = burgova(x_20(i,:),p_opt);
    a_2(i,:) = A(p_opt+1, 1:p_opt+1);
end
figure
    hold all
    plot(t_osa,a_2(:,1))
    plot(t_osa,a_2(:,2))
    plot(t_osa,a_2(:,3))
    xlim([tx max(t_osa)])
    hold off
    xlabel('t [s]'); ylabel('a[k]'); ylim([-2 1.25]) 
    title(['Promena estimiranih parametara a[k] za optimalan red p=',...
          num2str(p_opt)])
    legend('a[0]','a[1]','a[2]','Location','East')
    grid on