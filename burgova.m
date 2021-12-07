function [A, var] = burgova(x, maxp)

    N = length(x);
    var = zeros(1, maxp+1); 
    Kk = zeros(1, maxp);
    A = zeros(maxp + 1, maxp + 1); 

    eb = x(1:N-1);
    ef = x(2:N);

    % Inicijalizacija
    var(1) = sum(abs(x).^2)/N;
    A(:, 1) = 1;

    for k =1:maxp
     % Koeficijent refleksije
     Kk(k) = (-2*sum(ef.*conj(eb)))./sum(abs(ef).^2 + abs(eb).^2);

     % Greska predikcije
     var(k+1) = (1 - abs(Kk(k))^2)*var(k);

     % Parametri
     A(k+1, 2:k+1) = A(k, 2:k+1) + Kk(k)*conj(A(k, k:-1:1));

     ef_new = ef(2:end) + Kk(k)* eb(2:end);
     eb = eb(1:end-1) + Kk(k)* ef(1:end-1);
     ef = ef_new;
    end

end
