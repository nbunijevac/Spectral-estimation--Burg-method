function [rxx] = akf(x) 

    N = length(x);
    rxx = zeros(2*N-1,1);

    for k = N:2*N-1
      rxx(k) = sum(conj(x(1:2*N-k)).*x(1+k-N:2*N-k+k-N));
      rxx(k) = rxx(k)/(2*N-k); % delimo sa N - pomerena, 2*N-k - nepomerena
      rxx(2*N-k) = conj(rxx(k));
    end
    
end