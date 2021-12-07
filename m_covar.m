function [a , var] = m_covar(x,p)

  N = length(x);
  cxx = zeros(p+1,p+1);
 
  for j = 1:p+1
      for k = 1:p+1
         cxx(j,k) = sum(conj(x(p-j+2:N+1-j)).*x(p-k+2:N+1-k) + ...
                        conj(x(p-p+j:N-1-p+j).*x(p-p+k:N-1-p+k)))/(2*N-p);
      end
  end

  a = -(cxx(2:p+1,2:p+1)\cxx(2:p+1,1)).';
  var = cxx(1,1) + sum(a.*cxx(1,2:p+1));
  a = [1 a]; 
  
end