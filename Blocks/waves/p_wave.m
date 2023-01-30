function fs1 = p_wave(tp, td, tpq, tq, Ap)
N = .01;
dx = N;
Tinp = (tpq + td) - tq;
x = 0:0.01:Tinp;
L = Tinp;

y = - Ap/(tp/2)^2 * (x - (tp/2 +td)).^2 + Ap;
y(1:td*100) = 0*(1:(td*100));
y((td+tp)*100+1:Tinp*100+1) = 0*((td+tp)*100+1:Tinp*100+1);
%Fourier Series 
a0 = sum(y)*dx/L;
fs1 = 0;
for k=1:500
   
   A(k) = sum(y.*cos(pi*k*x/L))*dx/L ; 
   B(k) = sum(y.*sin(pi*k*x/L))*dx/L ; 
   fs1   = fs1 + A(k)*cos(k*x*pi/L) + B(k)*sin(k*x*pi/L); 

end
 fs1 = fs1 +a0/2;