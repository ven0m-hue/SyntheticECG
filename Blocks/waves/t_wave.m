function fs5 = t_wave(tt, tst, ts, td, At)
N = .01;
dx = N;
Tint = (tst + td) - ts;
x = 0:0.01:Tint;
L = Tint;

y = - At/(tt/2)^2 * (x - (tt/2 +ts)).^2 + At;
y(1:ts*100) = 0*(1:ts*100);
y((ts+tt)*100+1:Tint*100 +1) = 0*((ts+tt)*100+1:Tint*100 +1);
%Fourier Series 
a4 = sum(y)*dx/L;
fs5 = 0;

A4 = zeros;
x4 = zeros;

for k=1:500
 
   A4(k) = sum(y.*cos(pi*k*x/L))*dx/L ; 
   x4(k) = sum(y.*sin(pi*k*x/L))*dx/L ; 
   fs5   = fs5 + A4(k)*cos(k*x*pi/L) + x4(k)*sin(k*x*pi/L); 
 
end
 fs5 = fs5 +a4/2;

