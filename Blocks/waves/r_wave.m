function fs3 = r_wave(tr, Ar)
N = .01;
dx = N;
Tr = tr;
x = 0:.01:Tr;
L = Tr;
y = 0*x;
g1 = 2*Ar/tr * x;

g2 = -2*Ar/tr *x + 2*Ar;

y(1:tr/2*100) = g1(1:tr/2*100);
y(tr/2*100+1:tr*100) = g2(tr/2*100+1:tr*100);

 A2 = zeros;
 x2 = zeros;
 
 a2 = sum(y)*dx/L;
 fs3 = 0;
 for k=1:200
    
    A2(k) = sum(y.*cos(pi*k*x/L))*dx/L ; 
    x2(k) = sum(y.*sin(pi*k*x/L))*dx/L ; 
    fs3   = fs3 + A2(k)*cos(k*x*pi/L) + x2(k)*sin(k*x*pi/L); 
    
 end
 
  fs3 = fs3 +a2/2;