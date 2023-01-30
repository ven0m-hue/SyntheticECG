function fs2 = q_wave(tq, tr, Aq, Ar)
N = .01;
dx = N;
Tinq = tq;
x = 0:.01:Tinq;
y = 0*x;
L = Tinq;
mr = 2*Ar/tr;     % slope of R peak for correction
mx = (Ar+Aq)/mr;  % corrected slope by finding new value of x    
tdr = mx - tr/2;   % for mat point(x*100)
tqr = tq - tdr;% corrected x value for first equation

mq = Aq/tqr;    %corrected first slope

l1 = - mq * x ; %first slope

l2 = mr*x - mr*(tqr+tdr);     % second slope 

%Note!!
%(tqr+tdr) beacause, this equation is raising edge hence can't the shift is not sufficent according to the formula 
%there for geometrically shifted to the last point of the interval.

y(1:(tqr)*100) = l1(1:(tqr)*100);

y(tqr*100:(Tinq*100)+1) = l2(tqr*100:(Tinq*100)+1);

A1 = zeros;
x1 = zeros;
%Fourier Series
a1 = sum(y)*dx/L;
fs2 = 0;
for k=1:300
   %cc = jet(100);
   A1(k) = sum(y.*cos(pi*k*x/L))*dx/L ; 
   x1(k) = sum(y.*sin(pi*k*x/L))*dx/L ; 
   fs2   = fs2 + A1(k)*cos(k*x*pi/L) + x1(k)*sin(k*x*pi/L); 
  % plot(x,fs2/10,'-k','linewidth',1); 
   %pause(.1);
end

 fs2 = fs2 +a1/2;