function fs4 = s_wave(ts, tr, As, Ar)
N = .01;
dx = N;
Tins = ts;
x = 0:.01:Tins;
L = Tins;
mr = 2*Ar/tr;      % slope of R peak for correction
mx = (Ar+As)/mr;  % corrected slope by finding new value of x    

t1r = mx - tr/2;   % for mat point(x*100)
tsr = ts - t1r;% corrected x value for first equation

ms = As/tsr;   % corrected slope for second equation
%ms = (As)/(ts);   % wrong slope
y = 0*x;

a1 =  - mr*x  ;
a2 = ms * x - ms*(ts);

y(1:t1r*100) = a1(1:t1r*100);
y(t1r*100:Tins*100) = a2(t1r*100:Tins*100);

A3 = zeros;
x3 = zeros; 
a3 = sum(y)*dx/L;
fs4 = 0;
for k=1:100
   %cc = jet(100);
   A3(k) = sum(y.*cos(pi*k*x/L))*dx/L ; 
   x3(k) = sum(y.*sin(pi*k*x/L))*dx/L ; 
   fs4   = fs4 + A3(k)*cos(k*x*pi/L) + x3(k)*sin(k*x*pi/L); 
 %plot(x+80+20,fs4/10,'-k','linewidth',2.5);
 %pause(.1);  
end

 fs4 = fs4 +a3/2;