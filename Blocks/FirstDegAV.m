%Input patametrs of the ECG

td = 30;%10;%10;
tp =  6;%14;%15;
tq = 2;
tr = 2;
ts = 6;
tt = 16;

tst = 43;
tpq = 18;


Ap = 0.1;
Aq = 0.6;
Ar = .6;
As = 1;
At = .18 ;

%input to enter the BPM 
D = 93;


N = .01;
dx = N;
%b = 0:N:L;
%% P wave
%if tpq < tp+tq
  %  print("Error");
%else
Tinp = (tpq + td) - tq;
x = 0:0.01:Tinp;
L = Tinp;

y = - Ap/(tp/2)^2 * (x - (tp/2 +td)).^2 + Ap;
y(1:td*100) = 0*(1:(td*100));
y((td+tp)*100+1:Tinp*100+1) = 0*((td+tp)*100+1:Tinp*100+1);
plot(x,y,'-k','linewidth',1); %hold on ;
%Fourier Series 
a0 = sum(y)*dx/L;
fs1 = 0;

A = zeros;
B = zeros;
for k=1:500
   %cc = jet(20);
   A(k) = sum(y.*cos(pi*k*x/L))*dx/L ; 
   B(k) = sum(y.*sin(pi*k*x/L))*dx/L ; 
   fs1   = fs1 + A(k)*cos(k*x*pi/L) + B(k)*sin(k*x*pi/L); 
   %plot(b-10,fs1,'-k','color',cc(k,:),'linewidth',1.5); hold on;
   %pause(.1);
end
 fs1 = fs1 +a0/2;
plot(x,fs1,'-k','linewidth',1); 
hold on;

for k=D:D:1000
    plot(x+k,fs1,'-k','linewidth',1); 
end
%hold on;

%%
%%Q
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

y(1:(tqr)*100) = 0*(1:(tqr)*100);

y(tqr*100:(Tinq*100)+1) = 0*(tqr*100:(Tinq*100)+1);

%Fourier Series
a1 = sum(y)*dx/L;
fs2 = 0;

A1 = zeros; % Some pre-allocatation 
x1 = zeros;

for k=1:300
   %cc = jet(100);
   A1(k) = sum(y.*cos(pi*k*x/L))*dx/L ; 
   x1(k) = sum(y.*sin(pi*k*x/L))*dx/L ; 
   fs2   = fs2 + A1(k)*cos(k*x*pi/L) + x1(k)*sin(k*x*pi/L); 
  % plot(x,fs2/10,'-k','linewidth',1); 
   %pause(.1);
end

 fs2 = fs2 +a1/2;
 %plot(x,fs2/10,'-k','linewidth',2.5);

plot(x+Tinp, fs2, '-k','linewidth',1); 

%hold on
for k=D:D:1000
plot(x+k+Tinp,fs2,'-k','linewidth',1);   
end

%%
%qRs

Tr = tr;
x = 0:.01:Tr;
L = Tr;
y = 0*x;
g1 = 2*Ar/tr * x;

g2 = -2*Ar/tr *x + 2*Ar;
%g3 = 2*Ar/tr*x;
y(1:tr/2*100) = g1(1:tr/2*100);
y(tr/2*100+1:tr*100) = g2(tr/2*100+1:tr*100);

 a2 = sum(y)*dx/L;
 fs3 = 0;
 
 A2 = zeros; % Some pre-allocatation 
 x2 = zeros;
 for k=1:50
    cc = jet(200);
    A2(k) = sum(y.*cos(pi*k*x/L))*dx/L ; 
    x2(k) = sum(y.*sin(pi*k*x/L))*dx/L ; 
    fs3   = fs3 + A2(k)*cos(k*x*pi/L) + x2(k)*sin(k*x*pi/L); 
    %plot(x,fs3/10,'-k','linewidth',2); 
    %pause(.1);
 end
 
  fs3 = fs3 +a2/2;
  %plot(x,fs3/10,'-k','linewidth',2.5);

plot( x + Tinp + Tinq,fs3,'-k','linewidth',1); 
%axis([25 45 -10 135])
for k=D:D:1000
    plot(x+k+Tinp+Tinq,fs3,'-k','linewidth',1);   
end

%%
%%s
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

y(1:t1r*100) = a1(1:(t1r*100));
y(t1r*100:Tins*100) = a2(t1r*100:Tins*100);
plot(x + Tinp + Tinq + Tr ,y,'-k','linewidth',1);

a3 = sum(y)*dx/L;
fs4 = 0;

A3 = zeros; % Some pre-allocatation 
x3 = zeros;
for k=1:100
   cc = jet(100);
   A3(k) = sum(y.*cos(pi*k*x/L))*dx/L ; 
   x3(k) = sum(y.*sin(pi*k*x/L))*dx/L ; 
   fs4   = fs4 + A3(k)*cos(k*x*pi/L) + x3(k)*sin(k*x*pi/L); 
 %plot(x+80+20,fs4/10,'-k','linewidth',2.5);
 %pause(.1);  
end

 fs4 = fs4 +a3/2;
%plot(x+20,fs4/10,'-k','linewidth',2.5);
plot(x + Tinp + Tinq + Tr ,fs4,'-k','linewidth',1);

for k=D:D:1000
    plot(x+k+Tinp+Tinq+Tr,fs4,'-k','linewidth',1);   
end


%%
%T
%if tst < tt+ts
 %   print("Error");
%else
Tint = (tst + td) - ts;
x = 0:0.01:Tint;
L = Tint;
y = - At/(tt/2)^2 * (x - (tt/2 +ts)).^2 + At;
y(1:ts*100) = 0*(1:ts*100);
y((ts+tt)*100+1:Tint*100 +1) = 0*((ts+tt)*100+1:Tint*100 +1);

a4 = sum(y)*dx/L;
fs5 = 0;

A4 = zeros; % Some pre-allocatation 
x4 = zeros;

for k=1:500
   %cc = jet(80);
   A4(k) = sum(y.*cos(pi*k*x/L))*dx/L ; 
   x4(k) = sum(y.*sin(pi*k*x/L))*dx/L ; 
   fs5   = fs5 + A4(k)*cos(k*x*pi/L) + x4(k)*sin(k*x*pi/L); 
   %plot(x,fs,'linewidth',2);
   %pause(.1);
end
 fs5 = fs5 +a4/2;

plot(x + Tinp + Tinq + Tr +Tins, fs5,'-k', 'linewidth',1);   %hold off;

for k=D:D:1000
    plot(x+k+Tinp+Tinq+Tr+Tins,fs5,'-k','linewidth',1);   
end
hold off;

F = [fs1 fs2 fs3 fs4 fs5];


%%
axis([0 850 -2 1])
%shg
title('1st Degree AV Block');
xlabel('interval of 10 seconds x10^-2');
%xticks([1 2 3 4 5 6 7 8 9 10]);
ylabel('Apmlitude in mV');
%yticks([-5 0 5 10 15 20])
%axis([0 1002 -50 200])
grid on 
%AX.XMinorGrid = 'on';
shg