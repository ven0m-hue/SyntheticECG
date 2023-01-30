MIT License

Copyright (c) 2023 ven0m-hue

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.


%% The follwoing is my algorithm 
% Given a heart beat, it can determine it's rr interval and hence its period
% Given amplitude dueration of the each interval, compute the equation of
% each wave.
% compute its fourier transforms of each 
% plot the waveforms 
% Also do the same for 16 different abnormalities, by varying the
% parameters as studied from the abnormalities.
%% Eqautions 
% P wave
 % f(x) = a(p)*(x-(tp+td))^2 + Ap 
 % where, a(p) = Ap/(tp/2)^2
 %         tp is duration of p wave, Ap is amplitude
% similarly for T wave 
% Q wave
 % f(x)q = a(q)*x + Bq
 % where, a(q) = 2*Aq/tq , Bq = a(q)*tq
 % where  Aq is the amplitude, tq is duration of q wave
% S wave 
 %similarly for s wave 
%qRs wave
 % f(x)qrs =  mqrs*x - mqrs*tqrs  for x<0
 % f(x)qrs = -mqrs*x + mqrs*tqrs  for x>0
 % mqrs = (2*Aqrs/tqrs)
 % Aqrs = amplitude, tqrs = durtation of the qrs interval
 
%% 
%% just a prototype (doesn't follow the algorithm)
%%initials
L = 10;
N = .1;
dx = N;
b = 0:N:L;
%%p
f1 = -.4*(b-5).^2 + 2.5;
f1(1:25) = 0*(1:25);
f1(3*25+1:101) = 0*(3*25+1:101);
f1_ = [f1 f1];
%plot(b-10,f1,'-k','linewidth',2.5); hold on;
axis([0 10 -4 4])
%%  fourier series!
%P wave 
% init 
A = zeros;
B = zeros;
a0 = sum(f1)*dx/L;
fs1 = 0;
for k=1:100
   %cc = jet(80);
   A(k) = sum(f1.*cos(pi*k*b/L))*dx/L ; 
   B(k) = sum(f1.*sin(pi*k*b/L))*dx/L ; 
   fs1   = fs1 + A(k)*cos(k*b*pi/L) + B(k)*sin(k*b*pi/L); 
   %plot(x,fs,'linewidth',2);
   %pause(.1);
end
 fs1 = fs1 +a0/2;
 plot(b-10,fs1,'-k','linewidth',2.5); hold on;
 plot(b-10+80,fs1,'-k','linewidth',2.5); 
%axis([0 10 -1 10]);
%shg
%%
%%Q wave
% init 
A1 = zeros;
B1 = zeros;

b = 0:.01:10;    
z = 0*b;
g2 = -.2*b+1 ;
z(1:500) = 0*(1:500); 
%g2 = -3*a + 60;
z(501:1000) = g2(501:1000);
%plot(b,z,'-m','linewidth',2.5); hold on;
%axis([-100 100 -2 18])
%%
a1 = sum(z)*dx/L;
fs2 = 0;
for k=1:1000
   %cc = jet(80);
   A1(k) = sum(z.*cos(pi*k*b/L))*dx/L ; 
   B1(k) = sum(z.*sin(pi*k*b/L))*dx/L ; 
   fs2   = fs2 + A1(k)*cos(k*b*pi/L) + B1(k)*sin(k*b*pi/L); 
   %plot(x,fs,'linewidth',2);
   %pause(.1);
end
 fs2 = fs2 +a1/2;
 plot(b,fs2/10,'-k','linewidth',2.5);
 plot(b+80,fs2/10,'-k','linewidth',2.5);

%%

%% qRs
% init 
A2 = zeros;
B2 = zeros;

b = 10:.01:20;    
z = 0*b;
g1 = 3*b - 30;
z(1:500) = g1(1:500);
g2 = -3*b + 60;
z(501:1000) = g2(501:1000);
%plot(b,z,'-k','linewidth',2.5);
%axis([-10 100 -2 18])
%%
a2 = sum(z)*dx/L;
fs3 = 0;
for k=1:200
   %cc = jet(80);
   A2(k) = sum(z.*cos(pi*k*b/L))*dx/L ; 
   B2(k) = sum(z.*sin(pi*k*b/L))*dx/L ; 
   fs3   = fs3 + A2(k)*cos(k*b*pi/L) + B2(k)*sin(k*b*pi/L); 
   %plot(x,fs,'linewidth',2);
   %pause(.1);
end
 fs3 = fs3 +a2/2;
 plot(b,fs3/10,'-k','linewidth',2.5);
plot(b+80,fs3/10,'-k','linewidth',2.5);
%%
%%S
% init 
A3 = zeros;
B3 = zeros;

dx = .1;
L = 10;
b = 0:.01:10;
z = 1*b;
g1 = .4*b -2;
%z(1) = ones(1);
z(2:500) = g1(2:500); 
%g2 = -3*a + 60;
z(501:1001) = 0*(501:1001);
plot(b+20,z,'-k','linewidth',2.5); %hold off;
%%
a3 = sum(z)*dx/L;
fs4 = 0;
for k=1:1000
   %cc = jet(80);
   A3(k) = sum(z.*cos(pi*k*b/L))*dx/L ; 
   B3(k) = sum(z.*sin(pi*k*b/L))*dx/L ; 
   fs4   = fs4 + A3(k)*cos(k*b*pi/L) + B3(k)*sin(k*b*pi/L); 
   %plot(x,fs,'linewidth',2);
   %pause(.1);
end
 fs4 = fs4 +a3/2;
 plot(b+20,fs4/10,'-k','linewidth',2.5);
 plot(b+80+20,fs4/10,'-k','linewidth',2.5);
%%
%%T
% init 
A4 = zeros;
B4 = zeros;

l = 40;
dx = .1;
b = 30:.1:70;
f2 = -.025*(b-50).^2 + 2.5;
f2(1:100) = 0*(1:100);
f2(300:401) = 0*(300:401);
%plot(b,f2,'-k','linewidth',2.5);hold off;shg
a4 = sum(f2)*dx/l;
fs5 = 0;
for k=1:200
   %cc = jet(80);
   A4(k) = sum(f2.*cos(pi*k*b/l))*dx/l ; 
   B4(k) = sum(f2.*sin(pi*k*b/l))*dx/l ; 
   fs5   = fs5 + A4(k)*cos(k*b*pi/l) + B4(k)*sin(k*b*pi/l); 
   %plot(x,fs,'linewidth',2);
   %pause(.1);
end
 fs5 = fs5 +a4/2;
 plot(b,fs5,'-k','linewidth',2.5); 
 plot(b+80,fs5,'-k','linewidth',2.5);hold off;
 
 %%
 %AXIS
    axis([-20,160 -5 20]);
    shg
    





