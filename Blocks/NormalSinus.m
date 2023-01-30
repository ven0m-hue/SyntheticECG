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

%% P wave
%x = 0:0.01:Tinp;
p_wave(tp, td, tpq, tq, Ap);
plot(x,fs1,'-k','linewidth',1); 
hold on;

for k=D:D:1000
    plot(x+k,fs1,'-k','linewidth',1); 
end
%hold on;

%%
%%Q
q_wave(tq, tr, Aq, Ar);
plot(x+Tinp, fs2, '-k','linewidth',1); 

for k=D:D:1000
plot(x+k+Tinp,fs2,'-k','linewidth',1);   
end

%%
%qRs
r_wave(tr, Ar)
plot( x + Tinp + Tinq,fs3,'-k','linewidth',1); 

for k=D:D:1000
    plot(x+k+Tinp+Tinq,fs3,'-k','linewidth',1);   
end

%%
%%s
s_wave(ts, tr, As, Ar);
plot(x + Tinp + Tinq + Tr ,fs4,'-k','linewidth',1);

for k=D:D:1000
    plot(x+k+Tinp+Tinq+Tr,fs4,'-k','linewidth',1);   
end


%%
%T
t_wave(tt, tst, ts, td, At);
plot(x + Tinp + Tinq + Tr +Tins, fs5,'-k', 'linewidth',1);   %hold off;

for k=D:D:1000
    plot(x+k+Tinp+Tinq+Tr+Tins,fs5,'-k','linewidth',1);   
end
hold off;

F = [fs1 fs2 fs3 fs4 fs5];


%% Plots
axis([0 850 -2 1])
title('1st Degree AV Block');
xlabel('interval of 10 seconds x10^-2');
ylabel('Apmlitude in mV');
%axis([0 1002 -50 200])
grid on 
%AX.XMinorGrid = 'on';
shg