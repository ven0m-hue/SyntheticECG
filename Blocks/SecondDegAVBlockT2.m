%Input patametrs of the ECG
%block_type = input("Select the block type, if av typ1 :1 else if av type2: 2");
% for bunch of if statements 
flag_sb = 0;
flag_sa = 0;
flag_av1 = 0;
flag_av2 = 0;
flag_av = 0;
flag_RBB = 0;
flag_LBB = 0;

% inputs 
td = 8;%22;%39.5;
tp = 4;%10;%15;9
tq = 1;
tr = 3;%6;
ts = 2;
tt = 6;

tst = 12;%40;
tpq = 9;%18;


Ap = 0.3;
Aq = .08;
Ar = 1.65;
As = 0.26;
At = 0.31;

%input to enter the BPM 
D  = 34;%99.5;   %always add the td part again 
D3 = 34;%99.5;
N = .01;  %sampling frequency 100hz
dx = N;
count = 0;


%% P wave
Tinp = (tpq + td) - tq;
x = 0:0.01:Tinp;
fs1 = p_wave(tp, td, tpq, tq, Ap);  %function call for p wave

plot(x,fs1,'-k','linewidth',1); % first beat 
hold on;

temp1 = tpq;
temp2 = td;
temp3 = td/2;
D7 = D3*3;
% second eat 

while(D < 1000)
    while(td > 0)

        tpq = tpq + temp3;     %for now hardcoded  
        td =  td  - temp3;     %later decide the rate of change in the pq intervel      
        Tinp = (tpq + td) - tq; 
        x = 0:0.01:Tinp;
        fs1 = p_wave(tp, td, tpq, tq, Ap);
        plot(x+D,fs1,'-k','linewidth',1); % second beat onwards, till it hits the drop beat
        D = D + D3;
        count = count +1;
        D1 = D;
        
    end
        D5 = D1+D3;
      for I = D1:D3:D5  
        tpq = temp1;
        td  = temp2;
        if 
        fs1 = p_wave(tp, td, tpq, tq, 0);    %for sinus block Ap != 0 for AV type 2,1--> Ap = 0
        plot(x+D1,fs1,'-k','linewidth',1);
        fs1 = p_wave(tp, td, tpq, tq, Ap);   % for sinus arrest Ap = 0
        plot(x+D5,fs1,'-k','linewidth',1);   
      end
      D = D3+D5;
end
 


%     D6 = D5 + D3;
%     while(td > 0)
%     
%     tpq = tpq + temp3;     %for now hardcoded  
%     td =  td  - temp3;     %later decide the rate of change in the pq intervel      
%     Tinp = (tpq + td) - tq; 
%     x = 0:0.01:Tinp;
%     
%     fs1 = p_wave(tp, td, tpq, tq, Ap);
%     plot(x+D6,fs1,'-k','linewidth',1); % second beat onwards, till it hits the drop beat
%     D6 = D6 + D3;
%     count = count +1;
%     D2 = D6;
%     end
%    
%      tpq = temp1;
%      td  = temp2;
%      fs1 = p_wave(tp, td, tpq, tq, Ap);
%      plot(x+D2,fs1,'-k','linewidth',1);
%  
%%
%%Q
Tinq = tq;
x = 0:.01:tq;
fs2 = q_wave(tq, tr, Aq, Ar); %function call for q wave
plot(x+Tinp, fs2, '-k','linewidth',1); 

%hold on
Dr = 3;
for I=D3:D3:1000
     if I == D3 * Dr
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
        y(1:(tqr)*100) = 0*(1:(tqr)*100);
        y(tqr*100:(Tinq*100)+1) = 0*(tqr*100:(Tinq*100)+1);
     
        %Fourier Series
        a1 = sum(y)*dx/L;
        fs2 = 0;
        for k=1:300
           %cc = jet(100);
           A1(k) = sum(y.*cos(pi*k*x/L))*dx/L ; 
           x1(k) = sum(y.*sin(pi*k*x/L))*dx/L ; 
           fs2   = fs2 + A1(k)*cos(k*x*pi/L) + x1(k)*sin(k*x*pi/L); 
        end
     
         fs2 = fs2 +a1/2;
         plot(x+I+Tinp,fs2,'-k','linewidth',1);
         Dr = Dr + 4;
     else 
        fs2 = q_wave(tq, tr, Aq, Ar); %function call for q wave 
        plot(x+I+Tinp,fs2,'-k','linewidth',1);
     end
end

%%
%qRs
Tr = tr;
x = 0:.01:Tr;
fs3 = r_wave(tr, Ar);                                                       %function call for r wave
plot( x + Tinp + Tinq,fs3,'-k','linewidth',1); 

Dr = 3;  % reset the buffer
for I=D3:D3:1000
     if I == D3 * Dr 
        N = .01;
        dx = N;
        Tr = tr;
        x = 0:.01:Tr;
        L = Tr;
        y = 0*x;
        y(1:tr/2*100) = 0*(1:tr/2*100);
        y(tr/2*100+1:tr*100) = 0*(tr/2*100+1:tr*100);

         a2 = sum(y)*dx/L;
         fs3 = 0;
         for k=1:200
            A2(k) = sum(y.*cos(pi*k*x/L))*dx/L ; 
            x2(k) = sum(y.*sin(pi*k*x/L))*dx/L ; 
            fs3   = fs3 + A2(k)*cos(k*x*pi/L) + x2(k)*sin(k*x*pi/L);
         end

          fs3 = fs3 +a2/2;
         plot(x+I+Tinp+Tinq,fs3,'-k','linewidth',1); hold on;
         Dr = Dr + 4;
     else
        fs3 = r_wave(tr, Ar);                                               %function call for r wave 
        plot(x+I+Tinp+Tinq,fs3,'-k','linewidth',1);
     end
end

%%
%%s
Tins = ts;
x = 0:.01:Tins;
L = Tins;                                                                   %Period of S wave
fs4 = s_wave(ts, tr, As, Ar);                                               %function call for s wave
plot(x + Tinp + Tinq + Tr ,fs4,'-k','linewidth',1);
Dr = 3;
for I=D3:D3:1000
    if I == D3 * Dr
        Tins = ts;
        x = 0:.01:Tins;
        L = Tins;
        mr = 2*Ar/tr;                                                       % slope of R peak for correction
        mx = (Ar+As)/mr;                                                    % corrected slope by finding new value of x    
        t1r = mx - tr/2;                                                    % for mat point(x*100)
        y = 0*x;
        y(1:t1r*100) = 0*(1:t1r*100);
        y(t1r*100:Tins*100) = 0*(t1r*100:Tins*100);
        
        %Fourier Series
        
        a3 = sum(y)*dx/L;
        fs4 = 0;
        for k=1:100
           A3(k) = sum(y.*cos(pi*k*x/L))*dx/L ; 
           x3(k) = sum(y.*sin(pi*k*x/L))*dx/L ; 
           fs4   = fs4 + A3(k)*cos(k*x*pi/L) + x3(k)*sin(k*x*pi/L); 
        end
        fs4 = fs4 +a3/2;
        plot(x +I+ Tinp + Tinq + Tr ,fs4,'-k','linewidth',1); hold on;
        Dr = Dr + 4;
    else
        fs4 = s_wave(ts, tr, As, Ar);
        plot(x+I+Tinp+Tinq+Tr,fs4,'-k','linewidth',1);   
    end
end


%%
%T
td =  td  - temp3;                                                          % to remove the overlaping
Tint = (tst + td) - ts;
x = 0:0.01:Tint;

fs5 = t_wave(tt, tst, ts, td, At);                                          %function call for t wave
plot(x + Tinp + Tinq + Tr +Tins, fs5,'-k', 'linewidth',1);   

td =  td  - temp3/2;                                                        %invoked only for changing PR interval                   %To remove the overlaping.
Tint = (tst + td) - ts;                                                     %room for improvement, had no choice than seperately calling,
x = 0:0.01:Tint;                                                            %the function again 
fs5 = t_wave(tt, tst, ts, td, 0);                                           %function call for t wave
plot(x+D3 + Tinp + Tinq + Tr +Tins, fs5,'-k', 'linewidth',1);   
Dr = 3;
for I=D3+D3:D3:1000
    if I == D3 * Dr
        fs5 = t_wave(tt, tst, ts, td, 0);                                           % dropped beat
        plot(x+I+Tinp+Tinq+Tr+Tins,fs5,'-k','linewidth',1);
        Dr = Dr + 4;
    else
        fs5 = t_wave(tt, tst, ts, td, At);
        plot(x+I+Tinp+Tinq+Tr+Tins,fs5,'-k','linewidth',1);
    end
    
end
hold off;

F = [fs1 fs2 fs3 fs4 fs5];                                                  % entire signal in a single matrix


%%
axis([0 320 -1 2])

title('Second Degree AV Block, Type II');
xlabel('interval of 10 seconds x10^-2');
%xticks([1 2 3 4 5 6 7 8 9 10]);
ylabel('Apmlitude in mV');
%yticks([-5 0 5 10 15 20])
%axis([0 1002 -50 200])
grid on 
%AX.XMinorGrid = 'on';
shg