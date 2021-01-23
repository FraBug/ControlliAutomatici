%% Project of Automatic Control 
%%%%%%%%
%
%   Project:
%       Type 1.a
%   Team members:
%       Francesco Ciampana, Alessio Troffei, Vladyslav Tymofieiev,
%       Francesco Scavello
%
%%%%%%%%

% Project specs:
% - zero steady state error with a step reference signal of w(t)=W1(t)
% - Mf>45°
% - S_%<1%
% - T_a1< 0,3[s]

%System parameters definition
b = 0.3; a = 0.3; F_v = -8.7; 
T_a1 = 0.3;
S_p = 0.01; % Sovraelongazione percentuale
w_n_min=12 * 10^4;
A_n = 30; % abbattimento
%Matrixs
A = [0 , 1;
    a/4*sqrt(2)*F_v , -b];

B = [ 0 ;
      a ];
  
C = [1 , 0];
D = 0;

%Transfer function definition
s = tf("s");
[Num,Den] = ss2tf(A,B,C,D); %converts a state-space representation of a system into an equivalent transfer function.
G=tf(Num,Den);
zpk(G);

%Definition of bode frequency range
w_plot_min=10^(-2);%vedere e nel caso modificare
w_plot_max=10^6;%vedere e nel caso modificare

% Constraints variables
%xi=sqrt((log(S_p)^2)/(pi^2+log(S_p)^2)); %= circa 0.83;
xi=0.85;
Mf=xi*100;

w_c_min=460/(Mf*T_a1);
w_c_max=w_n_min; % [rad/s]
Xi_dB = 20*log10(1/A_n); %abbattimento di 30 volte

%Regolatore statico
Rs = 1/s;

%Ge = G* Rs

Ge = G * Rs;
zpk(Ge);

%REGOLATORE DINAMICO


%omega_cStar=(w_c_min+50);
%phi_star= (Mf+10) - 180 -(-270);display(phi_star);%angolo desiderato
%M_star = 90;
%tau= (M_star - cosd(phi_star))/(omega_cStar * sind(phi_star));
tau1 = 1;
tau2 = 100;
a1 = 10^(-4);
a2 = 10^(-5);
Rd = (1+ tau1*s) * (1+tau2*s)/ ((1+a2*tau2*s)*(1 + a1*tau1*s));

L = Rd*Ge;


%------------------------------------------
[MagG,phaseG,wG]=bode(Ge,{ w_plot_min, w_plot_max });
[MagL,phaseL,wL]=bode(L,{ w_plot_min, w_plot_max });

patch([w_n_min,w_plot_max,w_plot_max,w_n_min],[Xi_dB,Xi_dB,200,200],[1 0.5 0],'FaceAlpha',0.4,'EdgeAlpha',0);
patch([w_c_min,w_plot_min,w_plot_min,w_c_min],[-300,-300,0,0],[1 1 0],'FaceAlpha',0.4,'EdgeAlpha',0);
leg1 = legend('1) disturbo misura', '2) tempo di assest.', '3)Ge');
hold on,  zoom on;

margin(MagG, phaseG, wG);
margin(MagL, phaseL, wL);

leg2 = legend('G', 'L');

patch([w_c_min,w_c_max,w_c_max,w_c_min],[-180+Mf,-180+Mf,-280,-280],'r','FaceAlpha',0.35,'EdgeAlpha',0);


return 
