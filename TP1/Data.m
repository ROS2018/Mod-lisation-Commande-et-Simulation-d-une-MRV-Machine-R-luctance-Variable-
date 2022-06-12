clc,clear;
M1 = 30; M2 = 250; k= 2*10^4; k1 = 2*10^5; f = 1000;
A = [ 0 1 0 0
      -(k1+k)/M1 -f/M1 k/M1 f/M1
      0            0    0     1
      k/M2        f/M2   -k/M2  -f/M2];
B = [0; -1/M1; 0; 1/M2];
Bp = [0; k1/M1; 0 ; 1/M1];
CX= eye(4);
C = [-1 0 1 0];
D=[0;0;0;0]; 
sys = ss (A,B,CX,D);

t=0:0.001:5;
u=0*t;
X0= [0.0022, -0.1146, 0.0421, -0.0138]';

%lsim(sys,u,t,X0)

% La commandabilité:
CC = [B A*B A*A*B A*A*A*B];
rank(CC);

%%
%1
R = 10^(-7) ; %Qx = 10*eye(4);
Qy = 10 ;
Qx = C'*Qy*C;
N=0;
K = lqr(A,B,Qx,R,N);

%

%% Eestimateur d'etat optimal, Filtre de Kalman:
L4=A(4,:);
sys = ss (A,B,L4,0);
%==========================
Qx4 = 10 ;
Q2 = C'*Qy*C +  L4'*Qx4*L4;
%===============================
W = 0.01; V = 0.01;  Nn=0;
[Sys_ob,Kf] = kalman(sys,W,V,Nn);Kf

