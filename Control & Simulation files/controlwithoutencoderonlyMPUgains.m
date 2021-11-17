clear all

% Mb= 0.450 ;%kg Mass of robot
% Mw= 0.033;% kg Mass of the wheels,
% Jb= 0.0015;% kgm^2;% Moment of inertia about CoM
% r = 0.0335 ;%m radius of wheels
% Jw= (3.703*10^-3) ;  %add the motor J*n^2 %kgm^2 Moment of inertia for the wheels
% L = 0.068 ;%m Distance from wheel axle to CoM
% Ke= 0.477/48 ;%the real number is divided by 48        %0.477Vs/rad EMF constant
% Km= 0.313;
% 8/48 ;        %0.32Nm/A Motor constant
% R = 5.6/48;  %18 ;%Ohm Motor resistance
% b = 0.002 ;%Nms/rad Viscous friction constant
% g = 9.81 ;%m/s^2 Gravitational constant
Mb= 0.462-0.033*2 ;%kg Mass of robot
Mw= 0.033;% kg Mass of the wheels,
Jb= 0.0015;% kgm^2;% Moment of inertia about CoM
r = 0.0335 ;%m radius of wheels
Jw= 8.99*10^-4;%%(1.8517*10^-5) ;  %add the motor J*n^2 %kgm^2 Moment of inertia for the wheels
L = 0.068 ;%m Distance from wheel axle to CoM
Ke= 0.006 ;%the real number is divided by 48        %0.477Vs/rad EMF constant
Km= 0.0018 ;   %/48     %0.32Nm/A Motor constant
R = 5.6/48;  %18 ;%Ohm Motor resistance
b = -3.15*10^-4;%Nms/rad Viscous friction constant
g = 9.81 ;%m/s^2 Gravitational constant
z1=R*(2*(Jb*Jw+Jw*L^2*Mb+Jb*Mw*r^2+L^2*Mb*Mw*r^2)+Jb*Mb*r^2);

alp =(2*(R*b-Ke*Km)*(Mb*L^2+Mb*r*L+Jb))/z1;
beta= (-L^2*Mb^2*g*r^2)/(Jb*(2*Jw+Mb*r^2+2*Mw*r^2)+2*Jw*L^2*Mb+2*L^2*Mb*Mw*r^2);
gama=(-2*(R*b-Ke*Km)*(2*Jw+Mb*r^2+2*Mw*r^2+L*Mb*r))/(z1*r);
S=(L*Mb*g*(2*Jw+Mb*r^2+2*Mw*r^2))/(2*Jb*Jw+2*Jw*L^2*Mb+Jb*Mb*r^2+2*Jb*Mw*r^2+2*L^2*Mb*Mw*r^2);
eps=(Km*r)/(R*b-Ke*Km);

A=[0  1    0      0      0  ;
   0 alp beta   -r*alp   0  ;
   0  0    0      1      0  ;
   0 gama  S    -r*gama  0  ;
   0 0     1      0      0  ];

B=[0 alp*eps  0  gama*eps 0]';
C=[1 0 0 0 0 ;
   0 1 0 0 0 ;
   0 0 1 0 0 ;     
   0 0 0 1 0 ;
   0 0 0 0 1 ]; 
D=[0;0;0;0;0];
%increase the weigth of the corresponding state variable to better response
 Q = C'*C;
 Q(1,1) =0.1 ; %this is the weight of postion
 Q(2,2)=0.1;

 Q(3,3) =60000;   %weight of theta
 Q(4,4)=500;
 Q(5,5)=600;   %weight of position integrator
       %weight of theta integrator
 RR = 1;  %weight of the input gain
 K = lqr(A,B,Q,RR);
sys2=ss(A,B,C,D);
%%
%%K=zeros(1,5);


%K=[0 0 -200 -25 -0.0000];
Nbar=K(3); %constant gain called compensator makes response faster 
%K=zeros(1,6);
%K(6)=100;

Ac = (A-B*K);
Bc = Nbar*B+[0;0;0;0;-1];
Cc = C;
Dc = D;

states = {'x' 'x_dot' 'phi' 'phi_dot' 'theta integrator'};
inputs = {'r'};
outputs = {'x';'xdot'; 'phi';'phidot';'theta integral'};

sys_cl = ss(Ac,Bc,Cc,Dc,'statename',states,'inputname',inputs,'outputname',outputs);

t = 0:0.001:50;
set_x=0;%in meters % the set point on the position that we need the robot to go to it
x_0=[0 0 0.01 0 0] ;%intial state  %note theta in radian
rr =set_x*ones(size(t));
[y,t,x]=lsim(sys_cl,rr,t,x_0);
gain=Nbar*set_x-K*y';
gain=gain';
torque=(Km/R)*gain-((Ke*Km)/R)*((y(:,2)./r)-y(:,4));

%m=(y(:,2)./r)-y(:,4); %relative speed RBM
figure(14)
plot(t,y(:,1),'r',t,y(:,3),'g');
%axis([0, 5, -5, 5]);

title('Graph of position and  theta')
legend('y = postion in meter','y = theta in rad')
xlabel('time in seconds')
figure(2);
plot(t,gain,'b',t,torque*10.197,'y');
title('Graph of voltage and torque')
legend('y = voltage','y = torque in kg.cm')
xlabel('time in seconds')
fprintf('max gain=%f volt - min gain=%f volt\nmax torque=%f kg.cm -min torque= %fkg.cm\n',max(gain),min(gain),max(torque)*10,min(torque)*10);
