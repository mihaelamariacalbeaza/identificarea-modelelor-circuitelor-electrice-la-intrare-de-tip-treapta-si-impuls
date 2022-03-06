%% dreapta

clc; close all; clear all;
data_project_dreapta = readmatrix('C:\Users\myche\OneDrive\Desktop\scope_227.csv','Range','A3:C1002');
t_drepta=data_project_dreapta(:,1);
u_dreapta=data_project_dreapta(:,2);
y_treapta=data_project_dreapta(:,3);

figure,
plot(t_drepta,[u_dreapta,y_treapta])
%title('Suprapunerea dintre semnalul masurat si semnalul simulat cu modelul dedus.')
title('Reprezentare y stationar si y0')
grid


xlabel('T(sec.)');
ylabel('u(V)\y(V)');
hold on


i1=453;
i2=488;
i3=924; 
i4=944;


yst=mean(y_treapta(i3:i4));
y0=mean(y_treapta(i1:i2));
ust=mean(u_dreapta(i3:i4));
u0=mean(u_dreapta(i1:i2));
plot(t_drepta,yst.*ones(1,length(t_drepta)),'y')
plot(t_drepta,u0.*ones(1,length(t_drepta)),'g')
legend('u','y','y_s_t_a_t_i_o_n_a_r','y_0')


k=(yst-y0)/(ust-u0)
sigma=(max(y_treapta)-yst)/(yst-y0) %suprareglajul
zita=-log(sigma)/sqrt(pi^2+log(sigma).^2)

i5=497;
i6=666;
Tosc=2*(t_drepta(i6)-t_drepta(i5))
wosc=2*pi/Tosc
wn=wosc/sqrt(1-zita^2)
A=[0 1; -wn^2 -2*zita*wn]
B=[0;k*wn^2]
C=[1 0]
D=[0]
ysim=lsim(A,B,C,D,u_dreapta,t_drepta,[y_treapta(1),0]);
hold on
title('Suprapunerea dintre semnalul masurat si semnalul simulat cu modelul dedus.')
plot(t_drepta,ysim)
legend('u','y','y_s_t_a_t_i_o_n_a_r','y_0','semnalul simulat cu modelul ales')

J=sqrt((1/length(y_treapta))*sum((y_treapta-ysim).^2))%eroarea medie patratice
eMIN=norm(y_treapta-ysim)/norm(y_treapta-min(y_treapta)) %eroarea medie patratica normalizata(%)

data_project_impuls = readmatrix('C:\Users\myche\OneDrive\Desktop\scope_228.csv','Range','A3:C1002');
t_nou_impuls=data_project_impuls(:,1);
u_nou_impuls=data_project_impuls(:,2);
y_nou_impuls=data_project_impuls(:,3); 
%simularea modelului de la treapta cu date de la impuls
ysim_impuls_nou=lsim(A,B,C,D,u_nou_impuls,t_nou_impuls,[y_nou_impuls(1) 0]);
figure;
plot(t_nou_impuls,[u_nou_impuls, y_nou_impuls, ysim_impuls_nou])
title('simularea modelului de la treapta cu date de la impuls');
grid

xlabel('T(sec.)');
ylabel('u(v)\y(v)');
legend('u','y', 'treapta cu date de la impuls')

J_nou=sqrt((1/length(y_nou_impuls))*sum((y_nou_impuls-ysim_impuls_nou).^2))%eroarea medie patratice
eMIN_nou=norm(y_nou_impuls-ysim_impuls_nou)/norm(y_nou_impuls-min(y_nou_impuls)) %eroarea medie patratica normalizata(%)
