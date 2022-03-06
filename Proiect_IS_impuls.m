%% impuls

% H=K*wn^2/(s^2+2*tita*wn*s+wn^2);
clc; close all; clear all;
data_project_impuls = readmatrix('C:\Users\myche\OneDrive\Desktop\scope_228.csv','Range','A3:C1002');
t_impuls=data_project_impuls(:,1);
u_impuls=data_project_impuls(:,2);
y_impuls=data_project_impuls(:,3);


figure,
plot(t_impuls,[u_impuls,y_impuls])
title('Suprapunerea dintre semnalul masurat si semnalul simulat cu modelul dedus.')


xlabel('T(sec.)');
ylabel('u(V)\y(V)');
grid;
hold on

% i1=465;
% i2=475;
% i3=505; 
% i4=668;
% i5=835;

i1=937;
i2=950;
i3=505; 
i4=668;
i5=835;


ust=mean(u_impuls(i1:i2));
yst=mean(y_impuls(i1:i2));

plot(t_impuls,yst.*ones(1,length(t_impuls)),'g')
k=yst/ust

% determinarea suprareglajului
A_plus = (sum(abs(y_impuls(i3:i4)-yst)*(t_impuls(2)-t_impuls(1))))
A_minus = (sum(abs(y_impuls(i4:i5)-yst)*(t_impuls(2)-t_impuls(1))))
suprareglaj=A_minus/A_plus
factor_amortizare=-log(suprareglaj)/sqrt(pi^2+(log(suprareglaj))^2) %tita
tosc=2*(t_impuls(i5)-t_impuls(i4))


wosc=2*pi/tosc
wn=wosc/sqrt(1-factor_amortizare^2)
A=[0 1; -wn^2 -2*factor_amortizare*wn]
B=[0;k*wn^2]
C=[1 0]
D=[0]
ysim=lsim(A,B,C,D,u_impuls,t_impuls,[y_impuls(1),0]);
hold on

plot(t_impuls,ysim)
legend('u','y','y_s_t_a_t_i_o_n_a_r','semnalul simulat cu modelul ales')
grid
J=sqrt((1/length(y_impuls))*sum((y_impuls-ysim).^2))%eroarea medie patratice
eMIN=norm(y_impuls-ysim)/norm(y_impuls-min(y_impuls)) %eroarea medie patratica normalizata(%)


%simulate model impuls cu date de la treapra
data_project_dreapta = readmatrix('C:\Users\myche\OneDrive\Desktop\scope_227.csv','Range','A3:C1002');
t_nou_treapta=data_project_dreapta(:,1);
u_nou_treapta=data_project_dreapta(:,2);
y_nou_treapta=data_project_dreapta(:,3); 
ysim_treapta_nou=lsim(A,B,C,D,u_nou_treapta,t_nou_treapta, [y_nou_treapta(1) 0]);
figure,
plot(t_nou_treapta,[u_nou_treapta, y_nou_treapta, ysim_treapta_nou]);
title('simularea modelului de la impuls cu date de la treapta');
grid

xlabel('T(sec.)');
ylabel('u(v)\y(v)');
legend('u','y', 'impuls cu date de la treapta')


J_nou=sqrt((1/length(y_nou_treapta))*sum((y_nou_treapta-ysim_treapta_nou).^2))%eroarea medie patratice
eMIN_nou=norm(y_nou_treapta-ysim_treapta_nou)/norm(y_nou_treapta-min(y_nou_treapta)) %eroarea medie patratica normalizata(%)

