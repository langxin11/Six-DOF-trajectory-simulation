clc;clear;
%1:速度-v 2:弹道倾角-theta 3:弹道偏角-fi_v 
%4:w_x-omega_x1 5:w_y-omega_y1 6:w_x-omega_z1
%7:x 8 :y 9:z
%10: 俯仰角-Theta 11:偏航角-fi 12: 倾斜角-gama
%13：质量 m 速度倾斜角-gama_v
global Beta Alpha Gama_v global_delta global_dydt start_y
Beta=[0,0];Alpha=[0,0];Gama_v=[0,0];
global_delta=[0,0,0];global_dydt=zeros(1,16);
start_y=false;
in_0=[20;0.1*pi;0;0;0;0;0;20;0;0.1*pi;0;0;53.38;6000;20;0];
[t,out]=ode45(@missile,[0 15],in_0);
figure;
plot3(out(:,7),out(:,9),out(:,8),'k',LineWidth=2);
xlabel('x');ylabel("y");
grid on;
hold on;
plot3(out(:,14),out(:,16),out(:,15),'b',LineWidth=2);
% plot(out(:,7),out(:,8),'k',LineWidth=2);
% xlabel('x');ylabel("y");
% grid on;
% hold on;
% plot(out(:,14),out(:,15),'b',LineWidth=2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function  dydt=missile(t,in)
v=in(1);
theta=in(2);
fi_v=in(3);
omega_x1=in(4);
omega_y1=in(5);
omega_z1=in(6);
x=in(7);
y=in(8);
z=in(9);
Theta=in(10);
fi=in(11);
gama=in(12);
m=in(13);
%目标
xt=in(14);
yt=in(15);
zt=in(16);

global Beta Alpha Gama_v global_delta global_dydt start_y
beta=asin(cos(theta)*(cos(gama)*(sin(fi-fi_v)+sin(Theta)*sin(gama)*cos(fi-fi_v)))...
    -sin(theta)*cos(Theta)*sin(gama));
beta=real(beta);
Beta=[Beta;t,beta*180/pi];

alpha=asin((cos(theta)*(sin(Theta)*cos(gama)*cos(fi-fi_v)-...
    sin(gama)*sin(fi-fi_v))-sin(theta)*cos(Theta)*cos(gama))/cos(beta));
alpha=real(alpha);
Alpha=[Alpha;t,alpha*180/pi];
gama_v=asin((cos(alpha)*sin(beta)*sin(Theta)-...
    sin(alpha)*sin(beta)*cos(gama)*cos(Theta)+...
    cos(beta)*sin(gama)*cos(Theta))/cos(theta));
gama_v=real(gama_v);
Gama_v=[Gama_v;t,gama_v*180/pi];
g=9.81;
density=1.225;S=0.0227;L=1.8;Ma=v/343.13;l=0.5;

%弹目相对运动
Vt=1; theta_t=135*pi/180;  fi_vt=0;
%目标的质点弹道
dxt=Vt*cos(theta_t)*cos(fi_vt);
dyt=Vt*sin(theta_t);
dzt=-Vt*cos(theta_t)*sin(fi_vt);
Xs=[xt-x,yt-y,zt-z];Xs1=[xt-x,0,zt-z];
r=norm(Xs);
if r<2000||start_y
    start_y=true;
end
if xt>x
    if yt>y
        qy=acos((Xs*Xs1')/norm(Xs)/norm(Xs1));
    else
        qy=-acos((Xs*Xs1')/norm(Xs)/norm(Xs1));
    end
else
    if yt>y
        qy=pi-acos((Xs*Xs1')/norm(Xs)/norm(Xs1));
    else
        qy=acos((Xs*Xs1')/norm(Xs)/norm(Xs1))-pi;
    end
end
if xt>x
    if zt>z
        qz=-acos((Xs*[1,0,0]')/norm(Xs));
    else
        qz=acos((Xs*[1,0,0]')/norm(Xs));
    end
else
    if zt>z
        qz=acos((Xs*[1,0,0]')/norm(Xs))-pi;
    else
        qz=pi-acos((Xs*[1,0,0]')/norm(Xs));
    end
end

dr=Vt*(cos(theta_t)*cos(qy)*cos(fi_vt-qz)+sin(theta_t)*sin(qy))-...
    v*(cos(theta)*cos(qy)*cos(fi_v-qz)+sin(theta)*sin(qy));
dqy=(Vt*(sin(theta_t)*cos(qy)-cos(theta_t)*sin(qy)*cos(fi_vt-qz))-...
    v*(sin(theta)*cos(qy)-cos(theta)*sin(qy)*cos(fi_v-qz)))/r;
dqz=(Vt*cos(theta_t)*sin(fi_vt-qz)-v*cos(theta)*sin(fi_v-qz))/r/cos(qy);

%比例导引率
N=4;
ayc=N*v*dqy*cos(qz-fi_v);
ayc=max(min(ayc,40*g),-40*g);
azc=-N*v*dqz*cos(theta)+N*v*dqy*sin(qz-fi_v)*sin(theta);
azc=max(min(azc,40*g),-40*g);
ay=v*global_dydt(2);
az=-v*cos(theta)*global_dydt(3);
%系数插值
missile_areo=coe_interp_class(t,alpha,beta,Ma);
P=missile_areo.Thrust;
cx=missile_areo.Cx;
cy=missile_areo.Cy;
cz=missile_areo.Cz;
Xg=missile_areo.Xg;
mx_beta=missile_areo.Mx_beta;
mx_wx=missile_areo.Mx_wx;
mx_deltax=missile_areo.Mx_deltax;
mz_alpha0=missile_areo.Mz_alpha0;
my_beta0=missile_areo.My_beta0;
my_omegay=missile_areo.My_omegay;
mz_omegaz=missile_areo.Mz_omegaz;
mz_deltaz=missile_areo.Mz_deltaz;
J_z1=missile_areo.Jz;
mc=missile_areo.Mc;
cx_beta=missile_areo.Cx_beta;

X=cx*(0.5*density*v^2)*S;
Y=cy*(0.5*density*v^2)*S;
Z=cz*(0.5*density*v^2)*S;
mz_alpha=mz_alpha0+cy*(Xg-0.9381)/L;
J_x1=0.83;J_y1=J_z1;
%%动力系数

dynamic_coe=dynamic_coe_class(mz_omegaz,0.5*density*v^2,S,L,l,J_z1,J_y1,v,mz_alpha,...
                mz_deltaz,cy/alpha,P,theta,mx_deltax,mx_wx,m,J_x1,mx_beta,...
                my_beta0-cz*(Xg-0.9381)/L,mz_deltaz,my_omegay,Theta,cx_beta);
[B4,P4,E4]=dynamic_coe.Tran(alpha,Theta);
% if t>8.015
%       t=t;
% end
[delta_x,delta_y,delta_z]=control(dynamic_coe,t,global_dydt,in,ayc,azc,ay,az);

global_delta=[global_delta;delta_x,delta_y,delta_z];
M_x1=mx_beta*beta*(0.5*density*v^2)*S*l...
    +mx_wx*omega_x1*l/(2*v)*(0.5*density*v^2)*S*l...
    +mx_deltax*delta_x*(0.5*density*v^2)*S*l;
M_y1=(my_beta0-cz*(Xg-0.9381)/L)*(0.5*density*v^2)*S*L...
    +my_omegay*omega_y1*L/v*(0.5*density*v^2)*S*L...
    +mz_deltaz*delta_y*(0.5*density*v^2)*S*L;
M_z1=(mz_alpha0+cy*(Xg-0.9381)/L)*(0.5*density*v^2)*S*L...
    +mz_omegaz*omega_z1*L/v*(0.5*density*v^2)*S*L...
    +mz_deltaz*delta_z*(0.5*density*v^2)*S*L;

%%%%%%%%%微分方程
dvdt=(P*cos(alpha)*cos(beta)-X-m*g*sin(theta))/m;
dthetadt=(P*(sin(alpha)*cos(gama_v)+cos(alpha)*sin(beta)*sin(gama_v))...
    +Y*cos(gama_v)-Z*sin(gama_v)-m*g*cos(theta))/(m*v);
dfi_vdt=(P*(sin(alpha)*sin(gama_v)-cos(alpha)*sin(beta)*cos(gama_v))...
    +Y*sin(gama_v)+Z*cos(gama_v))/(-m*v*cos(theta));
domega_x1=(M_x1-(J_z1-J_y1)*omega_z1*omega_y1)/J_x1;
domega_y1=(M_y1-(J_x1-J_z1)*omega_z1*omega_x1)/J_y1;
domega_z1=(M_z1-(J_y1-J_x1)*omega_x1*omega_y1)/J_z1;
dx=v*cos(theta)*cos(fi_v);
dy=v*sin(theta);
dz=-v*cos(theta)*sin(fi_v);
dTheta=omega_y1*sin(gama)+omega_z1*cos(gama);
dfi=(omega_y1*cos(gama)-omega_z1*sin(gama))/cos(Theta);
dgama=omega_x1-tan(Theta)*(omega_y1*cos(gama)-omega_z1*sin(gama));
dm=-mc;

dydt=[dvdt;dthetadt;dfi_vdt;domega_x1;domega_y1;domega_z1;dx;dy;dz;...
    dTheta;dfi;dgama;dm;dxt;dyt;dzt];
global_dydt=dydt;
end
