clc;clear;
%1:速度-v 2:弹道倾角-theta 3:弹道偏角-gama_v 
%4:w_x-omega_x1 5:w_y-omega_y1 6:w_x-omega_z1
%7:x 8 :y 9:z
%10: 俯仰角-Theta 11:偏航角-fi 12: 倾斜角-gama
%13：质量 m
%速度倾斜角-gama_v
global Beta Alpha Gama_v
Beta=[];
Alpha=[];
Gama_v=[];
in_0=[20;0.1*pi;0;0;0;0;0;20;0;0.1*pi;0;0;53.38;0];
[t,out]=ode45(@missile,[0 12],in_0);
figure;
%plot3(out(:,7),out(:,9),out(:,8),'k',LineWidth=2);
plot(out(:,7),out(:,8),'k',LineWidth=2);
xlabel('x');ylabel("y");
grid on;
%figure;
%plot(t,out(:,1));

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
dTheta_dt=in(14);
global Beta Alpha Gama_v
beta=asin(cos(theta)*(cos(gama)*(sin(fi-fi_v)+sin(Theta)*sin(gama)*cos(fi-fi_v)))...
    -sin(theta)*cos(Theta)*sin(gama));
Beta=[Beta;t,beta];

alpha=asin((cos(theta)*(sin(Theta)*cos(gama)*cos(fi-fi_v)-...
    sin(gama)*sin(fi-fi_v))-sin(theta)*cos(Theta)*cos(gama))/cos(beta));
Alpha=[Alpha;t,alpha];
gama_v=asin((cos(alpha)*sin(beta)*sin(Theta)-...
    sin(alpha)*sin(beta)*cos(gama)*cos(Theta)+...
    cos(beta)*sin(gama)*cos(Theta))/cos(theta));
Gama_v=[Gama_v;t,gama_v];
g=9.81;
density=1.225;S=0.0227;L=1.8;Ma=v/343.13;l=0.5;
%系数插值
missile_areo=coe_interp_class(t,alpha,beta,Ma);
P=missile_areo.Thrust;
cx=missile_areo.Cx;
cy=missile_areo.Cy;
cz=missile_areo.Cz;
Xg=missile_areo.Xg;
mx_beta=missile_areo.Mx_beta;
mx_wx=missile_areo.Mx_wx;
mx_detlax=missile_areo.Mx_detlax;
mz_alpha0=missile_areo.Mz_alpha0;
mz_beta0=missile_areo.Mz_beta0;
mz_omgeaz=missile_areo.Mz_omgeaz;
mz_detaz=missile_areo.Mz_detaz;
J_z1=missile_areo.Jz;
mc=missile_areo.Mc;

X=cx*(0.5*density*v^2)*S;
Y=cy*(0.5*density*v^2)*S;
Z=cz*(0.5*density*v^2)*S;
mz_alpha=mz_alpha0+cy*(Xg-0.9381)/L;
J_x1=0.83;J_y1=J_z1;
%%动力系数
dynamic_coe=dynamic_coe_class(mz_omgeaz,0.5*density*v^2,S,L,J_z1,v,mz_alpha,...
mz_detaz,cy,P,0,theta,mx_detlax,mx_wx,m,J_x1);
deta_x=0;deta_y=0;deta_z=0;
[deta_x,deta_y,deta_z]=control(dynamic_coe,t,Theta,dTheta_dt);

M_x1=mx_beta*beta*(0.5*density*v^2)*S*l...
    +mx_wx*omega_x1*l/(2*v)*(0.5*density*v^2)*S*l...
    +mx_detlax*deta_x*(0.5*density*v^2)*S*l;

M_y1=(mz_beta0-cz*(Xg-0.9381)/L)*(0.5*density*v^2)*S*L...
    +mz_omgeaz*omega_y1*L/v*(0.5*density*v^2)*S*L...
    +mz_detaz*deta_y*(0.5*density*v^2)*S*L;

M_z1=(mz_alpha0+cy*(Xg-0.9381)/L)*(0.5*density*v^2)*S*L...
    +mz_omgeaz*omega_z1*L/v*(0.5*density*v^2)*S*L...
    +mz_detaz*deta_z*(0.5*density*v^2)*S*L;

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
dTheta2_dt2=-dTheta_dt+dthetadt;
dydt=[dvdt;dthetadt;dfi_vdt;domega_x1;domega_y1;domega_z1;dx;dy;dz;...
    dTheta;dfi;dgama;dm;dTheta2_dt2];
% if any(isnan(dydt))
%     error('Array contains NaN, stopping execution.');
% end
end
