clc;clear;
%1:速度-v 2:弹道倾角-theta 3:弹道偏角-gama_v 
%4:w_x-omega_x1 5:w_y-omega_y1 6:w_x-omega_z1
%7:x 8 :y 9:z
%10: 俯仰角-Theta 11:偏航角-fi 12: 倾斜角-gama
%13：质量 m
%速度倾斜角-gama_v
in_0=[20;0.1*pi;0;0;0;0;0;20;0;0.1*pi;0;0;53.38];
[t,out]=ode45(@missile,[0 30],in_0);
figure;
plot3(out(:,7),out(:,9),out(:,8),'k',LineWidth=2);
xlabel('x');ylabel("z");zlabel('y');
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
beta=asin(cos(theta)*(cos(gama)*(sin(fi-fi_v)+sin(Theta)*sin(gama)*cos(fi-fi_v)))...
    -sin(theta)*cos(Theta)*sin(gama));
alpha=asin((cos(theta)*(sin(Theta)*cos(gama)*cos(fi-fi_v)-...
    sin(gama)*sin(fi-fi_v))-sin(theta)*cos(Theta)*cos(gama))/cos(beta));
gama_v=asin((cos(alpha)*sin(beta)*sin(Theta)-...
    sin(alpha)*sin(beta)*cos(gama)*cos(Theta)+...
    cos(beta)*sin(gama)*cos(Theta))/cos(theta));
g=9.81;
density=1.225;S=0.0227;L=1.8;Ma=v/343.13;l=0.5;
P=thrust(t);
X=cx(alpha,Ma)*(0.5*density*v^2)*S;
Y=cy(alpha,Ma)*(0.5*density*v^2)*S;
dbstop if error
Z=-cy(beta,Ma)*(0.5*density*v^2)*S;
Xg=f_xg(t);
deta_x=0;
M_x1=mx_beta(abs(alpha),Ma)*beta*(0.5*density*v^2)*S*l...
    +mx_wx(Ma)*omega_x1*l/(2*v)*(0.5*density*v^2)*S*l...
    +mx_detlax(Ma)*deta_x*(0.5*density*v^2)*S*l;
deta_y=0;
% if t>0.1&&t<0.2
%     deta_y=0.1;
% end
M_y1=(mz_alpha0(beta,Ma)-cy(beta,Ma)*(Xg-0.9831)/L)*(0.5*density*v^2)*S*L...
    +mz_omgeaz(abs(beta),Ma,Xg)*omega_y1*L/v*(0.5*density*v^2)*S*L...
    +mz_detaz(Ma)*deta_y*(0.5*density*v^2)*S*L;
deta_z=0;
M_z1=(mz_alpha0(alpha,Ma)+cy(alpha,Ma)*(Xg-0.9831)/L)*(0.5*density*v^2)*S*L...
    +mz_omgeaz(abs(alpha),Ma,Xg)*omega_z1*L/v*(0.5*density*v^2)*S*L...
    +mz_detaz(Ma)*deta_z*(0.5*density*v^2)*S*L;
J_x1=0.83;J_y1=Jz(t);J_z1=Jz(t);
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
mc=m_c(t);
dm=-mc;

dydt=[dvdt;dthetadt;dfi_vdt;domega_x1;domega_y1;domega_z1;dx;dy;dz;...
    dTheta;dfi;dgama;dm];
if any(isnan(dydt))
    1;
end
end
%%
function z=cy(alpha,Ma)%升力系数插值
[X,Y]=meshgrid(0:2:10,0.1:0.1:0.9);
Z=[.0000	.6430	1.4758	2.2870	3.0713	3.8463;
.0000	.6454	1.4807	2.2942	3.0814	3.8598;
.0000	.6480	1.4858	2.3014	3.0915	3.8731;
.0000	.6512	1.4923	2.3107	3.1039	3.8891;
.0000	.6554	1.5007	2.3227	3.1197	3.9092;
.0000	.6617	1.5134	2.3409	3.1436	3.9401;
.0000	.6698	1.5304	2.3661	3.1775	3.9835;
.0000	.6792	1.5501	2.3950	3.2162	4.0323;
.0000	.6933	1.5935	2.4706	3.3273	4.1790;
];
if Ma<0.1
    Ma=0.1;
end
if Ma>0.9
    Ma=0.9;
end
if alpha>=0
    z=interp2(X,Y,Z,alpha,Ma,"linear");
else
    z=-interp2(X,Y,Z,-alpha,Ma,"linear");
end
end
%% 
function z=cx(alpha,Ma)%阻力系数插值
[X,Y]=meshgrid(0:2:10,0.1:0.1:0.9);
Z=[.4177	.4404	.5219	.6603	.8534	1.1023;
.3858	.4086	.4903	.6290	.8226	1.0723;
.3779	.4007	.4827	.6218	.8160	1.0666;
.3785	.4015	.4838	.6234	.8184	1.0700;
.3787	.4018	.4846	.6249	.8209	1.0738;
.3829	.4062	.4897	.6310	.8284	1.0835;
.3855	.4091	.4934	.6363	.8358	1.0938;
.4082	.4321	.5175	.6621	.8641	1.1254;
.4947	.5192	.6073	.7571	.9672	1.2392;
];
if Ma<0.1
    Ma=0.1;
end
if Ma>0.9
    Ma=0.9;
end
z=interp2(X,Y,Z,abs(alpha),Ma,"linear");
end
%% 
function y=thrust(t)%推力系数插值
T1=[.000	.15	.49	2.11];
P1=9.81*[331.2	614.3	505.4	607.8];
T2=[2.27	3.53	8.78	25.45	42.80	43.68  44.08];
P2=9.81*[486.5	439.7	420.1	410.0	408.0	407.9  22.2];
if t<=2.11
y=interp1(T1,P1,t,"linear",0);
elseif t<=2.27&&t>2.11
    y=486.5;
else
    y=interp1(T2,P2,t,"linear",0);
end
end
%% 
function y=f_xg(t)%导弹重心插值
T=[0.0	2.0	2.4	10.0	18.0	26.0	32.0	38.0	42.0	44.0];
Y=[.9381	.9095	.9091	.9026	.8969	.8928	.8907	.8896	.8895	.8896];
y=interp1(T,Y,t,"linear",0.8896);
end
%% 
function y=m_c(t)%质量流量插值
T=[0.	2.1	2.105	44.1	44.105	100];
dM=[2.362	2.362	0.21059	0.21059	0.	0.];
y=interp1(T,dM,t,"linear",0);
end
%% 静稳定力矩系数插值
function z=mz_alpha0(alpha,Ma)
[X,Y]=meshgrid(0:2:10,0.1:0.1:0.9);
Z=[0.0000  	-0.0104  	-0.0341  	-0.0564  	-0.0771  	-0.0985;
0.0000  	-0.0104  	-0.0341  	-0.0564  	-0.0770  	-0.0983;
0.0000  	-0.0104  	-0.0341  	-0.0564  	-0.0769  	-0.0982;
0.0000  	-0.0105  	-0.0342  	-0.0564  	-0.0768  	-0.0979;
0.0000  	-0.0104  	-0.0339  	-0.0560  	-0.0761  	-0.0969;
0.0000  	-0.0093  	-0.0314  	-0.0521  	-0.0708  	-0.0903;
0.0000  	-0.0080  	-0.0286  	-0.0477  	-0.0650  	-0.0829;
0.0000  	-0.0065  	-0.0252  	-0.0425  	-0.0578  	-0.0739;
0.0000  	-0.0053  	-0.0229  	-0.0391  	-0.0538  	-0.0693;
];
if Ma<0.1
    Ma=0.1;
end
if Ma>0.9
    Ma=0.9;
end
if alpha>=0
    z=interp2(X,Y,Z,alpha,Ma,"linear");
else
    z=-interp2(X,Y,Z,-alpha,Ma,"linear");  
end
end
%% 操纵力矩系数导数插值
function d=mz_detaz(ma)
X=[0.3 0.4 0.5 0.6 0.8];
Y=[-2.929	-3.022	-3.095	-3.2	-3.424];
if ma<0.3
    ma=0.3;
end
if ma>0.8
    ma=0.8;
end
d=interp1(X,Y,ma,"linear");
end
%% 
function d=mz_omgeaz(absalpha,ma,Xg)%阻尼力矩导数插值
[X,Y]=meshgrid(0:2:10,0.1:0.1:0.9);
Z_1=[-0.4686  	-0.4829  	-0.4982  	-0.5130  	-0.5272  	-0.5409;
-0.4707  	-0.4850  	-0.5003  	-0.5150  	-0.5292  	-0.5429;
-0.4744  	-0.4886  	-0.5039  	-0.5186  	-0.5327  	-0.5464;
-0.4797  	-0.4939  	-0.5090  	-0.5237  	-0.5378  	-0.5514;
-0.4882  	-0.5022  	-0.5173  	-0.5318  	-0.5458  	-0.5593;
-0.5089  	-0.5227  	-0.5376  	-0.5520  	-0.5658  	-0.5791;
-0.5366  	-0.5502  	-0.5649  	-0.5790  	-0.5927  	-0.6058;
-0.5738  	-0.5871  	-0.6014  	-0.6153  	-0.6287  	-0.6415;
-0.6272  	-0.6407  	-0.6553  	-0.6694  	-0.6830  	-0.6960;
];
Z_2=[-0.6179  	-0.6384  	-0.6600  	-0.6805  	-0.6999  	-0.7182;
-0.6207  	-0.6410  	-0.6626  	-0.6830  	-0.7024  	-0.7207;
-0.6253  	-0.6455  	-0.6670  	-0.6874  	-0.7067  	-0.7249;
-0.6319  	-0.6521  	-0.6734  	-0.6937  	-0.7129  	-0.7310;
-0.6424  	-0.6624  	-0.6835  	-0.7036  	-0.7226  	-0.7406;
-0.6669  	-0.6866  	-0.7074  	-0.7272  	-0.7459  	-0.7636;
-0.6997  	-0.7190  	-0.7395  	-0.7589  	-0.7774  	-0.7948;
-0.7435  	-0.7624  	-0.7824  	-0.8014  	-0.8194  	-0.8365;
-0.8069  	-0.8266  	-0.8474  	-0.8672  	-0.8859  	-0.9035;
];
if ma<0.1
    ma=0.1;
end
if ma>0.9
    ma=0.9;
end
d=(interp2(X,Y,Z_1,absalpha,ma,"linear")*(Xg-0.8896)+interp2(X,Y,Z_2,absalpha,ma,"linear")*(0.9381-Xg))/(0.938-0.8896);
end
%% 
function y=Jz(t)%转动惯量
T=[.0	2.0	2.4	6.4	10.4	14.4	18.4	22.4	26.4	30.4	34.0	38.4	42.4	44.0];
Y=[8.35	7.88	7.86	7.81	7.78	7.75	7.73	7.71	7.70	7.70	7.69	7.69	7.69	7.69];
y=interp1(T,Y,t,'linear',7.69);
end
%%滚动阻尼力矩系数导数
function y=mx_wx(ma)
X=[0.4:0.1:0.8];
Y=[-4.9881 -5.1414 -5.3397 -5.6206 -6.1041];
if ma<0.4
    ma=0.4;
end
if ma>0.8
    ma=0.8;
end
y=interp1(X,Y,ma,"linear");
end
%%滚转操纵力矩系数导数
function y=mx_detlax(ma)
X=[0.1:0.1:0.8];
Y=[-1.1708 -1.1786 -1.1918 -1.2110 -1.2428 -1.2848 -1.3436 -1.4273];
if ma<0.1
    ma=0.1;
end
if ma>0.8
    ma=0.8;
end
y=interp1(X,Y,ma,"linear");
end
%%横向静稳定力矩系数导数
function y=mx_beta(alpha,Ma)
[X,Y]=meshgrid(0:2:10,0.4:0.1:0.8);
if Ma<0.4
    Ma=0.4;
end
if Ma>0.8
    Ma=0.8;
end
Z=[0 -0.0452 -0.0904 -0.1356 -0.1808 -0.2260;
    0 -0.0479 -0.0958 -0.1437 -0.1916 -0.2395;
    0 -0.0516 -0.1033 -0.1549 -0.2066 -0.2582;
    0 -0.0566 -0.1133 -0.1699 -0.1929 -0.2266;
    0 -0.0643 -0.1286 -0.1929 -0.2572 -0.3214];
y=interp2(X,Y,Z,alpha,Ma,'linear');
end