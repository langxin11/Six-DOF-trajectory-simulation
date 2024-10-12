function [delta_x,delta_y,delta_z]=control(dynamic_coe,t,global_dydt,in,ayc,azc,ay,az)
    delta_x=0;delta_y=0;delta_z=0;
     
    K_M =dynamic_coe.K_M; 
    T_M=dynamic_coe.T_M;
    xi_M=dynamic_coe.xi_M;
    T_1=dynamic_coe.T_1;
    b17=dynamic_coe.b17; 
    b11=dynamic_coe.b11;
    delta_z=fun(K_M,T_M,xi_M,T_1,t,in,global_dydt,ayc,ay,dynamic_coe);
    delta_x=fun_x(b11,b17,in(12),global_dydt(12));
    delta_y=fun_y(azc,global_dydt(3),az);
    delta_y=0;
    delta_z=max(min(delta_z,10/180*pi),-10/180*pi);
    delta_y=max(min(delta_y,10/180*pi),-10/180*pi);
   
    delta_x=max(min(delta_x,10/180*pi),-10/180*pi);
end
%%%倾斜自动驾驶仪
function delta_x=fun_x(b11,b17,gama,dgama)
%     w_c=10*pi;xi=0.7;K_ACT=1;k_ac=1;
%     K_A=w_c^2/(b17*K_ACT*k_ac);
%     K_g=(2*xi*w_c+b11)/b17*K_ACT;
%     delta_x=((0-k_ac*gama)*K_A-dgama*K_g)*K_ACT;
    %%平飞段
    K_A=-1.5422;K_g=-0.06;
    K_ACT=1;k_ac=1;
    delta_x=((0-k_ac*gama)*K_A-dgama*K_g)*K_ACT;
end
%%设定姿态自动驾驶仪
function delta_z=fun(K_M,T_M,xi_M,T_1,t,in,dydt,ayc,ay,dynamic_coe)
    tau1=T_M/5;xi=0.7;
    delta_z=0;
    global start_y
    K_Theta=1;
    tau2=(-tau1^2/T_1+2*xi*tau1-T_1)/(2*xi*tau1/T_1-2*xi_M*tau1^2/(T_1*T_M)-1+tau1^2/T_M^2);
    K_A=T_M^2/(tau1^2*tau2*K_M*K_Theta);
    K_dTheta=((2*xi*tau1*tau2+tau1^2)*K_A*K_M*K_Theta-...
    2*xi_M*T_M)/(K_A*K_M*T_1);
    
    if t>1&&t<8.015
        Theta=in(10);dTheta=dydt(10);
        %deta_z=-K_A*(ast_Theta(t)-K_dTheta*dTheta-K_Theta*Theta);
        delta_z=-0.042*(ast_Theta(t)-0.5184*dTheta-1*Theta);
    end
    if (t>=8.015)&&(start_y==false)
        v=in(1);H=in(8);dH=dydt(8);
        Theta=in(10);dTheta=dydt(10);
        exp_Theta=0.1347*(304-H)+0.068984*(0-dH);%高度自动驾驶仪得到预期俯仰角
        delta_z=-0.042*(exp_Theta-0.5184*dTheta-1*Theta);
        
    elseif start_y==true
        %设定过载自动驾驶仪
        k_s=K;
%         k1=0.0000125;
%         k2=-0.0000001;
        dTheta=dydt(10);
        delta_z=(k1*(ayc-ay)-k2*dTheta)*1;
    end
end
function delta_y=fun_y(azc,dfi_v,az)
%偏航过载驾驶仪
global start_y
k1=-0.000007;
k2=-0.00001;
k3=-0.000007;
delta_y=-1*(k1*azc-k2*dfi_v-k3*az);
if start_y==true
    delta_y=-1*(k1*azc-k2*dfi_v-k3*az);
end
end
function out=ast_Theta(t)
    t0=1;t1=5.4;t2=8.1;
    out=0;
    if t>=t0&&t<t1
        out=18*pi/180;
    elseif t>=t1&&t<t2
        out=18*exp((t1-t))*pi/180;
    else
        out=0.1*pi/180;
    end
end