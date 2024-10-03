function [deta_x,deta_y,deta_z]=control(dynamic_coe,t,Theta,dTheta)
    deta_x=0;deta_y=0;deta_z=0;
    a22 =dynamic_coe.a22;  
    a24 =dynamic_coe.a24; 
    a241=dynamic_coe.a241; 
    a25 =dynamic_coe.a25; 
    a34 =dynamic_coe.a34; 
    a35=dynamic_coe.a35; 
    a33=dynamic_coe.a33; 
    K_M =dynamic_coe.K_M; 
    T_M=dynamic_coe.T_M;
    xi_M=dynamic_coe.xi_M;
    T_1=dynamic_coe.T_1;
    b17=dynamic_coe.b17; 
    b11=dynamic_coe.b11;
    deta_z=fun(K_M,T_M,xi_M,T_1,t,Theta,dTheta);
end
function deta_z=fun(K_M,T_M,xi_M,T_1,t,Theta,dTheta)
    %%设定姿态自动驾驶仪
    a=1;tau1=T_M/0.6;xi=0.7;
    deta_z=0;
    K_Theta=1;
    tau2=(-tau1^2/T_1+2*xi*tau1-T_1)/(2*xi*tau1/T_1-2*xi_M*tau1^2/(T_1*T_M)-1+tau1^2/T_M^2);
    K_A=T_M^2/(tau1^2*tau2*K_M*K_Theta);
    K_dTheta=((2*xi*tau1*tau2+tau1^2)*K_A*K_M*K_Theta-...
    2*xi_M*T_M)/(K_A*K_M*T_1);
    if t>1&&t<8
    deta_z=K_A*(ast_Theta(t)-K_dTheta*dTheta-K_Theta*Theta);
    end
end
function out=ast_Theta(t)
    t0=1;t1=4.85;t2=7.22;t3=8;
    out=0;
    if t>=t0&&t<t1
        out=16;
    elseif t>=t1&&t<t2
        out=16*exp((t1-t)/0.6);
    elseif t>=t2&&t<t3
        out=0;
    end    
end