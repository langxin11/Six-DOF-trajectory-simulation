classdef dynamic_coe_class
    %动力系数类
    properties(Access=public)
        a22;a24;a241;a25;a33;a34;a35;
        K_M;T_1;T_M;xi_M;
        b17;b11;b14;b18;
        b22;b24;b25;b28;
        b32;b34;b36;b38
    end

    methods
        function obj = dynamic_coe_class(mz_omegaz,q,S,L,l,J_z,J_y,v,mz_alpha,...
            mz_delta,cy_alpha,P,theta,mx_deltax,mx_omegax,m,J_x,mx_beta,...
            my_beta,my_deltay,my_omegay,Theta,Cx_beta)
            g = 9.8;
            obj.a22 = mz_omegaz*q*S*L^2/J_z/v;
            obj.a24 = mz_alpha*q*S*L/J_z;
            obj.a241 = 0;
            obj.a25 = mz_delta*q*S*L/J_z;
            obj.a34 = (cy_alpha*q*S+P)/m/v;
            obj.a35 = (mz_delta*q*S/m/v/100);
            obj.a33 = g*sin(theta)/v;
            obj.K_M = (-obj.a25*obj.a34+obj.a35*obj.a24)/(obj.a22*obj.a34+obj.a24);
            obj.T_M = 1/(-obj.a22*obj.a34-obj.a24)^0.5;
            obj.xi_M = (-obj.a22-obj.a241+obj.a34)/2/(-obj.a22*obj.a34-obj.a24)^0.5;
            obj.T_1 = (-obj.a35*obj.a241+obj.a25)/(obj.a25*obj.a34-obj.a35*obj.a24);
            obj.b17 = mx_deltax*q*S*l/J_x;
            obj.b11 = mx_omegax*q*S*l^2/J_x/v;
            obj.b14 = mx_beta*q*S*l/J_x;
            obj.b18=1/J_x;
           
            obj.b22=my_omegay*q*S*L/J_y;
            obj.b24=my_beta*q*S*L/J_y;
           
            obj.b25=my_deltay*q*S*L/J_y;
            obj.b28=1/J_y;
            obj.b32=-cos(theta)/cos(Theta);
            
            obj.b34=(P-Cx_beta*q*S)/(m*v);
            obj.b36=-g*cos(Theta)/v;
            obj.b38=-1/m/v;
        end
        function [B,P,D,E]=Tran(obj,alpha,Theta)
            B(1)=-obj.b25;
            B(2)=-obj.b25*(obj.b34-obj.a33-obj.b11);
            B(3)=-obj.b25*(-obj.b11*(obj.b34-obj.a33)-alpha*obj.b14);
            B(4)=-obj.b36*obj.b25*obj.b14;
            P(1)=-obj.a33-obj.b22+obj.b34-obj.b11;
            P(2)=-obj.b22*obj.b34+obj.a33*obj.b22+obj.b22*obj.b11-obj.b34*obj.b11+...
                obj.b11*obj.a33+obj.b24*obj.b32+obj.b24*tan(Theta)-obj.b14;
            P(3)=alpha*(obj.b22*obj.b14-obj.b24*tan(Theta)*obj.b11)-obj.b36*(obj.b24*tan(Theta)-obj.b14)...
                +obj.b11*obj.b34*obj.b22-obj.a33*obj.b22-obj.b24*obj.b32*obj.b11;
            P(4)=-obj.b36*(obj.b22*obj.b14-obj.b24*obj.b11*tan(Theta));
            E(1)=0;
            E(2)=obj.b25*tan(Theta);
            E(3)=obj.b25*(-tan(Theta)*(obj.b34-obj.a33)+tan(Theta)*obj.b11);
            E(4)=obj.b25*((obj.b34)*tan(Theta)*obj.b11-obj.b32*obj.b14);
            D(1)=0;
            D(2)=-obj.b25*(-alpha*tan(Theta)-obj.b32);
            D(3)=-obj.b25*(obj.b36*tan(Theta)-obj.b11*(-alpha*tan(Theta)-obj.b32));
            D(4)=obj.b25*obj.b36*tan(Theta)*obj.b11;
        end   
    end
end