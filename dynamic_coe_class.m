classdef dynamic_coe_class
    %动力系数类

    properties(Access=public)
        a22;a24;a241;a25;a33;a34;a35;
        K_M;T_1;T_M;xi_M;
        b17;b11;
    end

    methods
        function obj = dynamic_coe_class(mz_omega,q,S,L,J_z,v,mz_alpha,mz_delta,cy,P,Y_deltaz,sita,mx_deltax,mx_omegax,m,J_x)
            g = 9.8;
            obj.a22 = mz_omega*q*S*L^2/J_z/v;
            obj.a24 = mz_alpha*q*S*L/J_z;
            obj.a241 = 0;
            obj.a25 = mz_delta*q*S*L/J_z;
            obj.a34 = (cy*q*S+P)/m/v;
            obj.a35 = (Y_deltaz/m/v);
            obj.a33 = g*sin(sita)/v;
            obj.K_M = (-obj.a25*obj.a34+obj.a35*obj.a24)/(obj.a22*obj.a34+obj.a24);
            obj.T_M = 1/(-obj.a22*obj.a34-obj.a24)^0.5;
            obj.xi_M = (-obj.a22-obj.a241+obj.a34)/2/(-obj.a22*obj.a34-obj.a24)^0.5;
            obj.T_1 = (-obj.a35*obj.a241+obj.a25)/(obj.a25*obj.a34-obj.a35*obj.a24);
            obj.b17 = mx_deltax*q*S*L/J_x;
            obj.b11 = mx_omegax*q*S*L/J_x;
        end
        
        
    end
end