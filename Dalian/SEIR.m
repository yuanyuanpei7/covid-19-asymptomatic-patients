classdef SEIR
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        t;
        p; %pupulation
        ini;%ini=[S0 E0 I0 R0]
    end
    
    methods
        function obj = SEIR(t,p,ini)
            %UNTITLED Construct an instance of this class
            %   Detailed explanation goes here
            obj.t = t;
            obj.p = p;
            obj.ini=ini;
   
        end
        
      
        function dp = dynamS(obj,tx,x,para)  % tx is time
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            decays=para(10);

            c1=para(1);
            c2=para(2);
            Pi=para(3);
            beta=para(4)*exp(-decays*tx);
            alpha1=para(5);
            alpha2=para(6);
            gama1=para(7);
            gama2=para(8);
            N=para(9);
            
            S=N-x(2)-x(3)-x(4)-x(5); 

            %S=x(1);
            E=x(2);
            I=x(3);
            A=x(4);
            R=x(5);
            
            %N=S+E+I+R;
            
            dp=zeros(5,1);
            %%% dS/dt =  -c1*beta*S*Pi*E/N-c2*beta*S*(1-Pi)*E/N
            dp(1) = -c1*beta*S*Pi*E/N-c2*beta*S*(1-Pi)*E/N;     
            %%%% dE/dt = c1*beta*S*Pi*E/N+c2*beta*S*(1-Pi)*E/N-alpha1*Pi*E-alpha2*(1-Pi)*E;
            dp(2) =c1*beta*S*Pi*E/N+c2*beta*S*(1-Pi)*E/N-alpha1*Pi*E-alpha2*(1-Pi)*E;   
            %%%%  dI/dt = alpha1*Pi*E - gama1*I;
            dp(3) = alpha1*Pi*E - gama1*I;  
            %%%% dA/dt = alpha2*(1-Pi)*E-gama2*A;
            dp(4) = alpha2*(1-Pi)*E-gama2*A; 
            %%%% dR/dt = gama*I; 
            dp(5) = gama1*I+ gama2*A;               
            
            dp;   
        end
        
        function dy=loss(obj,para)

            SEIRin=SEIR(obj.t,obj.p,obj.ini);  %%%嵌套对类进行对象化
            [t,p]=ode45(@SEIRin.dynamS,obj.t,obj.ini,[],para); %%%对 对像化的类进行积分
            Iactural=obj.p(:,1);
            dyI=p(:,3)-Iactural;
            dyA=p(:,4)-obj.p(:,3);
            
            dyR=p(:,5)-obj.p(:,2);
            %dy=5*dyI'*dyI;
            %dy=5*dyI'*dyI+dyR'*dyR;
            %dy=dyI'*dyI+dyA'*dyA+dyR'*dyR;
            dy=dyI'*dyI+dyA'*dyA;%%%损失函数：（真实值-预测值）的平方求和
        end
        
    end
end

