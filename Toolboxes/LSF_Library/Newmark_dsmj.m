function [Y1 Y2 Y3]=Newmark_dsmj(K,M,C,F,dt,Y0,beta,gama)
%created at 2018.11.6
% inpit  K,M,C;F,dt; 初始状态Y0=[x0 v0] 算法参数 gama，beta  
% output 结构的位移、速度、加速度反应+结构的时程曲线
% %%%%%%%%%%%%%%%%%%%%%%%算法参数选择%%%%%%%%%%%%%%%%%%%%%%%%%%%
% gama=1/2； beta=1/4  即为平均加速度度算法
% gama=1/2； beta=1/2  即为中心差分算法
% gama=1/2； beta=1/6  即为线性加速度算法
% gama=1/2； beta=1/8  即为变加速度度算法
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 仅可计算结构线弹性状态
% 输入与输出格式：以列向量排列；F(i,t)为在第t时刻，第i自由度所受荷载
%%
% 初始化
    n=length(M); %结构自由度
    m=length(F); %时长步数
    Y1=zeros(n,m);Y2=zeros(n,m);Y3=zeros(n,m);
    Y1(:,1)=Y0(:,1);Y2(:,1)=Y0(:,2);Y3(:,1)=M\((F(1)-C*Y2(:,1)-K*Y1(:,1)));
%判断算法是否满足无条件稳定
     if gama<0.5||0.25*(0.5+gama)^2-beta<0
         disp('算法非无条件稳定');
         temp=input('输入0继续执行，非0值退出');
         if temp==0
             return
         end
     end
%积分常数
    A1=1/(beta*dt^2);A2=gama/(beta*dt);A3=1/(beta*dt);
    A4=1/(2*beta)-1;A5=dt/2*(gama/beta-2);A6=(gama/beta-1);
    B1=A4*M+A5*C;B2=A3*M+A6*C;B3=A1*M+A2*C;
%计算有效刚度矩阵
    Ke=K+A1*M+A2*C;
%循环
    for i=2:m
       %  计算t+dt时刻有效载荷向量(行)
       Fe=F(:,i)+B1*Y3(:,i-1)+B2*Y2(:,i-1)+B3*Y1(:,i-1);
       %   计算t+dt时刻结构位移
       Y1(:,i)=Ke\Fe;
       %   计算t+dt时刻结构的速度，加速度
       Y3(:,i)=A1*(Y1(:,i)-Y1(:,i-1))-A3*Y2(:,i-1)-A4*Y3(:,i-1);
       Y2(:,i)=Y2(:,i-1)+(1-gama)*dt*Y3(:,i-1)+gama*dt*Y3(:,i);
    end
return