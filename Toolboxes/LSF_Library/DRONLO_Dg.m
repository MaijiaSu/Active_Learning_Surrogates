function Dg = DRONLO_Dg(u)
% dynamic response of a non-linear oscillator
% input m*n ��mΪ���ݵ�����nΪx��ά��
% fun_par =[] ������ 
% output yΪ������

% ����׼��̬����uת��

    
 [Ns,Ndim] = size(u);
  X = u;  
  m=X(1);c1=X(2);c2=X(3);
  r=X(4);F1=X(5);t=X(6);
  w0=sqrt((c1+c2)./m);
  
  % ����1: Dg/Dm
  Dg(1,1) = 2*F1/(c1+c2)*cos(w0*t/2)*1/2*sqrt(m/(c1+c2))*(c1+c2)/m^2;
  
  % ����2: Dg/Dc1
  Dg(1,2) = 4*F1/m/w0^3*0.5*sqrt(m/(c1+c2))*sin(w0*t/2) -...
               2*F1/m/w0^2*cos(w0*t/2)*1/2*1/2*sqrt(m/(c1+c2));
           
  % ����3: Dg/Dc2
  Dg(1,3) = 4*F1/m/w0^3*0.5*sqrt(m/(c1+c2))*sin(w0*t/2) -...
               2*F1/m/w0^2*cos(w0*t/2)*1/2*1/2*sqrt(m/(c1+c2));          
  
  % ����4: Dg/Dr
  Dg(1,4) = 3;
  
  % ����5: Dg/DF1
  Dg(1,5) = -2/m/w0^2*sin(w0*t/2);
 
  % ����5: Dg/Dt1
  Dg(1,6) = -2*F1/m/w0^2*cos(w0*t/2)*w0/2;
  
  
%   g=3*r-abs(2*F1./m./w0.^2.*sin(w0.^1.*t1/2));
% ���ܺ������о���ֵ���ж��Ƿ�Ҫ�ı����
  S = 2*F1./m./w0.^2.*sin(w0.^1.*t/2); 
  if S<0
      Dg(1,[1,2,3,5,6]) = Dg(1,[1,2,3,5,6])*-1;
  end
  
  Dg = Dg(:);
return