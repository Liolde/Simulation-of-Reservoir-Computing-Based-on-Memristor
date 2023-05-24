function [x,t]=MackeyGlass(N,tau)
% Mackey-Glass混沌延迟微分方程 
%
% N为输出点数，tau为延迟时间
% x为序列返回值，t为时间返回值，h为采样间隔

t=zeros(N,1);
x=zeros(N,1);  
x(1)=1.2; t(1)=0; 
a=0.2;
b=0.1;
c=10;
h=5;
for k=1:N-1
  t(k+1)=t(k)+h; 
  if t(k)<tau
      k1=-b*x(k); 
      k2=-b*(x(k)+h*k1/2); 
      k3=-b*(x(k)+k2*h/2); 
      k4=-b*(x(k)+k3*h);
      x(k+1)=x(k)+(k1+2*k2+2*k3+k4)*h/6; 
  else 
      n=floor((t(k)-tau-t(1))/h+1); 
      k1=Df(x(n))-b*x(k); 
      k2=Df(x(n))-b*(x(k)+h*k1/2); 
      k3=Df(x(n))-b*(x(k)+2*k2*h/2); 
      k4=Df(x(n))-b*(x(k)+k3*h); 
      x(k+1)=x(k)+(k1+2*k2+2*k3+k4)*h/6; 
  end 
end
end

function y=Df(x)
a=0.2;
c=10;
y=a*x/(1+x^c);
end
