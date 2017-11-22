function y=alpha_channel(N,alpha,gama,beta,a)
%alpha=1.5; 
%gama=1; 
%beta=0;  
%a=0; 
%N=2000;
%N=length(signal);

v=-0.5*pi+pi*rand(1,N);       
w=exprnd(1,1,N);               
b=sin(alpha*v)./(cos(v).^(1/alpha));       
c=(cos((1-alpha)*v)./w).^((1-alpha)/alpha);   
x=b.*c;      
y=x;

