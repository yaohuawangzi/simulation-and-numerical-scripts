function y=Add_alpha(n,x,snr,N)

Ex=var(x);   
y=(20*10^(snr/10)/Ex)^0.5;  