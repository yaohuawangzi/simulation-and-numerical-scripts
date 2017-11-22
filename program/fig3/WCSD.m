% **********   WCSD ***********************
clc
clear
disp(['Ref: J. Lee and J. C. Principe, ...'])
disp(['Correntropy-based spectrum sensing for wireless microphones in ...'])
disp(['man-made noise environments in alpha noise'])
disp([' alpha£¬WCSD simulation'])
disp([])

tic
%% 
%*********************  **************************
delta = sqrt(0.5);
Gsnrdb = [-10 : 5:  0];
epsilon = 0;

%***************************************************
t = 0:0.25:64;
nSample = length(t);         % samples in signal
N = 5;

%****************** alpha noise   *******************
alpha = 1.5;
gama = 1;
beta = 0;
a=0;

%***************** thresh ******************************
thresh_max = 1/(sqrt(2*pi)) ;
thresh = (0:0.001: thresh_max )   ;

%%
%******************  calculation  *********************************

L = numel(Gsnrdb);

nloop = 1000;
hWait = waitbar(0,'please wait...');
for ii = 1:L
    d = 0;
    d_2 = 0 ;
    d_3 = 0 ;
    d_4 = 0 ;
    for jj = 1: nloop
                
        %**********  transmit signal   *******************************************
%         infoSignal = exp(-t.^2)' ; 
        infoSignal = cos(2 * pi * t /10)' ;
        txSignal = infoSignal ;
        
        %***************  attn ***********************
        attn_2 = sqrt( 10.^(Gsnrdb(ii)/10) * nSample / sum(abs(txSignal).^2));  
        
        %********** receive signal   *******************************************
        rxSignal = txSignal  *  attn_2;       % H1
        rxSignal_2 = 0 ;            % H0
        
%         %***************  attn of noise *********************
%         spow = sum(txSignal.^2);%/nSample;
%         attn = spow * 10 .^( -snrdb(ii)/10);
%         attn = sqrt(attn);                      %  attn of noise
        
        %***************** calculate noise ************************************
        RX =  zeros(nSample,N); 
        RX_2 =  zeros(nSample,N);
        for k = 1:1:N
            %*************  noise ************************ 
            noise = alpha_channel(nSample,alpha,gama,beta,a);  % alpha
            %*************  Signal + noise  ***********************
            RX(:,k) = rxSignal + noise' ;  % H1 
            RX_2(:,k) = rxSignal_2 + noise' ;  % H0
        end %for k

       %% ********* calculate  WCSD correntropy  *********************            
        %*************  noise ************************
        noise_2 = alpha_channel(nSample,alpha,gama,beta,a);  % alpha
        %*************  Signal + noise ***********************
        RX_3 = rxSignal + noise_2' ;  % H1
        RX_4 = rxSignal_2 + noise_2' ;  % H0        
       
        delta_2 = linspace(0.2,2,10)';
        %************ J ************************    
        for k = 1 : length(delta_2)  
            J1_3 = sum( exp(-1/(2 * delta_2(k).^2) * ( rxSignal - RX_3 ).^2 )) ;
            esti_J_3(k) = 1/( sqrt(2*pi) * nSample) * J1_3;
            
            J1_4 = sum( exp(-1/(2 * delta_2(k).^2) * ( rxSignal - RX_4  ).^2 )) ;
            esti_J_4(k) = 1/( sqrt(2*pi) * nSample) * J1_4;  
        end  % end of k = 1 : length(delta_2)
  

        % *****************************************   
        d_3 = d_3 + ( 1/length(delta_2) * sum(esti_J_3) > thresh);   %  H1      
        d_4 = d_4 + ( 1/length(delta_2) * sum(esti_J_4) > thresh);   %  H0    
                
    end %end of for  jj = 1: nloop

    pd_2(:,ii) = d_3;  
    pf_2(:,ii) = d_4;
    waitbar(ii/L,hWait);   
end   %end of for  ii = 1:L
pd_2 = pd_2/nloop;   % avg over 1000 simulation 
pf_2 = pf_2/nloop;   % avg over 1000 simulation

close(hWait);
toc

%% 
%***************  output result  ***************
figure(1);
xlabel('False alarm probability,{\it P_f}');
ylabel('Detection probability,{\it P_d }');
plot(pf_2(:,1),pd_2(:,1),'s:b',  ...
     pf_2(:,2),pd_2(:,2),'+:b', ...
     pf_2(:,3),pd_2(:,3),'*:b', 'LineWidth',1.5);

leg = legend('WCSD,{\it GSNR }= -15dB','WCSD,{\it GSNR }= -10dB','WCSD,{\it GSNR }= -5dB');
set(leg,'Location','SouthEast');
set(gcf, 'position', [400 400 400 300]);
grid on


hold on

%**********************  end of file  ***************************
 



