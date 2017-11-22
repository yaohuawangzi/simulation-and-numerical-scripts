
% *********************** MaxCorrentyopy ******************  
% 
% ROC

clc
clear
disp(['The simulation of MaxCorrentyopy  in alpha noise '])
disp([])
disp([])

tic
%% 
%**********************************************
delta = sqrt(0.5);  %% 
Gsnrdb = [ -10 ];
epsilon = 0;

%*********************************************
t = 0:0.25:64;
nSample = length(t);         % samples in signal
N = 5;

%****************** alpha noise**************
alpha = 1.5 ;
gama = 1;
beta = 0;
a=0;

%***************** thresh ******************************
thresh_max = 1/(sqrt(2*pi)*delta)  ;
thresh = (0:0.0005: thresh_max)   ;

%%
%******************  calculation  *********************************

L = numel(Gsnrdb);

nloop = 1000;
hWait = waitbar(0,'please wait...');
for ii = 1:L
    d = 0;
    d_2 = 0 ;
    for jj = 1: nloop
                
        %**********  transmit signal   *******************************************
%         infoSignal = exp(-t.^2)' ;   
%         txSignal = infoSignal  * Gsnrdb(ii) ;
        infoSignal = cos(2 * pi * t /10)' ;
        txSignal = infoSignal ;
        
        %***************  attn ***********************
        attn_2 = sqrt( 10.^(Gsnrdb(ii)/10) * nSample / sum(abs(txSignal).^2)); 
                
        %**********  receive signal   *******************************************
        rxSignal = txSignal  *  attn_2;       % H1
        rxSignal_2 = 0 ;            % H0
        
%         %***************  attn of noise *********************
%         spow = sum(txSignal.^2);%/nSample;
%         attn = spow * 10 .^( -snrdb(ii)/10);
%         attn = sqrt(attn);                      %    attn of noise
        
        %***************** calculate noise ************************************
        RX =  zeros(nSample,N); 
        RX_2 =  zeros(nSample,N);
        for k = 1:1:N
            %*************   noise  ************************
            noise = alpha_channel(nSample,alpha,gama,beta,a) ;  % alpha
            %*************  Signal + noise  ***********************
            RX(:,k) = rxSignal + noise' ;  % H1
            RX_2(:,k) = rxSignal_2 + noise' ;  % H0
            
            %*****************************************************
            RX(1:nSample/3,k) = noise(1:nSample/3)';  %30%          
%             RX(1:nSample*2/3,k) = noise(1:nSample*2/3)';  %50%
        end %for k
              
        P1 = -1* ones(1,nSample);  %P �ĳ�ֵ
        P = diag(-P1);
        P_2 = diag(-P1);                
       %% ********* calculate correntropy in H1 *************************************
        for  i = 1:1:1000   %1000 
            %************  kxi **************
            A = RX' * P * RX;            
            kxi = inv(A) * ( RX' * P * txSignal - epsilon/2)  ;  
            
            %************ J *****************
            J1 = sum( exp(-1/(2 * delta.^2) * ( txSignal - RX * kxi ).^2 )) - epsilon * sum( abs (kxi) );
            esti_J(i) = 1/( sqrt(2*pi) * delta ) * J1;      
            
            %**************  H1  **********************            
            if  ( i > 1 & ( esti_J(i) - esti_J(i-1) <= 0.0001 ) ) 
                esti_J = esti_J(i)/nSample;
                break;    
            else
            end
            %************ P *****************
            P = -exp(-(txSignal - RX * kxi).^2 /(2*delta.^2));    
            P = diag( -P ) ;    
                       
        end
        
        % ************************************   
        d = d + ( esti_J > thresh);         %  H1
        
       %% ********* calculate correntropy in H0 *************************************
        for  j = 1:1:1000   %1000 
            %************ kxi **************
            A_2 = RX_2' * P_2 * RX_2;
            kxi_2 = inv(A_2) * ( RX_2' * P * txSignal - epsilon/2) ;  
           
            %*********** J *****************
            J1_2 = sum( exp(-1/(2 * delta.^2) * ( txSignal - RX_2 * kxi_2 ).^2 )) - epsilon * sum( abs (kxi_2));
            esti_J_2(j) = 1/( sqrt(2*pi) * delta ) * J1_2;  
            
            %**************  H0  **********************
            if  ( j > 1 & (abs( esti_J_2(j-1) - esti_J_2(j)) <= 2 )  )
                esti_J_2 = esti_J_2(j)/nSample;
                break;      
            elseif  ( esti_J_2(j) < 5 )
                esti_J_2 = esti_J_2(j)/nSample;
                break;
            else
            end
            
            %***********P *****************
            P_2 = -exp(-(txSignal - RX_2 * kxi_2).^2 /(2*delta.^2)) ;    
            P_2 = diag( -P_2 ) ;    

        end
        
        % ***********************************   
        d_2 = d_2 + ( esti_J_2 > thresh);   %  H0        
    end %end of for  jj = 1: nloop
    pd(:,ii) = d;  
    pf(:,ii) = d_2;
    waitbar(ii/L,hWait);   
end   %end of for  ii = 1:L
pd = pd/nloop;   % avg over 1000 simulation 
pf = pf/nloop;   % avg over 1000 simulation

close(hWait);
toc
%% 
%***************  output result  ***************
hold on
figure(1);
plot(pf(:,1),pd(:,1),'o-b','LineWidth',1.5); %  S + * o
xlabel('False alarm probability,{\it P_f}��{\it GSNR} = -10dB');
ylabel('Detection probability,{\it P_d }');
leg = legend('100%','70%','40%');

set(leg,'Location','SouthEast');
set(gcf, 'position', [400 400 400 300]);
grid on
hold on


%**********************  end of file  ***************************
 



