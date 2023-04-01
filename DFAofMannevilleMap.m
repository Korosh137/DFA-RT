%%% Korosh Mahmoodi 03/20/2023
%%% If you apply this code, pleas cite the corresponding paper: 
%%% [1] Mahmoodi K, Kerick SE, Grigolini P, West BJ, & Franaszczuk PJ, Temporal Complexity Measure of Reaction Time Series: Operational vs. Chronological Time

tic
clc;
clear all;
close all;

                          PLOT = 1 ; %% if == 1, plots the DFA graphs
                          Ens = 1 ; %%% Number of DFA
                           if Ens>1
                                       PLOT = 0 ;  
                           end
                           
                            NumberofData = 1*10^5 ; %% Number of time durations
                            Time = 10^7 ;  %% Length of the time series          
                              

                              alpha = zeros(40 , Ens) ;
for rg = 1 : Ens  %% ens
    rg
                                          
                                         miu = zeros(40 , 1) ;
                                         alphaAve  = zeros(40 , 1) ;
                           
%%% Manneville Trajactories 

for aa = 1 :  40 
     aa
      miu(aa) =    1 + aa * 0.05 ;
 
      if miu(aa) > 2
                     Time = 10^6 ;            
      end
      
      
                                   g0 = 1/(miu(aa) -1) ;  
                                   TT = 1 ; 
Start = 2 ;
  
AA = 0 ;
LL = 0 ;
while LL < Time && AA <= NumberofData
AA = AA + 1 ;

    r = rand;
    g0 = 1/(miu(aa) -1) ;  
    S = TT * (-1 + 1/(r.^g0));
    Tau(AA) = round(S) ;
    
         while Tau(AA) == 0     ||    Tau(AA) >   10^6 
            
            r = rand;
            S = TT * (-1 + 1/(r.^g0));
            Tau(AA) = round(S) ;
            
         end
    LL = LL +   Tau(AA) ;  
end

                                                  
 %%%%% END Manneville 
                           
 
                              RG =  cumsum(Tau);  
                                TimeStep = RG(length(RG)) ; 

                                          a = 2 ;
                                          k0 =  length(Tau)  ;
                                   
                                     Event = zeros(k0, 1);        
                                     Event(1) = 0 ;
                                          
                                          for hh = 1 : k0 
                                              
                                          s = Tau(hh) ;
                                          
       %%%%%%%%%%%%%%%%%%%% Creating new noise from RT time series                                  
r = rand; 
if r < 0.5
    SIGN = 1 ;
else
    SIGN = -1 ;
end 

                                          for j = a : a+s

                                          Event(j) = SIGN ; 

                                     if j >= TimeStep 
                                       break 
                                     end
                                          
                                          end
                              
                          
                                     if j >= TimeStep 
                                       break 
                                     end
                                          a = j + 1 ;
                                          
                                          end
                                          
                                XDFA = Event ; % The secondary time series created by reaction times        
                                pts =  2000:50:5000 ; 
  
               %%%%%%%%%%%  DFA of order 1.
               
plot_fun = @(xp,A,ord) polyval(A,log(xp));

Interval = floor(length(XDFA)/100) ;
[A,F] = DFA_fun(XDFA, pts, 1);

if PLOT == 1
figure;
scatter(log(pts),log(F))

hold on
x =  pts ;
plot(log(x),plot_fun(x,A),'--')
xlabel('log_{10} W'), ylabel('log_{10} F(W)');
legend(['\alpha = ' num2str( sprintf('%.3f',A(1) ))],'Location','northwest');

hold off

end

alpha(aa, rg) =   A(1)   ;

alphaAve(aa) = alphaAve(aa) + A(1)   ;

end
end  %% Ensemble

alphaAve = alphaAve / Ens ;


plot(miu, alphaAve) ;

toc



%%%%%%%%%%%%%%%  The below function is from:
%%%  https://www.mathworks.com/matlabcentral/fileexchange/67889-detrended-fluctuation-analysis-dfa
function[A,F] = DFA_fun(data,pts,order)
% -----------------------------------------------------
% DESCRIPTION:
% Function for the DFA analysis.
% INPUTS: 
% data: a one-dimensional data vector.
% pts: sizes of the windows/bins at which to evaluate the fluctuation
% order: (optional) order of the polynomial for the local trend correction.
% if not specified, order == 1;
% OUTPUTS: 
% A: a 2x1 vector. A(1) is the scaling coefficient "alpha",
% A(2) the intercept of the log-log regression, useful for plotting (see examples).
% F: A vector of size Nx1 containing the fluctuations corresponding to the
% windows specified in entries in pts.
% -----------------------------------------------------
% Checking the inputs
if nargin == 2
   order = 1; 
end
sz = size(data);
if sz(1)< sz(2)
    data = data';
end
exit = 0;
if min(pts) == order+1
        disp(['WARNING: The smallest window size is ' num2str(min(pts)) '. DFA order is ' num2str(order) '.'])
        disp('This severly affects the estimate of the scaling coefficient')
        disp('(If order == [] (so 1), the corresponding fluctuation is zero.)')
elseif min(pts) < (order+1)
        disp(['ERROR: The smallest window size is ' num2str(min(pts)) '. DFA order is ' num2str(order) ':'])
        disp(['Aborting. The smallest window size should be of ' num2str(order+1) ' points at least.'])
        exit = 1;
end
if exit == 1
    return
end
% DFA
npts = numel(pts);
F = zeros(npts,1);
N = length(data);
for h = 1:npts
    
    w = pts(h);
    
    n = floor(N/w);
    Nfloor = n*pts(h);
    D = data(1:Nfloor);
    
    y  = cumsum(D-mean(D));
    
    bin = 0:w:(Nfloor-1);
    vec = 1:w;
    
    coeff = arrayfun(@(j) polyfit(vec',y(bin(j) + vec),order),1:n,'uni',0);
    y_hat = cell2mat(cellfun(@(y) polyval(y,vec),coeff,'uni',0));
    F(h)  = mean((y - y_hat').^2)^0.5;
    
end
A = polyfit(log(pts),log(F)',1);
end

