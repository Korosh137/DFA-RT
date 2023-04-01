%%% Korosh Mahmoodi 03/20/2023
%%% If you apply this code, pleas cite the corresponding paper: 
%%% [1] Mahmoodi K, Kerick SE, Grigolini P, West BJ, & Franaszczuk PJ, Temporal Complexity Measure of Reaction Time Series: Operational vs. Chronological Time

tic

clc
clear all;
close all;
                                Ens = 1 ; %%% Number of DFA. We used 100 
                                PLOT = 1 ; %% if == 1, plots the DFA graphs
                                if Ens>1
                                       PLOT = 0 ;  
                                end
                                 Maxplot = 20 ; %%% Maximum number of DFA graphs

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Import data
%%%% From the same repository in which this code is, download the complete_behavior_SK.csv file and import it as "X"

X = readtable('C:\Users\koros\Downloads\DFA_Crucial\complete_behavior_SK.csv');

data0 = table2cell(X); % convert table to cell structure

LL = length(data0) ;  %%% Length of data

ses = [data0{:,4}]; % session number

bloc = [data0{:,5}] ; % block number 

shot = [data0{:,10}] ; % target fired upon

hit = [data0{:,11}] ; % target hit

traw =  [data0{:,9}]; % raw ITIs

rtraw =  [data0{:,12}]; % raw RTs

id =  [data0{:,2}]; % subject IDs

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% End of Import data 

rt0Data = zeros(LL , 1);
rtraw(isnan(rtraw))=13 ;
Matt = zeros(32, 4) ; 
row =  1 ;

LIMITDATA = 70 ; %% datapoints valuable in block of 90 trials

    countplot = 0 ;
for par = 1 : 32
      par
      for ss = 1 :   6

                ID = par + 19 ;
                SIGN = 1 ;

ii = 0;
for tt = 1 : length(traw) 
if id(tt) == ID   && strcmp(data0(tt,6),'H') && ses(tt)== ss % && BLOCK(tt) ==bl 
    % For High Stress condition change 'L' to 'H'
    ii = ii +1 ;
    rt0Data(ii) = rtraw(tt) ;
    
if  strcmp(data0(tt,6) , "H") 
    Stress = 2;
end
% 
if  strcmp(data0(tt,6) , "L") 
        Stress = -2;
end


    
end
end
       
%%%%% missing data

rtData = zeros(ii, 1);
for jj = 1: ii
rtData(jj) = rt0Data(jj); 
end

                               Tau  = round(rtData ) ;
                               Tau = Tau(Tau~=13);
                               Tau = Tau(Tau ~=0);
                                      Tau = Tau(Tau >=100);
                                      Tau = Tau(Tau <=1500);

                             
for yu = 1 : Ens  %%% ensemble 
                                alpha = 0 ;
                                Stress = 0 ;
if length(Tau) > LIMITDATA
    
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
                                
%If you want to use DFA for original Reaction time series replace the last two commands with:
%                                 XDFA = Tau ; 
%                                 pts =  10:10:100 ; 
%%%% and Set Ens = 1;

                                
          
plot_fun = @(xp,A,ord) polyval(A,log(xp));
[A,F] = DFA_fun(XDFA,pts , 1);

alpha = A(1) ;

if PLOT == 1  && countplot < Maxplot
    countplot = countplot + 1 ;
figure

scatter(log(pts),log(F))

hold on
x = pts ; 
plot(log(x),plot_fun(x,A),'--')
xlabel('log_{10} W'), ylabel('log_{10} F(W)');
legend(['\alpha = ' num2str( sprintf('%.3f',alpha ))],'Location','northwest');

hold off
end
                            
end  %%% if data exist in the selected

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% MATRIX
Matt(row, 1) = par + 19 ;
Matt(row, 2) = ss ;
Matt(row, 3) = Stress ;
Matt(row, 4) = Matt(row, 4) + alpha ;

end %%% end ens
row = row + 1 ;

      end    %%% end session
end  %%% end participant

 Matt(:, 4) = Matt(:, 4) / Ens ;


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

