clc
clear all;
close all;

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

EC = zeros(32, 1) ; % Error of Commission
    
for par = 1 : 32 % participants
% for ss = 1 : 6 % sessions
%     for bl = 1 : 4 % blocks
    
                Partic = par + 19 ;

ErrorCommision= 0;
CorrectOmission=0;
for tt = 1 : length(rtraw) 
    
    % Add the limits in which you want to evaluate the Error of Commision
if id(tt) == Partic && strcmp(data0(tt,6),'H') %% && SESS(tt)== ss && BLOCK(tt) ==bl  %...
    % For Low Stress condition change 'H' to 'L'
    
if strcmp(data0(tt,7) ,'TargetUpFriendly')  
    
    if shot(tt) == 1  &&    rtraw(tt) >= 100   &&   rtraw(tt) < 1500
        ErrorCommision = ErrorCommision +1 ;
    end
    if shot(tt) == 0 
          CorrectOmission = CorrectOmission +1 ;
    end
    
end


end
end
      

%     end  %% end block
    
   
% end     %%% end session

ErrorCommision = ErrorCommision /( ErrorCommision + CorrectOmission) ;
EC(par, 1) =   ErrorCommision      ;


end  %%% end participant



toc