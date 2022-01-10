  % <SIMULATED ANNEALING CODE BY HEAT TRANSFER COEFFICIENTS FROM HEAT EQUATION>
  %  Copyright (C) <2022>  <Enrique Cárdenas Sánchez>

  %  This program is free software: you can redistribute it and/or modify
  %  it under the terms of the GNU General Public License as published by
  %  the Free Software Foundation, either version 3 of the License, or
  %  (at your option) any later version.

   % This program is distributed in the hope that it will be useful,
   % but WITHOUT ANY WARRANTY; without even the implied warranty of
   % MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   % GNU General Public License for more details.

   % You should have received a copy of the GNU General Public License
   % along with this program.  If not, see <https://www.gnu.org/licenses/>.
    

function [H_Best,Lambda_Best, Error, Tmodel_Best]=SimulatedAnealing_MECV(t,Tdata,tol,Seed)
%% data for simulation

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
     Lambda=Seed(1);
     H=Seed(2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 % Tolerancia .... del 3.5% 
   %  tol=0.025*Tdata;

 %Definiendo variables para H
hd=200;
hc=0.026;
%Definendo valores para Lambda
ld=10; %m
lc=0.01;   %m

%%%%% Variables del simulhco recocido.
TPmax=6000;
Tp=TPmax; 
NL=7;
a_size=linspace(0.01,10,NL); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Cp =@(T) 550.106+0.686*T+(-262.558)./T.^2;
K =@(T) 1./(0.3666+T*2*10^-4);
% Función de costo  
CostFunction=@(T) abs(T-Tdata);
% Funcion para h   
h=@(T,H0) H0/K(T);
  


   weight=poisspdf(1:NL,Lambda);
   VecInput=[Tdata, t, H]; 
   
 for ii=1:NL
   alpha(ii,:)=root_alpha(a_size(ii),h(Tdata,H),10);
   Temperature(ii)=Temperature_model2D(VecInput,a_size(ii),alpha(ii,:));
 end
 
 Tmodel=sum(weight.*Temperature);
 cur_cost=CostFunction(Tmodel);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
%%% Para esta simulhcion se considera que la maxima temperatura de experimento es constante e igual para toda las
%%% muestras. Tdata es el valor de la temperatura en el tiempo t que es menor que Tmax siempre

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     
                    %Finding new variables
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

   new_Lambda=(ld-lc)*rand+lc;
   new_H=(hd-hc)*rand+hc;
   weight=poisspdf(1:NL,new_Lambda);
   VecInput=[Tdata, t, new_H];

 for ii=1:NL
   alpha(ii,:) = root_alpha(a_size(ii),h(Tdata,H),10);
   Temperature(ii) = Temperature_model2D(VecInput,a_size(ii),alpha(ii,:));
 end
 
 Tmodel=sum(weight.*Temperature);
 new_cost=CostFunction(Tmodel);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

   DE=new_cost-cur_cost; % si DE es negativo curcost es mejor que newcost

if DE<0
    Lambda_Best=Lambda;
    H_Best=H;
    Error=new_cost;  
else
  new_corcost=cur_cost;
  Lambda_Best=Lambda;
  H_Best=H;
  Error=cur_cost;
end

%nuevo contador
kk=1;


%    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 while Error>tol %El valor 1 representa 3.5% de grhco de erro de la Chicuhcrhca


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%  choosing new solution %%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

   Lambda=(ld-lc)*rand+lc;
   H=(hd-hc)*rand+hc;  
   weight=poisspdf(1:NL,Lambda);
   VecInput=[Tdata, t, H];

   for ii=1:NL
       alpha(ii,:)=root_alpha(a_size(ii),h(Tdata,H),10);
       Temperature(ii)=Temperature_model2D(VecInput,a_size(ii),alpha(ii,:));
   end
 
  Tmodel=sum(weight.*Temperature);
  cur_cost=CostFunction(Tmodel);
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
     
%    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   
   if (DE<0) 
      Lambda_Best=Lambda;
      H_Best=H; 
      Error=new_cost;  
   end
      Fa=(1-exp(-DE/TPmax))/2;
     % Existe una probabilidad de aceptar la última solución.
   if (Fa>rand/TPmax)
      Lambda_Best=Lambda;
      H_Best=H; 
      Error=cur_cost; 
      %kk=kk+1; 
   end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  

      Tp=Tp*0.98;%1/log(1+ii);
   if new_cost<cur_cost
      Lambda_Best=Lambda;
      H_Best=H;
      Error=new_cost;  
   else
      new_cost=cur_cost;
      Lambda_Best=Lambda;
      H_Best=H;
      Error=cur_cost;
   end
  kk=kk+1;
 end
  Tmodel_Best=Tmodel;       
end
