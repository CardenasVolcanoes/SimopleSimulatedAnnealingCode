  % <Solve No-linear Equation by Newton-Rapson Method>
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
   
% Encuentra las n raices para a y h
 function alpha = root_alpha (a,h,n)
aa=0.95*pi/a;
alpha=linspace(1,n);
funy =@(x)  a*x*cot(a*x)+a*h-1;
  kk=1;
     for N=1:n
         x0=N*aa;
        [x, ~] = fzero(funy,x0);
        alpha(kk)=x; 
%        ffval(kk)=fval;
        kk=kk+1;
     end   
 end
      
