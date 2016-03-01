function [c,rk,q,l,g,e,dep,vent,etime_fv,o,oh,m,o2]=diffr(c,rk,q,l,g,e,dep,vent,etime_fv,o,oh,m,o2,varargin);

% evaluation of the generation (q), loss (l) and net growth (g)
% where g = dc/dt for each of the n species given in concentration
% array c.  rk is the array of reaction rates, and emissions, ventilation,
% and deposition given in this subroutine


 persistent dadj i ; 

 if isempty(dadj), dadj=0; end;
 if isempty(i), i=0; end;

 % use pssa to find radical concentrations

 o = rk(1).*c(1)./(rk(2).*m.*o2);
 oh = rk(7).*c(9).*c(2)./(rk(4).*c(4)+rk(5).*c(5)+rk(10).*c(1));

 % determine source rates for each species

 q(1) = e(1) + rk(3).*c(2).*c(3) + rk(7).*c(9).*c(2) + rk(8).*c(10).*c(2) + rk(9).*c(8).*c(2) + rk(12).*c(7);
 q(2) = e(2) + rk(1).*c(1);
 q(3) = rk(2).*o.*o2.*m;
 q(4) = e(4);
 q(5) = e(5) + rk(8).*c(10).*c(2);
 q(6) = rk(10).*oh.*c(1);
 q(7) = rk(11).*c(8).*c(1);
 q(8) = rk(5).*c(5).*oh + rk(12).*c(7);
 q(9) = rk(6).*c(5) + rk(8).*c(10).*c(2);
 q(10) = rk(4).*c(4).*oh + rk(6).*c(5) + rk(9).*c(8).*c(2);

 % at nighttime decrease deposition under stagnant conditions
 %  2000 hours to 0800 hours

 if(etime_fv < 8.0 || etime_fv > 20.0);
  dadj = 0.1;
 else;
  dadj = 1.0;
 end;

 % determine loss rate for each species

 l(1) = rk(1) + rk(10).*oh + rk(11).*c(8) + vent + dadj.*dep(1);
 l(2) = rk(3).*c(3) + rk(7).*c(9) + rk(8).*c(10) + rk(9).*c(8) + vent + dadj.*dep(2);
 l(3) = rk(3).*c(2) + vent + dadj.*dep(3);
 l(4) = rk(4).*oh + vent + dadj.*dep(4);
 l(5) = rk(5).*oh + rk(6) + vent + dadj.*dep(5);
 l(6) = vent + dadj.*dep(6);
 l(7) = rk(12) + vent + dadj.*dep(7);
 l(8) = rk(9).*c(2) + rk(11).*c(1) + vent + dadj.*dep(8);
 l(9) = rk(7).*c(2) + vent + dadj.*dep(9);
 l(10) = rk(8).*c(2) + vent + dadj.*dep(10);

 % determine net loss rates

 for i = 1: 10;

  g(i) = q(i) - l(i).*c(i);

 end;

 return;
end %subroutine diffr