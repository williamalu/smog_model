function [n,c,rk,tin,tout,e,dep,vent,etime_fv,o,oh,m,o2]=hybrid(n,c,rk,tin,tout,e,dep,vent,etime_fv,o,oh,m,o2,varargin);
%%
% subroutine for integration of stiff set of coupled differential equations
% using a hybrid integration technique
%
% stiffness: an ordinary differential equation problem is stiff if the 
% solution being sought is varying slowly, but there are nearby solutions 
% that vary rapidly, so the numerical method must take small steps to 
% obtain satisfactory results.
%
% the general for of the equations is assumed to be:
%
%   dc/dt = g(c) = q(c) - l(c)*c
%
% where q(c) is formation and emission and l(c) is first order loss
%
% must call a subroutine diffr(c,rk,q,l,g) to provide the g, q, and l
% functions from the concentrations (c) and the rate constants (rk)
% time step over that interval (tin to tout) should be small so that the
% rate constants are constant (i.e. photolysis does not change).
%
% in calling hybrid, n = number of equations
%                    c = vector of concentrations
%                        at time tin on input
%                        at time tout on output
%                    rk = rate constants
%                    tin = time at beginning of integration
%                    tout = time at end of integration

%%
% declare time values

 persistent c2 c3 convrg dt dtmin dtovr2 e1 e2 e3 e4 e5 e6 fascon firstCall g g2 ifast itmax j k l l2 q q2 remg stiff sumtau tau2 tnow xxx y2 y3 ; if isempty(firstCall),firstCall=1;end; 

 format_92=[ '\n ' ,' stepsize below minimum allowed in hybrid solution'];

 if isempty(remg), remg([1:1])=true; end;
 if isempty(tnow), tnow=0; end;
 if isempty(dtovr2), dtovr2=0; end;
 if isempty(xxx), xxx=0; end;
 if isempty(dt), dt=0; end;

 % declare work variables
 if isempty(y2), y2=0; end;
 if isempty(y3), y3=0; end;
 if isempty(tau2), tau2=0; end;
 if isempty(sumtau), sumtau=0; end;

 % declare integration control parameters

 % declare concentration variables
 if isempty(c2), c2=zeros(1,10); end;
 if isempty(c3), c3=zeros(1,10); end;

 % rate constants

 % formation rate (rxn + emission) for each species
 if isempty(q), q=zeros(1,10); end;
 if isempty(q2), q2=zeros(1,10); end;

 % sink (first order rate constants) for each species
 if isempty(l), l=zeros(1,10); end;
 if isempty(l2), l2=zeros(1,10); end;

 % net formation (g = q - l*c)
 if isempty(g), g=zeros(1,10); end;
 if isempty(g2), g2=zeros(1,10); end;

 % common blocks

 % flags for stiff species, fast convergence and convergence
 if isempty(stiff), stiff=false(1,10); end;
 if isempty(fascon), fascon=false; end;
 if isempty(convrg), convrg=false; end;
 if isempty(j), j=0; end;
 if isempty(k), k=0; end;

 % *time constants to set step size*  DO NOT CHANGE
 % e1 = scaling factor to determine initial step size
 % e2 = determines if individual ode is stiff or normal
 %      (if normal -> EULER'S METHOD, if stiff -> ASYMPTOTIC)
 % e3 = convergence criteria
 % e4 = if convergence is not achieved after a few steps, reduce
 %      the step size by this factor
 % e5 = if convergece is reached rapidly, increase the step size by
 %      10% using this factor
 % e6 = for species with low concentrations, automatically assume
 %      convergence.  If concentration below e6 -> assume convergence
 %
 % dtmin = minimum time step allowed before the integration halts with error
 % ifast = number of corrections for fast convergence (if <= ifast then
 %         consider fast convergence
 % itmax = if integration not completed by itmax, then reduce time step and
 %         try again

 if firstCall,   e1=[0.001];  end;
 if firstCall,  e2=[1.0];  end;
 if firstCall,  e3=[0.001];  end;
 if firstCall,  e4=[0.5];  end;
 if firstCall,  e5=[1.1];  end;
 if firstCall,  e6=[1.0e-9];  end;
 if firstCall,  dtmin=[1.0e-5];  end;
 if firstCall,  ifast=[2];  end;
 if firstCall,  itmax=[4];  end;
 firstCall=0;

 % determine q, l, and g by subroutine

 while (1);
  if(remg(1));
   [c,rk,q,l,g,e,dep,vent,etime_fv,o,oh,m,o2]=diffr(c,rk,q,l,g,e,dep,vent,etime_fv,o,oh,m,o2);

   % set the time step with a first estimate
   dt = 1.0;
   xxx = 1.0;
   for j = 1: n;
    if(g(j) ~= 0.0);
     xxx = abs(c(j)./g(j));
    end;
    if(xxx > 0.0);
     dt = min(dt,xxx);
    end;
   end;

   dt = e1.*dt;
   if(dt < dtmin);
    dt = dtmin;
   end;
   tnow = tin;

   % check step size will not put integration past tout
  end;
  remg(1)=true;
  if(tnow+dt > tout);
   dt = tout - tnow;
  end;
  dtovr2 = dt.*0.5;

  % check for stiffness
  for j = 1: n;
   stiff(j) =(dt.*l(j)) > e2;
  end;

  % predict concentration at t+dt
  for j = 1: n;
   if(stiff(j));
    tau2 = 2.0./l(j);
    c2(j) =(c(j).*(tau2-dt)+tau2.*dt.*q(j))./(tau2+dt);
   else;
    c2(j) = c(j) + g(j).*dt;
   end;

   % reset any negative concentrations
   if(c2(j) < 0.0);
    c2(j) = 0.0;
   end;
  end;

  % iteratively correct the prediction for t+dt
  k = 1;
  fascon = true;
  while (1);

   % compute q, l and g using lastest concentration at t+dt
   [c2,rk,q2,l2,g2,e,dep,vent,etime_fv,o,oh,m,o2]=diffr(c2,rk,q2,l2,g2,e,dep,vent,etime_fv,o,oh,m,o2);
   % corrected predictions for c(t+dt)
   for j = 1: n;
    if(stiff(j));
     sumtau = 1.0./l(j) + 1.0./l2(j);
     c3(j) =(dtovr2.*sumtau.*(q(j)+q2(j))+c(j).*(sumtau-dt))./(sumtau+dt);
    else;
     c3(j) = c(j) +(g(j)+g2(j)).*dtovr2;
    end;

    % reset any negative concentrations
    if(c3(j) < 0.0);
     c3(j) = 0.0;
    end;
   end;

   % test for convergence
   convrg = true;
   for j = 1: n;
    y2 = c2(j);
    y3 = c3(j);
    if(y3 > e6 && abs(y3-y2) > e3.*min(y2,y3));
     convrg = false;
     tempBreak=1;break;
    end;
   end;

   if( ~ convrg);

    % if it has converged, assign concentration values for next step
    % if too many iterations, decrease step size

    if(k < itmax);
     for j = 1: n;
      c2(j) = c3(j);
     end;
     fascon = k < ifast;
     k = fix(k + 1);
    else;
     % if iteration did not converge in itmax steps, reduce step size

     dt = e4.*dt;
     if(dt >= dtmin);
      remg(1)=false;
      tempBreak=1;break;
     end;
     [writeErrFlag]=writeFmt(1,[format_92]);
     warning(['stop encountered in original fortran code  ',char(10),';']);
     return
    end;

   else;

    % load c values at time d+dt into c array and assign q, l, g values
    % for next step

    for j = 1: n;
     c(j) = c3(j);
    end;

    tnow = tnow + dt;

    % if iteration converged in less than ifast tries, increase step size

    if(fascon);
     dt = e5.*dt;
    end;

    if(tnow >= tout);
     return;
    else;
     for j = 1: n;
      q(j) = q2(j);
      l(j) = l2(j);
      g(j) = g2(j);
     end;
     remg(1)=false;
     tempBreak=1;break;
    end;
   end;
  end;
  if(~(remg(1)));
   continue;
  end;
  tempBreak=1;break;
 end;

end %subroutine hybrid