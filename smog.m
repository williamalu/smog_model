function smog(varargin)

%  program smog
%
%  simulation of photochemical reactions of vocs, nox and hcho
%  and the formation of ozone
%
%  matt fraser, scott hersey
%  april 2002, feb 2016
%
%  coded in fortran (f77)
%
%
% photolysis rates
% deposition rates
% emission rates
% ventilation rate
% chemical reaction rate constants
% contant value for o2,m

%% Definitions

clear all; %clear functions;

 global unit2fid;  if ~isempty(unit2fid), unit2fid=[]; end
 persistent c dep e etime_fv firstCall i iday j1 j2 m n o o2 oh p_no2 p_rcho rk tin tinc tout treport trinc vent ; if isempty(firstCall),firstCall=1;end; 

 format_96=[ '\n ' ,'Emissions:',repmat(['%1x','%9.3e'] ,1,10)];
 format_97=[ '\n ' ,'end day','%2d', '\n ' ];
 format_98=['%1x','%4.1f','%5x',repmat(['%1x','%9.3e'] ,1,14)];
 format_99=['time(h)    NO2       NO        O3        RH        ''RCHO      HNO3      PAN       RCOO2     HO2       ''RO2       OOHj(NO2)    j(RCHO)'];

 if isempty(o2), o2=0; end;
 if isempty(m), m=0; end;
 % concentration of accumulation species
 % day of simulation
 if isempty(iday), iday=0; end;
 % time markers
 if isempty(etime_fv), etime_fv=0; end;
 if isempty(tinc), tinc=0; end;
 % time markers
 if isempty(treport), treport=0; end;
 if isempty(trinc), trinc=0; end;
 % interpolation indexes
 if isempty(j1), j1=0; end;
 if isempty(j2), j2=0; end;
 % pssa values for o and oh radicals
 if isempty(o), o=0; end;
 if isempty(oh), oh=0; end;
 % number of accumulation species
 % index
 if isempty(i), i=0; end;
 % times for integration
 if isempty(tin), tin=0; end;
 if isempty(tout), tout=0; end;




 %  ARRAY ASSIGNMENTS FOR SPECIES
 %  1 = no2
 %  2 = no
 %  3 = o3
 %  4 = rh
 %  5 = rcho
 %  6 = hno3
 %  7 = pan
 %  8 = rcoo2
 %  9 = ho2
 % 10 = ro2
 %  o and oh determined by pssa

 %  photolysis rates per minute from tabulated values for NO2 and RCHO. One rate for each
 %  hour of day. base case assumes first hour with light is 5 am, first hour of dark 8 pm, peak
 %  photolysis at 12 noon
 if firstCall,   p_no2=[0.0,0.0,0.0,0.0,0.010,0.136,0.292,0.392,0.464,0.503,0.523,0.529,0.519,0.493,0.450,0.379,0.260,0.104,0.0,0.0,0.0,0.0,0.0,0.0];  end;

 if firstCall,   p_rcho=[0.0,0.0,0.0,0.0,0.0,0.21e-3,0.67e-3,1.20e-3,1.64e-3,1.96e-3,2.14e-3,2.20e-3,2.13e-3,1.92e-3,1.59e-3,1.12e-3,0.60e-3,0.16e-3,0.0,0.0,0.0,0.0,0.0,0.0];  end;

 % reaction rate constants (units of ppm min). rate constant given for each
 % reaction of importance. rates 1 and 6 initially set at zero, as they are photolysis-driven. a later if
 % function re-sets these to appropriate values when the sun rises.
 if firstCall,   rk=[0.0,2.183e-5,26.59,3.775e3,2.341e4,0.0,1.214e4,1.127e4,3.8e3,1.613e4,2.07e3,2.143e-2];  end;

 % concentrations of species that are constant in ppm
 o2 = 2.1e5;
 m = 1.0e6;

 % initial concentrations of species in ppm (arbitrary)
 if firstCall,   c=[0.010,0.010,0.010,0.050,0.010,0.001,0.001,0.0,0.0,0.0];  end;

 % deposition rates dependant on deposition velocity
 % and thus mixing height and varies by how sticky compound is. you may
 % change these if you need to add a loss term to account for large amounts
 % of aerosol surface area in your city.
 %
 % no2 = 0.42 cm/s = 1.1% per hour
 % o3 = 2.5 cm/s = 6.4% per hour
 % rcho = 0.42 cm/s = 1.1% per hour
 % hno3 = 2.5 cm/s = 6.4% per hour
 % pan = 1.7 cm/s = 4.4% per hour

 % data for deposition given in fraction of concentration per minute

 if firstCall,   dep=[0.18e-3,0.0,1.1e-3,0.0,0.18e-3,1.1e-3,0.73e-3,0.0,0.0,0.0];  end;

 % emission rates given in ppm/min
 % base case emissions
 % you will assuredly change these for your city!

 if firstCall,   e=[6.2e-6,55.8e-6,0.0,125.0e-6,3.5e-6,0.0,0.0,0.0,0.0,0.0];  end;

 % loss rate through ventilation in fraction per minute
 % this is a windspeed term. higher windspeeds result in greater
 % ventilation rates. if you have a city with an inversion, you may want to
 % decrease this. in a city with high windspeeds, consider increasing it.

 if firstCall,   vent=[0.0007];  end;

 % number of accumulation species

 if firstCall,   n=[10];  end;
 firstCall=0;

 % time factors in hours
 iday = 1;
 etime_fv = 0.0;
 tinc = 0.1;
 treport = 0.3;
 trinc = 0.3;

 %% begin building output 
 %
 % output file = smog.out
 %
 thismlfid=fopen(strtrim('smog.out'),'w+');
 unit2fid=[unit2fid;1001,thismlfid];
 [writeErrFlag]=writeFmt(1001,[format_99]);
 while (1);


  % set photolysis rates by interpolation for reactions that are
  % photolysis-driven

  j1 = fix(fix(etime_fv));
  j2 = fix(j1 + 1);
  if(j1 == 0);
   rk(1) = 0.0;
   rk(6) = 0.0;
  else;
   rk(1) =(etime_fv-1.0.*j1).*p_no2(j2) +(1.0.*j2-etime_fv).*p_no2(j1);
   rk(6) =(etime_fv-1.0.*j1).*p_rcho(j2) +(1.0.*j2-etime_fv).*p_rcho(j1);
  end;

  % calculate a step for solution
  % convert from hours (etime) to min (tin)
  tin = 60.0.*etime_fv;
  tout = 60.0.*(etime_fv+tinc);
  [n,c,rk,tin,tout,e,dep,vent,etime_fv,o,oh,m,o2]=hybrid(n,c,rk,tin,tout,e,dep,vent,etime_fv,o,oh,m,o2);

  etime_fv = tout./60.0;

  % update time and write to output if needed
  if(etime_fv >= treport);
   treport = treport + trinc;
   [writeErrFlag]=writeFmt(1001,[format_98],'etime_fv',{'c(i)','i','1','1','10'},'o','oh','rk(1)','rk(6)');
  end;

  if(etime_fv >= 24.0);
   [writeErrFlag]=writeFmt(1,[format_97],'iday');
   [writeErrFlag]=writeFmt(1001,[format_97],'iday');
   iday = fix(iday + 1);
   etime_fv = etime_fv - 24.0;
   treport = treport - 24.0;
  end;

  if(iday >= 8);
   [writeErrFlag]=writeFmt(1001,[format_96],{'e(i)','i','1','1','10'});
   warning(['stop encountered in original fortran code  ',char(10),';']);
   return

  end;
 end;

end %program smog