% MS_VRH - n phase Voigt-Reuss-Hill average of elasticity.
%
% // Part of MSAT - The Matlab Seismic Anisotropy Toolkit //
%
% Given n elasticity matricies, densities and phase volume fractions, 
%     calculate the avarage of the Voigt upper and Reuss lower bound of the 
%     composite material. 
%
%  %  [Cave,rhave]=MS_VRH(VF, C, rh, ...)
%
% Usage: 
%     [Cave,rhave]=MS_VRH(VF, C, rh)                    
%          The volume fractions and densities (VF and rh) can be provided as
%          vectors of common length n and the elasticities as a (6,6,n)
%          matrix. 
%
%     [Cave,rhave]=MS_VRH(VF, C1, rh1, C2, rh2, C3, rh3, ...)                    
%          Alternitivly, VF can be a vector of volume fractions and the 
%          elasticity matricies and densities given as sequence of arguments
%          as (6,6) matricies and scalars. There should be a C,rho argument 
%          pair for each entry in the VF vector.
%
%     [Cave,rhave,Cvoigt,Creuss]=MS_VRH(VF, C, rh, ...)                    
%          For both input argument forms, it is possible to access the upper
%          and lower bound on the composite elasticity by adding output 
%          arguments.
%
% Notes:
%     The sum of the elements of VF is normalised to 1 before averaging.
%
% See also: MS_POLYAVERAGE

% Copyright (c) 2011, James Wookey and Andrew Walker
% Copyright (c) 2006, James Wookey
% All rights reserved.
% 
% Redistribution and use in source and binary forms, 
% with or without modification, are permitted provided 
% that the following conditions are met:
% 
%    * Redistributions of source code must retain the 
%      above copyright notice, this list of conditions 
%      and the following disclaimer.
%    * Redistributions in binary form must reproduce 
%      the above copyright notice, this list of conditions 
%      and the following disclaimer in the documentation 
%      and/or other materials provided with the distribution.
%    * Neither the name of the University of Bristol nor the names 
%      of its contributors may be used to endorse or promote 
%      products derived from this software without specific 
%      prior written permission.
% 
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS 
% AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED 
% WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED 
% WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A 
% PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL 
% THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY 
% DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR 
% CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, 
% PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF 
% USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER 
% CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN 
% CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE 
% OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS 
% SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

function [Cave,rhave,voigt_ave,reuss_ave]=MS_VRH(VF,C,r,varargin)


   optargs = size(varargin,2);
   if (optargs == 0)
       % VF, C and r should be arrays of size (n), (6,6,n) and (n),
       % respectivly. Check this and proceed.
       np = length(VF) ;% Number of phases
       sc = size(C);
       assert(length(sc) == 3, 'MS:VRH:args', ... 
           'For three argument form, C must be size (6,6,np)') 
       assert(sc(1) == 6, 'MS:VRH:args', ... 
           'For three argument form, C must be size (*6*,6,np)')
       assert(sc(2) == 6, 'MS:VRH:args', ... 
           'For three argument form, C must be size (6,*6*,np)')
       assert(sc(3) == np, 'MS:VRH:args', ... 
           'For three argument form, C must be size (6,6,*np*)')
       assert(length(r) == np, 'MS:VRH:args', ... 
           'For three argument form, r must be length np')
       Cs = C;
       rs = r;
   elseif (optargs > 0)
       % The old interface - VF should be a size (n) and we should see 2n-2
       % varargs in pairs of size (6,6) for Cs and scalars for the
       % densitys. Check this here then put the values into optargs == 0
       % format.
       np = length(VF); % Number of phases
       assert((np*2)-2  == (length(varargin)), 'MS:VRH:args', ...
           'For many argument form, there must be as many C, rho argument pairs as elements in VF')
       Cs = zeros(6,6,np);
       rs = zeros(1,np);
       Cs(:,:,1) = C(:,:);
       rs(1) = r;
       for i = 2:np
           Cs(:,:,i) = varargin{(i-2).*2+1};
           rs(i) = varargin{(i-2).*2+2} ;
       end
   end
  
%--normalise VF
   VF = VF ./ sum(VF) ;

   voigt_ave = zeros(6,6) ;
   reuss_ave = zeros(6,6) ;
   rhave = 0 ;

   for i=1:np
      voigt_ave = voigt_ave + Cs(:,:,i) .* VF(i) ;      % stiffness average
      reuss_ave = reuss_ave + inv(Cs(:,:,i)) .* VF(i) ; % compliance average
      rhave = rhave + rs(i) .* VF(i) ;
   end
   reuss_ave = inv(reuss_ave);
   Cave = (voigt_ave + reuss_ave) ./ 2  ;

return
