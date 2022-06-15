% MS_PHASEVELS - Wave velocities in anisotropic media.
%
% // Part of MSAT - The Matlab Seismic Anisotropy Toolkit //
%
% Calculate the phase velocity details for an elsticity matrix. 
%
%  [ pol, avs, vs1, vs2, vp, ...] = MS_phasevels( C, rh, inc, azi )
%
% Usage: 
%     [ pol, avs, vs1, vs2, vp ] = MS_phasevels( C, rh, inc, azi )                    
%         Calculate phase velocities from elasticity matrix C (in GPa) and
%         density rh (in kg/m^3) for a propogation direction defined by
%         an inclination and azimuth (both in degrees, see below). Output 
%         details are given below.
%
%     [ pol, avs, vs1, vs2, vp, SF, SS ] = MS_phasevels( C, rh, inc, azi )                    
%         Additionally output fast and slow S-wave polarisation in vector
%         form.
%
% Notes:
%     Azi is defined as the angle in degrees from the +ve 1-axis in x1-x2 
%     plane with +ve being clockwise when looking at origin from the
%     3-axis. Inc is defined as the angle in degrees from the x1-x2 plane
%     towards x3 with zero being in the x1-x2 plane. Inc and azi may be
%     scalars, or vectors of the same size. Outputs are:
%
%       'pol' = angle in plane normal to raypath of FSW                           
%              (deg, zero is x3 direction, +ve c'wise looking along             
%              raypath at origin)  
%       'avs' = shear-wave anisotropy
%       'vs1' = fast shear-wave velocity (km/s)
%       'vs2' = slow shear-wave velocity (km/s)
%       'vp'  = P-wave velocity (km/s)
%
%     and all are vectors of length equal to the input inc and azi vectors.
%     In the case of no S-wave splitting (vs1 and vs2 are equal to within
%     eps^1/2) pol is set to NaN. Optional outputs SF and SS are arrays of
%     size (length(inc),3), with each row corresponding to a polarisation
%     vector. This implementation is based on EMATRIX6 by D. Mainprice. 
%     Re-coded in MATLAB by James Wookey but now avoids transforming
%     from the 6x6 elasticity matrix into the 3x3x3x3 tensor form using the
%     method outlined in Winterstein (1990).
% 
% Reference: Mainprice D. (1990). An efficient
%            FORTRAN program to calculate seismic anisotropy from
%            the lattice preferred orientation of minerals.
%            Computers & Gesosciences, vol16, pp385-393.
%            Winterstein D. F. (1990). Velocity anisotropy terminology 
%            for geophysicists. Geophysics, vol 55, pp1070-1088. 

% Copyright (c) 2011-2012, James Wookey and Andrew Walker
% Copyright (c) 2007-2011, James Wookey
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

function [pol,avs,vs1,vs2,vp, S1P, S2P] = MS_phasevels(C,rh,inc,azi)

      if (length(inc)~=length(azi))
			error('MS:ListsMustMatch', ...
                'AZI and INC must be scalars or vectors of the same dimension');
      end
      
%  ** convert inc, azi to column vectors if necessary      
      inc = reshape(inc,length(inc),1) ;
      azi = reshape(azi,length(azi),1) ;

      isotol = sqrt(eps); % Mbars

      % Check that C is valid (if check not suppressed)
      MS_checkC(C);
      
      %  ** convert GPa to MB file units (Mbars), density to g/cc

      C(:,:) = C(:,:) * 0.01 ;
      rh = rh ./ 1e3 ;
      
		avs = zeros(size(azi)) ;
		vp = zeros(size(azi)) ;
		vs1 = zeros(size(azi)) ;
		vs2 = zeros(size(azi)) ;
		pol = zeros(size(azi)) ;
		S1 = zeros(length(azi),3) ;
		S1P = zeros(length(azi),3) ;
      S2P = zeros(length(azi),3) ;
        
%   ** Handle isotropic case quickly
     if isIsotropic(C, isotol)
         vp(:) = sqrt(( ((1.0/3.0)*(C(1,1)+2*C(1,2)))+ ...
                        ((4.0/3.0)*C(4,4)) )/rh)*10.0;
         vs1(:) = sqrt(C(4,4)/rh)*10.0; % Factor of 10 converts from
         vs2 = vs1;                     % Mbar to Pa.
         avs(:) = 0.0;
         pol(:) = NaN; % Both waves have same velocity... meaningless.
         S1P(:) = NaN;
         return
     end

%	** start looping
	for ipair = 1:length(inc)
		cazi = azi(ipair) ;
		cinc = inc(ipair) ;

%  ** create the cartesian vector
		XI = cart2(cinc,cazi) ;

%  ** compute phase velocities		
		[V,EIGVEC]=velo(XI,rh,C) ;
		
%  ** pull out the eigenvectors
      P  = EIGVEC(:,1) ;
      S1 = EIGVEC(:,2) ;

  		if ~isreal(S1)
            error_str = ['The S1 polarisation vector is not real!\n\n'...
                sprintf('inc = %f, azi = %f\n\n',cinc,cazi) ...
                sprintf('C = %f %f %f %f %f %f\n',C(1:6,1)) ...
                sprintf('    %f %f %f %f %f %f\n',C(1:6,2)) ...
                sprintf('    %f %f %f %f %f %f\n',C(1:6,3)) ...
                sprintf('    %f %f %f %f %f %f\n',C(1:6,4)) ...
                sprintf('    %f %f %f %f %f %f\n',C(1:6,5)) ...
                sprintf('    %f %f %f %f %f %f\n\n',C(1:6,6)) ...
                sprintf('S1 = %f %f %f\n',S1)];
  			error('MS:PHASEVELS:vectornotreal', error_str) ;
  		end
      S2 = EIGVEC(:,3) ;

%  ** calculate projection onto propagation plane      
      S1N = V_cross(XI,S1) ;
      S1P(ipair,:) = V_cross(XI,S1N);
      S2N = V_cross(XI,S2) ;
      S2P(ipair,:) = V_cross(XI,S2N);

%  ** rotate into y-z plane to calculate angles
%     (use functions optimised for the two needed 
%      rotations, see below).
      [S1PR] = V_rot_gam(S1P(ipair,:),cazi) ;
	   [S1PRR] = V_rot_bet(S1PR,cinc) ;


      
	   ph = atan2(S1PRR(2),S1PRR(3)) .* 180/pi ;

%  ** transform angle to between -90 and 90
      if (ph < -90.), ph = ph + 180.;end
      if (ph >  90.), ph = ph - 180.;end

%	** calculate some useful values
      dVS =  (V(2)-V(3)) ;
      VSmean = (V(2)+V(3))/2.0 ;

      avs(ipair) = 100.0*(dVS/VSmean) ;
      vp(ipair) =  V(1) ;
      vs1(ipair) = V(2) ;
      vs2(ipair) = V(3) ;
		
      pol(ipair) = ph ;
	end % ipair = 1:length(inc_in)

   % If any directions have zero avs (within machine accuracy)
   % set pol to NaN - array wise:
   isiso = real(avs > sqrt(eps)) ; % list of 1.0 and 0.0.
   pol = pol .* (isiso./isiso) ; % times by 1.0 or NaN. 

   S1P(:,1) = S1P(:,1) .* (isiso./isiso);
   S1P(:,2) = S1P(:,2) .* (isiso./isiso);
   S1P(:,3) = S1P(:,3) .* (isiso./isiso);
    
return
%=======================================================================================  

function [c] = V_cross(a, b)
    % This is ~10 times quicker than the Matlab cross() function
    % for me (AMW). We assume the arguments are both 3-vectors and
    % avoid the checks and reshaping needed for the more general case.
    % (According to the prfiler, this moves cross from the most expensive
    % child function costing ~50% of the time to the third most expensive 
    % child costing ~10% of the time).

    c = [a(2).*b(3)-a(3).*b(2)
         a(3).*b(1)-a(1).*b(3)
         a(1).*b(2)-a(2).*b(1)];
    
return


%=======================================================================================
% Rather useing a general case rotation about three angles we use two 
% special case rotations around the c and b axes. This saves time as
% we don't need to convert 0 degrees to radians six times per call to 
% phasevels, or do two double matrix multiplications. This saves ~60%
% of the time spent in vector rotation (and, with the cross thing above)
% reduces the 10000 evaulations needed for MS_anisotropy lma from ~19
% secs to ~10 secs.

function [VR] = V_rot_gam(V,gam)

%  Make rotation matrix
g = gam * pi/180. ;

RR = [ cos(g) sin(g) 0 ; -sin(g) cos(g) 0 ; 0 0 1 ] ;
VR = V * RR ;
 
return

function [VR] = V_rot_bet(V,bet)

%  Make rotation matrix
b = bet * pi/180. ;

RR = [ cos(b) 0 -sin(b) ; 0 1 0 ; sin(b) 0 cos(b) ] ;

VR = V * RR ;
 
return

%=======================================================================================  

%=======================================================================================  
	function [X] = cart2(inc,azm)
%=======================================================================================  
%c convert from spherical to cartesian co-ordinates
%c north x=100  west y=010 up z=001
%c irev=+1 positive vector x
%c irev=-1 negative vector x
% NB: pre-converting azm and inc to radians and using
%     cos and sin directly (instead to cosd and sind) 
%     is ~10x faster making the function ~4x faster.
    azmr = azm.*(pi/180.0);
    incr = inc.*(pi/180.0);
    caz=cos(azmr)  ;
    saz=sin(azmr)  ;
    cinc=cos(incr) ;
    sinc=sin(incr) ;
    X=[caz*cinc -saz*cinc sinc] ;
%c normalise to direction cosines
    r=sqrt(X(1)*X(1)+X(2)*X(2)+X(3)*X(3)) ;
   
	X = X./r ;
   return
%=======================================================================================  

%=======================================================================================  
	function [V,EIGVEC]=velo(X,rh,C)
%=======================================================================================  
% PHASE-VELOCITY SURFACES IN AN ANISOTROPIC MEDIUM
% revised April 1991
%     X(3) - DIRECTION OF INTEREST
%     RHO - DENSITY
%     V - PHASE VELOCITIES (1,2,3= P,S,SS)
%     EIGVEC(3,3) - eigenvectors stored by columns
%
% Translated to MATLAB by James Wookey         
% Revised, Andrew Walker 2012. 

% Form the 3x3 Christoffel tensor without converting the elasticity 
% from 6x6 to 3x3x3x3 form using the formula from page 1076 of 
% Winterstein 1990. This is > twice as fast as the quickest way I
% have found going via the full tensor form.
        gamma = [X(1) 0.0  0.0  0.0  X(3) X(2) ; ...
                 0.0  X(2) 0.0  X(3) 0.0  X(1) ; ...
                 0.0  0.0  X(3) X(2) X(1) 0.0 ];
        T = gamma * C * gamma';
         
% determine the eigenvalues of symmetric tij
        [EIVEC EIVAL] = eig(T) ;

% calculate velocities and sort
% note that we could get a significant speedup if
% we could avoid the sort - LAPACK usually does sort
% in order of incresing eigenvectors but I've found 
% cases where this does not happen (and we swap P and 
% S wave velocities). Note that MATLAB does not "guarantee
% that the eignevalues are not returned in sorted order"
% http://www.mathworks.com/matlabcentral/newsreader/view_original/394371 
% so we have to sort here...
		V_RAW = (sqrt([EIVAL(1,1) EIVAL(2,2) EIVAL(3,3)]./rh))*10. ;
		[V IND] = sort(V_RAW,2,'descend') ;
		EIGVEC = EIVEC ; % for dimensioning
		for i=1:3
			EIGVEC(:,i) = EIVEC(:,IND(i)) ;
		end

      return
%=======================================================================================  

function [ l ] = isIsotropic( C, tol )
    
    % Are we isotropic - assume matrix is symmetrical at this point.
    l = (abs(C(1,1)-C(2,2)) < tol) & (abs(C(1,1)-C(3,3)) < tol) & ...
        (abs(C(1,2)-C(1,3)) < tol) & (abs(C(1,2)-C(2,3)) < tol) & ...
        (abs(C(4,4)-C(5,5)) < tol) & (abs(C(4,4)-C(6,6)) < tol) & ...
        (abs(C(1,4)) < tol) & (abs(C(1,5)) < tol) & (abs(C(1,6)) < tol) & ...
        (abs(C(2,4)) < tol) & (abs(C(2,5)) < tol) & (abs(C(2,6)) < tol) & ...
        (abs(C(3,4)) < tol) & (abs(C(3,5)) < tol) & (abs(C(3,6)) < tol) & ...
        (abs(C(4,5)) < tol) & (abs(C(4,6)) < tol) & (abs(C(5,6)) < tol) & ...
        (((C(1,1)-C(1,2))/2.0)-C(4,4) < tol);
        
return
