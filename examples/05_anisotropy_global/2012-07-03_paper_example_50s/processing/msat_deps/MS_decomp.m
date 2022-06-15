% MS_DECOMP - Browaeys and Chevrot decomposition of the elasticity matrix.
%
% // Part of MSAT - The Matlab Seismic Anisotropy Toolkit //
%
%  Apply a decomposition of the elasticity tensor C, after: 
%     Browaeys and Chevrot (GJI, v159, 667-678, 2004)
%  
%  [Ciso] = MS_decomp(C)
%     Isotropic projection of the elastic tensor.
%
%  [Ciso,Chex] = MS_decomp(C)
%     Isotropic, and hexagonal parts of the elastic tensor
%
%  [Ciso,Chex,Ctet,Cort,Cmon,Ctri] = MS_decomp(C)
%     All parts of the elastic tensor
%     
%
% Notes:
%     Output matricies are the partial components of the input elasticity
%     matrix which maximise the norm of the high symmetry parts. This 
%     function assumes that C is in its optimal orientation for 
%     decomposition (use MS_axes to make this so).
%
% References:
%     Browaeys, J. T. and S. Chevrot (2004) Decomposition of the elastic
%         tensor and geophysical applications. Geophysical Journal 
%         international v159, 667-678.

%
% See also: MS_NORMS, MS_AXES

% Copyright (c) 2011, James Wookey and Andrew Walker
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

function [ varargout ] = MS_decomp( C )

   i=nargout ;
   
   if (nargout==6), i=5;, end
   
   for i=1:5
	   [X]=C2X(C) ;
	   M=Projector(i) ;
	   XH = M*X ;
	   CH = X2C(XH) ;
	   varargout{i} = CH ;
      C=C-CH ;
   end

   if (nargout==6), varargout{6} = C;, end
   

return

function M=Projector(order)
switch order 
case 1 % isotropic
   M = zeros(21,21) ;
   M(1:9,1:9) = [ ...
      3/15 3/15 3/15 sqrt(2)/15 sqrt(2)/15 sqrt(2)/15 2/15 2/15 2/15 ; ... 
      3/15 3/15 3/15 sqrt(2)/15 sqrt(2)/15 sqrt(2)/15 2/15 2/15 2/15 ; ... 
      3/15 3/15 3/15 sqrt(2)/15 sqrt(2)/15 sqrt(2)/15 2/15 2/15 2/15 ; ... 
      sqrt(2)/15 sqrt(2)/15 sqrt(2)/15 4/15 4/15 4/15 -sqrt(2)/15 -sqrt(2)/15 -sqrt(2)/15 ; ... 
      sqrt(2)/15 sqrt(2)/15 sqrt(2)/15 4/15 4/15 4/15 -sqrt(2)/15 -sqrt(2)/15 -sqrt(2)/15 ; ... 
      sqrt(2)/15 sqrt(2)/15 sqrt(2)/15 4/15 4/15 4/15 -sqrt(2)/15 -sqrt(2)/15 -sqrt(2)/15 ; ... 
      2/15 2/15 2/15 -sqrt(2)/15 -sqrt(2)/15 -sqrt(2)/15 1/5 1/5 1/5 ; ... 
      2/15 2/15 2/15 -sqrt(2)/15 -sqrt(2)/15 -sqrt(2)/15 1/5 1/5 1/5 ; ... 
      2/15 2/15 2/15 -sqrt(2)/15 -sqrt(2)/15 -sqrt(2)/15 1/5 1/5 1/5 ; ... 
   ] ;
case 2 % hexagonal
   M = zeros(21,21) ;
   M(1:9,1:9) = [...
      3/8 3/8 0 0 0 1./(4.*sqrt(2)) 0 0 1/4                            ; ...
      3/8 3/8 0 0 0 1./(4.*sqrt(2)) 0 0 1/4                            ; ...
      0 0 1 0 0 0 0 0 0                                                ; ...
      0 0 0 1/2 1/2 0 0 0 0                                            ; ...
      0 0 0 1/2 1/2 0 0 0 0                                            ; ...
      1./(4.*sqrt(2)) 1./(4.*sqrt(2)) 0 0 0 3/4 0 0 -1./(2.*sqrt(2))   ; ...
      0 0 0 0 0 0 1/2 1/2 0                                            ; ...
      0 0 0 0 0 0 1/2 1/2 0                                            ; ...
      1/4 1/4 0 0 0 -1./(2.*sqrt(2)) 0 0 1/2                           ; ...
   ] ;
case 3 % tetragonal
   M = zeros(21,21) ;
   M(1:9,1:9) = [...
      1/2 1/2 0 0 0 0 0 0 0  ; ...
      1/2 1/2 0 0 0 0 0 0 0  ; ...
      0 0 1 0 0 0 0 0 0      ; ...
      0 0 0 1/2 1/2 0 0 0 0  ; ...
      0 0 0 1/2 1/2 0 0 0 0  ; ...
      0 0 0 0 0 1 0 0 0      ; ...
      0 0 0 0 0 0 1/2 1/2 0  ; ...
      0 0 0 0 0 0 1/2 1/2 0  ; ...
      0 0 0 0 0 0 0 0 1      ; ...
   ] ;
case 4 % orthorhombic
   M = zeros(21,21) ;
   for jj=1:9
      M(jj,jj)=1;
   end
case 5 % monoclinic
   M = eye(21,21) ;
   for jj=[10, 11, 13, 14, 16, 17, 19, 20]
      M(jj,jj)=0;
   end
otherwise
   error('Unsupported symmetry class')
end
return


function [X]=C2X(C)
%  after Browaeys and Chevrot (GJI, 2004)
	X = zeros(21,1) ;
	
	X(1)  = C(1,1) ;
	X(2)  = C(2,2) ;
	X(3)  = C(3,3) ;
	X(4)  = sqrt(2).*C(2,3) ;
	X(5)  = sqrt(2).*C(1,3) ;
	X(6)  = sqrt(2).*C(1,2) ;
	X(7)  = 2.*C(4,4) ;
	X(8)  = 2.*C(5,5) ;
	X(9)  = 2.*C(6,6) ;
	X(10) = 2.*C(1,4) ;
	X(11) = 2.*C(2,5) ;
	X(12) = 2.*C(3,6) ;
	X(13) = 2.*C(3,4) ;
	X(14) = 2.*C(1,5) ;
	X(15) = 2.*C(2,6) ;
	X(16) = 2.*C(2,4) ;
	X(17) = 2.*C(3,5) ;
	X(18) = 2.*C(1,6) ;
	X(19) = 2.*sqrt(2).*C(5,6) ;
	X(20) = 2.*sqrt(2).*C(4,6) ;
	X(21) = 2.*sqrt(2).*C(4,5) ;
	
return

function [C]=X2C(X)
%  after Browaeys and Chevrot (GJI, 2004)
	C = zeros(6,6) ;
	
	C(1,1) = X(1);
	C(2,2) = X(2);
	C(3,3) = X(3);
	C(2,3) = 1./(sqrt(2)).*X(4);
	C(1,3) = 1./(sqrt(2)).*X(5);
	C(1,2) = 1./(sqrt(2)).*X(6);
	C(4,4) = 1./(2).*X(7);
	C(5,5) = 1./(2).*X(8);
	C(6,6) = 1./(2).*X(9);
	C(1,4) = 1./(2).*X(10);
	C(2,5) = 1./(2).*X(11);
	C(3,6) = 1./(2).*X(12);
	C(3,4) = 1./(2).*X(13);
	C(1,5) = 1./(2).*X(14);
	C(2,6) = 1./(2).*X(15);
	C(2,4) = 1./(2).*X(16);
	C(3,5) = 1./(2).*X(17);
	C(1,6) = 1./(2).*X(18);
	C(5,6) = 1./(2.*sqrt(2)).*X(19);
	C(4,6) = 1./(2.*sqrt(2)).*X(20);
	C(4,5) = 1./(2.*sqrt(2)).*X(21);
	
	for i=1:6
		for j=i:6
			C(j,i) = C(i,j) ;
		end
	end		
	
return
