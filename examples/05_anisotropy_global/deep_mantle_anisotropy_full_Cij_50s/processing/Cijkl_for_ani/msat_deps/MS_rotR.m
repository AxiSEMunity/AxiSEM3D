% MS_ROTR - Script to rotate a set of elastic constants by a rotation matrix
%
% // Part of MSAT - The Matlab Seismic Anisotropy Toolkit //
%
% Rotate an elasticity matrices using a rotation matrix
%
%  % [ CR ] = MS_rotR( C, R )
%
% Usage: 
%     For a three-by-three rotation matrix, R, rotate the elasticity matrix,
%     C to give a rotated matrix CR. 
%
% Notes:
%    The rotation is performed without transforming the elasticity 
%    matrix to the full tensor form following the method described in
%    Bowers. This eleminates eight nested loops and replaces them with pure
%    matrix-matrix operations, which is (~30 times) faster in Matlab. 
%    Unlike the other MSAT rotation functions, C and R cannot be lists but
%    must be 6x6 and 3x3 matricies, respectivly. Furthermore, the
%    corectness of the input arguments are not checked. Users are encoraged
%    to make use of MS_rot3 or MS_rotEuler for most rotation operations -
%    these make use of this function internally.
%
% References:
%    Bowers 'Applied Mechanics of Solids', Chapter 3
%
% See also: MS_ROT3 MS_ROTEULER

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

function [CR] = MS_rotR(C,R)

% form the K matrix (based on Bowers 'Applied Mechanics of Solids', Chapter 3)
K1 = [ R(1,1).^2 R(1,2).^2 R(1,3).^2 ; ...
       R(2,1).^2 R(2,2).^2 R(2,3).^2 ; ...
       R(3,1).^2 R(3,2).^2 R(3,3).^2 ] ;

K2 = [ R(1,2).*R(1,3) R(1,3).*R(1,1) R(1,1).*R(1,2) ; ...
       R(2,2).*R(2,3) R(2,3).*R(2,1) R(2,1).*R(2,2) ; ...
       R(3,2).*R(3,3) R(3,3).*R(3,1) R(3,1).*R(3,2) ] ;

K3 = [ R(2,1).*R(3,1) R(2,2).*R(3,2) R(2,3).*R(3,3) ; ...
       R(3,1).*R(1,1) R(3,2).*R(1,2) R(3,3).*R(1,3) ; ...
       R(1,1).*R(2,1) R(1,2).*R(2,2) R(1,3).*R(2,3) ] ;

K4 = [ R(2,2).*R(3,3)+R(2,3).*R(3,2) ...
                     R(2,3).*R(3,1)+R(2,1).*R(3,3) ...
                                 R(2,1).*R(3,2)+R(2,2).*R(3,1) ; ...
       R(3,2).*R(1,3)+R(3,3).*R(1,2) ...
                     R(3,3).*R(1,1)+R(3,1).*R(1,3) ...
                                 R(3,1).*R(1,2)+R(3,2).*R(1,1) ; ...      
       R(1,2).*R(2,3)+R(1,3).*R(2,2) ...
                     R(1,3).*R(2,1)+R(1,1).*R(2,3) ...
                                 R(1,1).*R(2,2)+R(1,2).*R(2,1)] ; 

K = [ K1  2.*K2 ; ...
      K3   K4   ] ;

CR = K * C * transpose(K) ;

return

