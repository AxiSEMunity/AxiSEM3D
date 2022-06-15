% MS_CHECKC - check consistency of a stiffness matrix against various criteria
%
% // Part of MSAT - The Matlab Seismic Anisotropy Toolkit //
%
% [ isgood ] = MS_checkC( C,... )
%
%	Inputs:
%     C = Stiffness tensor for checking
%
%  Outputs:
%     isgood = 1 for OK
%
%  Options:
%     Various checks can be disabled/enabled by specifying extra options:
%        'fast' - run a minimal set of checks (size and type only)
%        'warn' - warn only on more sophisticated errors than size/type
%        'nosymchk' - disable upper, lower symmetry checking
%        'nozerochk' - disable top left / trace zero checking 
%        'nopdefchk' - disable positive definiteness check 
% 
%     Various control parameters can be set by adding parameter,value pairs 
%     to the arguments
%        'thresh',x - set numerical threshold for zero to x (default = 1e-6)
%
%  Notes:
%  Currently, the following properties of the input elastic stiffness
%  matrix, C, are tested: (1) The maxtrix must be of rank 2 and size (6,6).
%  (2) The matrix must be symmetrical along the leading diagonal (i.e. 
%  the major symmetry of the elastic constants tensor, Cijkl == Cklij, must
%  be present). (3) The matrix must be positive-definite and be none zero
%  in the leading diagonal and top left corner.
%

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

function [isgood] = MS_checkC(C,varargin)

%  ** set defaults
      warn = 0 ;
      fast = 0 ;

      thresh = 1e-6 ;

      symchk = 1 ;
      zerochk = 1 ;
      pdefchk = 1 ;
      
%  ** process the optional arguments
      iarg = 1 ;
      while iarg<=(length(varargin))
         switch varargin{iarg}
            case 'fast'
               fast = 1 ;
               iarg = iarg + 1 ;
            case 'warn'
               warn = 1 ;
               iarg = iarg + 1 ;
            case 'nosymchk'
               symchk = 0 ;
               iarg = iarg + 1 ;
            case 'nozerochk'
               zerochk = 0 ;
               iarg = iarg + 1 ;
            case 'nopdefchk'
               pdefchk = 0 ;
               iarg = iarg + 1 ;
            case 'thresh'
               thresh = varargin{iarg+1} ;
               iarg = iarg + 2 ;
            otherwise 
               error('Unknown flag') ;   
         end   
      end 

      isgood = 1 ;

%  ** start with basic checks
      if ~isnumeric(C)
         error('MS:CHECKCNotNumeric', ...
	         'Elasticity matrix error: Appears not to be numeric.')
      end

      [nr nc] = size(C) ;
      if (nr~=6 | nc~=6)
         isgood = 0;
		   error('MS:CHECKCnot6x6',...
		      'Elasticity matrix error: Appears not to be a 6x6 matrix')
		end

%  ** if fast checking is selected, we're done here.       
      if fast, return, end

%  ** otherwise, perform more elaborate tests.
      
%  ** check symmetry    
      if symchk
         D=C-C';
         if (~isempty(find(abs(D)>thresh, 1)))
            isgood = 0;
            if warn
               warning('MS:CHECKCnotsym_W',...
                  'Elasticity matrix warning: non-symmetric')
            else
               error('MS:CHECKCnotsym',...
                  'Elasticity matrix error: non-symmetric')
            end   
         end   
      end

%  ** check that enough values are set non-zero
      if zerochk
         isopart = [C(1,1) C(1,2) C(1,3) C(2,2) C(2,3) ...
                    C(3,3) C(4,4) C(5,5) C(6,6)] ;
         if (~isempty(find(abs(isopart)<thresh, 1)))
            isgood = 0;
            if warn
               warning('MS:CHECKCbadzeros_W',...
                  'Elasticity matrix warning: zeros in trace or top left')
            else
               error('MS:CHECKCbadzeros',...
                  'Elasticity matrix error: zeros in trace or top left')
            end   
         end
      end

%  ** check that matrix is positive definite
      if pdefchk
         try
            [~]=chol(C) ;
         catch ME 
            % Only process the matlab posdef errors...
            if ~strcmp(ME.identifier, 'MATLAB:posdef')
                rethrow(ME)
            end
            isgood = 0;
            if warn
               warning('MS:CHECKCnotposdef_W',...
                  'Elasticity matrix warning: not positive definite')
            else
               error('MS:CHECKCnotposdef',...
                  'Elasticity matrix error: not positive definite')
            end
         end               
      end
return
%=======================================================================================  
