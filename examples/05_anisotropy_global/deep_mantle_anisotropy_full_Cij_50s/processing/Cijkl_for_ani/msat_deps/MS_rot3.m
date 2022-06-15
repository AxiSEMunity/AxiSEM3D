% MS_ROT3 - Elasticity matrix rotation.
%
% // Part of MSAT - The Matlab Seismic Anisotropy Toolkit //
%
% Rotates an elasticity matrix around the three axes. 
%
% Usage: 
%     [CR] = MS_rot3(C,alp,bet,gam)                    
%         C: input 6x6 elasticity matrix 
%         alp: clockwise rotation about 1-axis, looking at origin
%         bet: clockwise rotation about 2-axis, looking at origin
%         gam: clockwise rotation about 3-axis, looking at origin
%
%     [CR] = MS_rot3(...,'order',V)                    
%         order: V is a 3 element vector containing the order in which to apply the 3 rotations.
%                Default is [1 2 3].
%
% Notes:
%     Angles are given in degrees and correspond to yaw, -dip and aximuth,
%     respectvly. The rotations are applied in order, ie: alpha, then beta
%     then gamma (by default). The variables C, alp, bet and gam can be 
%     arrays or scalars but (1) the angles must all be the
%     same length and (2) either C or all the angles must be scalars unless
%     they are the same length.
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

function [CR] = MS_rot3(C,alp,bet,gam,varargin)

    orderV = [1 2 3] ;

    %  ** process the optional arguments
    iarg = 1 ;
    while iarg <= (length(varargin))
       switch lower(varargin{iarg})
          case 'order' 
             orderV = varargin{iarg+1} ;
             if (length(find(orderV==1))+length(find(orderV==2))+length(find(orderV==3)))~=3
                error('MS:ROT3:BadOrder',...
                'Order vector must be 3-element, containing 1,2 and 3') ;   
             end
             if (length(orderV)~=3)
                error('MS:ROT3:BadOrder',...
                'Order vector must be 3-element, containing 1,2 and 3') ;  
             end                  
             iarg = iarg + 2 ;
           otherwise 
             error('MS:ROT3:UnknownOption',...
                ['Unknown option: ' varargin{iarg}]) ;   
       end   
    end

    % How many angles - are lists correct length.
    numrs = size(alp, 1);
    assert(numrs == size(bet, 1), 'MS:ListsMustMatch',...
        'The number of alp angles must be the same as the number of bet angles') 
    assert(numrs == size(gam, 1), 'MS:ListsMustMatch',...
        'The number of alp angles must be the same as the number of gam angles')
    

    a = alp * pi/180. ;
    b = bet * pi/180. ;
    g = gam * pi/180. ;

    numdimCs = ndims(C);
    if ((numrs == 1) && (numdimCs == 2))
        % Two matrix scalar case
        CR = rot_3_scalar(C, a, b, g, orderV);  
    elseif ((numrs > 1) && (numdimCs == 2))
        % Many rotations and one matrix case
        CR = zeros(6,6,numrs);
        for i = 1:numrs
            CR(:,:,i) = rot_3_scalar(C, a(i), b(i), g(i), orderV);
        end
    elseif ((numrs == 1) && (numdimCs == 3))
        % Many matrix and one rotation case
        CR = zeros(6,6,size(C, 3));
        for i = 1:size(C, 3)
            CR(:,:,i) = rot_3_scalar(C(:,:,i), a, b, b, orderV);
        end
    elseif ((numrs > 1) && (numdimCs == 3))
        % List of rotatiosn and matrices case
        assert((numrs == size(C, 3)), 'MS:ListsMustMatch', ...
            'The length of the lists of rotations and matrices must match');
        CR = zeros(6,6,size(C, 3));
        for i = 1:size(C,3)
            CR(:,:,i) = rot_3_scalar(C(:,:,i), a(i), b(i), g(i), orderV);
        end
    end
    
end

function [CR] = rot_3_scalar(C, a, b, g, orderV)

    try 
        MS_checkC(C) ;
    catch ME
        error('MS:ROT3BadInputMatrix', ...
           ['Invalid input elasticity matrix: ' ME.message]) ; 
    end

    R = zeros(3,3,3) ;

    R(1,:,:) = [ 1 0 0 ; 0 cos(a) sin(a) ; 0 -sin(a) cos(a) ] ;
    R(2,:,:) = [ cos(b) 0 -sin(b) ; 0 1 0 ; sin(b) 0 cos(b) ] ;
    R(3,:,:) = [ cos(g) sin(g) 0 ; -sin(g) cos(g) 0 ; 0 0 1 ] ;

    RR =  squeeze(R(orderV(3),:,:)) * ...
             squeeze(R(orderV(2),:,:)) * ...
                squeeze(R(orderV(1),:,:));
 
    %  Delegate the rotation
    CR = MS_rotR(C, RR) ;

end

