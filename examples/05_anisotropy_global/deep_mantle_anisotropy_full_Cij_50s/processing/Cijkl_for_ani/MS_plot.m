% MS_PLOT - Plot phasevels/anisotropy on pole figures.
%
% // Part of MSAT - The Matlab Seismic Anisotropy Toolkit //
%
% Given an elasticity matrix and density, produce pole figures showing 
%     the P- and S-wave anisotrpy. 
%
%  % MS_plot(C, rh, ...)
%
% Usage: 
%     MS_plot(C, rh)                    
%         Produce three pole figures showing P-wave velocity, S-wave
%         anisotropy and fast S-wave polarisation direction
%
%     MS_plot(..., 'fontsize', f)                    
%          Set minimum fontsize in plots to f.
%
%     MS_plot(..., 'wtitle', S)                    
%          Set the (window) title to be string S. If a minus sign is prepended
%          to this title, the figure number is also suppressed (the minus sign
%          is not shown). 
%
%     MS_plot(..., 'cmap', CM)                    
%          Redefine the colormap. CM can either be a (nx3) matrix containing a 
%          colormap, or a string describing a function to generate such a 
%          matrix (such as the built-in MATLAB colormap functions). E.g.:
%             MS_plot(C,rh,'cmap','cool(64)') - uses the MATLAB function cool
%                to generate a cyan-to-purple colourmap.
%          The default is 'jet', reversed so blue is fast/high, as is 
%          conventional for seismic velocity colorscales. 
%
%     MS_plot(..., 'reverse')
%          Unreverse the sense of the colour map. Default (no argument) is 
%          for blue to be fast/high as is conventional for seismic velocity 
%          colorscales. Setting this option results in red being fast/high.
%
%     MS_plot(..., 'pcontours', pcvect)
%     MS_plot(..., 'avscontours', ascvect)
%     MS_plot(..., 'scontours', scvect)
%          Set the contour levels for the P-wave velocity, S-wave
%          anisotropy and (optional) S-wave velocity plots. The values of 
%          pcvect, ascvect and scvect can be scalars indicating the number 
%          of contour lines to use or vectors, where each value represents 
%          a different contour line. Thus, the minimum, maximum values and 
%          spacing of contours can easily be set using the syntax 
%          [minval:spacing:maxval] for pcvect ascvect or scvect.
%
%     MS_plot(..., 'limitsonpol')
%          Include markers for the maximum and minimum values of S-wave
%          anisotropy on the plot of fast S-wave polarisation direction.
%          Not shown by default.
%
%     MS_plot(..., 'polsize' s1, s2, w1, w2)
%          Set the size of the markers used to indicate the fast S-wave
%          polarisation direction. There are two markers used a larger
%          white "background" marker with a length set by s1 (defaults to
%          0.18) and width set by w1 (defaults to 3.0) and an inner black
%          marker with width set by w2 (defaults to 2.0) and length set by
%          s2 (defaults to 0.18). A useful "neet" alternitive to the
%          default is s1 = 0.18, s2 = 0.16, w1 = 2.0 and w2 = 1.0.
%
%     MS_plot(..., 'plotmap', pmap)
%          Change the order, or content, of the subplots. This allows a
%          single pole figure (e.g. the P-wave velocity) to be drawn, or
%          the pole figures to be drawn in a different order. The argument
%          pmap must be a cell array of strings. Each string represents a
%          different pole figure and its location in the cell array
%          represents the location of the pole figure. Thus the default
%          {'VP', 'AVS', 'AVSPOL'} draws three pole figures in a row. If
%          the argument were {'VP'; 'AVS'} two pole figures would be drawn,
%          one above the other. Valid strings are 'VP", 'AVS' and 'AVSPOL'.
%          These are not case sensitive.
%
%     MS_plot(..., 'pdata', azimuth, inclination, vp)
%     MS_plot(..., 'sdata', azimuth, inclination, polarisation, avs)
%          Add data points to the P-wave velocity ('VP') or S-wave 
%          polarisation plot ('AVSPOL'), respectively. This can be used,
%          for example, to compare an elastic model with shear-wave
%          splitting measurements. In each case data points relate
%          to a ray propagating in a direction described by the azimuth 
%          and inclination. The point is coloured according to either 
%          the P-wave velocity, vp, or S-wave anisotropy, avs, in km/s
%          or %. In the case of S-wave anisotropy, the fast polarisation 
%          direction is also shown. Each argument (azimuth, inclination, 
%          vp, polarisation and avs) is an array and each must be the 
%          same length.
% 
%     MS_plot(..., 'quiet')
%          Don't write isotropic velocities to the terminal. 
%
% See also: MS_SPHERE, MS_PHASEVELS

% Copyright (c) 2011-2013 James Wookey and Andrew Walker
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


%=========================================================================
function MS_plot(C,rh,varargin)
%=========================================================================

%  ** Set defaults, these can be overriden in the function call
%  ** configure contouring options
      cvect = 10 ;     % number of contours
      VPcvect = cvect ;
      AVScvect = cvect ;
      VScvect = cvect ;
      
%  ** Put markers on SWS pol plot. 
      limitsonpol = 0;
   
%  ** Scaling for SWS plo plot.
      qwhite_scale = 0.18;
      qblack_scale = 0.18;
      qwhite_width = 3.0;
      qblack_width = 2.0;
      
%  ** configure font options
      fntsz = 12;

%  ** configure colormap options
      cmap = jet(64) ;
      icmapflip = 1 ; % reverse the sense of the colourscale
      
%  ** Write to matlab terminal?
      silentterm = 0 ;

%  ** default window title
      wtitle = 'MSAT polefigure.' ;
      
%  ** map of what we want to plot
      plotmap = {'VP', 'AVS', 'AVSPOL'};
      sdata_plot = 0; % no sdata overlay
      pdata_plot = 0; % no pdata overlay      
      
%  ** process the optional arguments
      iarg = 1 ;
      while iarg <= (length(varargin))
         switch lower(varargin{iarg})
            case 'reverse'
               icmapflip = 0 ;
               iarg = iarg + 1 ;
            case 'quiet'
               silentterm = 1;
               iarg = iarg + 1 ;
            case 'contours'
               % We keep this for backwards compat.
               % but note that it is not (and has never been)
               % documented.
               cvect = varargin{iarg+1} ;
               iarg = iarg + 2 ;
               VPcvect = cvect ;
               AVScvect = cvect ;
            case 'pcontours'
               VPcvect = varargin{iarg+1} ;
               iarg = iarg + 2 ;
            case 'avscontours'
               AVScvect = varargin{iarg+1} ;
               iarg = iarg + 2 ;
            case 'scontours'
               VScvect = varargin{iarg+1} ;
               iarg = iarg + 2 ;
            case 'limitsonpol'
               limitsonpol = 1;
               iarg = iarg + 1;
            case 'polsize'
               qwhite_scale = varargin{iarg+1};
               qblack_scale = varargin{iarg+2};
               qwhite_width = varargin{iarg+3};
               qblack_width = varargin{iarg+4};
               iarg = iarg + 5 ;
            case 'fontsize'
               fntsz = varargin{iarg+1} ;
               iarg = iarg + 2 ;
            case 'wtitle'
               wtitle = varargin{iarg+1} ;
               iarg = iarg + 2 ;
            case 'cmap'
               cmarg = varargin{iarg+1} ;
               if isstr(cmarg)
                  eval(['cmap = ' cmarg ';']) ;
               else
                  cmap = cmarg ;
               end
               iarg = iarg + 2 ;
            case 'plotmap'
               plotmap = varargin{iarg+1} ;
               iarg = iarg + 2 ;
            case 'sdata'
               sdata_azi = varargin{iarg+1};
               sdata_inc = varargin{iarg+2};
               sdata_pol = varargin{iarg+3};
               sdata_mag = varargin{iarg+4};
               sdata_plot = 1;
               iarg = iarg + 5;
            case 'pdata'
               pdata_azi = varargin{iarg+1};
               pdata_inc = varargin{iarg+2};
               pdata_mag = varargin{iarg+3};
               pdata_plot = 1;
               iarg = iarg + 4;
            otherwise 
               error(['Unknown option: ' varargin{iarg}]) ;   
         end   
      end 

      % What to plot - work out size and shape for later.
      assert(iscellstr(plotmap), 'MS:PLOT:BADPLOTMAP', ...
          'The plotmap must be a cell string');
      plotsize = size(plotmap);
      nrow = plotsize(1);
      ncol = plotsize(2);
      
      % check the inputs: C
      assert(MS_checkC(C)==1, 'MS:PLOT:badC', 'MS_checkC error MS_plot') ;

      % Check any data inputs
      if pdata_plot
          assert(length(pdata_azi)==length(pdata_inc),...
              'MS:PLOT:pdata_mismatch', ...
              'P-wave data arrays must be the same size.') ;
          assert(length(pdata_azi)==length(pdata_mag),...
              'MS:PLOT:pdata_mismatch', ...
              'P-wave data arrays must be the same size.') ;
      end
      if sdata_plot
          assert(length(sdata_azi)==length(sdata_inc),...
              'MS:PLOT:sdata_mismatch', ...
              'S-wave data arrays must be the same size.') ;
          assert(length(sdata_azi)==length(sdata_mag),...
              'MS:PLOT:sdata_mismatch', ...
              'S-wave data arrays must be the same size.') ;
          assert(length(sdata_azi)==length(sdata_pol),...
              'MS:PLOT:sdata_mismatch', ...
              'S-wave data arrays must be the same size.') ;
      end
%  ** buggy MATLAB contourf (version 7.1-7.3)
%
%     in these versions of MATLAB, the contourf routine doesn't seem able to
%     to cope with the spherical geometry surfaces. I therefore use 
%     contouf('v6',...) instead
%
      vstr = ver ;
      if vstr(1).Version >= 7.1 & vstr(1).Version <= 7.3
         buggyMATLAB = 1;
      else
         buggyMATLAB = 0;
      end   

%  ** set the viewing angle
      view_angle = [-90,90] ; % 1-North 2-West 3-Out of screen

      rad = pi./180 ;
      deg = 180./pi ;
      
      % Set up inc-az grids...
      [INC,AZ] = meshgrid([90:-6:0],[0:6:360]) ;
      
      % Invoke MS_phasevels to get wave velocities etc.
      [~,~,vs1,vs2,vp, S1P] = MS_phasevels(C,rh,...
        reshape(INC,61*16,1),reshape(AZ,61*16,1));
    
      % reverse so sph2cart() works properly
      AZ = -AZ;
      
      % Reshape results back to grids
      VS1 = reshape(vs1,61,16);
      VS2 = reshape(vs2,61,16);
      VP =  reshape(vp,61,16);
      VS1_x = reshape(S1P(:,1),61,16);
      VS1_y = reshape(S1P(:,2),61,16);
      VS1_z = reshape(S1P(:,3),61,16);
        
      % aVS data
      dVS =  (VS1-VS2) ;
      VSmean = (VS1+VS2)./2.0 ;
      AVS = 100.0*(dVS./VSmean) ;
      
%  ** output average velocities
      VPiso=mean(mean(VP)) ;
      VSiso=mean([mean(mean(VS1)) mean(mean(VS2))]) ;

      if ~silentterm
          fprintf('Isotropic average velocities: VP=%f, VS=%f\n',...
              VPiso,VSiso) ;
      end
%  ** Prepare window
      win_width = 466.667.*ncol;
      win_height = 400.*nrow;
      if strcmp(wtitle(1),'-')
         figure('Position',[1 1 win_width win_height],'name', ...
                wtitle(2:end),'NumberTitle','off') ;
      else
         figure('Position',[1 1 win_width win_height],'name',wtitle) ;
      end   
      
%  ** Setup a seismic colourmap (i.e. red->green->blue)
      if icmapflip
         cmap = flipud(cmap) ;
      end   
%  ** generate X/Y matrices for plotting
      [X,Y,Z] = sph2cart(AZ.*rad,INC.*rad,ones(size(AZ))) ;

      k = 0;
      for j = 1:nrow;
          for i = 1:ncol
              k = k + 1;
              subplot(nrow,ncol,k)
              
              switch lower(plotmap{j,i})
                  case 'vp'
                      % VP velocity plot
                      contour_pole(X, Y, VP, view_angle, VPcvect, ...
                          cmap, fntsz, buggyMATLAB, 'V_P (km/s)')
                      
                      [VPmin, VPmax] = max_min_pole(AZ, INC, VP);
                      
                      if pdata_plot
                          add_data(pdata_azi, pdata_inc, 0, pdata_mag, 0);
                      end
                      
                      %  ** add some information to the plot
                      VPlabel1 = sprintf(['Min. V_P =%6.2f, max.' ...
                          ' V_P =%6.2f'],VPmin,VPmax) ;
                      text(-1.15,0.8,VPlabel1,'FontSize',fntsz, ...
                          'FontWeight','bold') ;
                      VPmean = (VPmax+VPmin)./2.0 ;
                      VPani = (VPmax-VPmin)/VPmean .* 100 ;
                      VPlabel2 = sprintf('Anisotropy =%6.1f%%',VPani) ;
                      text(-1.3,0.8,VPlabel2,'FontSize',fntsz, ...
                          'FontWeight','bold') ;
                      
                  case 'avs'
                      % dVS anisotropy plot
                      contour_pole(X, Y, AVS, view_angle, AVScvect,...
                          cmap, fntsz, ...
                          buggyMATLAB, 'dV_S (%)')
                      
                      [AVSmin, AVSmax] = max_min_pole(AZ, INC, AVS);
                      
                      %  ** add some information to the plot
                      AVSlabel1 = sprintf(['V_S anisotropy, min' ...
                          ' =%6.2f, max =%6.2f'],AVSmin,AVSmax) ;
                      text(-1.15,1.0,AVSlabel1,'FontSize',fntsz,...
                          'FontWeight','bold') ;
                      
                  case 'vs1'
                      % VS1 anisotropy plot
                      contour_pole(X, Y, VS1, view_angle, VScvect,...
                          cmap, fntsz, ...
                          buggyMATLAB, 'V_{S1} (km/s)')
                      
                      [VS1min, VS1max] = max_min_pole(AZ, INC, VS1);
                      
                      %  ** add some information to the plot
                      VS1label1 = sprintf(['Min. V_{S1} =%6.2f, max.' ...
                          ' V_{S1} =%6.2f'],VS1min,VS1max) ;
                      text(-1.15,0.8,VS1label1,'FontSize',fntsz, ...
                          'FontWeight','bold') ;
                      VS1mean = (VS1max+VS1min)./2.0 ;
                      VS1ani = (VS1max-VS1min)/VS1mean .* 100 ;
                      VS1label2 = sprintf('Anisotropy =%6.1f%%',VS1ani) ;
                      text(-1.3,0.8,VS1label2,'FontSize',fntsz, ...
                          'FontWeight','bold') ;
                      
                  case 'vs2'
                      % VS2 anisotropy plot
                      contour_pole(X, Y, VS2, view_angle, VScvect,...
                          cmap, fntsz, ...
                          buggyMATLAB, 'V_{S2} (km/s)')
                      
                      [VS2min, VS2max] = max_min_pole(AZ, INC, VS1);
                      
                      %  ** add some information to the plot
                      VS2label1 = sprintf(['Min. V_{S2} =%6.2f, max.' ...
                          ' V_{S2} =%6.2f'],VS1min,VS1max) ;
                      text(-1.15,0.8,VS2label1,'FontSize',fntsz, ...
                          'FontWeight','bold') ;
                      VS2mean = (VS2max+VS2min)./2.0 ;
                      VS2ani = (VS2max-VS2min)/VS2mean .* 100 ;
                      VS2label2 = sprintf('Anisotropy =%6.1f%%',VS2ani) ;
                      text(-1.3,0.8,VS2label2,'FontSize',fntsz, ...
                          'FontWeight','bold') ;
                      
                  case 'avspol'
                      % Fast shear-wave polarisation plot
                      contour_pole(X, Y, AVS, view_angle, AVScvect, cmap, ...
                          fntsz, buggyMATLAB, 'Fast-shear polarisation')
                      
                      if (limitsonpol)
                          [~, ~] = max_min_pole(AZ, INC, AVS);
                      end
                      
                      pol_pole(VS1_x,VS1_y,VS1_z, X, Y, Z, AZ, INC, ...
                          qwhite_scale, qblack_scale, qwhite_width, ...
                          qblack_width);
                      
                      if sdata_plot
                          add_data(sdata_azi, sdata_inc, sdata_pol, ...
                              sdata_mag, 1);
                      end
                      
                  otherwise
                      error('MS:PLOT:BADPLOTMAP', ['An elelment of the '...
                          'plotmap was not recognised']);
              end
          end
      end
      
end
%===============================================================================

%===============================================================================
   function [AN,BN,CN] = vnormalise2(A,B,C)
%===============================================================================
%
%     normalise a 3 vector
%
      VMAG = sqrt( A.^2 + B.^2 + C.^2 )  ; 

      AN = A./VMAG ;
      BN = B./VMAG ;
      CN = C./VMAG ;
      
   end
%===============================================================================

%===============================================================================
   function [xr,yr,zr] = rotate_pm_vector(x,y,z,az,in)
%===============================================================================
%           
%     Rotate the particle motion vector so propagation direction is vertical
%     (to plot in a 2D way)
%
      rad = pi./180;
      V=[x,y,z]' ;
      
%  ** build up rotation matrices      
      gamma = az .* rad ;
      
%  ** first rotation (about Z)
      R1 = [  cos(gamma), sin(gamma), 0 ;...
             -sin(gamma), cos(gamma), 0 ;...
                 0      ,    0      , 1 ] ;
                    
%  ** second rotation (about Y)
      beta = (90-in).*rad ;
      R2 = [  cos(beta), 0 , -sin(beta) ;...
                 0     , 1 ,     0      ;...
              sin(beta), 0 ,  cos(beta) ] ;

%  ** third rotation (about Z, reversing first rotation)
      gamma = -az .* rad ;
      R3 = [  cos(gamma), sin(gamma), 0 ;...
             -sin(gamma), cos(gamma), 0 ;...
                 0      ,    0      , 1 ] ;

%  ** apply rotation
      VR = R3 * R2 * (R1 * V) ;
      xr = VR(1) ; 
      yr = VR(2) ;
      zr = VR(3) ;

   end
%===============================================================================

%===============================================================================
   function [imin,jmin,Zmin,imax,jmax,Zmax] = minmax2d(Z)
%===============================================================================
%           
%     Find the maximum and minimum values in a (2d) matrix and return their
%     values and indices
%
      [A,I] = min(Z) ; [Zmin,II] = min(A) ; imin = I(II) ; jmin = II ;
      [A,I] = max(Z) ; [Zmax,II] = max(A) ; imax = I(II) ; jmax = II ;

   end
%===============================================================================

function add_data(data_azi, data_inc, data_pol, data_mag, with_pol)

    % Get data points as XYZ
    % reverse so sph2cart() works properly
    data_azi = -data_azi;
    rad = pi./180 ;
    % Data points
    [X,Y,Z] = sph2cart(data_azi.*rad, data_inc.*rad, ones(size(data_azi)));
    
    if with_pol
        % Vectors in ray frame:
        nsdata = length(data_pol);
        vec = [ zeros(nsdata,1), sind(data_pol)', cosd(data_pol)'];
        V_x = zeros(nsdata,1);
        V_y = zeros(nsdata,1);
        V_z = zeros(nsdata,1);
        % Rotate into fame used to describe ray direction
        for i = 1:nsdata
            vecout = V_rot_gam(V_rot_bet(vec(i,:),-data_inc(i)),data_azi(i))';
            V_x(i) = vecout(1);
            V_y(i) = vecout(2);
            V_z(i) = vecout(3);
        end
    
        %  ** transform vectors
        [V_x,V_y,V_z] = vnormalise2(V_x,V_y,V_z) ;
        [XN,YN,ZN] = vnormalise2(X,Y,Z) ;

        A = zeros(3,nsdata) ;
        B = zeros(3,nsdata) ;
     
        A(1,:) = XN ;
        A(2,:) = YN ;
        A(3,:) = ZN ;
                  
        B(1,:) =  V_x ;
        B(2,:) =  V_y ;
        B(3,:) =  V_z ;
     
        C=cross(A,B) ;
        D=cross(A,C) ;
        
        VR_x = D(1,:) ;
        VR_y = D(2,:) ;
        VR_z = D(3,:) ;
      
        [VR_x,VR_y,VR_z] = vnormalise2(VR_x,VR_y,VR_z) ;
    
        %  ** rotate the particle motion vector so the normal to sphere is vertical
        VR_xR = zeros(1,nsdata);
        VR_yR = zeros(1,nsdata);
        VR_zR = zeros(1,nsdata);
        for ip = 1:nsdata
             [VR_xR(ip),VR_yR(ip),VR_zR(ip)] = ...
                        rotate_pm_vector(...
             VR_x(ip),VR_y(ip),VR_z(ip),...
             data_azi(ip),data_inc(ip));          
        end 
        
        h=quiver(X,Y,VR_xR,VR_yR,0.1,'w.') ;
        set(h,'LineWidth',3.0) ;

        h=quiver(X,Y,-VR_xR,-VR_yR,0.1,'w.') ;
        set(h,'LineWidth',3.0) ;
      
        %pol_pole(V_X', V_Y', V_Z', X, Y, Z, data_azi, data_inc, ...
        %      0.20, 0.20, 5.0, 3.0)

    end
    % Pivot arrays to prevent data_mag looking like a colour.
    scatter(X',Y',25,data_mag','.')
    scatter(X,Y,26,'wo')
    
end
      
function contour_pole(X, Y, vals, view_angle, cvect, cmap, fntsz, ...
    buggyMATLAB, titletext)

     if (min(min(vals)) ~= max(max(vals)))
         if buggyMATLAB
            [h1,h2]=contourf('v6',X,Y,vals,cvect) ;
            for j=1:length(h2)
               set(h2(j),'LineStyle','none')
            end
         else
            contourf(X,Y,vals,cvect,'LineStyle','none') ;
            if (~isscalar(cvect))
                % Limit the contours to the input vector.
                caxis([min(cvect) max(cvect)]);
            end
         end
      else
         surf(X,Y,zeros(size(vals)),vals) ;
         shading flat ;
      end   
      view(view_angle)
      colormap(cmap) ;
      daspect([1 1 1]);
      colorbar('FontSize',fntsz) ;
      axis off
      
      title(titletext,'FontSize',fntsz+4,'FontWeight','bold') ;

      hold on

end

function [VALmin, VALmax] = max_min_pole(AZ, INC, VAL)

      rad = pi./180 ;

      % get min and max values and positions      
      [imin,jmin,VALmin,imax,jmax,VALmax] = minmax2d(VAL) ;
      
      [VALmin_x,VALmin_y]=sph2cart(AZ(imin,jmin)*rad,INC(imin,jmin)*rad,1) ;   
      [VALmax_x,VALmax_y]=sph2cart(AZ(imax,jmax)*rad,INC(imax,jmax)*rad,1) ;   

      % mark the max. min. values 
      h=plot(VALmax_x,VALmax_y,'ws') ;
      set(h,'MarkerFaceColor','black');
      h=plot(VALmin_x,VALmin_y,'wo') ;
      set(h,'MarkerFaceColor','black');
end

function pol_pole(V_x,V_y,V_z, X, Y, Z, AZ, INC, ...
          qwhite_scale, qblack_scale, qwhite_width, qblack_width)
      
%  ** transform vectors
      [V_x,V_y,V_z] = vnormalise2(V_x,V_y,V_z) ;
      [XN,YN,ZN] = vnormalise2(X,Y,Z) ;

      A = zeros(3,61,16) ;
      B = zeros(3,61,16) ;
     
      A(1,:,:) = XN ;
      A(2,:,:) = YN ;
      A(3,:,:) = ZN ;
                  
      B(1,:,:) =  V_x ;
      B(2,:,:) =  V_y ;
      B(3,:,:) =  V_z ;
     
      C=cross(A,B) ;
      D=cross(A,C) ;
      
      VR_x(:,:) = D(1,:,:) ;
      VR_y(:,:) = D(2,:,:) ;
      VR_z(:,:) = D(3,:,:) ;
      
      [VR_x,VR_y,VR_z] = vnormalise2(VR_x,VR_y,VR_z) ;
 
%  ** define subset of polarisations to plot      
      cl = [1,3,5,7,9,11,15] ;
      drw = [60 10 5 3 2 2 2] ;

      ii=0 ;
      ip=0 ;
      np=136 ;
      ind1 = zeros(1,np) ;
      ind2 = zeros(1,np) ;
      
      for iinc=cl
         ii=ii+1 ;
         for iaz=1:drw(ii):61
            ip = ip + 1;
            ind1(ip) = iaz ;
            ind2(ip) = iinc ;
         end
      end   
      
%  ** rotate the particle motion vector so the normal to sphere is vertical
      for ip = 1:np
         [VR_xR(ind1(ip),ind2(ip)),VR_yR(ind1(ip),...
                  ind2(ip)),VR_zR(ind1(ip),ind2(ip))] = ...
                    rotate_pm_vector(...
         VR_x(ind1(ip),ind2(ip)),VR_y(ind1(ip),...
                ind2(ip)),VR_z(ind1(ip),ind2(ip)),...
         AZ(ind1(ip),ind2(ip)),INC(ind1(ip),ind2(ip)));          
      end   

%  ** form the subsets
      X2 = zeros(1,np) ;
      Y2 = zeros(1,np) ;
      Z2 = zeros(1,np) ;
      U2 = zeros(1,np) ;
      V2 = zeros(1,np) ;
      W2 = zeros(1,np) ;

      for i=1:length(ind1) ;
         X2(i)=X(ind1(i),ind2(i)) ;
         Y2(i)=Y(ind1(i),ind2(i)) ;
         Z2(i)=Z(ind1(i),ind2(i)) ;
         U2(i)=VR_xR(ind1(i),ind2(i)) ;
         V2(i)=VR_yR(ind1(i),ind2(i)) ;
         W2(i)=VR_zR(ind1(i),ind2(i)) ;
      end      
      
%  ** plot the vectors (changed to 2D to allow plotting with contourf surface)     
%      h=quiver3(X2,Y2,Z2,U2,V2,W2,0.18,'k.') ;      
      h=quiver(X2,Y2,U2,V2,qwhite_scale,'w.') ;
      set(h,'LineWidth',qwhite_width) ;

      h=quiver(X2,Y2,-U2,-V2,qwhite_scale,'w.') ;
      set(h,'LineWidth',qwhite_width) ;

      h=quiver(X2,Y2,U2,V2,qblack_scale,'k.') ;
      set(h,'LineWidth',qblack_width) ;

      h=quiver(X2,Y2,-U2,-V2,qblack_scale,'k.') ;
      set(h,'LineWidth',qblack_width) ;
end

function [VR] = V_rot_gam(V,gam)

    %  Make rotation matrix
    g = gam * pi/180. ;

    RR = [ cos(g) sin(g) 0 ; -sin(g) cos(g) 0 ; 0 0 1 ] ;
    VR = V * RR ;
 
end

function [VR] = V_rot_bet(V,bet)

    %  Make rotation matrix
    b = bet * pi/180. ;

    RR = [ cos(b) 0 -sin(b) ; 0 1 0 ; sin(b) 0 cos(b) ] ;

    VR = V * RR ;
 
end