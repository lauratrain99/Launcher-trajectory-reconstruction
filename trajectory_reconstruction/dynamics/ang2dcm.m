function [eulermat] = ang2dcm(phi,theta,psi,ijk)
% ang2dcm() computes the direction cosine matrix from given Euler angles
%
%   [eulermat] = ang2dcm(phi,theta,psi,ijk) calculates the direction 
%   cosine matrix defined by successive rotations by angles phi, theta
%   and psi respectively, about coordinate axes i,j,k. Euler angle
%   rotations means that the first rotation by angle "phi" is
%   about the "i" axis, the second rotation by angle "theta" is
%   about the "j" axis, and the the third rotation by angle "psi"
%   is about the "k" axis.
%
% Syntax: c = ang2dcm(a1, a2, a3, ijk)
%
% input:    phi	   : first euler angle [rad]
%           theta  : second euler angle [rad]
%           psi    : third euler angle [rad]
% option:   ijk	   : rotation type sequence
%	
% output:   c      : euler angle rotation matrix



  
ijk='123';  % default ijk = "123"

% convert euler angles: phi, theta and psi from [deg] to [rad]
%torad = pi/180;
%phi   = phi*torad;
%theta = theta*torad;
%psi   = psi*torad;

% sines and cosines of euler angles
cphi   = cos(phi);
sphi   = sin(phi);
ctheta = cos(theta);
stheta = sin(theta);
cpsi   = cos(psi);
spsi   = sin(psi);

% limit for assigning elements to zero
eps    = 1.0e-9;

switch ijk

% euler angle representation: 1-2-3
case '123'
    eulermat = [ cpsi*ctheta, cpsi*stheta*sphi+spsi*cphi,-cpsi*stheta*cphi+spsi*sphi;...
                -spsi*ctheta,-spsi*stheta*sphi+cpsi*cphi, spsi*stheta*cphi+cpsi*sphi;...
        		      stheta,     -ctheta*sphi          ,       ctheta*cphi];

% euler angle representation: 1-3-2
case '132'
    eulermat = [ cpsi*ctheta, cpsi*stheta*cphi+spsi*sphi, cpsi*stheta*sphi-spsi*cphi;...
     		 -stheta     ,      ctheta*cphi		 ,      ctheta*sphi          ;...
		  spsi*ctheta,spsi*stheta*cphi-cpsi*sphi, spsi*stheta*sphi+cpsi*cphi];

% euler angle representation: 2-3-1
case '231'
    eulermat = [      ctheta*cphi  	      ,      stheta,     -ctheta*sphi;...
     		-cpsi*stheta*cphi + sphi*spsi , cpsi*ctheta, cpsi*stheta*sphi+spsi*cphi;...
		  spsi*stheta*cphi + cpsi*sphi,-spsi*ctheta,-spsi*stheta*sphi+cpsi*cphi];
-cpsi*stheta*cphi + sphi*spsi
% euler angle representation: 2-1-3
case  '213'
    eulermat = [ cpsi*cphi+spsi*stheta*sphi, spsi*ctheta,-cpsi*sphi+spsi*stheta*cphi;...
		 -spsi*cphi+cpsi*stheta*sphi, cpsi*ctheta, spsi*sphi+cpsi*stheta*cphi;...
		                 ctheta*sphi,     -stheta,                ctheta*cphi];

% euler angle representation: 3-1-2
case '312'
    eulermat = [ cpsi*cphi-spsi*stheta*sphi, cpsi*sphi+spsi*stheta*cphi,-spsi*ctheta;...
	                -ctheta*sphi,                ctheta*cphi,      stheta;...
		  spsi*cphi+cpsi*stheta*sphi, spsi*sphi-cpsi*stheta*cphi, cpsi*ctheta];

% euler angle representation: 3-2-1
case '321'
    eulermat = [             ctheta*cphi,                ctheta*sphi,     -stheta;...
           	 -cpsi*sphi+spsi*stheta*cphi, cpsi*cphi+spsi*stheta*sphi, spsi*ctheta;...
    		  spsi*sphi+cpsi*stheta*cphi,-spsi*cphi+cpsi*stheta*sphi, cpsi*ctheta];

% euler angle representation: 1-2-1
case '121'
    eulermat = [      ctheta,                stheta*sphi,               -stheta*cphi;...
		  spsi*stheta, cpsi*cphi-spsi*ctheta*sphi, cpsi*sphi+spsi*ctheta*cphi;...
		  cpsi*stheta,-spsi*cphi-cpsi*ctheta*sphi,-spsi*sphi+cpsi*ctheta*cphi];

% euler angle representation: 1-3-1
case '131'
    eulermat = [   ctheta,      stheta*cphi          ,      stheta*sphi          ;...
     		 -cpsi*stheta, cpsi*ctheta*cphi-spsi*sphi, cpsi*ctheta*sphi+spsi*cphi;...
     		  spsi*stheta,-spsi*ctheta*cphi-cpsi*sphi,-spsi*ctheta*sphi+cpsi*cphi];

% euler angle representation: 2-1-2
case '212'
     eulermat = [ cpsi*cphi-spsi*ctheta*sphi, spsi*stheta,-cpsi*sphi-spsi*ctheta*cphi;...
		                 stheta*sphi,      ctheta,                stheta*cphi;...
		  spsi*cphi+cpsi*ctheta*sphi,-cpsi*stheta,-spsi*sphi+cpsi*ctheta*cphi];

% euler angle representation: 2-3-2
case '232'
     eulermat = [ cpsi*ctheta*cphi-spsi*sphi,cpsi*stheta,-cpsi*ctheta*sphi-spsi*cphi;...
		      -stheta*cphi          ,      ctheta,      stheta*sphi          ;...
		  spsi*ctheta*cphi+cpsi*sphi, spsi*stheta,-spsi*ctheta*sphi+cpsi*cphi];

% euler angle representation: 3-1-3
case '313'
     eulermat = [ cpsi*cphi-spsi*ctheta*sphi, cpsi*sphi+spsi*ctheta*cphi, spsi*stheta;...
		 -spsi*cphi-cpsi*ctheta*sphi,-spsi*sphi+cpsi*ctheta*cphi, cpsi*stheta;...
		                 stheta*sphi,               -stheta*cphi,      ctheta];

% euler angle representation: 3-2-3
case '323'
     eulermat = [ cpsi*ctheta*cphi-spsi*sphi, cpsi*ctheta*sphi+spsi*cphi,-cpsi*stheta;...
		 -spsi*ctheta*cphi-cpsi*sphi,-spsi*ctheta*sphi+cpsi*cphi, spsi*stheta;...
		       stheta*cphi          ,      stheta*sphi          ,      ctheta];

otherwise
    error(['euler2dcm: unknown rotation type sequency: ''' ijk '''']);
      
end

dum           = find(abs(eulermat) < eps);
eulermat(dum) = 0.0;

return
% end of function
