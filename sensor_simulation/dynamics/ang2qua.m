function [qua]=ang2qua(ang)
% ang2qua() converts Euler(3,2,1) angles into a quaternion notation
%
%   [qua]=ang2qua(ang) converts Euler(3,2,1) angles into a quaternion 
%   notation (cos term in the first element). It can handle single vector
%   transformation as well as a series of vectors.
%
% Syntax:[qua]=angqua(ang)
%
% Input:  vector with the 3 angles ang(3,:)
%
% Output: vector with the quaternion qua(4,:)
%

% input checking
[m,n] = size(ang);
if ( m ~= 3 )
     error('Input must be a vector or matrix with 3 rows')
end

for i=1:n
   a1=ang(3,i)/2; % the first rotational angle
   a2=ang(2,i)/2;
   a3=ang(1,i)/2;
   c1=cos(a1);
   c2=cos(a2);
   c3=cos(a3);
   s1=sin(a1);
   s2=sin(a2);
   s3=sin(a3);
% to save multiplications
   d1=c1*c2;
   d2=s2*s3;
   d3=s2*c3;
   d4=s1*c2;
   qua(1,i)=d1*c3 + s1*d2;
   qua(2,i)=d1*s3 - s1*d3;
   qua(3,i)=c1*d3 + d4*s3;
   qua(4,i)=d4*c3 - c1*d2;
end

return
%endfunction
