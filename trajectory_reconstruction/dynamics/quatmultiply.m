function qout = quatmultiply( q, varargin )
%  QUATMULTIPLY Calculate the product of two quaternions.
%   N = QUATMULTIPLY( Q, R ) calculates the quaternion product, N, for two
%   given quaternions, Q and R.  Inputs Q and R can be either M-by-4 matrices 
%   containing M quaternions, or a single 1-by-4 quaternion.  N returns an 
%   M-by-4 matrix of quaternion products.  Each element of Q and R must be a
%   real number.  Additionally, Q and R have their scalar number as the first 
%   column.
%
%   Examples:
%
%   Determine the product of two 1-by-4 quaternions:
%      q = [1 0 1 0];
%      r = [1 0.5 0.5 0.75];
%      mult = quatmultiply(q, r)
%
%   Determine the product of a 1-by-4 quaternion with itself:
%      q = [1 0 1 0];
%      mult = quatmultiply(q)
%
%   Determine the product of 1-by-4 and 2-by-4 quaternions:
%      q = [1 0 1 0];
%      r = [1 0.5 0.5 0.75; 2 1 0.1 0.1];
%      mult = quatmultiply(q, r)
%
%   Note: Quaternion multiplication is not commutative.
%#eml
error(nargchk(1, 2, nargin,'struct'));

if any(~isreal(q(:)))
    error('First input elements are not real.');
end

if (size(q,2) ~= 4)
    error('First input dimension is not M-by-4.');
end

if nargin == 1
    r = q;
else
    r = varargin{1};
    if any(~isreal(r(:)))
        error('Second input elements are not real.');
    end
    if (size(r,2) ~= 4)
        error('Second input dimension is not M-by-4.');
    end
    if (size(r,1) ~= size(q,1) && ~( size(r,1) == 1 || size(q,1) == 1))
         error('Number of input rows are neither equal nor one.');
    end
end

% Calculate vector portion of quaternion product
% vec = s1*v2 + s2*v1 + cross(v1,v2)
vec = [q(:,1).*r(:,2) q(:,1).*r(:,3) q(:,1).*r(:,4)] + ...
         [r(:,1).*q(:,2) r(:,1).*q(:,3) r(:,1).*q(:,4)]+...
         [ q(:,3).*r(:,4)-q(:,4).*r(:,3) ...
           q(:,4).*r(:,2)-q(:,2).*r(:,4) ...
           q(:,2).*r(:,3)-q(:,3).*r(:,2)];

% Calculate scalar portion of quaternion product
% scalar = s1*s2 - dot(v1,v2)
scalar = q(:,1).*r(:,1) - q(:,2).*r(:,2) - ...
             q(:,3).*r(:,3) - q(:,4).*r(:,4);
    
qout = [scalar  vec];
       
