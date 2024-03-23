function qout = quatconj( qin )
%  QUATCONJ Calculate the conjugate of a quaternion.
%   N = QUATCONJ( Q ) calculates the conjugate, N, for a given quaternion, Q.  
%   Input Q is an M-by-4 matrix containing M quaternions.  N returns an 
%   M-by-4 matrix of conjugates.  Each element of Q must be a real number.  
%   Additionally, Q has its scalar number as the first column.
%
%   Examples:
%
%   Determine the conjugate of q = [1 0 1 0]:
%      conj = quatconj([1 0 1 0])
%
%#eml
if any(~isreal(qin(:)))
  error('Input elements are not real.');
end

if (size(qin,2) ~= 4)
  error('Input dimension is not M-by-4.');
end

qout = [ qin(:,1)  -qin(:,2:4) ];
