function [out] = qua2dcm(in)
% qua2dcm() converts a quaternion into the equivalent Direction Cosine Matrix
%


out  = zeros(3,3);
[m,n]=size(in);
if m*n ~= 4
    error('qua2dcm: The input must be q quaternion vector (4x1)');
end

q12=in(2)*in(2);
q22=in(3)*in(3);
q32=in(4)*in(4);

q1_2=in(2)*in(3);
q0_3=in(1)*in(4);
q1_3=in(2)*in(4);
q0_2=in(1)*in(3);
q2_3=in(3)*in(4);
q0_1=in(1)*in(2);

%building rotation matrix row by row

out(1,1)=1-2*(q22+q32);
out(1,2)=2*(q1_2+q0_3);
out(1,3)=2*(q1_3-q0_2);
out(2,1)=2*(q1_2-q0_3);
out(2,2)=1-2*(q12+q32);
out(2,3)=2*(q2_3+q0_1);
out(3,1)=2*(q1_3+q0_2);
out(3,2)=2*(q2_3-q0_1);
out(3,3)=1-2*(q12+q22);

return 
% endfunction
