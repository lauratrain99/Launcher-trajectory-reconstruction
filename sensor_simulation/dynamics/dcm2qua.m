function [q]=dcm2qua(c)
% dcm2qua() converts a Direction Cosine Matrix into the equivalent quaternion
% 
%   [q]=dcm2qua(c) performs a conversion from a Direction Cosine Matrix C_b/a from
%   frame A to frame B, into the equivalent quaternion q_b/a.
%
% Syntax:   [q] = dcm2qua(c)
%
% Input:    c: direction cosine matrix c(3x3)
%
% Output:   q: quaternion vector(4x1)

[m,n] = size(c);

if m*n ~= 9
   error('dcm2qua: the input matrix must be 3x3');
end

out     = zeros(4,1);
out(1)  = 0.5*sqrt(1+c(1,1)+c(2,2)+c(3,3));

if out(1) ~= 0
    d      = 1/(4*out(1));
    out(2) = d*(c(2,3)-c(3,2));
    out(3) = d*(c(3,1)-c(1,3));
    out(4) = d*(c(1,2)-c(2,1));
else
    out(2) = 0.5*sqrt(1+c(1,1)-c(2,2)-c(3,3));
    if out(2) ~= 0
        out(3) = c(2,1)/(2*out(2));
        out(4) = c(3,1)/(2*out(2));
    else
        out(3) = 0.5*sqrt(1-c(1,1)+c(2,2)-c(3,3));
        if out(3) ~= 0
            d = 0.5*sqrt(1-c(1,1)-c(2,2)+c(3,3));
            if (c(2,3) < 0) d = -d; end
            out(4) = d;
        else
            out(4) = 1;
        end
        
    end
end

q = out(:); % to output a column vector

end
%endfunction
