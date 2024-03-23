function [qdot] = qderivative(t,q,tw,w)
    qdot = zeros(4,1);

    w1 = interp1(tw,w(:,1)',t);
    w2 = interp1(tw,w(:,2)',t);
    w3 = interp1(tw,w(:,3)',t);
    
    
    omega_mat = [0, -w1, -w2, -w3;
             w1, 0, w3, -w2;
             w2, -w3, 0, w1;
             w3, w2, -w1, 0];

    qdot_ = 1/2 * omega_mat * q;
    qdot(1) = qdot_(1);
    qdot(2) = qdot_(2);
    qdot(3) = qdot_(3);
    qdot(4) = qdot_(4);
    
    
end