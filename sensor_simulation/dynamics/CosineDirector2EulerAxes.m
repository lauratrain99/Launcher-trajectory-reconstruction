function [eHat,PHI] = CosineDirector2EulerAxes(C_Matrix)


PHI = acos(1/2*(trace(C_Matrix)-1));
eHat(1) = 1/(2*sin(PHI))*(C_Matrix(2,3)-C_Matrix(3,2));
eHat(2) = 1/(2*sin(PHI))*(C_Matrix(3,1)-C_Matrix(1,3));
eHat(3) = 1/(2*sin(PHI))*(C_Matrix(1,2)-C_Matrix(2,1));


end
