function C = CosineDirectorMatrix(uVect,eastVect)
global r0Vect

bxVect = uVect/norm(uVect);
byVect = eastVect/norm(eastVect);
bzVect = cross(bxVect,byVect);

IxVect = [1 0 0];
IyVect = [0 1 0];
IzVect = [0 0 1];

C(1,1) = dot(bxVect,IxVect);
C(1,2) = dot(bxVect,IyVect);
C(1,3) = dot(bxVect,IzVect);
C(2,1) = dot(byVect,IxVect);
C(2,2) = dot(byVect,IyVect);
C(2,3) = dot(byVect,IzVect);
C(3,1) = dot(bzVect,IxVect);
C(3,2) = dot(bzVect,IyVect);
C(3,3) = dot(bzVect,IzVect);