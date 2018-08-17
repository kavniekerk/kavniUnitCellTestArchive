charLen = 1;
Point(1) = {0,0,0,charLen};

Point(2) = {1,0,0,charLen};
Point(3) = {0,1,0,charLen};
Point(4) = {-1,0,0,charLen};
Point(5) = {0,-1,0,charLen};

Circle(1) = {2,1,3};
Circle(2) = {3,1,4};
Circle(3) = {4,1,5};
Circle(4) = {5,1,2};

Line Loop(5) = {1,2,3,4};
Plane Surface(6) = {5};

extr1[] = Extrude {0,0,20} {
    Surface{6};
};
extr2[] = Extrude { {0,-1,0}, // direction of rotation axis
                     {2,0,20}, // a point on the rotation axis
                     -Pi/2 } { // the rotation angle
    Surface{extr1[0]};
};
Extrude {10,0,0} {
   Surface{extr2[0]};
}