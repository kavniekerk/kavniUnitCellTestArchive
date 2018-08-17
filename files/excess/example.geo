Mesh.Smoothing = 10;
//Mesh.SecondOrderLinear = 1;
//Mesh.SubdivisionAlgorithm = 0;
//Mesh.RemeshAlgorithm = 0;
//Mesh.RemeshParametrization = 1;
//Mesh.RecombineAll = 1;
Mesh.Algorithm3D = 6;
Mesh.ElementOrder = 2;
Mesh.Optimize = 1;
//Mesh 3;
Geometry.CopyMeshingMethod=1;

rad1=1.0;
midRad1=Sqrt(0.5);
rad2=0.5;
midRad2=Sqrt(0.125);

transLine = 10;
transCirc = 2;

cl1 = 1e+22;
Point(1) = {0, 0, 0};
Point(2) = {0, rad1, 0};
Point(3) = {midRad1, midRad1, 0};
Point(4) = {0, rad2, 0};
Point(5) = {midRad2, midRad2, 0};
Circle(1) = {4, 1, 5};
Circle(2) = {2, 1, 3};
Line(3) = {5, 3};
Line(4) = {4, 2};
Line Loop(6) = {4, 2, -3, -1};
Ruled Surface(6) = {6};

Transfinite Line {1,2} = transCirc Using Progression 0.5;
Transfinite Line {3,4} = transLine Using Progression 0.5;
Transfinite Surface{6};
Recombine Surface{6};
Physical Surface("arc") = {6};

Extrude {0, 0, 0.2} {
   Surface{6}; Layers{ 2 }; Recombine;
}