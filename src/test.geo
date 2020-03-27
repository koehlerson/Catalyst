//+
SetFactory("OpenCASCADE");
Box(1) = {0, 0, 0, 1, 1, 1};
//+
Sphere(2) = {0.5, 0.5, 0.5, 0.25, -Pi/2, Pi/2, 2*Pi};
//+
BooleanDifference(3) = { Volume{1}; Delete; }{ Volume{2}; Delete;};
//+
Physical Surface(1) = {7};
//+
Physical Surface("Sphere") = {7};
//+
Physical Volume("Domain") = {3};
//+
Physical Surface("Sphere") += {7};
//+
Physical Surface("Boundary") = {4, 3, 1, 2, 5, 6};
//+
Physical Volume("Domain") += {3};
//+
Physical Volume("Domain") += {3};
//+
Physical Surface("Sphere") += {7};
//+
Physical Surface("Boundary") += {4, 3, 6, 2, 5, 1};
//+
Physical Surface("Sphere") += {7};
//+
Physical Surface("Sphere") += {7};
//+
Physical Surface("Boundary") += {4, 6, 5, 2, 3, 1};
//+
Physical Volume("Domain") += {3};
//+
Physical Volume("Domain") += {3};
//+
Physical Surface("Boundary") += {1, 4, 3, 5, 6, 2};
