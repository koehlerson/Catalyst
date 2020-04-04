// Gmsh project created on Fri Apr  3 12:15:54 2020
SetFactory("OpenCASCADE");
//+
Sphere(1) = {0, 0, 0, 0.25, -Pi/2, Pi/2, 2*Pi};
//+
Physical Volume("Domain") = {1};
//+
Physical Surface("Boundary") = {1};
//+
Physical Surface("Boundary") += {1};
//+
Physical Volume("Domain") += {1};
