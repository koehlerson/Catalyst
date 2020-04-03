// Gmsh project created on Fri Apr  3 12:15:54 2020
SetFactory("OpenCASCADE");
//+
Sphere(1) = {0, 0, 0, 2.5e-7, -Pi/2, Pi/2, 2*Pi};
//+
Physical Surface("Boundary") = {1};
//+
Physical Volume("Domain") = {1};
