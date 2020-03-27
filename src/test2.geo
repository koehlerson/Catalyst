// Gmsh project created on Mon Mar 23 13:47:35 2020
SetFactory("OpenCASCADE");
//+
Box(1) = {0, 0, 0, 1, 1, 1};
//+
Physical Surface("Boundary") = {1, 4, 6, 2, 3, 5};
//+
Physical Volume("Domain") = {1};
