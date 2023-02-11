//+
Point(1) = {-1, -1, 0, 1.0};
//+
Point(2) = {1, -1, 0, 1.0};
//+
Point(3) = {1, 1, 0, 1.0};
//+
Point(4) = {-1, 1, 0, 1.0};
//+
Point(5) = {-0.5, 0.0, 0, 1.0};
//+
Point(6) = {-0.45, 0.0, 0, 1.0};
//+
Point(7) = {-0.5, 0.05, 0, 1.0};
//+
Point(8) = {-0.55, 0.0, 0, 1.0};
//+
Point(9) = {-0.5, -0.05, 0, 1.0};
//+
Line(1) = {1, 2};
//+
Line(2) = {2, 3};
//+
Line(3) = {3, 4};
//+
Line(4) = {4, 1};
//+
Circle(5) = {6, 5, 7};
//+
Circle(6) = {7, 5, 8};
//+
Circle(7) = {8, 5, 9};
//+
Circle(8) = {9, 5, 6};
//+
Curve Loop(1) = {6, 7, 8, 5};
//+
Curve Loop(2) = {4, 1, 2, 3};
//+
Plane Surface(1) = {1, 2};
//+
Physical Curve("WALL", 9) = {3, 1};
//+
Physical Curve("BODY", 10) = {5, 6, 7, 8};
//+
Physical Curve("INLET", 11) = {4};
//+
Physical Curve("OUTLET", 12) = {2};
