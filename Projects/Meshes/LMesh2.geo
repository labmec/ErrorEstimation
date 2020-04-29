//+
Point(1) = {0, 0, 0, 1.0};
//+
Point(2) = {1, 0, 0, 1.0};
//+
Point(3) = {1, 1, 0, 1.0};
//+
Point(4) = {-1, 1, 0, 1.0};
//+
Point(5) = {-1, -1, 0, 1.0};
//+
Point(6) = {0, -1, 0, 1.0};
//+
Line(1) = {1, 2};
//+
Line(2) = {2, 3};
//+
Line(3) = {3, 4};
//+
Line(4) = {4, 5};
//+
Line(5) = {5, 6};
//+
Line(6) = {6, 1};
//+
Curve Loop(1) = {5, 6, 1, 2, 3, 4};
//+
Plane Surface(1) = {1};
//+
Physical Surface("domain") = {1};
//+
Physical Curve("dirichlet") = {3, 2, 1, 6, 5, 4};//+
Point(7) = {0, 0.5, 0, 1.0};
//+
Recursive Delete {
  Point{7}; 
}
//+
Point(7) = {0, 1, 0, 1.0};
//+
Point(8) = {0, 1,0, 1.0};
//+
Point(9) = {-1, 0, 0, 1.0};
//+
Line(7) = {4, 1};
//+
Line(8) = {7, 9};
//+
Line(9) = {2, 7};
//+
Line(10) = {3, 1};
//+
Line(11) = {9, 6};
//+
Line(12) = {1, 5};
//+
Point(10) = {0, 0.5, 0, 1.0};
//+
Point(11) = {-1, 0.5, -0, 1.0};
//+
Point(12) = {-0.5, 1, 0, 1.0};
//+
Point(13) = {-0.5, 1, 0, 1.0};
//+
Point(14) = {0.5, 1, 0, 1.0};
//+
Point(15) = {0.5, -0, 0, 1.0};
//+
Point(16) = {1, 0.5, 0, 1.0};
//+
Point(17) = {-1, -0.5, -0, 1.0};
//+
Point(18) = {-0.5, -0.5, -0, 1.0};
//+
Point(19) = {-0.5, -0, -0, 1.0};
//+
Point(20) = {-0, -0.5, -0, 1.0};
//+
Point(21) = {-0.5, -1, -0, 1.0};
//+
Line(13) = {11, 10};
//+
Line(14) = {10, 16};
//+
Line(15) = {14, 15};
//+
Line(16) = {12, 19};
//+
Line(17) = {17, 20};
//+
Line(18) = {21, 21};
//+
Line(19) = {19, 21};
//+
Line(20) = {9, 1};
//+
Line(21) = {1, 7};