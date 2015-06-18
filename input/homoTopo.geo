cl = 4;

Point(1) = { 0,    0,   0, cl };
Point(2) = { 1000, 0,   0, cl };
Point(3) = { 0,    500, 0, cl };
Point(4) = { 1000, 500, 0, cl };

Point(5) = { 500, 530, 0, cl };
Point(6) = { 420, 510, 0, cl };
Point(7) = { 580, 510, 0, cl };

Line(1) = { 1, 2 };
Line(2) = { 1, 3 };
Line(3) = { 2, 4 };

Spline(4) = { 3, 6, 5, 7, 4 };

Line Loop(5) = { 1, 3, -4, -2 };

Plane Surface(1) = { 5 };

Physical Surface(1) = { 1 };

left   = 11;
right  = 12;
bottom = 13;
top    = 14;
Physical Line(left)   = { 2 };
Physical Line(right)  = { 3 };
Physical Line(bottom) = { 1 };
Physical Line(top)    = { 4 };

