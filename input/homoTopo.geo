cl = 10;

Point(1) = { 0,    0,   0, cl };
Point(2) = { 1000, 0,   0, cl };
Point(3) = { 0,    500, 0, cl };
Point(4) = { 1000, 500, 0, cl };

Point(5) = { 500, 530, 0, cl };
Point(6) = { 450, 510, 0, cl };
Point(7) = { 550, 510, 0, cl };

Line(1) = { 1, 2 };
Line(2) = { 1, 3 };
Line(3) = { 2, 4 };

BSpline(4) = { 3, 6, 5, 7, 4 };

Line Loop(5) = { 1, 3, -4, -2 };

Plane Surface(1) = { 5 };

Physical Surface(1) = { 1 };

