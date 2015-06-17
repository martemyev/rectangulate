// Gmsh project created on Mon Oct 13 22:08:00 2014
lc_1 = 2;
lc_2 = 2;
lc_3 = 2;
lc_4 = 2;
lc_5 = 2;
lc_6 = 2;
lc_7 = 2;
lc_8 = 2;
lc_9 = 2;
lc_10= 2;

Point(1) = {1000, 0, 0, lc_1};
Point(2) = {0, 0, 0, lc_1};
Point(3) = {0, -500, 0, lc_7};
Point(4) = {1000, -500, 0, lc_7};

// topo
Point(5) = {750, 10, 0, lc_1};
Point(6) = {500, 40, 0, lc_1};
Point(7) = {250, 10, 0, lc_1};

Point(8) = {1000, -20, 0, lc_1};
Point(9) = {750, -30, 0, lc_1};
Point(10) = {500, -10, 0, lc_1};
Point(11) = {250, -20, 0, lc_1};
Point(12) = {0, -20, 0, lc_1};

Point(13) = {1000, -40, 0, lc_2};
Point(14) = {750, -70, 0, lc_2};
Point(15) = {500, -100, 0, lc_2};
Point(16) = {250, -30, 0, lc_2};
Point(17) = {0, -40, 0, lc_2};

Point(18) = {1000, -100, 0, lc_4};
Point(19) = {500, -100, 0, lc_4};
Point(20) = {0, -100, 0, lc_5};

BSpline(1) = {1, 5, 6, 7, 2};
BSpline(2) = {8, 9, 10, 11, 12};
BSpline(3) = {13, 14, 15, 16, 17};
Line(4) = {18, 19};
Line(5) = {19, 20};

// bedrock
Point(21) = {1000, -200, 0, lc_4};
Point(22) = {800, -150, 0, lc_4};
Point(23) = {600, -170, 0, lc_4};

Point(24) = {1000, -300, 0, lc_5};
Point(25) = {700, -350, 0, lc_5};
Point(26) = {500, -230, 0, lc_5};
Point(27) = {250, -220, 0, lc_5};
Point(28) = {0, -200, 0, lc_5};

Point(29) = {1000, -400, 0, lc_6};
Point(30) = {750, -450, 0, lc_6};
Point(31) = {500, -350, 0, lc_6};
Point(32) = {250, -300, 0, lc_6};
Point(33) = {0, -300, 0, lc_6};

BSpline(6) = {21, 22, 23, 19};
BSpline(7) = {24, 25, 26, 27, 28};
BSpline(8) = {29, 30, 31, 32, 33};

// Left, Bottom, Right
Line(9) = {2, 12};
Line(10) = {12, 17};
Line(11) = {17, 20};
Line(12) = {20, 28};
Line(13) = {28, 33};
Line(14) = {33, 3};

Line(15) = {3, 4};

Line(16) = {4, 29};
Line(17) = {29, 24};
Line(18) = {24, 21};
Line(19) = {21, 18};
Line(20) = {18, 13};
Line(21) = {13, 8};
Line(22) = {8, 1};

left=11;
right=12;
bottom=13;
top=14;

Physical Line(top) = {1};
Physical Line(left) = {9, 10, 11, 12, 13, 14};
Physical Line(bottom) = {15};
Physical Line(right) = {16, 17, 18, 19, 20, 21, 22};

// karst collapse structure 1
Point(34) = {485, -30, 0, lc_8};
Point(35) = {480, -25, 0, lc_8};
Point(36) = {460, -30, 0, lc_8};
Point(37) = {440, -25, 0, lc_8};
Point(38) = {435, -30, 0, lc_8};
Point(39) = {445, -30, 0, lc_8}; //center
Point(40) = {440, -38.66, 0, lc_8};
Point(41) = {460, -45, 0, lc_8};
Point(42) = {470, -38.66, 0, lc_8};
Point(43) = {475, -30, 0, lc_8}; //center

BSpline(23) = {34, 35, 36, 37, 38, 40, 41, 42, 34};

// karst collapse structure 2
Point(44) = {635, -50, 0, lc_8};
Point(45) = {630, -45, 0, lc_8};
Point(46) = {610, -50, 0, lc_8};
Point(47) = {590, -45, 0, lc_8};
Point(48) = {585, -50, 0, lc_8};
Point(49) = {595, -50, 0, lc_8}; //center
Point(50) = {590, -58.66, 0, lc_8};
Point(51) = {610, -65, 0, lc_8};
Point(52) = {630, -58.66, 0, lc_8};
Point(53) = {625, -50, 0, lc_8}; //center

BSpline(24) = {44, 45, 46, 47, 48, 50, 51, 52, 44};

// karst collapse structure 3
Point(54) = {310, -367, 0, lc_9};
Point(55) = {295, -350, 0, lc_9};
Point(56) = {240, -365, 0, lc_9};
Point(57) = {165, -350, 0, lc_9};
Point(58) = {145, -365, 0, lc_9};
Point(59) = {210, -360, 0, lc_9}; 
Point(60) = {187.5, -398.97, 0, lc_9};
Point(61) = {240, -400, 0, lc_9};
Point(62) = {293.5, -398.97, 0, lc_9};
Point(63) = {270, -360, 0, lc_9};

BSpline(25) = {54, 55, 56, 57, 58, 60, 61, 62, 54};

// karst collapse structure 4
Point(64) = {670, -440, 0, lc_10};
Point(65) = {660, -435, 0, lc_10};
Point(66) = {610, -445, 0, lc_10};
Point(67) = {560, -435, 0, lc_10};
Point(68) = {550, -440, 0, lc_10};
Point(69) = {570, -440, 0, lc_10}; 
Point(70) = {560, -457.32, 0, lc_10};
Point(71) = {610, -470, 0, lc_10};
Point(72) = {660, -457.32, 0, lc_10};
Point(73) = {650, -440, 0, lc_10};

BSpline(26) = {64, 65, 66, 67, 68, 70, 71, 72, 64};

// karst cavities
Point(74) = {543, -38, 0, lc_8};
Point(75) = {535, -36, 0, lc_8};
Point(76) = {533, -38, 0, lc_8};
Point(77) = {535, -45, 0, lc_8};
Point(78) = {545, -50, 0, lc_8};
Point(79) = {550, -47, 0, lc_8};
Point(80) = {554, -40, 0, lc_8};
Point(81) = {550, -37, 0, lc_8};
BSpline(27) = {74, 75, 76, 77, 78, 79, 80, 81, 74};

Point(82) = {345, -32, 0, lc_8};
Point(83) = {340, -34, 0, lc_8};
Point(84) = {328, -28, 0, lc_8};
Point(85) = {320, -30, 0, lc_8};
Point(86) = {319, -32, 0, lc_8};
Point(87) = {320, -34, 0, lc_8};
Point(88) = {329, -38, 0, lc_8};
Point(89) = {340, -40, 0, lc_8};
BSpline(28) = {82, 83, 84, 85, 86, 87, 88, 89, 82};

Point(90) = {400, -43, 0, lc_8};
Point(91) = {395, -41, 0, lc_8};
Point(92) = {381, -43, 0, lc_8};
Point(93) = {379, -44, 0, lc_8};
Point(94) = {379, -47, 0, lc_8};
Point(95) = {380, -49, 0, lc_8};
Point(96) = {395, -45, 0, lc_8};
Point(97) = {400, -47, 0, lc_8};
BSpline(29) = {90, 91, 92, 93, 94, 95, 96, 97, 90};

Point(98)  = {500, -58, 0, lc_8};
Point(99)  = {495, -60, 0, lc_8};
Point(100) = {487, -58, 0, lc_8};
Point(101) = {485, -61, 0, lc_8};
Point(102) = {483, -62, 0, lc_8};
Point(103) = {488, -68, 0, lc_8};
Point(104) = {495, -65, 0, lc_8};
Point(105) = {500, -62, 0, lc_8};
BSpline(30) = {98, 99, 100, 101, 102, 103, 104, 105, 98};

Point(108) = {682, -38, 0, lc_8};
Point(109) = {677, -42, 0, lc_8};
Point(110) = {665, -40, 0, lc_8};
Point(111) = {663, -47, 0, lc_8};
Point(112) = {668, -48, 0, lc_8};
Point(113) = {672, -50, 0, lc_8};
Point(114) = {685, -45, 0, lc_8};
Point(115) = {684, -40, 0, lc_8};
BSpline(31) = {108, 109, 110, 111, 112, 113, 114, 115, 108};

// surface
Line Loop(32) = {1, 9, -2, 22};
Plane Surface(33) = {32};
Line Loop(34) = {2, 10, -3, 21};
Line Loop(35) = {28};
Line Loop(36) = {29};
Line Loop(37) = {23};
Line Loop(38) = {30};
Line Loop(39) = {27};
Line Loop(40) = {24};
Line Loop(41) = {31};
Plane Surface(42) = {34, 35, 36, 37, 38, 39, 40, 41};
Line Loop(43) = {3, 11, -5, -4, 20};
Plane Surface(44) = {43};
Line Loop(45) = {4, -6, 19};
Plane Surface(46) = {45};
Line Loop(47) = {5, 12, -7, 18, 6};
Plane Surface(48) = {47};
Line Loop(49) = {7, 13, -8, 17};
Plane Surface(50) = {49};
Line Loop(51) = {8, 14, 15, 16};
Line Loop(52) = {25};
Line Loop(53) = {26};
Plane Surface(54) = {51, 52, 53};

Plane Surface(55) = {35}; //small left
Plane Surface(56) = {36};
Plane Surface(57) = {37};
Plane Surface(58) = {38};
Plane Surface(59) = {39};
Plane Surface(60) = {40};
Plane Surface(61) = {41}; //small right

Plane Surface(62) = {52}; //big1
Plane Surface(63) = {53}; //big2

// material
Physical Surface(1) = {33};
Physical Surface(2) = {42};
Physical Surface(3) = {44};
Physical Surface(4) = {46};
Physical Surface(5) = {48};
Physical Surface(6) = {50};
Physical Surface(7) = {54};

Physical Surface(8) = {55, 56, 57, 58, 59, 60, 61}; //water
Physical Surface(9) = {62};
Physical Surface(10)= {63};

Mesh.SubdivisionAlgorithm = 0; // 0 - triangles, 1 - quadrilaterals, 2 - hexahedrals

source=newp; Point(source) = { 250, 12, 0, lc_1 };
Point{source} In Surface{33};

source=newp; Point(source) = { 500, 12, 0, lc_1 };
Point{source} In Surface{33};

source=newp; Point(source) = { 750, 12, 0, lc_1 };
Point{source} In Surface{33};

