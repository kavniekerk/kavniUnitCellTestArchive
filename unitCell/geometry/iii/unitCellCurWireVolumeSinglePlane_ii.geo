
mesh_level = 0.025;
a_x1 = 0.05;

mwf = 1;
r_w = 0.0025;
p_0 = 0.025;
p_x1 = 0.025 - ( 0 * r_w/mwf );
R_x1 = ( p_x1 * p_x1 + r_w * r_w )/( 2 * r_w );
alpha = Asin( (p_x1 / R_x1) );

lcWireMesh = 0.001;

// Wire 1a0

p0_1a0 = newp; Point(p0_1a0) = {2*p_x1,	2*p_x1,       -1*r_w+mesh_level, 	lcWireMesh};			// centre circle
p1_1a0 = newp; Point(p1_1a0) = {2*p_x1,	2*p_x1,       -2*r_w+mesh_level, 	lcWireMesh};			// bottom circle
p2_1a0 = newp; Point(p2_1a0) = {2*p_x1,	2*p_x1 + r_w, -1*r_w+mesh_level, 	lcWireMesh};			// right circle
p3_1a0 = newp; Point(p3_1a0) = {2*p_x1,	2*p_x1, 	   0*r_w+mesh_level,    lcWireMesh};			// top circle
p4_1a0 = newp; Point(p4_1a0) = {2*p_x1,	2*p_x1 - r_w, -1*r_w+mesh_level, 	lcWireMesh};			// left circle

lcl1_1a0 = newc; Circle(lcl1_1a0) = {p1_1a0, p0_1a0, p2_1a0};
// Transfinite Curve {lcl1_1a0} = lcl1_1a0;
lcl2_1a0 = newc; Circle(lcl2_1a0) = {p2_1a0, p0_1a0, p3_1a0};
// Transfinite Curve {lcl2_1a0} = lcl2_1a0;
lcl3_1a0 = newc; Circle(lcl3_1a0) = {p3_1a0, p0_1a0, p4_1a0};
// Transfinite Curve {lcl3_1a0} = lcl3_1a0;
lcl4_1a0 = newc; Circle(lcl4_1a0) = {p4_1a0, p0_1a0, p1_1a0};
// Transfinite Curve {lcl4_1a0} = lcl4_1a0;

llps1_1a0 = newll; Line Loop(llps1_1a0) = {lcl1_1a0, lcl2_1a0, lcl3_1a0, lcl4_1a0};
sps1_1a0 = news;  Plane Surface(sps1_1a0) = {llps1_1a0};

// Wire 1a1

p0_1a1 = newp; Point(p0_1a1) = {2*p_x1 + a_x1/2,	2*p_x1,       -1*r_w+mesh_level, 	lcWireMesh};			// centre circle
p1_1a1 = newp; Point(p1_1a1) = {2*p_x1 + a_x1/2,	2*p_x1,       -2*r_w+mesh_level, 	lcWireMesh};			// bottom circle
p2_1a1 = newp; Point(p2_1a1) = {2*p_x1 + a_x1/2,	2*p_x1 + r_w, -1*r_w+mesh_level, 	lcWireMesh};			// right circle
p3_1a1 = newp; Point(p3_1a1) = {2*p_x1 + a_x1/2,	2*p_x1, 	   0*r_w+mesh_level,    lcWireMesh};			// top circle
p4_1a1 = newp; Point(p4_1a1) = {2*p_x1 + a_x1/2,	2*p_x1 - r_w, -1*r_w+mesh_level, 	lcWireMesh};			// left circle

lcl1_1a1 = newc; Circle(lcl1_1a1) = {p1_1a1, p0_1a1, p2_1a1};
// Transfinite Curve {lcl1_1a1} = lcl1_1a1;
lcl2_1a1 = newc; Circle(lcl2_1a1) = {p2_1a1, p0_1a1, p3_1a1};  
// Transfinite Curve {lcl2_1a1} = lcl2_1a1;
lcl3_1a1 = newc; Circle(lcl3_1a1) = {p3_1a1, p0_1a1, p4_1a1};  
// Transfinite Curve {lcl3_1a1} = lcl3_1a1;
lcl4_1a1 = newc; Circle(lcl4_1a1) = {p4_1a1, p0_1a1, p1_1a1};
// Transfinite Curve {lcl4_1a1} = lcl4_1a1;

// llps1_1a1 = newll; Line Loop(llps1_1a1) = {lcl1_1a1, lcl2_1a1, lcl3_1a1, lcl4_1a1};
// sps1_1a1 = news;  Plane Surface(sps1_1a1) = {llps1_1a1};

lsl1_1a1 = newl; Line(lsl1_1a1) = {p1_1a0, p1_1a1};
// Transfinite Line {lsl1_1a1} = lsl1_1a1;
lsl2_1a1 = newl; Line(lsl2_1a1) = {p2_1a0, p2_1a1};  
// Transfinite Line {lsl2_1a1} = lsl2_1a1;
lsl3_1a1 = newl; Line(lsl3_1a1) = {p3_1a0, p3_1a1};  
// Transfinite Line {lsl3_1a1} = lsl3_1a1;
lsl4_1a1 = newl; Line(lsl4_1a1) = {p4_1a0, p4_1a1};
// Transfinite Line {lsl4_1a1} = lsl4_1a1;

llcs1_1a1 = newll; Line Loop(llcs1_1a1) = {lsl1_1a1, lcl1_1a1, -lsl2_1a1, -lcl1_1a0};
llcs2_1a1 = newll; Line Loop(llcs2_1a1) = {lsl2_1a1, lcl2_1a1, -lsl3_1a1, -lcl2_1a0};
llcs3_1a1 = newll; Line Loop(llcs3_1a1) = {lsl3_1a1, lcl3_1a1, -lsl4_1a1, -lcl3_1a0};
llcs4_1a1 = newll; Line Loop(llcs4_1a1) = {lsl4_1a1, lcl4_1a1, -lsl1_1a1, -lcl4_1a0};

scs1_1a1 = news;  Ruled Surface(scs1_1a1) = {llcs1_1a1};
scs2_1a1 = news;  Ruled Surface(scs2_1a1) = {llcs2_1a1};
scs3_1a1 = news;  Ruled Surface(scs3_1a1) = {llcs3_1a1};
scs4_1a1 = news;  Ruled Surface(scs4_1a1) = {llcs4_1a1};

slcs_1a1[] = {};
slcs_1a1[] += { scs1_1a1, scs2_1a1, scs3_1a1, scs4_1a1 };

// Wire 1a2

p0_1a2 = newp; Point(p0_1a2) = {2*p_x1 + a_x1,	2*p_x1,       -1*r_w+mesh_level, 	lcWireMesh};			// centre circle
p1_1a2 = newp; Point(p1_1a2) = {2*p_x1 + a_x1,	2*p_x1,       -2*r_w+mesh_level, 	lcWireMesh};			// bottom circle
p2_1a2 = newp; Point(p2_1a2) = {2*p_x1 + a_x1,	2*p_x1 + r_w, -1*r_w+mesh_level, 	lcWireMesh};			// right circle
p3_1a2 = newp; Point(p3_1a2) = {2*p_x1 + a_x1,	2*p_x1, 	   0*r_w+mesh_level,    lcWireMesh};			// top circle
p4_1a2 = newp; Point(p4_1a2) = {2*p_x1 + a_x1,	2*p_x1 - r_w, -1*r_w+mesh_level, 	lcWireMesh};			// left circle

lcl1_1a2 = newc; Circle(lcl1_1a2) = {p1_1a2, p0_1a2, p2_1a2};
// Transfinite Curve {lcl1_1a2} = lcl1_1a2;
lcl2_1a2 = newc; Circle(lcl2_1a2) = {p2_1a2, p0_1a2, p3_1a2};  
// Transfinite Curve {lcl2_1a2} = lcl2_1a2;
lcl3_1a2 = newc; Circle(lcl3_1a2) = {p3_1a2, p0_1a2, p4_1a2};  
// Transfinite Curve {lcl3_1a2} = lcl3_1a2;
lcl4_1a2 = newc; Circle(lcl4_1a2) = {p4_1a2, p0_1a2, p1_1a2};
// Transfinite Curve {lcl4_1a2} = lcl4_1a2;

llps1_1a2 = newll; Line Loop(llps1_1a2) = {lcl1_1a2, lcl2_1a2, lcl3_1a2, lcl4_1a2};
sps1_1a2 = news;  Plane Surface(sps1_1a2) = {llps1_1a2};

lsl1_1a2 = newl; Line(lsl1_1a2) = {p1_1a1, p1_1a2};
// Transfinite Line {lsl1_1a2} = lsl1_1a2;
lsl2_1a2 = newl; Line(lsl2_1a2) = {p2_1a1, p2_1a2};
// Transfinite Line {lsl2_1a2} = lsl2_1a2;
lsl3_1a2 = newl; Line(lsl3_1a2) = {p3_1a1, p3_1a2};
// Transfinite Line {lsl3_1a2} = lsl3_1a2;
lsl4_1a2 = newl; Line(lsl4_1a2) = {p4_1a1, p4_1a2};
// Transfinite Line {lsl4_1a2} = lsl4_1a2;

llcs1_1a2 = newll; Line Loop(llcs1_1a2) = {lsl1_1a2, lcl1_1a2, -lsl2_1a2, -lcl1_1a1};
llcs2_1a2 = newll; Line Loop(llcs2_1a2) = {lsl2_1a2, lcl2_1a2, -lsl3_1a2, -lcl2_1a1};
llcs3_1a2 = newll; Line Loop(llcs3_1a2) = {lsl3_1a2, lcl3_1a2, -lsl4_1a2, -lcl3_1a1};
llcs4_1a2 = newll; Line Loop(llcs4_1a2) = {lsl4_1a2, lcl4_1a2, -lsl1_1a2, -lcl4_1a1};

scs1_1a2 = news;  Ruled Surface(scs1_1a2) = {llcs1_1a2};
scs2_1a2 = news;  Ruled Surface(scs2_1a2) = {llcs2_1a2};
scs3_1a2 = news;  Ruled Surface(scs3_1a2) = {llcs3_1a2};
scs4_1a2 = news;  Ruled Surface(scs4_1a2) = {llcs4_1a2};

slcs_1a2[] = {};
slcs_1a2[] += { scs1_1a2, scs2_1a2, scs3_1a2, scs4_1a2 };

sps1_1a0[] = sps1_1a0;
sps1_1a2[] = sps1_1a2;

sl_wire_exterior_surface_1a[] = newreg; Surface Loop(sl_wire_exterior_surface_1a) = { sps1_1a0[], slcs_1a1[], slcs_1a2[], sps1_1a2[] };
vol_1a_wire = newreg; Volume(vol_1a_wire) = sl_wire_exterior_surface_1a[];
physvol_1a_wire = newreg; Physical Volume(physvol_1a_wire) = vol_1a_wire;
physsurf_1a_wire = newreg; Physical Surface(physsurf_1a_wire) = { sps1_1a0[], slcs_1a1[], slcs_1a2[], sps1_1a2[] };





