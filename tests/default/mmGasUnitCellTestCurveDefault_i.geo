a = 0.045;                // the "pitch", or distance between GEM pillars, in mm

mwf = 1;                                      // mesh_window_factor

mm = 1;                                       // geometrical scaling
r_w = 0.010 * mm;                             // R of Wiremesh, in microns
hp_0 = 0.020;                                    // half pitch of the window, in mm
hp = 0.020 * mm - 0*r_w/mwf * mm;                // half pitch of the window, in microns

p = hp_0;                                       // half pitch of the window, in mm

R = (p * p + r_w * r_w)/( (2 * r_w) );            // R
alpha = Asin((p/R));                          // angle in radians

mesh_level = 0.000;                            // mesh level, in mm
mesh_window = 0.020;                           // mesh window, in mm

sp_fac = p*0.00;

x1_sp_wind_fac = p*0.00;
x2_sp_wind_fac = p*0.00;
y1_sp_wind_fac = p*0.00;
y2_sp_wind_fac = p*0.00;

Rtp = R - R*0.00;
Rtn = R + R*0.00;

h_f = 0*r_w;                                  // Heightfactor

geo_f_x = 1;
geo_f_y = 1;

m = 0;
n = 0;

// Characteristic lengths

  lcCopperPlateBdry = 0.0001;
  lcExtElectrodeBdry = 0.0001;
  LcWiremesh = 0.0001; 


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////// Face 1b - half wire (y - z) extrude in x direction - Corner 1 to Corner 2
// Wire 1b1

p0_1b = newp; Point(p0_1b) = {-p,-p,r_w+mesh_level*mm, LcWiremesh * mm};              // centre circle
p1b_1_1 = newp; Point(p1b_1_1) = {-p,-p,0+mesh_level*mm, LcWiremesh * mm};            // bottom circle
p2_1b = newp; Point(p2_1b) = {-p,-p+r_w,r_w+mesh_level*mm, LcWiremesh * mm};          // right circle
p1b_3_1 = newp; Point(p1b_3_1) = {-p,-p,2*r_w+mesh_level*mm, LcWiremesh * mm};        // top circle
// p4_1b = newp; Point(p4_1b) = {-p,-p-r_w,r_w+mesh_level*mm, LcWiremesh * mm};       // left circle

l1_1b = newl; Circle(l1_1b) = {p1b_1_1, p0_1b, p2_1b};
l2_1b = newl; Circle(l2_1b) = {p2_1b, p0_1b, p1b_3_1};
l2_1bs = newl; Line(l2_1bs) = {p1b_1_1, p1b_3_1};

ll1_1b = newll; Line Loop(ll1_1b) = {l1_1b, l2_1b, -l2_1bs};

s_1_1b = news; Plane Surface(s_1_1b) = {ll1_1b};

sb_1_1[] = {};
tmpb_1_1a[] = {s_1_1b};

tmpb_1_1b[] = Extrude {{x1_sp_wind_fac,0,0},{0,1,0},{-p-sp_fac,-p,-Rtn+r_w}, alpha} {
  Surface{tmpb_1_1a[0]}; 
};

sb_1_1[] += tmpb_1_1b[{2:4}];


// Wire 1b2

sb_1_2[] = {};
tmpb_1_2a[] = {tmpb_1_1b[0]};

tmpb_1_2b[] = Extrude {{x1_sp_wind_fac,0,0},{0,-1,0},{p+sp_fac,-p,Rtp-r_w}, alpha} {
  Surface{tmpb_1_2a[0]}; 
};

sb_1_2[] += tmpb_1_2b[{2:4}];


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////// Face 2a - half wire (x - z) extrude in y direction - Corner 3 to Corner 2
// Wire 2a1

p0_2a = newp; Point(p0_2a) = {p,p,-r_w+mesh_level*mm, LcWiremesh * mm};               // centre circle
p2a_1_1 = newp; Point(p2a_1_1) = {p,p,-2*r_w+mesh_level*mm, LcWiremesh * mm};         // bottom circle
// p2_2a = newp; Point(p2_2a) = {p+r_w,p,-r_w+mesh_level*mm, LcWiremesh * mm};        // right circle
p2a_3_1 = newp; Point(p2a_3_1) = {p,p,0+mesh_level*mm, LcWiremesh * mm};              // top circle
p4_2a = newp; Point(p4_2a) = {p-r_w,p,-r_w+mesh_level*mm, LcWiremesh * mm};           // left circle

l2_2as = newl; Line(l2_2as) = {p2a_1_1, p2a_3_1};
l3_2a = newl; Circle(l3_2a) = {p2a_3_1, p0_2a, p4_2a};
l4_2a = newl; Circle(l4_2a) = {p4_2a, p0_2a, p2a_1_1};

ll1_2a = newll; Line Loop(ll1_2a) = {l3_2a, l4_2a, l2_2as};

s_1_2a = news; Plane Surface(s_1_2a) = {ll1_2a};

sa_2_1[] = {};
tmpa_2_1a[] = {s_1_2a};

tmpa_2_1b[] = Extrude {{0,y1_sp_wind_fac,0},{-1,0,0},{p,p+sp_fac,Rtp-r_w}, alpha} {
  Surface{tmpa_2_1a[0]};
};

sa_2_1[] += tmpa_2_1b[{2:4}];


// Wire 2a2

sa_2_2[] = {};
tmpa_2_2a[] = {tmpa_2_1b[0]};

tmpa_2_2b[] = Extrude {{0,y1_sp_wind_fac,0},{1,0,0},{p,-p-sp_fac,-Rtn+r_w}, alpha} {
  Surface{tmpa_2_2a[0]}; 
};

sa_2_2[] += tmpa_2_2b[{2:4}];


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////// Face 1a - half wire (y - z) extrude in x direction - Corner 3 to Corner 4
// Wire 1a1

p0_1a = newp; Point(p0_1a) = {p,p,r_w+mesh_level*mm, LcWiremesh * mm};                  // centre circle
p1a_1_1 = newp; Point(p1a_1_1) = {p,p,0+mesh_level*mm, LcWiremesh * mm};                // bottom circle
// p2_1a = newp; Point(p2_1a) = {p,p+r_w,r_w+mesh_level*mm, LcWiremesh * mm};           // right circle
p1a_3_1 = newp; Point(p1a_3_1) = {p,p,2*r_w+mesh_level*mm, LcWiremesh * mm};            // top circle
p4_1a = newp; Point(p4_1a) = {p,p-r_w,r_w+mesh_level*mm, LcWiremesh * mm};              // left circle

l2_1as = newl; Line(l2_1as) = {p1a_1_1, p1a_3_1};
l3_1a = newl; Circle(l3_1a) = {p1a_3_1, p0_1a, p4_1a};
l4_1a = newl; Circle(l4_1a) = {p4_1a, p0_1a, p1a_1_1};

ll1_1a = newll; Line Loop(ll1_1a) = {l3_1a, l4_1a, l2_1as};

s_1_1a = news; Plane Surface(s_1_1a) = {ll1_1a};

sa_1_1[] = {};
tmpa_1_1a[] = {};
tmpa_1_1a[] = {s_1_1a};

tmpa_1_1b[] = Extrude {{x2_sp_wind_fac,0,0},{0,-1,0},{p+sp_fac,p,-Rtn+r_w}, alpha} {
  Surface{tmpa_1_1a[0]};
};

sa_1_1[] += tmpa_1_1b[{2:4}];


// Wire 1a2

sa_1_2[] = {};
tmpa_1_2a[] = {tmpa_1_1b[0]};

tmpa_1_2b[] = Extrude {{x2_sp_wind_fac,0,0},{0,1,0},{-p-sp_fac,p,Rtp-r_w}, alpha} {
  Surface{tmpa_1_2a[0]};
};

sa_1_2[] += tmpa_1_2b[{2:4}];


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////// Face 2b - half wire (x - z) extrude in y direction - Corner 1 to Corner 4
// Wire 2b1

p0_2b = newp; Point(p0_2b) = {-p,-p,-r_w+mesh_level*mm, LcWiremesh * mm};                 // centre circle
p2b_1_1 = newp; Point(p2b_1_1) = {-p,-p,-2*r_w+mesh_level*mm, LcWiremesh * mm};           // bottom circle
p2_2b = newp; Point(p2_2b) = {-p+r_w,-p,-r_w+mesh_level*mm, LcWiremesh * mm};             // right circle
p2b_3_1 = newp; Point(p2b_3_1) = {-p,-p,0+mesh_level*mm, LcWiremesh * mm};                // top circle
// p4_2b = newp; Point(p4_2b) = {-p-r_w,-p,-r_w+mesh_level*mm, LcWiremesh * mm};          // left circle

l1_2b = newl; Circle(l1_2b) = {p2b_1_1, p0_2b, p2_2b};
l2_2b = newl; Circle(l2_2b) = {p2_2b, p0_2b, p2b_3_1};
l2_2bs = newl; Line(l2_2bs) = {p2b_1_1, p2b_3_1};

ll1_2b = newll; Line Loop(ll1_2b) = {l1_2b, l2_2b, -l2_2bs};

s_1_2b = news; Plane Surface(s_1_2b) = {ll1_2b};

sb_2_1[] = {};
tmpb_2_1a[] = {s_1_2b};

tmpb_2_1b[] = Extrude {{0,y2_sp_wind_fac,0},{1,0,0},{-p,-p-sp_fac,Rtp-r_w}, alpha} {
  Surface{tmpb_2_1a[0]};  
};

sb_2_1[] += tmpb_2_1b[{2:4}];

// Wire 2b2

sb_2_2[] = {};
tmpb_2_2a[] = {tmpb_2_1b[0]};

tmpb_2_2b[] = Extrude {{0,y2_sp_wind_fac,0},{-1,0,0},{-p,p+sp_fac,-Rtn+r_w}, alpha} {
  Surface{tmpb_2_2a[0]}; 
};

sb_2_2[] += tmpb_2_2b[{2:4}];


// SHELL

// --------------------------------------------------------------------------

// *******************************
// Corner 1
// *******************************
pc1_1 = newp; Point(pc1_1) = {geo_f_x*-mesh_window/2+geo_f_x*m*-mesh_window/2, geo_f_y*-mesh_window/2+geo_f_y*n*-mesh_window/2, 0,lcCopperPlateBdry};


// *******************************
// Corner 2
// *******************************
pc1_2 = newp; Point(pc1_2) = {geo_f_x*mesh_window/2+geo_f_x*m*mesh_window/2, geo_f_y*-mesh_window/2+geo_f_y*n*-mesh_window/2, 0,lcCopperPlateBdry};


// *******************************
// Corner 3
// *******************************
pc1_3 = newp; Point(pc1_3) = {geo_f_x*mesh_window/2+geo_f_x*m*mesh_window/2, geo_f_y*mesh_window/2+geo_f_y*n*mesh_window/2, 0,lcCopperPlateBdry};


// *******************************
// Corner 4
// *******************************
pc1_4 = newp; Point(pc1_4) = {geo_f_x*-mesh_window/2+geo_f_x*m*-mesh_window/2, geo_f_y*mesh_window/2+geo_f_y*n*mesh_window/2, 0,lcCopperPlateBdry};


// UPPER SQUARE

// *******************************
// Corner 1
// *******************************
ptpR1_1 = newp; Point(ptpR1_1) = {-p-sp_fac, -p-sp_fac, R-r_w, lcCopperPlateBdry};


// *******************************
// Corner 2
// *******************************
ptpR1_2 = newp; Point(ptpR1_2) = {p+sp_fac, -p-sp_fac, R-r_w, lcCopperPlateBdry};


// *******************************
// Corner 3
// *******************************
ptpR1_3 = newp; Point(ptpR1_3) = {-p-sp_fac, p+sp_fac, R-r_w, lcCopperPlateBdry};


// *******************************
// Corner 4
// *******************************
ptpR1_4 = newp; Point(ptpR1_4) = {p+sp_fac, p+sp_fac, R-r_w, lcCopperPlateBdry};

// UPPER SQUARE

// *******************************
// Corner 1
// *******************************
ptnR1_1 = newp; Point(ptnR1_1) = {-p-sp_fac, -p-sp_fac, -R+r_w, lcCopperPlateBdry};


// *******************************
// Corner 2
// *******************************
ptnR1_2 = newp; Point(ptnR1_2) = {p+sp_fac, -p-sp_fac, -R+r_w, lcCopperPlateBdry};


// *******************************
// Corner 3
// *******************************
ptnR1_3 = newp; Point(ptnR1_3) = {-p-sp_fac, p+sp_fac, -R+r_w, lcCopperPlateBdry};


// *******************************
// Corner 4
// *******************************
ptnR1_4 = newp; Point(ptnR1_4) = {p+sp_fac, p+sp_fac, -R+r_w, lcCopperPlateBdry};
