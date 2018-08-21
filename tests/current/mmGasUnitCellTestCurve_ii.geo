a = 0.040;                                        // the "pitch", or distance between GEM pillars, in mm

mwf = 1;                                          // mesh_window_factor

mm = 1;                                           // geometrical scaling
r_w = 0.009 * mm;                                 // R of Wiremesh, in microns
hp_0 = 0.0225;                                    // half pitch of the window, in mm
hp = 0.0225 * mm - 0*r_w/mwf * mm;                // half pitch of the window, in microns

p = hp_0;                                         // half pitch of the window, in mm

R = ( p * p + r_w * r_w )/( (2 * r_w^2) );          // R
alpha = Asin( (p/R) );                            // angle in radians

mesh_level = 0.000;                               // mesh level, in mm
mesh_window = 0.020;                              // mesh window, in mm

x1_sp_wind_fac = 1.00;
x2_sp_wind_fac = 1.00;
y1_sp_wind_fac = 1.00;
y2_sp_wind_fac = 1.00;

sp_fac1 = p*0.00;

x1_sp_wind_fac2 = p*0.00;
x2_sp_wind_fac2 = p*0.00;
y1_sp_wind_fac2 = p*0.00;
y2_sp_wind_fac2 = p*0.00;

Rtp = R - R*0.00;
Rtn = R + R*0.00;

h_f = 0*r_w;                                      // Heightfactor

geo_f_x = 1;
geo_f_y = 1;

m = 0;
n = 0;

// Characteristic lengths

  lcCopperPlateBdry = 0.005;
  lcExtElectrodeBdry = 0.005;
  LcWiremesh = 0.001;


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////// Face 2b - half wire (x - z) extrude in y direction - Corner 1 to Corner 4
// Wire 2b1

p0_2b = newp; Point(p0_2b) = {-p,-p,-r_w+mesh_level*mm, LcWiremesh * mm};                                               // centre circle
p2b_1_1 = newp; Point(p2b_1_1) = {-p,-p,-2*r_w+mesh_level*mm, LcWiremesh * mm};                                         // bottom circle
p2_2b = newp; Point(p2_2b) = {-p+r_w,-p,-r_w+mesh_level*mm, LcWiremesh * mm};                                           // right circle
p2b_3_1 = newp; Point(p2b_3_1) = {-p,-p,0+mesh_level*mm, LcWiremesh * mm};                                              // top circle
// p4_2b = newp; Point(p4_2b) = {-p-r_w,-p,-r_w+mesh_level*mm, LcWiremesh * mm};                                        // left circle

l1_2b = newl; Circle(l1_2b) = {p2b_1_1, p0_2b, p2_2b};
l2_2b = newl; Circle(l2_2b) = {p2_2b, p0_2b, p2b_3_1};
l2_2bs = newl; Line(l2_2bs) = {p2b_1_1, p2b_3_1};

ll1_2b = newll; Line Loop(ll1_2b) = {l1_2b, l2_2b, -l2_2bs};

s_1_2b = news; Plane Surface(s_1_2b) = {ll1_2b};

sb_2_1[] = {};
tmpb_2_1a[] = {s_1_2b};

/*
tmpb_2_1b[] = Extrude {0,y2_sp_wind_fac,0} {
  Surface{tmpb_2_1a[0]};
};

sb_2_1b[] = tmpb_2_1b[{2:4}];
*/

tmpb_2_1c[] = Extrude {{0,y2_sp_wind_fac2,0},{1,0,0},{-p,-p+1*sp_fac1,Rtp-r_w}, alpha} {
  Surface{tmpb_2_1a[0]};
};

sb_2_1c[] = tmpb_2_1c[{2:4}];

// Wire 2b2

sb_2_2[] = {};
tmpb_2_2a[] = {tmpb_2_1c[0]};

tmpb_2_2b[] = Extrude {{0,y2_sp_wind_fac2,0},{-1,0,0},{-p,p+1*sp_fac1,-Rtn+r_w}, alpha} {
  Surface{tmpb_2_2a[0]};
};

sb_2_2b[] = tmpb_2_2b[{2:4}];

/*
tmpb_2_2c[] = Extrude {0,y2_sp_wind_fac,0} {
  Surface{tmpb_2_2b[0]};
};

sb_2_2c[] = tmpb_2_2c[{2:4}];
*/

s_1_2b1[] = s_1_2b;
s_1_2b2[] = tmpb_2_2b[0];

sl_wire_exterior_surface_2b[] = newreg; Surface Loop(sl_wire_exterior_surface_2b) = { s_1_2b1[0], sb_2_1c[], sb_2_2b[], s_1_2b2[0] }; // sb_2_1b[], sb_2_2c[]
// vol_2b_wire1 = newreg; Volume(vol_2b_wire1) = { tmpb_2_1c[] }; // tmpb_2_1b[1], tmpb_2_2c[1]
// vol_2b_wire2 = newreg; Volume(vol_2b_wire2) = { tmpb_2_2b[] }; // tmpb_2_1b[1], tmpb_2_2c[1]
vol_2b_wire = newreg; Volume(vol_2b_wire) = sl_wire_exterior_surface_2b[];

physvol_2b_wire = newreg; Physical Volume(physvol_2b_wire) = vol_2b_wire;
// physvol_2b_wire1 = newreg; Physical Volume(physvol_2b_wire1) = vol_2b_wire1;
// physvol_2b_wire2 = newreg; Physical Volume(physvol_2b_wire2) = vol_2b_wire2;
physsurf_2b_wire = newreg; Physical Surface(physsurf_2b_wire) = { s_1_2b1[], sb_2_1c[], sb_2_2b[], s_1_2b2[] }; // sb_2_1b[], sb_2_2c[]


/*
// SHELL

// --------------------------------------------------------------------------

// *******************************
// Corner 1
// *******************************
pc1_1 = newp; Point(pc1_1) = {geo_f_x*-mesh_window/2+geo_f_x*m*-mesh_window/2+1*sp_fac1, geo_f_y*-mesh_window/2+geo_f_y*n*-mesh_window/2+1*sp_fac1, 0,lcCopperPlateBdry};


// *******************************
// Corner 2
// *******************************
pc1_2 = newp; Point(pc1_2) = {geo_f_x*mesh_window/2+geo_f_x*m*mesh_window/2+1*sp_fac1, geo_f_y*-mesh_window/2+geo_f_y*n*-mesh_window/2+1*sp_fac1, 0,lcCopperPlateBdry};


// *******************************
// Corner 3
// *******************************
pc1_3 = newp; Point(pc1_3) = {geo_f_x*mesh_window/2+geo_f_x*m*mesh_window/2+1*sp_fac1, geo_f_y*mesh_window/2+geo_f_y*n*mesh_window/2+1*sp_fac1, 0,lcCopperPlateBdry};


// *******************************
// Corner 4
// *******************************
pc1_4 = newp; Point(pc1_4) = {geo_f_x*-mesh_window/2+geo_f_x*m*-mesh_window/2+1*sp_fac1, geo_f_y*mesh_window/2+geo_f_y*n*mesh_window/2+1*sp_fac1, 0,lcCopperPlateBdry};


// UPPER SQUARE

// *******************************
// Corner 1
// *******************************
ptR1_0 = newp; Point(ptR1_0) = {-p+1*sp_fac1, -p+1*sp_fac1, R-r_w, lcCopperPlateBdry};


// *******************************
// Corner 1
// *******************************
ptpR1_1 = newp; Point(ptpR1_1) = {-p, -p+1*sp_fac1, R-r_w, lcCopperPlateBdry};


// *******************************
// Corner 2
// *******************************
ptpR1_2 = newp; Point(ptpR1_2) = {p+1*sp_fac1, -p, R-r_w, lcCopperPlateBdry};


// *******************************
// Corner 3
// *******************************
ptpR1_3 = newp; Point(ptpR1_3) = {p+1*sp_fac1, p+1*sp_fac1, R-r_w, lcCopperPlateBdry};


// *******************************
// Corner 4
// *******************************
ptpR1_4 = newp; Point(ptpR1_4) = {-p+1*sp_fac1, p+1*sp_fac1, R-r_w, lcCopperPlateBdry};

// UPPER SQUARE

// *******************************
// Corner 1
// *******************************
ptnR1_1 = newp; Point(ptnR1_1) = {-p+1*sp_fac1, -p, -R+r_w, lcCopperPlateBdry};


// *******************************
// Corner 2
// *******************************
ptnR1_2 = newp; Point(ptnR1_2) = {p+1*sp_fac1, -p+1*sp_fac1, -R+r_w, lcCopperPlateBdry};


// *******************************
// Corner 3
// *******************************
ptnR1_3 = newp; Point(ptnR1_3) = {p+1*sp_fac1, p+1*sp_fac1, -R+r_w, lcCopperPlateBdry};


// *******************************
// Corner 4
// *******************************
ptnR1_4 = newp; Point(ptnR1_4) = {-p, p+1*sp_fac1, -R+r_w, lcCopperPlateBdry};
*/