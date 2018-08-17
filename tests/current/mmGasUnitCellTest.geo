// *********************************************************************
// mmGasRotationTest.geo
//
// Description:
// Geometry file for a MM cell.
// This cell can be repeated any number of times within Garfield 
// to construct an arbitrarily large MM
//
// *********************************************************************  


tC = 0.0035;                                      // copper thickness, in mm
tD = 0.0035;                                      // dielectric thickness, in mm
lE = 0.40;                                        // distance from GEM plates to upper exterior electrode, in mm
lP = 0.02;                                        // distance from lower LEM plate to pad (readout) plane, in mm
a = 0.040;                                        // the "pitch", or distance between GEM pillars, in mm

mwf = 1;                                          // mesh_window_factor

mm = 1;                                           // geometrical scaling
r_w = 0.010 * mm;                                 // R of Wiremesh, in microns
hp_0 = 0.010;                                     // half pitch of the window, in mm
hp = 0.010 * mm - 0*r_w/mwf * mm;                 // half pitch of the window, in microns

p = hp_0;                                         // half pitch of the window, in mm

R = (p * p + r_w * r_w)/( (2 * r_w) );            // R
alpha = Asin((p/R));                              // angle in radians

mesh_level = 0.00;                               // mesh level, in mm
mesh_window = 0.020;                              // mesh window, in mm

sp_fac_r1 = 0.725;
sp_fac_r2 = 0.725;

x1_sp_wind_fac = p*0.25;
x2_sp_wind_fac = p*0.25;
y1_sp_wind_fac = p*0.25;
y2_sp_wind_fac = p*0.25;

sp_fac1 = p*1.000;
sp_fac2 = p*1.000;

sp_fac_str1 = p*1.00;
sp_fac_str2 = p*1.00;
sp_fac_str3 = p*0.00;
sp_fac_str4 = p*0.00;

x1_sp_wind_fac2 = p*0.00;
x2_sp_wind_fac2 = p*0.00;
y1_sp_wind_fac2 = p*0.00;
y2_sp_wind_fac2 = p*0.00;

Rtp = R + R*0.0;
Rtn = R - R*0.0;

h_f = 0*r_w;                                      // Heightfactor

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


p1b1_0 = newp; Point(p1b1_0) = {-p+p,-p+p,r_w+mesh_level*mm, LcWiremesh * mm};                                           // centre circle
p1b1_1 = newp; Point(p1b1_1) = {-p+p,-p+p,0+mesh_level*mm, LcWiremesh * mm};                                             // bottom circle
p1b1_2 = newp; Point(p1b1_2) = {-p+p,-p+p+r_w,r_w+mesh_level*mm, LcWiremesh * mm};                                       // right circle
p1b1_3 = newp; Point(p1b1_3) = {-p+p,-p+p,2*r_w+mesh_level*mm, LcWiremesh * mm};                                         // top circle
// p1b1_4 = newp; Point(p1b1_4) = {-p+p,-p+p-r_w,r_w+mesh_level*mm, LcWiremesh * mm};                                    // left circle

l1b1_1 = newl; Circle(l1b1_1) = {p1b1_1, p1b1_0, p1b1_2};
l1b1_2 = newl; Circle(l1b1_2) = {p1b1_2, p1b1_0, p1b1_3};
l1b1s_3 = newl; Line(l1b1s_3) = {p1b1_1, p1b1_3};

ll1b1_1 = newll; Line Loop(ll1b1_1) = {l1b1_1, l1b1_2, -l1b1s_3};

s1b1_1 = news; Plane Surface(s1b1_1) = {ll1b1_1};

tmpba_1_1[] = {s1b1_1};
sbaf_1_1_1[] = tmpba_1_1[0];


tmpbb_1_1[] = Extrude {1.325*x1_sp_wind_fac +0.00*sp_fac_str1,0,0} {
  Surface{tmpba_1_1[0]};
};

sbb_1_1[] = {};
sbb_1_1[] += tmpbb_1_1[{2:4}];


ptsbc_1_1 = newp; Point(ptsbc_1_1) = { -p+p +0.00*p +1.75*x1_sp_wind_fac +0.00*R*Cos(alpha/2) +0.00*sp_fac1 +0.00*sp_fac_str1, -p+p +0.00*p +0.00*sp_fac1, Rtn-r_w +1.00*mesh_level*mm -1.175*sp_fac_r1*r_w +0.00*R*Sin(alpha/2) +0.00*x1_sp_wind_fac +0.00*sp_fac_str2 };



tmpbc_1_1[] = Extrude {{x1_sp_wind_fac2,0,0}, {0,1,0}, { -p+p +0.00*p +1.75*x1_sp_wind_fac +0.00*R*Cos(alpha/2) +0.00*sp_fac1 +0.00*sp_fac_str1, -p+p +0.00*p +0.00*sp_fac1, Rtn-r_w +1.00*mesh_level*mm -1.175*sp_fac_r1*r_w +0.00*R*Sin(alpha/2) +0.00*x1_sp_wind_fac +0.00*sp_fac_str2 }, 0.20*alpha } {
  Surface{tmpbb_1_1[0]};
};


sbc_1_1[] = {};
sbc_1_1[] += tmpbc_1_1[{2:4}];



// tmpbc2_1_1[] = Extrude {0,0,-1.00*x1_sp_wind_fac +0.00*sp_fac_str2} {
//  Surface{tmpbd_1_1[0]};
// };
//
// sbc2_1_1[] = {};
// sbc2_1_1[] += tmpbc2_1_1[{2:4}];


ptsbd_1_1 = newp; Point(ptsbd_1_1) = { -p+p +0.00*p +2.10*x1_sp_wind_fac +0.00*R*Cos(alpha/2) +0.00*sp_fac1 +0.00*sp_fac_str1, -p+p +0.00*p +0.00*sp_fac1, Rtp-r_w +1.00*mesh_level*mm -0.825*sp_fac_r1*r_w +0.00*R*Sin(alpha/2) +0.00*x1_sp_wind_fac +0.00*sp_fac_str2 };



tmpbd_1_1[] = Extrude {{x1_sp_wind_fac2,0,0}, {0,1,0}, { -p+p +0.00*p +2.10*x1_sp_wind_fac +0.00*R*Cos(alpha/2) +0.00*sp_fac1 +0.00*sp_fac_str1, -p+p +0.00*p +0.00*sp_fac1, Rtp-r_w +1.00*mesh_level*mm -0.825*sp_fac_r1*r_w +0.00*R*Sin(alpha/2) +0.00*x1_sp_wind_fac +0.00*sp_fac_str2 }, 0.60*alpha} {
  Surface{tmpbc_1_1[0]};
};


sbd_1_1[] = {};
sbd_1_1[] += tmpbd_1_1[{2:4}];



// Wire 1a2


// tmpba_1_2[] = {tmpbd_1_1[0]};
// 
// tmpbb_1_2[] = Extrude {1.00*x1_sp_wind_fac +0.00*sp_fac_str1,0,0} {
//   Surface{tmpba_1_2[0]};
// };
// 
// sbb_1_2[] = {};
// sbb_1_2[] += tmpbb_1_2[{2:4}];


ptsbc_1_2 = newp; Point(ptsbc_1_2) = { p+p +2.00*p -2.10*x1_sp_wind_fac +0.00*R*Cos(alpha/2) +0.00*sp_fac1 +0.00*sp_fac_str1, -p+p +0.00*p +0.00*sp_fac1, -Rtp+r_w +1.00*mesh_level*mm +0.825*sp_fac_r1*r_w +0.00*R*Sin(alpha/2) +0.00*x1_sp_wind_fac +0.00*sp_fac_str2 };



tmpbc_1_2[] = Extrude {{x1_sp_wind_fac2,0,0}, {0,-1,0}, { p+p +2.00*p -2.10*x1_sp_wind_fac +0.00*R*Cos(alpha/2) +0.00*sp_fac1 +0.00*sp_fac_str1, -p+p +0.00*p +0.00*sp_fac1, -Rtp+r_w +1.00*mesh_level*mm +0.825*sp_fac_r1*r_w +0.00*R*Sin(alpha/2) +0.00*x1_sp_wind_fac +0.00*sp_fac_str2 }, 0.60*alpha} {
  Surface{tmpbd_1_1[0]};
};


sbc_1_2[] = {};
sbc_1_2[] += tmpbc_1_2[{2:4}];



// tmpbc2_1_2[] = Extrude {0,0,1.00*x1_sp_wind_fac +0.00*sp_fac_str2} {
//   Surface{tmpbc_1_2[0]};
// };
// 
// sbc2_1_2[] = {};
// sbc2_1_2[] += tmpbc2_1_2[{2:4}];


ptsbd_1_2 = newp; Point(ptsbd_1_2) = { p+p +2.00*p -1.75*x1_sp_wind_fac +0.00*R*Cos(alpha/2) +0.00*sp_fac1 +0.00*sp_fac_str1, -p+p +0.00*p +0.00*sp_fac1, -Rtp+r_w +1.00*mesh_level*mm +1.175*sp_fac_r1*r_w +0.00*R*Sin(alpha/2) +0.00*x1_sp_wind_fac +0.00*sp_fac_str2 };



tmpbd_1_2[] = Extrude {{x1_sp_wind_fac2,0,0}, {0,-1,0},{ p+p +2.00*p -1.75*x1_sp_wind_fac +0.00*R*Cos(alpha/2) +0.00*sp_fac1 +0.00*sp_fac_str1, -p+p +0.00*p +0.00*sp_fac1, -Rtp+r_w +1.00*mesh_level*mm +1.175*sp_fac_r1*r_w +0.00*R*Sin(alpha/2) +0.00*x1_sp_wind_fac +0.00*sp_fac_str2 }, 0.20*alpha} {
  Surface{tmpbc_1_2[0]};
};


sbd_1_2[] = {};
sbd_1_2[] += tmpbd_1_2[{2:4}];



tmpbe_1_2[] = Extrude {1.325*x1_sp_wind_fac,0,0} {
  Surface{tmpbd_1_2[0]};
};

sbe_1_2[] = {};
sbe_1_2[] += tmpbe_1_2[{2:4}];

sbef_2_1_2[] = {};
sbef_2_1_2[] = {tmpbe_1_2[0]};



/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////// Face 2a - half wire (x - z) extrude in y direction - Corner 3 to Corner 2
// Wire 2a1

p2a1_0 = newp; Point(p2a1_0) = {p+p+2*sp_fac1,p+p+2*sp_fac1,-r_w+mesh_level*mm, LcWiremesh * mm};                        // centre circle
p2a1_1 = newp; Point(p2a1_1) = {p+p+2*sp_fac1,p+p+2*sp_fac1,-2*r_w+mesh_level*mm, LcWiremesh * mm};                      // bottom circle
// p2a1_2 = newp; Point(p2a1_2) = {p+p+2*sp_fac1+r_w,p+p+2*sp_fac1,-r_w+mesh_level*mm, LcWiremesh * mm};                 // right circle
p2a1_3 = newp; Point(p2a1_3) = {p+p+2*sp_fac1,p+p+2*sp_fac1,0+mesh_level*mm, LcWiremesh * mm};                           // top circle
p2a1_4 = newp; Point(p2a1_4) = {p+p+2*sp_fac1-r_w,p+p+2*sp_fac1,-r_w+mesh_level*mm, LcWiremesh * mm};                    // left circle

l2a1s_1 = newl; Line(l2a1s_1) = {p2a1_1, p2a1_3};
l2a1_2 = newl; Circle(l2a1_2) = {p2a1_3, p2a1_0, p2a1_4};
l2a1_3 = newl; Circle(l2a1_3) = {p2a1_4, p2a1_0, p2a1_1};

ll2a1_1 = newll; Line Loop(ll2a1_1) = {l2a1_2, l2a1_3, l2a1s_1};

s2a1_1 = news; Plane Surface(s2a1_1) = {ll2a1_1};

tmpaa_2_1[] = {s2a1_1};
saaf_2_2_1[] = tmpaa_2_1[0];

tmpab_2_1[] = Extrude {0,-1.325*y1_sp_wind_fac +0.00*sp_fac_str1,0} {
  Surface{tmpaa_2_1[0]};
};

sab_2_1[] = {};
sab_2_1[] += tmpab_2_1[{2:4}];


ptsac_2_1 = newp; Point(ptsac_2_1) = { p+p +2.00*p +0.00*sp_fac1, p+p +2.00*p -1.75*y1_sp_wind_fac +0.00*R*Cos(alpha/2) +0.00*sp_fac1 +0.00*sp_fac_str1, -Rtn+r_w +1.00*mesh_level*mm +1.175*sp_fac_r1*r_w +0.00*R*Sin(alpha/2) +0.00*y1_sp_wind_fac +0.00*sp_fac_str2 };



tmpac_2_1[] = Extrude {{y1_sp_wind_fac2,0,0}, {-1,0,0}, { p+p +2.00*p +0.00*sp_fac1, p+p +2.00*p -1.75*y1_sp_wind_fac +0.00*R*Cos(alpha/2) +0.00*sp_fac1 +0.00*sp_fac_str1, -Rtn+r_w +1.00*mesh_level*mm +1.175*sp_fac_r1*r_w +0.00*R*Sin(alpha/2) +0.00*y1_sp_wind_fac +0.00*sp_fac_str2 }, 0.20*alpha } {
  Surface{tmpab_2_1[0]};
};

sac_2_1[] = {};
sac_2_1[] += tmpac_2_1[{2:4}];



// tmpac2_2_1[] = Extrude {0,0,-1.00*y1_sp_wind_fac +0.00*sp_fac_str2} {
//   Surface{tmpad_2_1[0]};
// };
// 
// sac2_2_1[] = {};
// sac2_2_1[] += tmpac2_2_1[{2:4}];


ptsad_2_1 = newp; Point(ptsad_2_1) = { p+p +2.00*p +0.00*sp_fac1, p+p +2.00*p -2.10*y1_sp_wind_fac +0.00*R*Cos(alpha/2) +0.00*sp_fac1 +0.00*sp_fac_str1, -Rtp+r_w +1.00*mesh_level*mm +0.825*sp_fac_r1*r_w +0.00*R*Sin(alpha/2) +0.00*y1_sp_wind_fac +0.00*sp_fac_str2 };



tmpad_2_1[] = Extrude {{y1_sp_wind_fac2,0,0}, {-1,0,0}, { p+p +2.00*p +0.00*sp_fac1, p+p +2.00*p -2.10*y1_sp_wind_fac +0.00*R*Cos(alpha/2) +0.00*sp_fac1 +0.00*sp_fac_str1, -Rtp+r_w +1.00*mesh_level*mm +0.825*sp_fac_r1*r_w +0.00*R*Sin(alpha/2) +0.00*y1_sp_wind_fac +0.00*sp_fac_str2 }, 0.60*alpha} {
  Surface{tmpac_2_1[0]};
};

sad_2_1[] = {};
sad_2_1[] += tmpad_2_1[{2:4}];



// Wire 2a2


// tmpaa_2_2[] = {tmpad_2_1[0]};
// 
// tmpab_2_2[] = Extrude {0,1.00*y1_sp_wind_fac +0.00*sp_fac_str1,0} {
//   Surface{tmpaa_2_2[0]};
// };
// 
// sab_2_2[] = {};
// sab_2_2[] += tmpab_2_2[{2:4}];


ptsac_2_2 = newp; Point(ptsac_2_2) = { p+p +2.00*p +0.00*sp_fac1, -p+p +2.10*y1_sp_wind_fac +0.00*R*Cos(alpha/2) +0.00*sp_fac1 +0.00*sp_fac_str1, Rtp-r_w +1.00*mesh_level*mm -0.825*sp_fac_r1*r_w +0.00*R*Sin(alpha/2) +0.00*y1_sp_wind_fac +0.00*sp_fac_str2 };



tmpac_2_2[] = Extrude {{x1_sp_wind_fac2,0,0}, {1,0,0}, { p+p +2.00*p +0.00*sp_fac1, -p+p +2.10*y1_sp_wind_fac +0.00*R*Cos(alpha/2) +0.00*sp_fac1 +0.00*sp_fac_str1, Rtp-r_w +1.00*mesh_level*mm -0.825*sp_fac_r1*r_w +0.00*R*Sin(alpha/2) +0.00*y1_sp_wind_fac +0.00*sp_fac_str2 }, 0.60*alpha} {
  Surface{tmpad_2_1[0]};
};

sac_2_2[] = {};
sac_2_2[] += tmpac_2_2[{2:4}];



// tmpac2_2_2[] = Extrude {0,0,1.00*y1_sp_wind_fac +0.00*sp_fac_str2} {
//   Surface{tmpac_2_2[0]};
// };
// 
// sac2_2_2[] = {};
// sac2_2_2[] += tmpac2_2_2[{2:4}];


ptsad_2_2 = newp; Point(ptsad_2_2) = { p+p +2.00*p +0.00*sp_fac1, -p+p +1.75*y1_sp_wind_fac +0.00*R*Cos(alpha/2) +0.00*sp_fac1 +0.00*sp_fac_str1, Rtp-r_w +1.00*mesh_level*mm -1.175*sp_fac_r1*r_w +0.00*R*Sin(alpha/2) +0.00*y1_sp_wind_fac +0.00*sp_fac_str2 };



tmpad_2_2[] = Extrude {{y1_sp_wind_fac2,0,0}, {1,0,0},{ p+p +2.00*p +0.00*sp_fac1, -p+p +1.75*y1_sp_wind_fac +0.00*R*Cos(alpha/2) +0.00*sp_fac1 +0.00*sp_fac_str1, Rtp-r_w +1.00*mesh_level*mm -1.175*sp_fac_r1*r_w +0.00*R*Sin(alpha/2) +0.00*y1_sp_wind_fac +0.00*sp_fac_str2 }, 0.20*alpha} {
  Surface{tmpac_2_2[0]};
};

sad_2_2[] = {};
sad_2_2[] += tmpad_2_2[{2:4}];



tmpae_2_2[] = Extrude {0,-1.325*y1_sp_wind_fac +0.00*sp_fac_str2,0} {
  Surface{tmpad_2_2[0]};
};

sae_2_2[] = {};
sae_2_2[] += tmpae_2_2[{2:4}];

saef_2_2_2[] = {};
saef_2_2_2[] = {tmpae_2_2[0]};



/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////// Face 1a - half wire (y - z) extrude in x direction - Corner 3 to Corner 4
// Wire 1a1

p1a1_0 = newp; Point(p1a1_0) = {p+p+2*sp_fac1,p+p+2*sp_fac1,r_w+mesh_level*mm, LcWiremesh * mm};                         // centre circle
p1a1_1 = newp; Point(p1a1_1) = {p+p+2*sp_fac1,p+p+2*sp_fac1,0+mesh_level*mm, LcWiremesh * mm};                           // bottom circle
// p1a1_2 = newp; Point(p1a1_2) = {p+p+2*sp_fac1,p+p+2*sp_fac1+r_w,r_w+mesh_level*mm, LcWiremesh * mm};                  // right circle
p1a1_3 = newp; Point(p1a1_3) = {p+p+2*sp_fac1,p+p+2*sp_fac1,2*r_w+mesh_level*mm, LcWiremesh * mm};                       // top circle
p1a1_4 = newp; Point(p1a1_4) = {p+p+2*sp_fac1,p+p+2*sp_fac1-r_w,r_w+mesh_level*mm, LcWiremesh * mm};                     // left circle

l1a1s_1 = newl; Line(l1a1s_1) = {p1a1_1, p1a1_3};
l1a1_2 = newl; Circle(l1a1_2) = {p1a1_3, p1a1_0, p1a1_4};
l1a1_3 = newl; Circle(l1a1_3) = {p1a1_4, p1a1_0, p1a1_1};

ll1a1_1 = newll; Line Loop(ll1a1_1) = {l1a1_2, l1a1_3, l1a1s_1};

s1a1_1 = news; Plane Surface(s1a1_1) = {ll1a1_1};

tmpaa_1_1[] = {s1a1_1};
saaf_1_1_1[] = tmpaa_1_1[0];

tmpab_1_1[] = Extrude {-1.325*x2_sp_wind_fac +0.00*sp_fac_str1,0,0} {
  Surface{tmpaa_1_1[0]};
};

sab_1_1[] = {};
sab_1_1[] += tmpab_1_1[{2:4}];


ptsac_1_1 = newp; Point(ptsac_1_1) = { p+p +2.00*p -1.75*x2_sp_wind_fac +0.00*R*Cos(alpha/2) +0.00*sp_fac1 +0.00*sp_fac_str1, p+p +2.00*p +0.00*sp_fac1, Rtn-r_w +1.00*mesh_level*mm -1.175*sp_fac_r2*r_w +0.00*R*Sin(alpha/2) +0.00*x2_sp_wind_fac +0.00*sp_fac_str2 };



tmpac_1_1[] = Extrude {{x2_sp_wind_fac2,0,0}, {0,-1,0}, { p+p +2.00*p -1.75*x2_sp_wind_fac +0.00*R*Cos(alpha/2) +0.00*sp_fac1 +0.00*sp_fac_str1, p+p +2.00*p +0.00*sp_fac1, Rtn-r_w +1.00*mesh_level*mm -1.175*sp_fac_r2*r_w +0.00*R*Sin(alpha/2) +0.00*x2_sp_wind_fac +0.00*sp_fac_str2 }, 0.20*alpha } {
  Surface{tmpab_1_1[0]};
};



sac_1_1[] = {};
sac_1_1[] += tmpac_1_1[{2:4}];

// tmpac2_1_1[] = Extrude {0,0,-1.00*x2_sp_wind_fac +0.00*sp_fac_str2} {
//   Surface{tmpad_1_1[0]};
// };
// 
// sac2_1_1[] = {};
// sac2_1_1[] += tmpac2_1_1[{2:4}];


ptsad_1_1 = newp; Point(ptsad_1_1) = { p+p +2.00*p -2.10*x2_sp_wind_fac +0.00*R*Cos(alpha/2) +0.00*sp_fac1 +0.00*sp_fac_str1, p+p +2.00*p +0.00*sp_fac1, Rtp-r_w +1.00*mesh_level*mm -0.825*sp_fac_r2*r_w +0.00*R*Sin(alpha/2) +0.00*x2_sp_wind_fac +0.00*sp_fac_str2 };



tmpad_1_1[] = Extrude {{x2_sp_wind_fac2,0,0}, {0,-1,0}, { p+p +2.00*p -2.10*x2_sp_wind_fac +0.00*R*Cos(alpha/2) +0.00*sp_fac1 +0.00*sp_fac_str1, p+p +2.00*p +0.00*sp_fac1, Rtp-r_w +1.00*mesh_level*mm -0.825*sp_fac_r2*r_w +0.00*R*Sin(alpha/2) +0.00*x2_sp_wind_fac +0.00*sp_fac_str2 }, 0.60*alpha} {
  Surface{tmpac_1_1[0]};
};



sad_1_1[] = {};
sad_1_1[] += tmpad_1_1[{2:4}];


// Wire 1a2


// tmpaa_1_2[] = {tmpad_1_1[0]};
// 
// tmpab_1_2[] = Extrude {1.00*x1_sp_wind_fac +0.00*sp_fac_str1,0,0} {
//   Surface{tmpaa_1_2[0]};
// };
// 
// sab_1_2[] = {};
// sab_1_2[] += tmpab_1_2[{2:4}];


ptsac_1_2 = newp; Point(ptsac_1_2) = { -p+p +2.10*x2_sp_wind_fac +0.00*R*Cos(alpha/2) +0.00*sp_fac1 +0.00*sp_fac_str1, p+p +2.00*p +0.00*sp_fac1, -Rtp+r_w +1.00*mesh_level*mm +0.825*sp_fac_r2*r_w +0.00*R*Sin(alpha/2) +0.00*x2_sp_wind_fac +0.00*sp_fac_str2 };



tmpac_1_2[] = Extrude {{x1_sp_wind_fac2,0,0}, {0,1,0}, { -p+p +2.10*x2_sp_wind_fac +0.00*R*Cos(alpha/2) +0.00*sp_fac1 +0.00*sp_fac_str1, p+p +2.00*p +0.00*sp_fac1, -Rtp+r_w +1.00*mesh_level*mm +0.825*sp_fac_r2*r_w +0.00*R*Sin(alpha/2) +0.00*x2_sp_wind_fac +0.00*sp_fac_str2 }, 0.60*alpha} {
  Surface{tmpad_1_1[0]};
};



sac_1_2[] = {};
sac_1_2[] += tmpac_1_2[{2:4}];


// tmpac2_1_2[] = Extrude {0,0,1.00*x1_sp_wind_fac +0.00*sp_fac_str2} {
//   Surface{tmpac_1_2[0]};
// };
// 
// sac2_1_2[] = {};
// sac2_1_2[] += tmpac2_1_2[{2:4}];


ptsad_1_2 = newp; Point(ptsad_1_2) = { -p+p +1.75*x2_sp_wind_fac +0.00*R*Cos(alpha/2) +0.00*sp_fac1 +0.00*sp_fac_str1, p+p +2.00*p +0.00*sp_fac1, -Rtp+r_w +1.00*mesh_level*mm +1.175*sp_fac_r2*r_w +0.00*R*Sin(alpha/2) +0.00*x2_sp_wind_fac +0.00*sp_fac_str2 };



tmpad_1_2[] = Extrude {{x2_sp_wind_fac2,0,0}, {0,1,0},{ -p+p +1.75*x2_sp_wind_fac +0.00*R*Cos(alpha/2) +0.00*sp_fac1 +0.00*sp_fac_str1, p+p +2.00*p +0.00*sp_fac1, -Rtp+r_w +1.00*mesh_level*mm +1.175*sp_fac_r2*r_w +0.00*R*Sin(alpha/2) +0.00*x2_sp_wind_fac +0.00*sp_fac_str2 }, 0.20*alpha} {
  Surface{tmpac_1_2[0]};
};



sad_1_2[] = {};
sad_1_2[] += tmpad_1_2[{2:4}];

tmpae_1_2[] = Extrude {-1.325*x2_sp_wind_fac +0.00*sp_fac_str2,0,0} {
  Surface{tmpad_1_2[0]};
};

sae_1_2[] = {};
sae_1_2[] += tmpae_1_2[{2:4}];

saef_1_1_2[] = {};
saef_1_1_2[] = {tmpae_1_2[0]};



/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////// Face 2b - half wire (x - z) extrude in y direction - Corner 1 to Corner 4
// Wire 2b1

p2b1_0 = newp; Point(p2b1_0) = {-p+p,-p+p,-r_w+mesh_level*mm, LcWiremesh * mm};                                          // centre circle
p2b1_1 = newp; Point(p2b1_1) = {-p+p,-p+p,-2*r_w+mesh_level*mm, LcWiremesh * mm};                                        // bottom circle
p2b1_2 = newp; Point(p2b1_2) = {-p+p+r_w,-p+p,-r_w+mesh_level*mm, LcWiremesh * mm};                                      // right circle
p2b1_3 = newp; Point(p2b1_3) = {-p+p,-p+p,0+mesh_level*mm, LcWiremesh * mm};                                             // top circle
// p2b1_4 = newp; Point(p2b1_4) = {-p+p-r_w,-p+p,-r_w+mesh_level*mm, LcWiremesh * mm};                                   // left circle

l2b1_1 = newl; Circle(l2b1_1) = {p2b1_1, p2b1_0, p2b1_2};
l2b1_2 = newl; Circle(l2b1_2) = {p2b1_2, p2b1_0, p2b1_3};
l2b1s_3 = newl; Line(l2b1s_3) = {p2b1_1, p2b1_3};

ll2b1_1 = newll; Line Loop(ll2b1_1) = {l2b1_1, l2b1_2, -l2b1s_3};

s2b1_1 = news; Plane Surface(s2b1_1) = {ll2b1_1};

tmpba_2_1[] = {s2b1_1};
saaf_2_2_1[] = tmpba_2_1[0];

tmpbb_2_1[] = Extrude {0,1.325*y2_sp_wind_fac +0.00*sp_fac_str1,0} {
  Surface{tmpba_2_1[0]};
};

sab_2_1[] = {};
sab_2_1[] += tmpbb_2_1[{2:4}];


ptsac_2_1 = newp; Point(ptsac_2_1) = { p+p +2.00*p +0.00*sp_fac1, -p+p +0.00*p +1.75*y2_sp_wind_fac +0.00*R*Cos(alpha/2) +0.00*sp_fac1 +0.00*sp_fac_str1, -Rtn+r_w +1.00*mesh_level*mm +1.175*sp_fac_r2*r_w +0.00*R*Sin(alpha/2) +0.00*y2_sp_wind_fac +0.00*sp_fac_str2 };



tmpbc_2_1[] = Extrude {{y2_sp_wind_fac2,0,0}, {1,0,0}, { p+p +2.00*p +0.00*sp_fac1, -p+p +0.00*p +1.75*y2_sp_wind_fac +0.00*R*Cos(alpha/2) +0.00*sp_fac1 +0.00*sp_fac_str1, -Rtn+r_w +1.00*mesh_level*mm +1.175*sp_fac_r2*r_w +0.00*R*Sin(alpha/2) +0.00*y2_sp_wind_fac +0.00*sp_fac_str2 }, 0.20*alpha } {
  Surface{tmpbb_2_1[0]};
};

sac_2_1[] = {};
sac_2_1[] += tmpbc_2_1[{2:4}];



// tmpbc2_2_1[] = Extrude {0,0,1.00*y2_sp_wind_fac +0.00*sp_fac_str2} {
//   Surface{tmpbd_2_1[0]};
// };
// 
// sac2_2_1[] = {};
// sac2_2_1[] += tmpbc2_2_1[{2:4}];


ptsad_2_1 = newp; Point(ptsad_2_1) = { p+p +2.00*p +0.00*sp_fac1, -p+p +0.00*p +2.10*y2_sp_wind_fac +0.00*R*Cos(alpha/2) +0.00*sp_fac1 +0.00*sp_fac_str1, -Rtp+r_w +1.00*mesh_level*mm +0.825*sp_fac_r2*r_w +0.00*R*Sin(alpha/2) +0.00*y2_sp_wind_fac +0.00*sp_fac_str2 };



tmpbd_2_1[] = Extrude {{y2_sp_wind_fac2,0,0}, {1,0,0}, { p+p +2.00*p +0.00*sp_fac1, -p+p +0.00*p +2.10*y2_sp_wind_fac +0.00*R*Cos(alpha/2) +0.00*sp_fac1 +0.00*sp_fac_str1, -Rtp+r_w +1.00*mesh_level*mm +0.825*sp_fac_r2*r_w +0.00*R*Sin(alpha/2) +0.00*y2_sp_wind_fac +0.00*sp_fac_str2 }, 0.60*alpha} {
  Surface{tmpbc_2_1[0]};
};

sad_2_1[] = {};
sad_2_1[] += tmpbd_2_1[{2:4}];



// Wire 2b2


// tmpba_2_2[] = {tmpbd_2_1[0]};
// 
// tmpbb_2_2[] = Extrude {0,1.00*y2_sp_wind_fac +0.00*sp_fac_str1,0} {
//   Surface{tmpba_2_2[0]};
// };
// 
// sab_2_2[] = {};
// sab_2_2[] += tmpbb_2_2[{2:4}];


ptsac_2_2 = newp; Point(ptsac_2_2) = { p+p +2.00*p +0.00*sp_fac1, p+p +2.00*p -2.10*y2_sp_wind_fac +0.00*R*Cos(alpha/2) +0.00*sp_fac1 +0.00*sp_fac_str1, Rtp-r_w +1.00*mesh_level*mm -0.825*sp_fac_r2*r_w +0.00*R*Sin(alpha/2) +0.00*y2_sp_wind_fac +0.00*sp_fac_str2 };



tmpbc_2_2[] = Extrude {{x1_sp_wind_fac2,0,0}, {-1,0,0}, { p+p +2.00*p +0.00*sp_fac1, p+p +2.00*p -2.10*y2_sp_wind_fac +0.00*R*Cos(alpha/2) +0.00*sp_fac1 +0.00*sp_fac_str1, Rtp-r_w +1.00*mesh_level*mm -0.825*sp_fac_r2*r_w +0.00*R*Sin(alpha/2) +0.00*y2_sp_wind_fac +0.00*sp_fac_str2 }, 0.60*alpha} {
  Surface{tmpbd_2_1[0]};
};

sac_2_2[] = {};
sac_2_2[] += tmpbc_2_2[{2:4}];



// tmpbc2_2_2[] = Extrude {0,0,1.00*y2_sp_wind_fac +0.00*sp_fac_str2} {
//   Surface{tmpbc_2_2[0]};
// };
// 
// sac2_2_2[] = {};
// sac2_2_2[] += tmpbc2_2_2[{2:4}];


ptsad_2_2 = newp; Point(ptsad_2_2) = { p+p +2.00*p +0.00*sp_fac1, p+p +2.00*p -1.75*y2_sp_wind_fac +0.00*R*Cos(alpha/2) +0.00*sp_fac1 +0.00*sp_fac_str1, Rtp-r_w +1.00*mesh_level*mm -1.175*sp_fac_r2*r_w +0.00*R*Sin(alpha/2) +0.00*y2_sp_wind_fac +0.00*sp_fac_str2 };



tmpbd_2_2[] = Extrude {{y2_sp_wind_fac2,0,0}, {-1,0,0},{ p+p +2.00*p +0.00*sp_fac1, p+p +2.00*p -1.75*y2_sp_wind_fac +0.00*R*Cos(alpha/2) +0.00*sp_fac1 +0.00*sp_fac_str1, Rtp-r_w +1.00*mesh_level*mm -1.175*sp_fac_r2*r_w +0.00*R*Sin(alpha/2) +0.00*y2_sp_wind_fac +0.00*sp_fac_str2 }, 0.20*alpha} {
  Surface{tmpbc_2_2[0]};
};

sad_2_2[] = {};
sad_2_2[] += tmpbd_2_2[{2:4}];



tmpbe_2_2[] = Extrude {0,1.325*y2_sp_wind_fac +0.00*sp_fac_str2,0} {
  Surface{tmpbd_2_2[0]};
};

sae_2_2[] = {};
sae_2_2[] += tmpbe_2_2[{2:4}];

saef_2_2_2[] = {};
saef_2_2_2[] = {tmpbe_2_2[0]};



/*
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////// Line Definition


// *******************************
//
// Comparative IF Statement
//
// *******************************

For q In {1:2}
  For r In {1:2}
    For s In {1:2}

    // *******************************
    //
    // Face physsurf_bdh_1_1 (Corner 1 - Corner 2) - sbb_1_1 / sbb_1_2
    //
    // *******************************

      If(q == 1 && r == 1) 

      ll_bdhbbt~{q}~{r}~{s}[] = {};
      ll_bdhbbb~{q}~{r}~{s}[] = {};

      pbdhbbt~{q}~{r}~{s}() = {};
      pbdhbbb~{q}~{r}~{s}() = {};

      ll_bdhbb~{q}~{r}~{s}[] = {};

        ll_bdhbb~{q}~{r}~{s}[] += Boundary{ Surface{sbb~{r}~{s}[2]}; };

      For t In {0:3}

        If(t == 0 && s == 2)
          ll_bdhbbc0~{q}~{r}~{s}[] = {};
          pbdhbbc0~{q}~{r}~{s}() = {};
          ll_bdhbbc~{q}~{r}~{s}[] = {};

          ll_bdhbbc~{q}~{r}~{s}[] += Boundary{ Surface{tmpbb~{r}~{s}[t]}; };
          ll_bdhbbc0~{q}~{r}~{s}[] = Abs(ll_bdhbbc~{q}~{r}~{s}[t]);
          pbdhbbc0~{q}~{r}~{s}() = PointsOf{ Line{Abs(ll_bdhbbc~{q}~{r}~{s}[t])}; };
        EndIf

        If(t == 1)
          ll_bdhbbb~{q}~{r}~{s}[] = Abs(ll_bdhbb~{q}~{r}~{s}[t]);
          pbdhbbb~{q}~{r}~{s}() = PointsOf{ Line{Abs(ll_bdhbb~{q}~{r}~{s}[t])}; };
        EndIf  

        If(t == 1 && s == 2)
          ll_bdhbbc1~{q}~{r}~{s}[] = {};
          pbdhbbc1~{q}~{r}~{s}() = {};

          ll_bdhbbc1~{q}~{r}~{s}[] = Abs(ll_bdhbbc~{q}~{r}~{s}[t]);
          pbdhbbc1~{q}~{r}~{s}() = PointsOf{ Line{Abs(ll_bdhbbc~{q}~{r}~{s}[t])}; };
        EndIf

        If(t == 2 && s == 2)
          ll_bdhbb1~{q}~{r}~{s}[] = {};
          pbdhbb1~{q}~{r}~{s}() = {};

          ll_bdhbb1~{q}~{r}~{s}[] = Abs(ll_bdhbb~{q}~{r}~{s}[t]);
          pbdhbb1~{q}~{r}~{s}() = PointsOf{ Line{Abs(ll_bdhbb~{q}~{r}~{s}[t])}; };
        EndIf

        If(t == 3)  
          ll_bdhbbt~{q}~{r}~{s}[] = Abs(ll_bdhbb~{q}~{r}~{s}[t]);
          pbdhbbt~{q}~{r}~{s}() = PointsOf{ Line{Abs(ll_bdhbb~{q}~{r}~{s}[t])}; };
        EndIf

      EndFor

      EndIf


    // *******************************
    //
    // Face physsurf_bdh_1_1 (Corner 1 - Corner 2) - sbc_1_1 / sbc_1_2
    //
    // *******************************

      If(q == 1 && r == 1) 

      ll_bdhbct~{q}~{r}~{s}[] = {};
      ll_bdhbbc~{q}~{r}~{s}[] = {};

      pbdhbct~{q}~{r}~{s}() = {};
      pbdhbbc~{q}~{r}~{s}() = {};

      ll_bdhbc~{q}~{r}~{s}[] = {};

        ll_bdhbc~{q}~{r}~{s}[] += Boundary{ Surface{sbc~{r}~{s}[2]}; };

      For t In {0:3}

        If(t == 0 && s == 2)
          ll_bdhbcc0~{q}~{r}~{s}[] = {};
          pbdhbcc0~{q}~{r}~{s}() = {};
          ll_bdhbcc~{q}~{r}~{s}[] = {};

          ll_bdhbcc~{q}~{r}~{s}[] += Boundary{ Surface{tmpbc~{r}~{s}[t]}; };
          ll_bdhbcc0~{q}~{r}~{s}[] = Abs(ll_bdhbcc~{q}~{r}~{s}[t]);
          pbdhbcc0~{q}~{r}~{s}() = PointsOf{ Line{Abs(ll_bdhbcc~{q}~{r}~{s}[t])}; };
        EndIf

        If(t == 1)
          ll_bdhbbc~{q}~{r}~{s}[] = Abs(ll_bdhbc~{q}~{r}~{s}[t]);
          pbdhbbc~{q}~{r}~{s}() = PointsOf{ Line{Abs(ll_bdhbc~{q}~{r}~{s}[t])}; };
        EndIf  

        If(t == 1 && s == 2)
          ll_bdhbcc1~{q}~{r}~{s}[] = {};
          pbdhbcc1~{q}~{r}~{s}() = {};

          ll_bdhbcc1~{q}~{r}~{s}[] = Abs(ll_bdhbcc~{q}~{r}~{s}[t]);
          pbdhbcc1~{q}~{r}~{s}() = PointsOf{ Line{Abs(ll_bdhbcc~{q}~{r}~{s}[t])}; };
        EndIf

        If(t == 2 && s == 2)
          ll_bdhbc1~{q}~{r}~{s}[] = {};
          pbdhbc1~{q}~{r}~{s}() = {};

          ll_bdhbc1~{q}~{r}~{s}[] = Abs(ll_bdhbc~{q}~{r}~{s}[t]);
          pbdhbc1~{q}~{r}~{s}() = PointsOf{ Line{Abs(ll_bdhbc~{q}~{r}~{s}[t])}; };
        EndIf

        If(t == 3)  
          ll_bdhbct~{q}~{r}~{s}[] = Abs(ll_bdhbc~{q}~{r}~{s}[t]);
          pbdhbct~{q}~{r}~{s}() = PointsOf{ Line{Abs(ll_bdhbc~{q}~{r}~{s}[t])}; };
        EndIf

      EndFor

      EndIf


    // *******************************
    //
    // Face physsurf_bdh_1_2 (Corner 2 - Corner 3) - sab_2_1 / sab_2_2
    //
    // *******************************

      If(q == 1 && r == 2)

      ll_bdhabt~{q}~{r}~{s}[] = {};
      ll_bdhabb~{q}~{r}~{s}[] = {}; 

      pbdhabt~{q}~{r}~{s}() = {};
      pbdhabb~{q}~{r}~{s}() = {};

      ll_bdhab~{q}~{r}~{s}[] = {};

        ll_bdhab~{q}~{r}~{s}[] += Boundary{ Surface{sab~{r}~{s}[2]}; };

      For t In {0:3}

        If(t == 0 && s == 2)
          ll_bdhabc0~{q}~{r}~{s}[] = {};
          pbdhabc0~{q}~{r}~{s}() = {};
          ll_bdhabc~{q}~{r}~{s}[] = {};

          ll_bdhabc~{q}~{r}~{s}[] += Boundary{ Surface{tmpab~{r}~{s}[t]}; };
          ll_bdhabc0~{q}~{r}~{s}[] = Abs(ll_bdhabc~{q}~{r}~{s}[t]);
          pbdhabc0~{q}~{r}~{s}() = PointsOf{ Line{Abs(ll_bdhabc~{q}~{r}~{s}[t])}; };
        EndIf

        If(t == 1)
          ll_bdhabb~{q}~{r}~{s}[] = Abs(ll_bdhab~{q}~{r}~{s}[t]);
          pbdhabb~{q}~{r}~{s}() = PointsOf{ Line{Abs(ll_bdhab~{q}~{r}~{s}[t])}; };
        EndIf

        If(t == 1 && s == 2)
          ll_bdhabc1~{q}~{r}~{s}[] = {};
          pbdhabc1~{q}~{r}~{s}() = {};

          ll_bdhabc1~{q}~{r}~{s}[] = Abs(ll_bdhabc~{q}~{r}~{s}[t]);
          pbdhabc1~{q}~{r}~{s}() = PointsOf{ Line{Abs(ll_bdhabc~{q}~{r}~{s}[t])}; };
        EndIf

        If(t == 2 && s == 2)
          ll_bdhab1~{q}~{r}~{s}[] = {}; 
          pbdhab1~{q}~{r}~{s}() = {};

          ll_bdhab1~{q}~{r}~{s}[] = Abs(ll_bdhab~{q}~{r}~{s}[t]);
          pbdhab1~{q}~{r}~{s}() = PointsOf{ Line{Abs(ll_bdhab~{q}~{r}~{s}[t])}; };
        EndIf

        If(t == 3)
          ll_bdhabt~{q}~{r}~{s}[] = Abs(ll_bdhab~{q}~{r}~{s}[t]);
          pbdhabt~{q}~{r}~{s}() = PointsOf{ Line{Abs(ll_bdhab~{q}~{r}~{s}[t])}; };
        EndIf

      EndFor

      EndIf


    // *******************************
    //
    // Face physsurf_bdh_1_2 (Corner 2 - Corner 3) - sac_2_1 / sac_2_2
    //
    // *******************************

      If(q == 1 && r == 2)

      ll_bdhact~{q}~{r}~{s}[] = {};
      ll_bdhacb~{q}~{r}~{s}[] = {}; 

      pbdhact~{q}~{r}~{s}() = {};
      pbdhacb~{q}~{r}~{s}() = {};

      ll_bdhac~{q}~{r}~{s}[] = {};

        ll_bdhac~{q}~{r}~{s}[] += Boundary{ Surface{sac~{r}~{s}[2]}; };

      For t In {0:3}

        If(t == 0 && s == 2)
          ll_bdhacc0~{q}~{r}~{s}[] = {};
          pbdhacc0~{q}~{r}~{s}() = {};
          ll_bdhacc~{q}~{r}~{s}[] = {};

          ll_bdhacc~{q}~{r}~{s}[] += Boundary{ Surface{tmpac~{r}~{s}[t]}; };
          ll_bdhacc0~{q}~{r}~{s}[] = Abs(ll_bdhacc~{q}~{r}~{s}[t]);
          pbdhacc0~{q}~{r}~{s}() = PointsOf{ Line{Abs(ll_bdhacc~{q}~{r}~{s}[t])}; };
        EndIf

        If(t == 1)
          ll_bdhacb~{q}~{r}~{s}[] = Abs(ll_bdhac~{q}~{r}~{s}[t]);
          pbdhacb~{q}~{r}~{s}() = PointsOf{ Line{Abs(ll_bdhac~{q}~{r}~{s}[t])}; };
        EndIf

        If(t == 1 && s == 2)
          ll_bdhacc1~{q}~{r}~{s}[] = {};
          pbdhacc1~{q}~{r}~{s}() = {};

          ll_bdhacc1~{q}~{r}~{s}[] = Abs(ll_bdhacc~{q}~{r}~{s}[t]);
          pbdhacc1~{q}~{r}~{s}() = PointsOf{ Line{Abs(ll_bdhacc~{q}~{r}~{s}[t])}; };
        EndIf

        If(t == 2 && s == 2)
          ll_bdhac1~{q}~{r}~{s}[] = {}; 
          pbdhac1~{q}~{r}~{s}() = {};

          ll_bdhac1~{q}~{r}~{s}[] = Abs(ll_bdhac~{q}~{r}~{s}[t]);
          pbdhac1~{q}~{r}~{s}() = PointsOf{ Line{Abs(ll_bdhac~{q}~{r}~{s}[t])}; };
        EndIf

        If(t == 3)
          ll_bdhact~{q}~{r}~{s}[] = Abs(ll_bdhac~{q}~{r}~{s}[t]);
          pbdhact~{q}~{r}~{s}() = PointsOf{ Line{Abs(ll_bdhac~{q}~{r}~{s}[t])}; };
        EndIf

      EndFor

      EndIf


    // *******************************
    //
    // Face physsurf_bdh_2_1 (Corner 3 - Corner 4) - sab_1_1[] / sab_1_2[]
    //
    // *******************************

      If(q == 2 && r == 1)

      ll_bdhabt~{q}~{r}~{s}[] = {};
      ll_bdhabb~{q}~{r}~{s}[] = {}; 

      pbdhabt~{q}~{r}~{s}() = {};
      pbdhabb~{q}~{r}~{s}() = {};

      ll_bdhab~{q}~{r}~{s}[] = {};

        ll_bdhab~{q}~{r}~{s}[] += Boundary{ Surface{sab~{r}~{s}[2]}; };

      For t In {0:3}

        If(t == 0 && s == 2)
          ll_bdhabc0~{q}~{r}~{s}[] = {};
          pbdhabc0~{q}~{r}~{s}() = {};
          ll_bdhabc~{q}~{r}~{s}[] = {};

          ll_bdhabc~{q}~{r}~{s}[] += Boundary{ Surface{tmpab~{r}~{s}[t]}; };
          ll_bdhabc0~{q}~{r}~{s}[] = Abs(ll_bdhabc~{q}~{r}~{s}[t]);
          pbdhabc0~{q}~{r}~{s}() = PointsOf{ Line{Abs(ll_bdhabc~{q}~{r}~{s}[t])}; };
        EndIf
      
        If(t == 3)
          ll_bdhabb~{q}~{r}~{s}[] = Abs(ll_bdhab~{q}~{r}~{s}[t]);
          pbdhabb~{q}~{r}~{s}() = PointsOf{ Line{Abs(ll_bdhab~{q}~{r}~{s}[t])}; };
        EndIf

        If(t == 2 && s == 2)
          ll_bdhab1~{q}~{r}~{s}[] = {};
          pbdhab1~{q}~{r}~{s}() = {};

          ll_bdhab1~{q}~{r}~{s}[] = Abs(ll_bdhab~{q}~{r}~{s}[t]);
          pbdhab1~{q}~{r}~{s}() = PointsOf{ Line{Abs(ll_bdhab~{q}~{r}~{s}[t])}; };
        EndIf

        If(t == 1)
          ll_bdhabt~{q}~{r}~{s}[] = Abs(ll_bdhab~{q}~{r}~{s}[t]);
          pbdhabt~{q}~{r}~{s}() = PointsOf{ Line{Abs(ll_bdhab~{q}~{r}~{s}[t])}; };
        EndIf

        If(t == 1 && s == 2)
          ll_bdhabc1~{q}~{r}~{s}[] = {};
          pbdhabc1~{q}~{r}~{s}() = {};

          ll_bdhabc1~{q}~{r}~{s}[] = Abs(ll_bdhabc~{q}~{r}~{s}[t]);
          pbdhabc1~{q}~{r}~{s}() = PointsOf{ Line{Abs(ll_bdhabc~{q}~{r}~{s}[t])}; };
        EndIf

      EndFor

      EndIf


    // *******************************
    //
    // Face physsurf_bdh_2_1 (Corner 3 - Corner 4) - sac_1_1[] / sac_1_2[]
    //
    // *******************************

      If(q == 2&& r == 1)

      ll_bdhact~{q}~{r}~{s}[] = {};
      ll_bdhacb~{q}~{r}~{s}[] = {}; 

      pbdhact~{q}~{r}~{s}() = {};
      pbdhacb~{q}~{r}~{s}() = {};

      ll_bdhac~{q}~{r}~{s}[] = {};

        ll_bdhac~{q}~{r}~{s}[] += Boundary{ Surface{sac~{r}~{s}[2]}; };

      For t In {0:3}

        If(t == 0 && s == 2)
          ll_bdhacc0~{q}~{r}~{s}[] = {};
          pbdhacc0~{q}~{r}~{s}() = {};
          ll_bdhacc~{q}~{r}~{s}[] = {};

          ll_bdhacc~{q}~{r}~{s}[] += Boundary{ Surface{tmpac~{r}~{s}[t]}; };
          ll_bdhacc0~{q}~{r}~{s}[] = Abs(ll_bdhacc~{q}~{r}~{s}[t]);
          pbdhacc0~{q}~{r}~{s}() = PointsOf{ Line{Abs(ll_bdhacc~{q}~{r}~{s}[t])}; };
        EndIf
      
        If(t == 3)
          ll_bdhacb~{q}~{r}~{s}[] = Abs(ll_bdhac~{q}~{r}~{s}[t]);
          pbdhacb~{q}~{r}~{s}() = PointsOf{ Line{Abs(ll_bdhac~{q}~{r}~{s}[t])}; };
        EndIf

        If(t == 2 && s == 2)
          ll_bdhac1~{q}~{r}~{s}[] = {};
          pbdhac1~{q}~{r}~{s}() = {};

          ll_bdhac1~{q}~{r}~{s}[] = Abs(ll_bdhac~{q}~{r}~{s}[t]);
          pbdhac1~{q}~{r}~{s}() = PointsOf{ Line{Abs(ll_bdhac~{q}~{r}~{s}[t])}; };
        EndIf

        If(t == 1)
          ll_bdhact~{q}~{r}~{s}[] = Abs(ll_bdhac~{q}~{r}~{s}[t]);
          pbdhact~{q}~{r}~{s}() = PointsOf{ Line{Abs(ll_bdhac~{q}~{r}~{s}[t])}; };
        EndIf

        If(t == 1 && s == 2)
          ll_bdhacc1~{q}~{r}~{s}[] = {};
          pbdhacc1~{q}~{r}~{s}() = {};

          ll_bdhacc1~{q}~{r}~{s}[] = Abs(ll_bdhacc~{q}~{r}~{s}[t]);
          pbdhacc1~{q}~{r}~{s}() = PointsOf{ Line{Abs(ll_bdhacc~{q}~{r}~{s}[t])}; };
        EndIf

      EndFor

      EndIf


    // *******************************
    //
    // Face physsurf_bdh_2_2 (Corner 4 - Corner 1) - sbb_2_1[] / sbb_2_2[]
    //
    // *******************************

      If(q == 2&& r == 2) 

      ll_bdhbbt~{q}~{r}~{s}[] = {};
      ll_bdhbbb~{q}~{r}~{s}[] = {};
      pbdhbbt~{q}~{r}~{s}() = {};
      pbdhbbb~{q}~{r}~{s}() = {};

      ll_bdhbb~{q}~{r}~{s}[] = {};

        ll_bdhbb~{q}~{r}~{s}[] += Boundary{ Surface{sbb~{r}~{s}[2]}; };

      For t In {0:3}
    
        If(t == 0 && s == 2)
          ll_bdhbbc0~{q}~{r}~{s}[] = {};
          pbdhbbc0~{q}~{r}~{s}() = {};
          ll_bdhbbc~{q}~{r}~{s}[] = {};

          ll_bdhbbc~{q}~{r}~{s}[] += Boundary{ Surface{tmpbb~{r}~{s}[t]}; };
          ll_bdhbbc0~{q}~{r}~{s}[] = Abs(ll_bdhbbc~{q}~{r}~{s}[t]);
          pbdhbbc0~{q}~{r}~{s}() = PointsOf{ Line{Abs(ll_bdhbbc~{q}~{r}~{s}[t])}; };
        EndIf

        If(t == 3)
          ll_bdhbbb~{q}~{r}~{s}[] = Abs(ll_bdhbb~{q}~{r}~{s}[t]);
          pbdhbbb~{q}~{r}~{s}() = PointsOf{ Line{Abs(ll_bdhbb~{q}~{r}~{s}[t])}; };
        EndIf

        If(t == 2 && s == 2)
          ll_bdhbb1~{q}~{r}~{s}[] = {};
          pbdhbb1~{q}~{r}~{s}() = {};

          ll_bdhbb1~{q}~{r}~{s}[] = Abs(ll_bdhbb~{q}~{r}~{s}[t]);
          pbdhbb1~{q}~{r}~{s}() = PointsOf{ Line{Abs(ll_bdhbb~{q}~{r}~{s}[t])}; };
        EndIf

        If(t == 1)  
          ll_bdhbbt~{q}~{r}~{s}[] = Abs(ll_bdhbb~{q}~{r}~{s}[t]);
          pbdhbbt~{q}~{r}~{s}() = PointsOf{ Line{Abs(ll_bdhbb~{q}~{r}~{s}[t])}; };
        EndIf

        If(t == 1 && s == 2)  
          ll_bdhbbc1~{q}~{r}~{s}[] = {};
          pbdhbbc1~{q}~{r}~{s}() = {};

          ll_bdhbbc1~{q}~{r}~{s}[] = Abs(ll_bdhbbc~{q}~{r}~{s}[t]);
          pbdhbbc1~{q}~{r}~{s}() = PointsOf{ Line{Abs(ll_bdhbbc~{q}~{r}~{s}[t])}; };
        EndIf

      EndFor

      EndIf


    // *******************************
    //
    // Face physsurf_bdh_2_2 (Corner 4 - Corner 1) - sbc_2_1[] / sbc_2_2[]
    //
    // *******************************

      If(q == 2&& r == 2) 

      ll_bdhbct~{q}~{r}~{s}[] = {};
      ll_bdhbcb~{q}~{r}~{s}[] = {};
      pbdhbct~{q}~{r}~{s}() = {};
      pbdhbcb~{q}~{r}~{s}() = {};

      ll_bdhbc~{q}~{r}~{s}[] = {};

        ll_bdhbc~{q}~{r}~{s}[] += Boundary{ Surface{sbc~{r}~{s}[2]}; };

      For t In {0:3}
    
        If(t == 0 && s == 2)
          ll_bdhbcc0~{q}~{r}~{s}[] = {};
          pbdhbcc0~{q}~{r}~{s}() = {};
          ll_bdhbcc~{q}~{r}~{s}[] = {};

          ll_bdhbcc~{q}~{r}~{s}[] += Boundary{ Surface{tmpbc~{r}~{s}[t]}; };
          ll_bdhbcc0~{q}~{r}~{s}[] = Abs(ll_bdhbcc~{q}~{r}~{s}[t]);
          pbdhbcc0~{q}~{r}~{s}() = PointsOf{ Line{Abs(ll_bdhbc~{q}~{r}~{s}[t])}; };
        EndIf

        If(t == 3)
          ll_bdhbcb~{q}~{r}~{s}[] = Abs(ll_bdhbc~{q}~{r}~{s}[t]);
          pbdhbcb~{q}~{r}~{s}() = PointsOf{ Line{Abs(ll_bdhbc~{q}~{r}~{s}[t])}; };
        EndIf

        If(t == 2 && s == 2)
          ll_bdhbc1~{q}~{r}~{s}[] = {};
          pbdhbc1~{q}~{r}~{s}() = {};

          ll_bdhbc1~{q}~{r}~{s}[] = Abs(ll_bdhbc~{q}~{r}~{s}[t]);
          pbdhbc1~{q}~{r}~{s}() = PointsOf{ Line{Abs(ll_bdhbc~{q}~{r}~{s}[t])}; };
        EndIf

        If(t == 1)  
          ll_bdhbct~{q}~{r}~{s}[] = Abs(ll_bdhbc~{q}~{r}~{s}[t]);
          pbdhbct~{q}~{r}~{s}() = PointsOf{ Line{Abs(ll_bdhbc~{q}~{r}~{s}[t])}; };
        EndIf

        If(t == 1 && s == 2)  
          ll_bdhbcc1~{q}~{r}~{s}[] = {};
          pbdhbcc1~{q}~{r}~{s}() = {};

          ll_bdhbcc1~{q}~{r}~{s}[] = Abs(ll_bdhbcc~{q}~{r}~{s}[t]);
          pbdhbcc1~{q}~{r}~{s}() = PointsOf{ Line{Abs(ll_bdhbcc~{q}~{r}~{s}[t])}; };
        EndIf

      EndFor

      EndIf

  EndFor
 EndFor
EndFor


For q In {2:2}
  For r In {1:2}
    For s In {2:2}

    // *******************************
    //
    // End faces - sacf2_1_2, sacf2_2_2, sbcf2_1_2, sbcf2_2_2
    //
    // *******************************

      If(q == 2 && s == 2) 

      ll_facesacb~{q}~{r}~{s}[] = {};

        ll_facesacb~{q}~{r}~{s}[] += Boundary{ Surface{sacf~{q}~{r}~{s}[]}; };

        ptsacfmin~{q}~{r}~{s}() = PointsOf{ Line{Abs(ll_facesacb~{q}~{r}~{s}[0])}; };
        psacfminchz~{q}~{r}~{s}[] = Point{ ptsacfmin~{q}~{r}~{s}[0] };

      For w In {1:1}

        pfacesacbc1chz~{q}~{r}~{s}[] = Point{ ptsacfmin~{q}~{r}~{s}[w] };

        If( pfacesacbc1chz~{q}~{r}~{s}[2] < psacfminchz~{q}~{r}~{s}[2] )

          psacfminz[] = pfacesacbc1chz~{q}~{r}~{s}[2];
          ptsacfminz~{q}~{r}~{s}() = ptsacfmin~{q}~{r}~{s}[w];

        Else

          psacfminz[]  = psacfminchz~{q}~{r}~{s}[2];
          ptsacfminz~{q}~{r}~{s}() = ptsacfmin~{q}~{r}~{s}[0];

        EndIf

      EndFor


      For t In {1:2}

        ptfacesacbc1~{q}~{r}~{s}() = PointsOf{ Line{Abs(ll_facesacb~{q}~{r}~{s}[t])}; };

        For l In {0:1}

          pfacesacbc1chz~{q}~{r}~{s}[] = Point{ ptfacesacbc1~{q}~{r}~{s}[l] };

          If( pfacesacbc1chz~{q}~{r}~{s}[2] < psacfminchz~{q}~{r}~{s}[2] )

            psacfminz[]  = pfacesacbc1chz~{q}~{r}~{s}[2];
            ptsacfminz~{q}~{r}~{s}() = ptfacesacbc1~{q}~{r}~{s}[l];

          EndIf

        EndFor

      EndFor


      ll_facesacb~{q}~{r}~{s}[] = {};

        ll_facesacb~{q}~{r}~{s}[] += Boundary{ Surface{sacf~{q}~{r}~{s}[]}; };

        ptsacfmin~{q}~{r}~{s}() = PointsOf{ Line{Abs(ll_facesacb~{q}~{r}~{s}[0])}; };
        psacfminchz~{q}~{r}~{s}[] = Point{ ptsacfmin~{q}~{r}~{s}[0] };

      For w In {1:1}

        pfacesacbc1chz~{q}~{r}~{s}[] = Point{ ptsacfmin~{q}~{r}~{s}[w] };

        If( pfacesacbc1chz~{q}~{r}~{s}[2] < psacfminchz~{q}~{r}~{s}[2] )

          psacfminz[]  = pfacesacbc1chz~{q}~{r}~{s}[2];
          ptsacfminz~{q}~{r}~{s}() = ptsacfmin~{q}~{r}~{s}[w];

        Else

          psacfminz[]  = psacfminchz~{q}~{r}~{s}[2];
          ptsacfminz~{q}~{r}~{s}() = ptsacfmin~{q}~{r}~{s}[0];

        EndIf

      EndFor

      
      For t In {1:2}

        ptfacesacbc1~{q}~{r}~{s}() = PointsOf{ Line{Abs(ll_facesacb~{q}~{r}~{s}[t])}; };

        For l In {0:1}

          pfacesacbc1chz~{q}~{r}~{s}[] = Point{ ptfacesacbc1~{q}~{r}~{s}[l] };
    
          If( pfacesacbc1chz~{q}~{r}~{s}[2] < psacfminchz~{q}~{r}~{s}[2] )

            psacfminz[]  = pfacesacbc1chz~{q}~{r}~{s}[2];
            ptsacfminz~{q}~{r}~{s}() = ptfacesacbc1~{q}~{r}~{s}[l];

          EndIf

        EndFor

      EndFor


      ll_facesacb~{q}~{r}~{s}[] = {};

        ll_facesacb~{q}~{r}~{s}[] += Boundary{ Surface{sacf~{q}~{r}~{s}[]}; };

        ptsacfmax~{q}~{r}~{s}() = PointsOf{ Line{Abs(ll_facesacb~{q}~{r}~{s}[0])}; };
        psacfmaxchz~{q}~{r}~{s}[] = Point{ ptsacfmax~{q}~{r}~{s}[0] };

      For w In {1:1}

        pfacesacbc1chz~{q}~{r}~{s}[] = Point{ ptsacfmax~{q}~{r}~{s}[w] };

        If( pfacesacbc1chz~{q}~{r}~{s}[2] > psacfmaxchz~{q}~{r}~{s}[2] )

          psacfmaxz~{q}~{r}~{s}[] = pfacesacbc1chz~{q}~{r}~{s}[2];
          ptsacfmaxz~{q}~{r}~{s}() = ptsacfmax~{q}~{r}~{s}[w];

        Else

          psacfmaxz~{q}~{r}~{s}[] = psacfmaxchz~{q}~{r}~{s}[2];
          ptsacfmaxz~{q}~{r}~{s}() = ptsacfmax~{q}~{r}~{s}[0];

        EndIf

      EndFor

      
      For t In {1:2}

        ptfacesacbc1~{q}~{r}~{s}() = PointsOf{ Line{Abs(ll_facesacb~{q}~{r}~{s}[t])}; };

        For l In {0:1}

          pfacesacbc1chz~{q}~{r}~{s}[] = Point{ ptfacesacbc1~{q}~{r}~{s}[l] };
    
          If( pfacesacbc1chz~{q}~{r}~{s}[2] > psacfmaxchz~{q}~{r}~{s}[2] )

            psacfmaxz~{q}~{r}~{s}[] = pfacesacbc1chz~{q}~{r}~{s}[2];
            ptsacfmaxz~{q}~{r}~{s}() = ptfacesacbc1~{q}~{r}~{s}[l];

          EndIf

        EndFor

      EndFor


      ll_facesacb~{q}~{r}~{s}[] = {};

        ll_facesacb~{q}~{r}~{s}[] += Boundary{ Surface{sacf~{q}~{r}~{s}[]}; };

        ptsacfmax~{q}~{r}~{s}() = PointsOf{ Line{Abs(ll_facesacb~{q}~{r}~{s}[0])}; };
        psacfmaxchz~{q}~{r}~{s}[] = Point{ ptsacfmax~{q}~{r}~{s}[0] };

      For w In {1:1}

        pfacesacbc1chz~{q}~{r}~{s}[] = Point{ ptsacfmax~{q}~{r}~{s}[w] };

        If( pfacesacbc1chz~{q}~{r}~{s}[2] > psacfmaxchz~{q}~{r}~{s}[2] )

          psacfmaxz~{q}~{r}~{s}[] = pfacesacbc1chz~{q}~{r}~{s}[2];
          ptsacfmaxz~{q}~{r}~{s}() = ptsacfmax~{q}~{r}~{s}[w];

        Else

          psacfmaxz~{q}~{r}~{s}[] = psacfmaxchz~{q}~{r}~{s}[2];
          ptsacfmaxz~{q}~{r}~{s}() = ptsacfmax~{q}~{r}~{s}[0];

        EndIf

      EndFor

      
      For t In {1:2}

        ptfacesacbc1~{q}~{r}~{s}() = PointsOf{ Line{Abs(ll_facesacb~{q}~{r}~{s}[t])}; };

        For l In {0:1}

          pfacesacbc1chz~{q}~{r}~{s}[] = Point{ ptfacesacbc1~{q}~{r}~{s}[l] };
    
          If( pfacesacbc1chz~{q}~{r}~{s}[2] > psacfmaxchz~{q}~{r}~{s}[2] )

            psacfmaxz~{q}~{r}~{s}[] = pfacesacbc1chz~{q}~{r}~{s}[2];
            ptsacfmaxz~{q}~{r}~{s}() = ptfacesacbc1~{q}~{r}~{s}[l];

          EndIf

        EndFor

      EndFor


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


      ll_facesbcb~{q}~{r}~{s}[] = {};

        ll_facesbcb~{q}~{r}~{s}[] += Boundary{ Surface{sbcf~{q}~{r}~{s}[]}; };

        ptsbcfmin~{q}~{r}~{s}() = PointsOf{ Line{Abs(ll_facesbcb~{q}~{r}~{s}[0])}; };
        psbcfminchz~{q}~{r}~{s}[] = Point{ ptsbcfmin~{q}~{r}~{s}[0] };

      For w In {1:1}

        pfacesbcbc1chz~{q}~{r}~{s}[] = Point{ ptsbcfmin~{q}~{r}~{s}[w] };

        If( pfacesbcbc1chz~{q}~{r}~{s}[2] < psbcfminchz~{q}~{r}~{s}[2] )

          psbcfminz~{q}~{r}~{s}[] = pfacesbcbc1chz~{q}~{r}~{s}[2];
          ptsbcfminz~{q}~{r}~{s}() = ptsbcfmin~{q}~{r}~{s}[w];

        Else

          psbcfminz~{q}~{r}~{s}[] = psbcfminchz~{q}~{r}~{s}[2];
          ptsbcfminz~{q}~{r}~{s}() = ptsbcfmin~{q}~{r}~{s}[0];

        EndIf

      EndFor

      
      For t In {1:2}

        ptfacesbcbc1~{q}~{r}~{s}() = PointsOf{ Line{Abs(ll_facesbcb~{q}~{r}~{s}[t])}; };

        For l In {0:1}

          pfacesbcbc1chz~{q}~{r}~{s}[] = Point{ ptfacesbcbc1~{q}~{r}~{s}[l] };
    
          If( pfacesbcbc1chz~{q}~{r}~{s}[2] < psbcfminchz~{q}~{r}~{s}[2] )

            psbcfminz~{q}~{r}~{s}[] = pfacesbcbc1chz~{q}~{r}~{s}[2];
            ptsbcfminz~{q}~{r}~{s}() = ptfacesbcbc1~{q}~{r}~{s}[l];

          EndIf

        EndFor

      EndFor


      ll_facesbcb~{q}~{r}~{s}[] = {};

        ll_facesbcb~{q}~{r}~{s}[] += Boundary{ Surface{sbcf~{q}~{r}~{s}[]}; };

        ptsbcfmin~{q}~{r}~{s}() = PointsOf{ Line{Abs(ll_facesbcb~{q}~{r}~{s}[0])}; };
        psbcfminchz~{q}~{r}~{s}[] = Point{ ptsbcfmin~{q}~{r}~{s}[0] };

      For w In {1:1}

        pfacesbcbc1chz~{q}~{r}~{s}[] = Point{ ptsbcfmin~{q}~{r}~{s}[w] };

        If( pfacesbcbc1chz~{q}~{r}~{s}[2] < psbcfminchz~{q}~{r}~{s}[2] )

          psbcfminz~{q}~{r}~{s}[] = pfacesbcbc1chz~{q}~{r}~{s}[2];
          ptsbcfminz~{q}~{r}~{s}() = ptsbcfmin~{q}~{r}~{s}[w];

        Else

          psbcfminz~{q}~{r}~{s}[] = psbcfminchz~{q}~{r}~{s}[2];
          ptsbcfminz~{q}~{r}~{s}() = ptsbcfmin~{q}~{r}~{s}[0];

        EndIf

      EndFor

    
      For t In {1:2}

        ptfacesbcbc1~{q}~{r}~{s}() = PointsOf{ Line{Abs(ll_facesbcb~{q}~{r}~{s}[t])}; };

        For l In {0:1}

          pfacesbcbc1chz~{q}~{r}~{s}[] = Point{ ptfacesbcbc1~{q}~{r}~{s}[l] };
    
          If( pfacesbcbc1chz~{q}~{r}~{s}[2] < psbcfminchz~{q}~{r}~{s}[2] )

            psbcfminz~{q}~{r}~{s}[] = pfacesbcbc1chz~{q}~{r}~{s}[2];
            ptsbcfminz~{q}~{r}~{s}() = ptfacesbcbc1~{q}~{r}~{s}[l];

          EndIf

        EndFor

      EndFor


      ll_facesbcb~{q}~{r}~{s}[] = {};

        ll_facesbcb~{q}~{r}~{s}[] += Boundary{ Surface{sbcf~{q}~{r}~{s}[]}; };

        ptsbcfmax~{q}~{r}~{s}() = PointsOf{ Line{Abs(ll_facesbcb~{q}~{r}~{s}[0])}; };
        psbcfmaxchz~{q}~{r}~{s}[] = Point{ ptsbcfmax~{q}~{r}~{s}[0] };

      For w In {1:1}

        pfacesbcbc1chz~{q}~{r}~{s}[] = Point{ ptsbcfmax~{q}~{r}~{s}[w] };

        If( pfacesbcbc1chz~{q}~{r}~{s}[2] > psbcfmaxchz~{q}~{r}~{s}[2] )

          psbcfmaxz~{q}~{r}~{s}[] = pfacesbcbc1chz~{q}~{r}~{s}[2];
          ptsbcfmaxz~{q}~{r}~{s}() = ptsbcfmax~{q}~{r}~{s}[w];

        Else

          psbcfmaxz~{q}~{r}~{s}[] = psbcfmaxchz~{q}~{r}~{s}[2];
          ptsbcfmaxz~{q}~{r}~{s}() = ptsbcfmax~{q}~{r}~{s}[0];

        EndIf

      EndFor

      
      For t In {1:2}

        ptfacesbcbc1~{q}~{r}~{s}() = PointsOf{ Line{Abs(ll_facesbcb~{q}~{r}~{s}[t])}; };

        For l In {0:1}

          pfacesbcbc1chz~{q}~{r}~{s}[] = Point{ ptfacesbcbc1~{q}~{r}~{s}[l] };
    
          If( pfacesbcbc1chz~{q}~{r}~{s}[2] > psbcfmaxchz~{q}~{r}~{s}[2] )

            psbcfmaxz~{q}~{r}~{s}[] = pfacesbcbc1chz~{q}~{r}~{s}[2];
            ptsbcfmaxz~{q}~{r}~{s}() = ptfacesbcbc1~{q}~{r}~{s}[l];

          EndIf

        EndFor

      EndFor


      ll_facesbcb~{q}~{r}~{s}[] = {};

        ll_facesbcb~{q}~{r}~{s}[] += Boundary{ Surface{sbcf~{q}~{r}~{s}[]}; };

        ptsbcfmax~{q}~{r}~{s}() = PointsOf{ Line{Abs(ll_facesbcb~{q}~{r}~{s}[0])}; };
        psbcfmaxchz~{q}~{r}~{s}[] = Point{ ptsbcfmax~{q}~{r}~{s}[0] };

      For w In {1:1}

        pfacesbcbc1chz~{q}~{r}~{s}[] = Point{ ptsbcfmax~{q}~{r}~{s}[w] };

        If( pfacesbcbc1chz~{q}~{r}~{s}[2] > psbcfmaxchz~{q}~{r}~{s}[2] )

          psbcfmaxz~{q}~{r}~{s}[] = pfacesbcbc1chz~{q}~{r}~{s}[2];
          ptsbcfmaxz~{q}~{r}~{s}() = ptsbcfmax~{q}~{r}~{s}[w];

        Else

          psbcfmaxz~{q}~{r}~{s}[] = psbcfmaxchz~{q}~{r}~{s}[2];
          ptsbcfmaxz~{q}~{r}~{s}() = ptsbcfmax~{q}~{r}~{s}[0];

        EndIf

      EndFor

      
      For t In {1:2}

        ptfacesbcbc1~{q}~{r}~{s}() = PointsOf{ Line{Abs(ll_facesbcb~{q}~{r}~{s}[t])}; };

        For l In {0:1}

          pfacesbcbc1chz~{q}~{r}~{s}[] = Point{ ptfacesbcbc1~{q}~{r}~{s}[l] };
    
          If( pfacesbcbc1chz~{q}~{r}~{s}[2] > psbcfmaxchz~{q}~{r}~{s}[2] )

            psbcfmaxz~{q}~{r}~{s}[] = pfacesbcbc1chz~{q}~{r}~{s}[2];
            ptsbcfmaxz~{q}~{r}~{s}() = ptfacesbcbc1~{q}~{r}~{s}[l];

          EndIf

        EndFor

      EndFor


  EndFor
 EndFor
EndFor
*/


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////// Shell Definition

/*
// SHELL

// --------------------------------------------------------------------------

// *******************************
//
// Corner 1
//
// *******************************
pc1_1 = newp; Point(pc1_1) = {geo_f_x*0+geo_f_x*m*a, geo_f_y*0+geo_f_y*n*a, tD/2,lcCopperPlateBdry};
pc2_1 = newp; Point(pc2_1) = {geo_f_x*0+geo_f_x*m*a, geo_f_y*0+geo_f_y*n*a, -1*tD/2,lcCopperPlateBdry};
pc3_1 = newp; Point(pc3_1) = {geo_f_x*0+geo_f_x*m*a, geo_f_y*0+geo_f_y*n*a, (2*tC+tD)/2,lcCopperPlateBdry};
pc4_1 = newp; Point(pc4_1) = {geo_f_x*0+geo_f_x*m*a, geo_f_y*0+geo_f_y*n*a, -1*(2*tC+tD)/2,lcCopperPlateBdry};

// *******************************
//
// Corner 2
//
// *******************************
pc1_2 = newp; Point(pc1_2) = {geo_f_x*a+geo_f_x*m*a, geo_f_y*0+geo_f_y*n*a, tD/2,lcCopperPlateBdry};
pc2_2 = newp; Point(pc2_2) = {geo_f_x*a+geo_f_x*m*a, geo_f_y*0+geo_f_y*n*a, -1*tD/2,lcCopperPlateBdry};
pc3_2 = newp; Point(pc3_2) = {geo_f_x*a+geo_f_x*m*a, geo_f_y*0+geo_f_y*n*a, (2*tC+tD)/2,lcCopperPlateBdry};
pc4_2 = newp; Point(pc4_2) = {geo_f_x*a+geo_f_x*m*a, geo_f_y*0+geo_f_y*n*a, -1*(2*tC+tD)/2,lcCopperPlateBdry};

// *******************************
//
// Corner 3
//
// *******************************
pc1_3 = newp; Point(pc1_3) = {geo_f_x*a+geo_f_x*m*a, geo_f_y*a+geo_f_y*n*a, tD/2,lcCopperPlateBdry};
pc2_3 = newp; Point(pc2_3) = {geo_f_x*a+geo_f_x*m*a, geo_f_y*a+geo_f_y*n*a, -1*tD/2,lcCopperPlateBdry};
pc3_3 = newp; Point(pc3_3) = {geo_f_x*a+geo_f_x*m*a, geo_f_y*a+geo_f_y*n*a, (2*tC+tD)/2,lcCopperPlateBdry};
pc4_3 = newp; Point(pc4_3) = {geo_f_x*a+geo_f_x*m*a, geo_f_y*a+geo_f_y*n*a, -1*(2*tC+tD)/2,lcCopperPlateBdry};

// *******************************
//
// Corner 4
//
// *******************************
pc1_4 = newp; Point(pc1_4) = {geo_f_x*0+geo_f_x*m*a, geo_f_y*a+geo_f_y*n*a, tD/2,lcCopperPlateBdry};
pc2_4 = newp; Point(pc2_4) = {geo_f_x*0+geo_f_x*m*a, geo_f_y*a+geo_f_y*n*a, -1*tD/2,lcCopperPlateBdry};
pc3_4 = newp; Point(pc3_4) = {geo_f_x*0+geo_f_x*m*a, geo_f_y*a+geo_f_y*n*a, (2*tC+tD)/2,lcCopperPlateBdry};
pc4_4 = newp; Point(pc4_4) = {geo_f_x*0+geo_f_x*m*a, geo_f_y*a+geo_f_y*n*a, -1*(2*tC+tD)/2,lcCopperPlateBdry};

// --------------------------------------------------------------------------

// *******************************************************
//
// Copper planes
//
// *******************************************************

// Points between two half pillars on upper LEM
ptmc_1 = newp; Point(ptmc_1) = {geo_f_x*a/2+geo_f_x*m*a, geo_f_y*0+geo_f_y*n*a, (2*tC+tD)/2, lcCopperPlateBdry};
ptmd_1 = newp; Point(ptmd_1) = {geo_f_x*a/2+geo_f_x*m*a, geo_f_y*0+geo_f_y*n*a, tD/2, lcCopperPlateBdry};

ptmc_2 = newp; Point(ptmc_2) = {geo_f_x*a+geo_f_x*m*a, geo_f_y*a/2+geo_f_y*n*a, (2*tC+tD)/2, lcCopperPlateBdry};
ptmd_2 = newp; Point(ptmd_2) = {geo_f_x*a+geo_f_x*m*a, geo_f_y*a/2+geo_f_y*n*a, tD/2, lcCopperPlateBdry};

ptmc_3 = newp; Point(ptmc_3) = {geo_f_x*a/2+geo_f_x*m*a, geo_f_y*a+geo_f_y*n*a, (2*tC+tD)/2, lcCopperPlateBdry};
ptmd_3 = newp; Point(ptmd_3) = {geo_f_x*a/2+geo_f_x*m*a, geo_f_y*a+geo_f_y*n*a, tD/2, lcCopperPlateBdry};

ptmc_4 = newp; Point(ptmc_4) = {geo_f_x*0+geo_f_x*m*a, geo_f_y*a/2+geo_f_y*n*a, (2*tC+tD)/2, lcCopperPlateBdry};
ptmd_4 = newp; Point(ptmd_4) = {geo_f_x*0+geo_f_x*m*a, geo_f_y*a/2+geo_f_y*n*a, tD/2, lcCopperPlateBdry};

// Top lower boundary
pcptl1 = newp; Point(pcptl1) = {geo_f_x*0+geo_f_x*m*a, geo_f_y*0+geo_f_y*n*a, tD/2,lcCopperPlateBdry};
pcptl2 = newp; Point(pcptl2) = {geo_f_x*a+geo_f_x*m*a, geo_f_y*0+geo_f_y*n*a, tD/2,lcCopperPlateBdry};
pcptl3 = newp; Point(pcptl3) = {geo_f_x*a+geo_f_x*m*a, geo_f_y*a+geo_f_y*n*a, tD/2,lcCopperPlateBdry};
pcptl4 = newp; Point(pcptl4) = {geo_f_x*0+geo_f_x*m*a, geo_f_y*a+geo_f_y*n*a, tD/2,lcCopperPlateBdry};

// Top upper boundary
pcptu1 = newp; Point(pcptu1) = {geo_f_x*0+geo_f_x*m*a, geo_f_y*0+geo_f_y*n*a, (2*tC+tD)/2,lcCopperPlateBdry};
pcptu2 = newp; Point(pcptu2) = {geo_f_x*a+geo_f_x*m*a, geo_f_y*0+geo_f_y*n*a, (2*tC+tD)/2,lcCopperPlateBdry};
pcptu3 = newp; Point(pcptu3) = {geo_f_x*a+geo_f_x*m*a, geo_f_y*a+geo_f_y*n*a, (2*tC+tD)/2,lcCopperPlateBdry};
pcptu4 = newp; Point(pcptu4) = {geo_f_x*0+geo_f_x*m*a, geo_f_y*a+geo_f_y*n*a, (2*tC+tD)/2,lcCopperPlateBdry};

// Border lines
// Upper boundary
lcptub1a = newl; Line(lcptub1a) = {pc3_1,ptmc_1};
lcptub1b = newl; Line(lcptub1b) = {ptmc_1,pc3_2};
lcptub2a = newl; Line(lcptub2a) = {pc3_2,ptmc_2};
lcptub2b = newl; Line(lcptub2b) = {ptmc_2,pc3_3};
lcptub3a = newl; Line(lcptub3a) = {pc3_3,ptmc_3};
lcptub3b = newl; Line(lcptub3b) = {ptmc_3,pc3_4};
lcptub4a = newl; Line(lcptub4a) = {pc3_4,ptmc_4};
lcptub4b = newl; Line(lcptub4b) = {ptmc_4,pc3_1};

// Lower boundary
lcptlb5a = newl; Line(lcptlb5a) = {pc1_1,ptmd_1};
lcptlb5b = newl; Line(lcptlb5b) = {ptmd_1,pc1_2};
lcptlb6a = newl; Line(lcptlb6a) = {pc1_2,ptmd_2};
lcptlb6b = newl; Line(lcptlb6b) = {ptmd_2,pc1_3};
lcptlb7a = newl; Line(lcptlb7a) = {pc1_3,ptmd_3};
lcptlb7b = newl; Line(lcptlb7b) = {ptmd_3,pc1_4};
lcptlb8a = newl; Line(lcptlb8a) = {pc1_4,ptmd_4};
lcptlb8b = newl; Line(lcptlb8b) = {ptmd_4,pc1_1};

// Connect the upper and lower points with lines to form the plate
lcptib9 = newl; Line(lcptib9) = {pc3_1, pc1_1};
lcptib10 = newl; Line(lcptib10) = {pc3_2, pc1_2};
lcptib11 = newl; Line(lcptib11) = {pc3_3, pc1_3};
lcptib12 = newl; Line(lcptib12) = {pc3_4, pc1_4};

// ---------------------------------------------

// Points between two half pillars on lower LEM
pbmd_1 = newp; Point(pbmd_1) = {geo_f_x*a/2+geo_f_x*m*a, geo_f_y*0+geo_f_y*n*a, -1*tD/2, lcCopperPlateBdry};
pbmc_1 = newp; Point(pbmc_1) = {geo_f_x*a/2+geo_f_x*m*a, geo_f_y*0+geo_f_y*n*a, -1*(2*tC+tD)/2, lcCopperPlateBdry};

pbmd_2 = newp; Point(pbmd_2) = {geo_f_x*a+geo_f_x*m*a, geo_f_y*a/2+geo_f_y*n*a, -1*tD/2, lcCopperPlateBdry};
pbmc_2 = newp; Point(pbmc_2) = {geo_f_x*a+geo_f_x*m*a, geo_f_y*a/2+geo_f_y*n*a, -1*(2*tC+tD)/2, lcCopperPlateBdry};

pbmd_3 = newp; Point(pbmd_3) = {geo_f_x*a/2+geo_f_x*m*a, geo_f_y*a+geo_f_y*n*a, -1*tD/2, lcCopperPlateBdry};
pbmc_3 = newp; Point(pbmc_3) = {geo_f_x*a/2+geo_f_x*m*a, geo_f_y*a+geo_f_y*n*a, -1*(2*tC+tD)/2, lcCopperPlateBdry};

pbmd_4 = newp; Point(pbmd_4) = {geo_f_x*0+geo_f_x*m*a, geo_f_y*a/2+geo_f_y*n*a, -1*tD/2, lcCopperPlateBdry};
pbmc_4 = newp; Point(pbmc_4) = {geo_f_x*0+geo_f_x*m*a, geo_f_y*a/2+geo_f_y*n*a, -1*(2*tC+tD)/2, lcCopperPlateBdry};

// Bottom lower boundary
pcpbl1 = newp; Point(pcpbl1) = {geo_f_x*0+geo_f_x*m*a, geo_f_y*0+geo_f_y*n*a, -1*(2*tC+tD)/2,lcCopperPlateBdry};
pcpbl2 = newp; Point(pcpbl2) = {geo_f_x*a+geo_f_x*m*a, geo_f_y*0+geo_f_y*n*a, -1*(2*tC+tD)/2,lcCopperPlateBdry};
pcpbl3 = newp; Point(pcpbl3) = {geo_f_x*a+geo_f_x*m*a, geo_f_y*a+geo_f_y*n*a, -1*(2*tC+tD)/2,lcCopperPlateBdry};
pcpbl4 = newp; Point(pcpbl4) = {geo_f_x*0+geo_f_x*m*a, geo_f_y*a+geo_f_y*n*a, -1*(2*tC+tD)/2,lcCopperPlateBdry};

// Bottom upper boundary
pcpbu1 = newp; Point(pcpbu1) = {geo_f_x*0+geo_f_x*m*a, geo_f_y*0+geo_f_y*n*a, -1*tD/2,lcCopperPlateBdry};
pcpbu2 = newp; Point(pcpbu2) = {geo_f_x*a+geo_f_x*m*a, geo_f_y*0+geo_f_y*n*a, -1*tD/2,lcCopperPlateBdry};
pcpbu3 = newp; Point(pcpbu3) = {geo_f_x*a+geo_f_x*m*a, geo_f_y*a+geo_f_y*n*a, -1*tD/2,lcCopperPlateBdry};
pcpbu4 = newp; Point(pcpbu4) = {geo_f_x*0+geo_f_x*m*a, geo_f_y*a+geo_f_y*n*a, -1*tD/2,lcCopperPlateBdry};

// Border lines
// Upper boundary
lcpbub1a = newl; Line(lcpbub1a) = {pc4_1,pbmc_1};
lcpbub1b = newl; Line(lcpbub1b) = {pbmc_1,pc4_2};
lcpbub2a = newl; Line(lcpbub2a) = {pc4_2,pbmc_2};
lcpbub2b = newl; Line(lcpbub2b) = {pbmc_2,pc4_3};
lcpbub3a = newl; Line(lcpbub3a) = {pc4_3,pbmc_3};
lcpbub3b = newl; Line(lcpbub3b) = {pbmc_3,pc4_4};
lcpbub4a = newl; Line(lcpbub4a) = {pc4_4,pbmc_4};
lcpbub4b = newl; Line(lcpbub4b) = {pbmc_4,pc4_1};

// Lower boundary
lcpblb5a = newl; Line(lcpblb5a) = {pc2_1,pbmd_1};
lcpblb5b = newl; Line(lcpblb5b) = {pbmd_1,pc2_2};
lcpblb6a = newl; Line(lcpblb6a) = {pc2_2,pbmd_2};
lcpblb6b = newl; Line(lcpblb6b) = {pbmd_2,pc2_3};
lcpblb7a = newl; Line(lcpblb7a) = {pc2_3,pbmd_3};
lcpblb7b = newl; Line(lcpblb7b) = {pbmd_3,pc2_4};
lcpblb8a = newl; Line(lcpblb8a) = {pc2_4,pbmd_4};
lcpblb8b = newl; Line(lcpblb8b) = {pbmd_4,pc2_1};

// Connect the upper and lower points with lines to form the plate
lcpbib9 = newl; Line(lcpbib9) = {pc4_1, pc2_1};
lcpbib10 = newl; Line(lcpbib10) = {pc4_2, pc2_2};
lcpbib11 = newl; Line(lcpbib11) = {pc4_3, pc2_3};
lcpbib12 = newl; Line(lcpbib12) = {pc4_4, pc2_4};

// Lines connecting the upper and lower level corners
lcorner1 = newl; Line(lcorner1) = {pc1_1, pc2_1};
lcorner2 = newl; Line(lcorner2) = {pc1_2, pc2_2};
lcorner3 = newl; Line(lcorner3) = {pc1_3, pc2_3};
lcorner4 = newl; Line(lcorner4) = {pc1_4, pc2_4};

// Lines splitting the LEM in half
lmid1_1 = newl; Line(lmid1_1) = {ptmc_1, ptmd_1};
lmid1_2 = newl; Line(lmid1_2) = {ptmd_1, pbmd_1};
lmid1_3 = newl; Line(lmid1_3) = {pbmd_1, pbmc_1};

lmid2_1 = newl; Line(lmid2_1) = {ptmc_2, ptmd_2};
lmid2_2 = newl; Line(lmid2_2) = {ptmd_2, pbmd_2};
lmid2_3 = newl; Line(lmid2_3) = {pbmd_2, pbmc_2};

lmid3_1 = newl; Line(lmid3_1) = {ptmc_3, ptmd_3};
lmid3_2 = newl; Line(lmid3_2) = {ptmd_3, pbmd_3};
lmid3_3 = newl; Line(lmid3_3) = {pbmd_3, pbmc_3};

lmid4_1 = newl; Line(lmid4_1) = {ptmc_4, ptmd_4};
lmid4_2 = newl; Line(lmid4_2) = {ptmd_4, pbmd_4};
lmid4_3 = newl; Line(lmid4_3) = {pbmd_4, pbmc_4};

// **********************************************
//
// External Electrodes
//
// **********************************************

// Top electrode
pexet1 = newp; Point(pexet1) = {geo_f_x*0+geo_f_x*m*a, geo_f_y*0+geo_f_y*n*a, (2*tC+tD)/2+lE,lcExtElectrodeBdry};
pexet2 = newp; Point(pexet2) = {geo_f_x*a/2+geo_f_x*m*a, geo_f_y*0+geo_f_y*n*a, (2*tC+tD)/2+lE,lcExtElectrodeBdry};
pexet3 = newp; Point(pexet3) = {geo_f_x*a+geo_f_x*m*a, geo_f_y*0+geo_f_y*n*a, (2*tC+tD)/2+lE,lcExtElectrodeBdry};
pexet4 = newp; Point(pexet4) = {geo_f_x*a+geo_f_x*m*a, geo_f_y*a+geo_f_y*n*a, (2*tC+tD)/2+lE,lcExtElectrodeBdry};
pexet5 = newp; Point(pexet5) = {geo_f_x*a/2+geo_f_x*m*a, geo_f_y*a+geo_f_y*n*a, (2*tC+tD)/2+lE,lcExtElectrodeBdry};
pexet6 = newp; Point(pexet6) = {geo_f_x*0+geo_f_x*m*a, geo_f_y*a+geo_f_y*n*a, (2*tC+tD)/2+lE,lcExtElectrodeBdry};

// Top electrode lines
lexet1 = newl; Line(lexet1) = {pexet1, pexet2};
lexet2 = newl; Line(lexet2) = {pexet2, pexet3};
lexet3 = newl; Line(lexet3) = {pexet3, pexet4};
lexet4 = newl; Line(lexet4) = {pexet4, pexet5};
lexet5 = newl; Line(lexet5) = {pexet5, pexet6};
lexet6 = newl; Line(lexet6) = {pexet6, pexet1};

// *************************************************
//
// Define surfaces
//
// *************************************************

// Copper plate surfaces

llcp_up_border1 = newreg; Line Loop(llcp_up_border1) = {lcptib9, lcptlb5a, lcptlb5b, -lcptib10, -lcptub1a, -lcptub1b};
pscp_up_border1 = newreg; Plane Surface(pscp_up_border1) = {llcp_up_border1};
llcp_up_border2 = newreg; Line Loop(llcp_up_border2) = {lcptib10, lcptlb6a, lcptlb6b, -lcptib11, -lcptub2a, -lcptub2b};

pscp_up_border2 = newreg; Plane Surface(pscp_up_border2) = {llcp_up_border2};
llcp_up_border3 = newreg; Line Loop(llcp_up_border3) = {lcptib11, lcptlb7a, lcptlb7b, -lcptib12, -lcptub3a, -lcptub3b};
pscp_up_border3 = newreg; Plane Surface(pscp_up_border3) = {llcp_up_border3};
llcp_up_border4 = newreg; Line Loop(llcp_up_border4) = {lcptib12, lcptlb8a, lcptlb8b, -lcptib9, -lcptub4a, -lcptub4b};
pscp_up_border4 = newreg; Plane Surface(pscp_up_border4) = {llcp_up_border4};

llcp_low_border1 = newreg; Line Loop(llcp_low_border1) = {lcpbib9, lcpblb5a, lcpblb5b, -lcpbib10, -lcpbub1a, -lcpbub1b};
pscp_low_border1 = newreg; Plane Surface(pscp_low_border1) = {llcp_low_border1}; 
llcp_low_border2 = newreg; Line Loop(llcp_low_border2) = {lcpbib10, lcpblb6a, lcpblb6b, -lcpbib11, -lcpbub2a, -lcpbub2b};
pscp_low_border2 = newreg; Plane Surface(pscp_low_border2) = {llcp_low_border2};
llcp_low_border3 = newreg; Line Loop(llcp_low_border3) = {lcpbib11, lcpblb7a, lcpblb7b, -lcpbib12, -lcpbub3a, -lcpbub3b};
pscp_low_border3 = newreg; Plane Surface(pscp_low_border3) = {llcp_low_border3};
llcp_low_border4 = newreg; Line Loop(llcp_low_border4) = {lcpbib12, lcpblb8a, lcpblb8b, -lcpbib9, -lcpbub4a, -lcpbub4b};
pscp_low_border4 = newreg; Plane Surface(pscp_low_border4) = {llcp_low_border4};

llcp_face1 = newreg; Line Loop(llcp_face1) = {lcptub1a, lcptub1b, lcptub2a, lcptub2b, lcptub3a, lcptub3b, lcptub4a, lcptub4b};
pscp_face1 = newreg; Plane Surface(pscp_face1) = {llcp_face1};

llcp_face3 = newreg; Line Loop(llcp_face3) = {lcpbub1a, lcpbub1b, lcpbub2a, lcpbub2b, lcpbub3a, lcpbub3b, lcpbub4a, lcpbub4b};
pscp_face2 = newreg; Plane Surface(pscp_face2) = {llcp_face3};

// Gas & Dielectric surfaces

ll_top_gas1 = newreg; Line Loop(ll_top_gas1) = {lcptlb5a, lcptlb5b, lcptlb6a, lcptlb6b, lcptlb7a, lcptlb7b, lcptlb8a, lcptlb8b};
ps_top_gas = newreg; Plane Surface(ps_top_gas) = {ll_top_gas1};

ll_bottom_gas1 = newreg; Line Loop(ll_bottom_gas1) = {lcpblb5a, lcpblb5b, lcpblb6a, lcpblb6b, lcpblb7a, lcpblb7b, lcpblb8a, lcpblb8b};
ps_bottom_gas = newreg; Plane Surface(ps_bottom_gas) = {ll_bottom_gas1};

ll_side_gas1a = newreg; Line Loop(ll_side_gas1a) = {lcptlb5a, lmid1_2, -lcpblb5a, -lcorner1};
ps_side_gas1a = newreg; Plane Surface(ps_side_gas1a) = {ll_side_gas1a};

ll_side_gas2a = newreg; Line Loop(ll_side_gas2a) = {lcptlb6a, lmid2_2, -lcpblb6a,  -lcorner2};
ps_side_gas2a = newreg; Plane Surface(ps_side_gas2a) = {ll_side_gas2a};

ll_side_gas3a = newreg; Line Loop(ll_side_gas3a) = {lcptlb7a, lmid3_2, -lcpblb7a, -lcorner3};
ps_side_gas3a = newreg; Plane Surface(ps_side_gas3a) = {ll_side_gas3a};

ll_side_gas4a = newreg; Line Loop(ll_side_gas4a) = {lcptlb8a, lmid4_2, -lcpblb8a, -lcorner4};
ps_side_gas4a = newreg; Plane Surface(ps_side_gas4a) = {ll_side_gas4a};

ll_side_gas1b = newreg; Line Loop(ll_side_gas1b) = {lcptlb5b, lcorner2, -lcpblb5b, -lmid1_2};
ps_side_gas1b = newreg; Plane Surface(ps_side_gas1b) = {ll_side_gas1b};

ll_side_gas2b = newreg; Line Loop(ll_side_gas2b) = {lcptlb6b, lcorner3, -lcpblb6b, -lmid2_2};
ps_side_gas2b = newreg; Plane Surface(ps_side_gas2b) = {ll_side_gas2b};

ll_side_gas3b = newreg; Line Loop(ll_side_gas3b) = {lcptlb7b, lcorner4, -lcpblb7b, -lmid3_2};
ps_side_gas3b = newreg; Plane Surface(ps_side_gas3b) = {ll_side_gas3b};

ll_side_gas4b = newreg; Line Loop(ll_side_gas4b) = {lcptlb8b, lcorner1, -lcpblb8b, -lmid4_2};
ps_side_gas4b = newreg; Plane Surface(ps_side_gas4b) = {ll_side_gas4b};

Coherence;
*/

/*
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////// Bounding and intersecting surfaces Definition

// SURFACES

// *******************************
//
// Face physsurf_bdh_1_1 (Corner 1 to Corner 2) - sbb_1_1 / sbb_1_2 / sbc_1_1 / sbc_1_2
//
// *******************************

// Face physsurf_bdh_1_1 (Corner 1 - Corner 2)

// psbcfmin_2_1_2, psacfmax_2_2_2


l1bdh_1_1_bsurft1 = newl; Line(l1bdh_1_1_bsurft1) = {pexet3, psacfmaxz_2_2_2[]};
l2bdh_1_1_bsurft1 = newl; Line(l2bdh_1_1_bsurft1) = {pexet1, pb2_1_1};
l1bdh_1_1_bsurfb1 = newl; Line(l1bdh_1_1_bsurfb1) = {pc3_2, psbcfminz_2_1_2[]};
l2bdh_1_1_bsurfb1 = newl; Line(l2bdh_1_1_bsurfb1) = {pc3_1, pb1_1_1};


llbdh_1_1_bsurf1t = newreg; Line Loop(llbdh_1_1_bsurf1t) = {lexet1, lexet2, l1bdh_1_1_bsurft1,  -lb2_2, -lb1_2, -ll_bdhbbt_1_1_1[0], -ll_bdhbbt_1_1_2[0], -l2bdh_1_1_bsurft1};
llbdh_1_1_bsurf1b = newreg; Line Loop(llbdh_1_1_bsurf1b) = {lcptub1a, lcptub1b, l1bdh_1_1_bsurfb1, ll_bdhabc0_1_2_2[0], ll_bdhabc1_1_2_2[0], -ll_bdhbbb_1_1_1[0], -ll_bdhbbb_1_1_2[0], -l2bdh_1_1_bsurfb1};

psbdh_1_1_bsurf1t = newreg; Plane Surface(psbdh_1_1_bsurf1t) = {llbdh_1_1_bsurf1t};
psbdh_1_1_bsurf1b = newreg; Plane Surface(psbdh_1_1_bsurf1b) = {llbdh_1_1_bsurf1b};

// *******************************
//
// Face physsurf_bdh_1_2 (Corner 2 to Corner 3) - sab_2_1 / sab_2_2 / sac_2_1 / sac_2_2
//
// *******************************

// Face physsurf_bdh_1_2 (Corner 2 - Corner 3)

// psaaf1_1_1max, psbcf2_1_2min

l1bdh_1_2_bsurft1 = newl; Line(l1bdh_1_2_bsurft1) = {pexet4, pa2_1_1};
l1bdh_1_2_bsurfb1 = newl; Line(l1bdh_1_2_bsurfb1) = {pc3_3, pa1_1_1};

llbdh_1_2_bsurf3t = newreg; Line Loop(llbdh_1_2_bsurf3t) = {lexet3, l1bdh_1_2_bsurft1, ll_bdhabt_1_2_2[0], ll_bdhabt_1_2_1[0], ll_bdhbbc0_1_1_2[0], ll_bdhbbc1_1_1_2[0], -l1bdh_1_1_bsurft1};
llbdh_1_2_bsurf3b = newreg; Line Loop(llbdh_1_2_bsurf3b) = {lcptub2a, lcptub2b, l1bdh_1_2_bsurfb1, -la4_1, -la3_1, ll_bdhabb_1_2_2[0], ll_bdhabb_1_2_1[0], -l1bdh_1_1_bsurfb1};

psbdh_1_2_bsurf3t = newreg; Plane Surface(psbdh_1_2_bsurf3t) = {llbdh_1_2_bsurf3t};
psbdh_1_2_bsurf3b = newreg; Plane Surface(psbdh_1_2_bsurf3b) = {llbdh_1_2_bsurf3b};

// *******************************
//
// Face physsurf_bdh_2_1 (Corner 3 to Corner 4) - sab_1_1[] / sab_1_2[] / sac_1_1[] / sac_1_2[] 
//
// *******************************

// Face physsurf_bdh_2_1 (Corner 3 to Corner 4)

// psaaf1_2_1min, psbcf2_2_2max

l1bdh_2_1_bsurft4 = newl; Line(l1bdh_2_1_bsurft4) = {pexet6, pbdhabt_2_1_2[1]};
l1bdh_2_1_bsurfb4 = newl; Line(l1bdh_2_1_bsurfb4) = {pc3_4, pbdhbbb_2_2_2[1]};

llbdh_2_1_bsurf4t = newreg; Line Loop(llbdh_2_1_bsurf4t) = {lexet4, lexet5, l1bdh_2_1_bsurft4, -ll_bdhabt_2_1_1[0], -ll_bdhabt_2_1_2[0], la3_2, la4_2, -l1bdh_1_2_bsurft1};
llbdh_2_1_bsurf4b = newreg; Line Loop(llbdh_2_1_bsurf4b) = {lcptub3a, lcptub3b, l1bdh_2_1_bsurfb4, -ll_bdhbbc1_2_2_2[0], -ll_bdhbbc0_2_2_2[0], -ll_bdhabb_2_1_1[0], -ll_bdhabb_2_1_2[0], -l1bdh_1_2_bsurfb1}; 

psbdh_2_1_bsurf4t = newreg; Plane Surface(psbdh_2_1_bsurf4t) = {llbdh_2_1_bsurf4t};
psbdh_2_1_bsurf4b = newreg; Plane Surface(psbdh_2_1_bsurf4b) = {llbdh_2_1_bsurf4b};

// *******************************
//
// Face physsurf_bdh_2_2 (Corner 4 to Corner 1) - sbb_2_1[] / sbb_2_2[] / sbc_2_1[] / sbc_2_2[]
//
// *******************************

// Face physsurf_bdh_2_2 (Corner 4 to Corner 1)

// psacf2_1_2min, psbaf1_1_1max

llbdh_2_2_bsurf6t = newreg; Line Loop(llbdh_2_2_bsurf6t) = {lexet6, l2bdh_1_1_bsurft1, ll_bdhbbt_2_2_2[0], ll_bdhbbt_2_2_1[0], -ll_bdhabc1_2_1_2[0], -ll_bdhabc0_2_1_2[0], -l1bdh_2_1_bsurft4};
llbdh_2_2_bsurf6b = newreg; Line Loop(llbdh_2_2_bsurf6b) = {lcptub4a, lcptub4b, l2bdh_1_1_bsurfb1, lb1_1, lb2_1, ll_bdhbbb_2_2_2[0], ll_bdhbbb_2_2_1[0], -l1bdh_2_1_bsurfb4};

psbdh_2_2_bsurf6t = newreg; Plane Surface(psbdh_2_2_bsurf6t) = {llbdh_2_2_bsurf6t};
psbdh_2_2_bsurf6b = newreg; Plane Surface(psbdh_2_2_bsurf6b) = {llbdh_2_2_bsurf6b};

// Bounding surfaces

ll_bsurf7 = newreg; Line Loop(ll_bsurf7) = {lexet1, lexet2, lexet3, lexet4, lexet5, lexet6};
ps_bsurf7 = newreg; Plane Surface(ps_bsurf7) = {ll_bsurf7};

// Volumes
sl_gas = newreg; Surface Loop(sl_gas) = {ps_side_gas1a, ps_side_gas2a, ps_side_gas3a, ps_side_gas4a, ps_side_gas1b, ps_side_gas2b, ps_side_gas3b, ps_side_gas4b, psbdh_1_1_bsurf1t, psbdh_1_1_bsurf1b, psbdh_1_2_bsurf3t, psbdh_1_2_bsurf3b, psbdh_2_1_bsurf4t, psbdh_2_1_bsurf4b, psbdh_2_2_bsurf6t, psbdh_2_2_bsurf6b, ps_bsurf7, pscp_up_border1, pscp_up_border2, pscp_up_border3, pscp_up_border4, ps_bottom_gas, sa_1_1[0], sa_1_1[1], sa_1_2[0], sa_1_2[1], sb_1_1[0], sb_1_1[1], sb_1_2[0], sb_1_2[1], sa_2_1[0], sa_2_1[1], sa_2_2[0], sa_2_2[1], sb_2_1[0], sb_2_1[1], sb_2_2[0], sb_2_2[1]};
vol_gas = newreg; Volume(vol_gas) = {sl_gas};

sl_lower_cp = newreg; Surface Loop(sl_lower_cp) = {pscp_face2, pscp_low_border1, pscp_low_border2, pscp_low_border3, pscp_low_border4, ps_bottom_gas};
vol_lower_cp = newreg; Volume(vol_lower_cp) = {sl_lower_cp};

sl_wire = newreg; Surface Loop(sl_wire) = {sbaf1_1_1[], sbb_1_1[], sbc_1_1[], sbb_1_2[], sbc_1_2[], sbcf2_1_2[], saaf1_2_1[], sab_2_1[], sac_2_1[], sab_2_2[], sac_2_2[], sacf2_2_2[], saaf1_1_1[], sab_1_1[], sac_1_1[], sab_1_2[], sac_1_2[], sacf2_1_2[], sbaf1_2_1[], sbb_2_1[], sbc_2_1[], sbb_2_2[], sbc_2_2[], sbcf2_2_2[]};
vol_wire = newreg; Volume(vol_wire) = {sl_wire};

// Physical volumes

// Volume 1
physvol_gas = newreg; Physical Volume(physvol_gas) = {vol_gas};

// Volume 2
physvol_lower_cp = newreg; Physical Volume(physvol_lower_cp) = {vol_lower_cp};

// Volume 3
physvol_wire = newreg; Physical Volume(physvol_wire) = {vol_wire};

// Physical surfaces

// Gas physical surface
// Bounding surface 1
physsurf_gas = newreg; Physical Surface(physsurf_gas) = {ps_side_gas1a, ps_side_gas2a, ps_side_gas3a, ps_side_gas4a, ps_side_gas1b, ps_side_gas2b, ps_side_gas3b, ps_side_gas4b, psbdh_1_1_bsurf1t, psbdh_1_1_bsurf1b, psbdh_1_2_bsurf3t, psbdh_1_2_bsurf3b, psbdh_2_1_bsurf4t, psbdh_2_1_bsurf4b, psbdh_2_2_bsurf6t, psbdh_2_2_bsurf6b, ps_bsurf7, pscp_up_border1, pscp_up_border2, pscp_up_border3, pscp_up_border4, ps_bottom_gas, sa_1_1[0], sa_1_1[1], sa_1_2[0], sa_1_2[1], sb_1_1[0], sb_1_1[1], sb_1_2[0], sb_1_2[1], sa_2_1[0], sa_2_1[1], sa_2_2[0], sa_2_2[1], sb_2_1[0], sb_2_1[1], sb_2_2[0], sb_2_2[1]};

// Recombine Surface{physsurf_gas};

// Surfaces to which voltages will be applied
// Bounding surface 2
physsurf_lower_cp = newreg; Physical Surface(physsurf_lower_cp) = {pscp_face2, pscp_low_border1, pscp_low_border2, pscp_low_border3, pscp_low_border4, ps_bottom_gas};

// Recombine Surface{physsurf_lower_cp};

// Bounding surface 3
physsurf_wire = newreg; Physical Surface(physsurf_wire) = {sbaf1_1_1[], sbb_1_1[], sbc_1_1[], sbb_1_2[], sbc_1_2[], sbcf2_1_2[], saaf1_2_1[], sab_2_1[], sac_2_1[], sab_2_2[], sac_2_2[], sacf2_2_2[], saaf1_1_1[], sab_1_1[], sac_1_1[], sab_1_2[], sac_1_2[], sacf2_1_2[], sbaf1_2_1[], sbb_2_1[], sbc_2_1[], sbb_2_2[], sbc_2_2[], sbcf2_2_2[]};

// Recombine Surface{physsurf_wire};

// Bounding surface 4
physsurf_upper_el = newreg; Physical Surface(physsurf_upper_el) = {ps_bsurf7};

// Recombine Surface{physsurf_gas};

// Surfaces for periodic boundary conditions

// Bounding surface 5
physsurf_bdh_1_1 = newreg; Physical Surface(physsurf_bdh_1_1) = {pscp_up_border1, ps_side_gas1a, ps_side_gas1b, psbdh_1_1_bsurf1t, psbdh_1_1_bsurf1b};

// Bounding surface 6
physsurf_bdh_1_2 = newreg; Physical Surface(physsurf_bdh_1_2) = {pscp_up_border2, ps_side_gas2a, ps_side_gas2b, psbdh_1_2_bsurf3t, psbdh_1_2_bsurf3b};

// Bounding surface 7
physsurf_bdh_2_1 = newreg; Physical Surface(physsurf_bdh_2_1) = {pscp_up_border3, ps_side_gas3a, ps_side_gas3b, psbdh_2_1_bsurf4t, psbdh_2_1_bsurf4b};

// Bounding surface 8
physsurf_bdh_2_2 = newreg; Physical Surface(physsurf_bdh_2_2) = {pscp_up_border4, ps_side_gas4a, ps_side_gas4b, psbdh_2_2_bsurf6t, psbdh_2_2_bsurf6b};
*/


/*
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////// Window Definition

// WINDOWS / PITCH

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

// End