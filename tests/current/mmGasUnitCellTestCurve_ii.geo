a = 0.040;                                       // the "pitch", or distance between GEM pillars, in mm

mwf = 1;                                         // meshWindow_factor

mm = 1;                                          // geometrical scaling
rW = 0.009 * mm;                                 // R of Wiremesh, in microns
hp0 = 0.0225;                                    // half pitch of the window, in mm
hp = 0.0225 * mm - 0*rW/mwf * mm;                // half pitch of the window, in microns

p = hp0;                                         // half pitch of the window, in mm

R = ( p * p + rW * rW )/( (2 * rW) );            // R
alpha = Asin( (p/R) );                           // angle in radians

meshLevel = 0.000;                               // mesh level, in mm
meshWindow = 0.020;                             // mesh window, in mm

x1SPWindFac = 1.00;
x2SPWindFac = 1.00;
y1SPWindFac = 1.00;
y2SPWindFac = 1.00;

spFac1 = p*0.00;
spFac2 = p*0.00;

x1SPWindFac2 = p*0.00;
x2SPWindFac2 = p*0.00;
y1SPWindFac2 = p*0.00;
y2SPWindFac2 = p*0.00;

Rtp = R - R*0.00;
Rtn = R + R*0.00;

h_f = 0*rW;                                      // Heightfactor

geofx = 1;
geofy = 1;

m = 0;
n = 0;

// Characteristic lengths

  lcCopperPlateBdry = 0.005;
  lcExtElectrodeBdry = 0.005;
  lcWireMesh = 0.001;


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////// Face 1a - half wire (y - z) extrude in x direction - Corner 3 to Corner 4
//___________________________________________________________________________________________________________________________________________________
// Wire 1a1

p1a1_0 = newp; Point(p1a1_0) = {p+spFac1,p+spFac1,rW+meshLevel*mm, lcWireMesh * mm};               // centre circle
p1a1_1 = newp; Point(p1a1_1) = {p+spFac1,p+spFac1,0+meshLevel*mm, lcWireMesh * mm};                // bottom circle
// p1a1_2 = newp; Point(p1a1_2) = {p+spFac1,p+spFac1+rW,rW+meshLevel*mm, lcWireMesh * mm};         // right circle
p1a1_3 = newp; Point(p1a1_3) = {p+spFac1,p+spFac1,2*rW+meshLevel*mm, lcWireMesh * mm};             // top circle
p1a1_4 = newp; Point(p1a1_4) = {p+spFac1,p+spFac1-rW,rW+meshLevel*mm, lcWireMesh * mm};            // left circle

l1a1s_2 = newl; Line(l1a1s_2) = {p1a1_1, p1a1_3};
l1a1_3 = newl; Circle(l1a1_3) = {p1a1_3, p1a1_0, p1a1_4};
l1a1_4 = newl; Circle(l1a1_4) = {p1a1_4, p1a1_0, p1a1_1};

ll1a1_1 = newll; Line Loop(ll1a1_1) = {l1a1_3, l1a1_4, l1a1s_2};

s1a1_1 = news; Plane Surface(s1a1_1) = {ll1a1_1};

tmpaa_1_1[] = {};
tmpaa_1_1[] = {s1a1_1};
saaf_1_1_1[] = s1a1_1;

tmpab_1_1[] = Extrude {{x2SPWindFac2,0,0},{0,-1,0},{p+1*spFac1, p+1*spFac1, -Rtn+rW}, alpha} {
  Surface{tmpaa_1_1[0]};
};

sab_1_1[] = {};
sab_1_1[] = tmpab_1_1[{2:4}];


//___________________________________________________________________________________________________________________________________________________
// Wire 1a2


tmpaa_1_2[] = {tmpab_1_1[0]};

tmpab_1_2[] = Extrude {{x2SPWindFac2,0,0},{0,1,0},{-p+1*spFac1, p+1*spFac1, Rtp-rW}, alpha} {
  Surface{tmpaa_1_2[0]};
};

sab_1_2[] = {};
sab_1_2[] = tmpab_1_2[{2:4}];

sabf_1_1_2[] = {};
sabf_1_1_2[] = tmpab_1_2[0];

sl_wire_exterior_surface_1a[] = newreg; Surface Loop(sl_wire_exterior_surface_1a) = { saaf_1_1_1[0], sab_1_1[], sab_1_2[], sabf_1_1_2[0] };
// vol_1a_wire = newv; Compound Volume(vol_1a_wire) = { tmpab_1_1[1], tmpab_1_2[1] };
vol_1a_wire = newreg; Volume(vol_1a_wire) = sl_wire_exterior_surface_1a[];

physvol_1a_wire = newreg; Physical Volume(physvol_1a_wire) = vol_1a_wire;
physsurf_1a_wire = newreg; Physical Surface(physsurf_1a_wire) = { saaf_1_1_1[0], sab_1_1[], sab_1_2[], sabf_1_1_2[0] };


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////// Face 1b - half wire (y - z) extrude in x direction - Corner 1 to Corner 2
//___________________________________________________________________________________________________________________________________________________
// Wire 1b1

p1b1_0 = newp; Point(p1b1_0) = {-p+spFac1,-p+spFac1,rW+meshLevel*mm, lcWireMesh * mm};             // centre circle
p1b1_1 = newp; Point(p1b1_1) = {-p+spFac1,-p+spFac1,0+meshLevel*mm, lcWireMesh * mm};              // bottom circle
p1b1_2 = newp; Point(p1b1_2) = {-p+spFac1,-p+spFac1+rW,rW+meshLevel*mm, lcWireMesh * mm};          // right circle
p1b1_3 = newp; Point(p1b1_3) = {-p+spFac1,-p+spFac1,2*rW+meshLevel*mm, lcWireMesh * mm};           // top circle
// p1b1_4 = newp; Point(p1b1_4) = {-p+spFac1,-p+spFac1-rW,rW+meshLevel*mm, lcWireMesh * mm};       // left circle

l1b1_1 = newc; Circle(l1b1_1) = {p1b1_1, p1b1_0, p1b1_2};
l1b1_2 = newc; Circle(l1b1_2) = {p1b1_2, p1b1_0, p1b1_3};
l1b1s_3 = newl; Line(l1b1s_3) = {p1b1_1, p1b1_3};

ll1b1_1 = newll; Line Loop(ll1b1_1) = {l1b1_1, l1b1_2, -l1b1s_3};

s1b1_1 = news; Plane Surface(s1b1_1) = {ll1b1_1};

tmpba_1_1[] = {s1b1_1};
sbaf_1_1_1[] = tmpba_1_1[0];

tmpbb_1_1[] = Extrude {{x1SPWindFac2,0,0},{0,1,0},{-p+spFac1, -p+spFac1, -Rtn+rW}, alpha} {
  Surface{tmpba_1_1[0]};
};

sbb_1_1[] = {};
sbb_1_1[] += tmpbb_1_1[{2:4}];


//___________________________________________________________________________________________________________________________________________________
// Wire 1b2


tmpba_1_2[] = {tmpbb_1_1[0]};

tmpbb_1_2[] = Extrude {{x1SPWindFac2,0,0},{0,-1,0},{p+spFac1, -p+spFac1, Rtp-rW}, alpha} {
  Surface{tmpba_1_2[0]};
};

sbb_1_2[] = {};
sbb_1_2[] = tmpbb_1_2[{2:4}];

sbbf_2_1_2[] = {};
sbbf_2_1_2[] = {tmpbb_1_2[0]};

sl_wire_exterior_surface_1b[] = newreg; Surface Loop(sl_wire_exterior_surface_1b) = { sbaf_1_1_1[0], sbb_1_1[], sbb_1_2[], sbbf_2_1_2[0] };
// vol_1b_wire = newreg; Compound Volume(vol_1b_wire) = { tmpbb_1_1[1], tmpbb_1_2[1] };
vol_1b_wire = newreg; Volume(vol_1b_wire) = sl_wire_exterior_surface_1b[];

physvol_1b_wire = newreg; Physical Volume(physvol_1b_wire) = vol_1b_wire;
physsurf_1b_wire = newreg; Physical Surface(physsurf_1b_wire) = { sbaf_1_1_1[0], sbb_1_1[], sbb_1_2[], sbbf_2_1_2[0] };


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////// Face 2a - half wire (x - z) extrude in y direction - Corner 3 to Corner 2
//___________________________________________________________________________________________________________________________________________________
// Wire 2a1

p2a1_0 = newp; Point(p2a1_0) = {p+spFac2,p+spFac2,-rW+meshLevel*mm, lcWireMesh * mm};              // centre circle
p2a1_1 = newp; Point(p2a1_1) = {p+spFac2,p+spFac2,-2*rW+meshLevel*mm, lcWireMesh * mm};            // bottom circle
// p2a1_2 = newp; Point(p2a1_2) = {p+spFac2+rW,p+spFac2,-rW+meshLevel*mm, lcWireMesh * mm};        // right circle
p2a1_3 = newp; Point(p2a1_3) = {p+spFac2,p+spFac2,0+meshLevel*mm, lcWireMesh * mm};                // top circle
p2a1_4 = newp; Point(p2a1_4) = {p+spFac2-rW,p+spFac2,-rW+meshLevel*mm, lcWireMesh * mm};           // left circle

l2a1s_1 = newl; Line(l2a1s_1) = {p2a1_1, p2a1_3};
l2a1_2 = newl; Circle(l2a1_2) = {p2a1_3, p2a1_0, p2a1_4};
l2a1_3 = newl; Circle(l2a1_3) = {p2a1_4, p2a1_0, p2a1_1};

ll2a1_1 = newll; Line Loop(ll2a1_1) = {l2a1_2, l2a1_3, l2a1s_1};

s2a1_1 = news; Plane Surface(s2a1_1) = {ll2a1_1};

tmpaa_2_1[] = {s2a1_1};
saaf_2_2_1[] = tmpaa_2_1[0];


tmpab_2_1[] = Extrude {{0,y1SPWindFac2,0},{-1,0,0},{p+spFac2, p+spFac2, Rtp-rW}, alpha} {
  Surface{tmpaa_2_1[0]};
};

sab_2_1[] = {};
sab_2_1[] += tmpab_2_1[{2:4}];


//___________________________________________________________________________________________________________________________________________________
// Wire 2a2


tmpaa_2_2[] = {tmpab_2_1[0]};

tmpab_2_2[] = Extrude {{0,y1SPWindFac2,0},{1,0,0},{p+spFac2, -p+spFac2, -Rtn+rW}, alpha} {
  Surface{tmpaa_2_2[0]};
};

sab_2_2[] = {};
sab_2_2[] = tmpab_2_2[{2:4}];

sabf_2_2_2[] = {};
sabf_2_2_2[] = {tmpab_2_2[0]};


sl_wire_exterior_surface_2a[] = newreg; Surface Loop(sl_wire_exterior_surface_2a) = { saaf_2_2_1[0], sab_2_1[], sab_2_2[], sabf_2_2_2[0] };
// vol_2a_wire1 = newreg; Volume(vol_2a_wire1) = { tmpab_2_1[1] };
// vol_2a_wire2 = newreg; Volume(vol_2a_wire2) = { tmpab_2_2[1] };
vol_2a_wire = newreg; Volume(vol_2a_wire) = sl_wire_exterior_surface_2a[];

physvol_2a_wire = newreg; Physical Volume(physvol_2a_wire) = vol_2a_wire;
// physvol_2a_wire1 = newreg; Physical Volume(physvol_2a_wire1) = vol_2a_wire1;
// physvol_2a_wire2 = newreg; Physical Volume(physvol_2a_wire2) = vol_2a_wire2;
physsurf_2a_wire = newreg; Physical Surface(physsurf_2a_wire) = { saaf_2_2_1[0], sab_2_1[], sab_2_2[], sabf_2_2_2[0] };


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////// Face 2b - half wire (x - z) extrude in y direction - Corner 1 to Corner 4
//___________________________________________________________________________________________________________________________________________________
// Wire 2b1


p2b1_0 = newp; Point(p2b1_0) = {-p+spFac2,-p+spFac2,-rW+meshLevel*mm, lcWireMesh * mm};        // centre circle
p2b1_1 = newp; Point(p2b1_1) = {-p+spFac2,-p+spFac2,-2*rW+meshLevel*mm, lcWireMesh * mm};      // bottom circle
p2b1_2 = newp; Point(p2b1_2) = {-p+spFac2+rW,-p+spFac2,-rW+meshLevel*mm, lcWireMesh * mm};     // right circle
p2b1_3 = newp; Point(p2b1_3) = {-p+spFac2,-p+spFac2,0+meshLevel*mm, lcWireMesh * mm};          // top circle
// p2b1_4 = newp; Point(p2b1_4) = {-p+spFac2-rW,-p+spFac2,-rW+meshLevel*mm, lcWireMesh * mm};  // left circle

l2b1_1 = newl; Circle(l2b1_1) = {p2b1_1, p2b1_0, p2b1_2};
l2b1_2 = newl; Circle(l2b1_2) = {p2b1_2, p2b1_0, p2b1_3};
l2b1s_3 = newl; Line(l2b1s_3) = {p2b1_1, p2b1_3};

ll2b1_1 = newll; Line Loop(ll2b1_1) = {l2b1_1, l2b1_2, -l2b1s_3};

s2b1_1 = news; Plane Surface(s2b1_1) = {ll2b1_1};

tmpba_2_1[] = {s2b1_1};
sbaf_2_2_1[] = tmpba_2_1[0];


tmpbb_2_1[] = Extrude {{0,y2SPWindFac2,0},{1,0,0},{-p+spFac2, -p+spFac2, Rtp-rW}, alpha} {
  Surface{tmpba_2_1[0]};
};

sbb_2_1[] = {};
sbb_2_1[] += tmpbb_2_1[{2:4}];


//___________________________________________________________________________________________________________________________________________________
// Wire 2b2


tmpba_2_2[] = {tmpbb_2_1[0]};

tmpbb_2_2[] = Extrude {{0,y2SPWindFac2,0},{-1,0,0},{-p+spFac2,p+spFac2,-Rtn+rW}, alpha} {
  Surface{tmpba_2_2[0]};
};

sbb_2_2[] = {};
sbb_2_2[] = tmpbb_2_2[{2:4}];


sbbf_2_2_2[] = {};
sbbf_2_2_2[] = {tmpbb_2_2[0]};


sl_wire_exterior_surface_2b[] = newreg; Surface Loop(sl_wire_exterior_surface_2b) = { sbaf_2_2_1[0], sbb_2_1[], sbb_2_2[], sbbf_2_2_2[0] };
// vol_2b_wire1 = newreg; Volume(vol_2b_wire1) = { tmpbb_2_1[1] };
// vol_2b_wire2 = newreg; Volume(vol_2b_wire2) = { tmpbb_2_2[1] };
vol_2b_wire = newreg; Volume(vol_2b_wire) = sl_wire_exterior_surface_2b[];

physvol_2b_wire = newreg; Physical Volume(physvol_2b_wire) = vol_2b_wire;
// physvol_2b_wire1 = newreg; Physical Volume(physvol_2b_wire1) = vol_2b_wire1;
// physvol_2b_wire2 = newreg; Physical Volume(physvol_2b_wire2) = vol_2b_wire2;
physsurf_2b_wire = newreg; Physical Surface(physsurf_2b_wire) = { sbaf_2_2_1[0], sbb_2_1[], sbb_2_2[], sbbf_2_2_2[0] };


/*
// SHELL

// --------------------------------------------------------------------------

// *******************************
// Corner 1
// *******************************
pc1_1 = newp; Point(pc1_1) = {geofx*-meshWindow/2+geofx*m*-meshWindow/2+1*spFac1, geofy*-meshWindow/2+geofy*n*-meshWindow/2+1*spFac1, 0,lcCopperPlateBdry};


// *******************************
// Corner 2
// *******************************
pc1_2 = newp; Point(pc1_2) = {geofx*meshWindow/2+geofx*m*meshWindow/2+1*spFac1, geofy*-meshWindow/2+geofy*n*-meshWindow/2+1*spFac1, 0,lcCopperPlateBdry};


// *******************************
// Corner 3
// *******************************
pc1_3 = newp; Point(pc1_3) = {geofx*meshWindow/2+geofx*m*meshWindow/2+1*spFac1, geofy*meshWindow/2+geofy*n*meshWindow/2+1*spFac1, 0,lcCopperPlateBdry};


// *******************************
// Corner 4
// *******************************
pc1_4 = newp; Point(pc1_4) = {geofx*-meshWindow/2+geofx*m*-meshWindow/2+1*spFac1, geofy*meshWindow/2+geofy*n*meshWindow/2+1*spFac1, 0,lcCopperPlateBdry};


// UPPER SQUARE

// *******************************
// Corner 1
// *******************************
ptR1_0 = newp; Point(ptR1_0) = {-p+1*spFac1, -p+1*spFac1, R-rW, lcCopperPlateBdry};


// *******************************
// Corner 1
// *******************************
ptpR1_1 = newp; Point(ptpR1_1) = {-p, -p+1*spFac1, R-rW, lcCopperPlateBdry};


// *******************************
// Corner 2
// *******************************
ptpR1_2 = newp; Point(ptpR1_2) = {p+1*spFac1, -p, R-rW, lcCopperPlateBdry};


// *******************************
// Corner 3
// *******************************
ptpR1_3 = newp; Point(ptpR1_3) = {p+1*spFac1, p+1*spFac1, R-rW, lcCopperPlateBdry};


// *******************************
// Corner 4
// *******************************
ptpR1_4 = newp; Point(ptpR1_4) = {-p+1*spFac1, p+1*spFac1, R-rW, lcCopperPlateBdry};

// UPPER SQUARE

// *******************************
// Corner 1
// *******************************
ptnR1_1 = newp; Point(ptnR1_1) = {-p+1*spFac1, -p, -R+rW, lcCopperPlateBdry};


// *******************************
// Corner 2
// *******************************
ptnR1_2 = newp; Point(ptnR1_2) = {p+1*spFac1, -p+1*spFac1, -R+rW, lcCopperPlateBdry};


// *******************************
// Corner 3
// *******************************
ptnR1_3 = newp; Point(ptnR1_3) = {p+1*spFac1, p+1*spFac1, -R+rW, lcCopperPlateBdry};


// *******************************
// Corner 4
// *******************************
ptnR1_4 = newp; Point(ptnR1_4) = {-p, p+1*spFac1, -R+rW, lcCopperPlateBdry};
*/