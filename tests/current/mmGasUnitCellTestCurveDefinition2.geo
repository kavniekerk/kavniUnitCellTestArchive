// *********************************************************************
// mmGasRotationTest.geo
//
// Description:
// Geometry file for a MM cell.
// This cell can be repeated any number of times within Garfield 
// to construct an arbitrarily large MM
//
// *********************************************************************  


tC = 0.0035;                                     // copper thickness, in mm
tD = 0.0035;                                     // dielectric thickness, in mm
lE = 0.40;                                       // distance from GEM plates to upper exterior electrode, in mm
lP = 0.02;                                       // distance from lower LEM plate to pad (readout) plane, in mm
a = 0.045;                                       // the "pitch", or distance between GEM pillars, in mm

mwf = 1;                                         // meshWindow_factor

mm = 1;                                          // geometrical scaling
rW = 0.009 * mm;                                 // R of Wiremesh, in microns
hp0 = 0.0225;                                    // half pitch of the window, in mm
hp = 0.010 * mm - 0*rW/mwf * mm;                 // half pitch of the window, in microns

p = hp0;                                         // half pitch of the window, in mm

R = (p * p + rW * rW)/( (2 * rW) );              // R
alpha = Asin( (p/R) );                           // angle in radians

meshLevel = 0.10;                                // mesh level, in mm
meshWindow = -0.0135;                            // mesh window, in mm

spFacR1 = 0.00;
spFacR2 = 0.00;

x1SPWindFac = p*1.00;
x2SPWindFac = p*1.00;
y1SPWindFac = p*1.00;
y2SPWindFac = p*1.00;

spFac1 = p*0.000;
spFac2 = p*0.000;

spFacStr1 = p*1.00;
spFacStr2 = p*1.00;
spFacStr3 = p*0.00;
spFacStr4 = p*0.00;

x1SPWindFac2 = p*0.00;
x2SPWindFac2 = p*0.00;
y1SPWindFac2 = p*0.00;
y2SPWindFac2 = p*0.00;

Rtp = R + R*0.0;
Rtn = R - R*0.0;

hF = 0*rW;                                       // height factor

geofx = 1;
geofy = 1;

m = 0;
n = 0;

// Characteristic lengths

  lcCopperPlateBdry = 0.0001;
  lcExtElectrodeBdry = 0.0001;
  lcWireMesh = 0.001;


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////// Face 1b - half wire (y - z) extrude in x direction - Corner 1 to Corner 2
//___________________________________________________________________________________________________________________________________________________
// Wire 1b1


p1b1_0 = newp; Point(p1b1_0) = {-p+p+spFac1,-p+p+spFac1,rW+meshLevel*mm, lcWireMesh * mm};                                  // centre circle
p1b1_1 = newp; Point(p1b1_1) = {-p+p+spFac1,-p+p+spFac1,0+meshLevel*mm, lcWireMesh * mm};                                   // bottom circle
p1b1_2 = newp; Point(p1b1_2) = {-p+p+spFac1,-p+p+spFac1+rW,rW+meshLevel*mm, lcWireMesh * mm};                               // right circle
p1b1_3 = newp; Point(p1b1_3) = {-p+p+spFac1,-p+p+spFac1,2*rW+meshLevel*mm, lcWireMesh * mm};                                // top circle
// p1b1_4 = newp; Point(p1b1_4) = {-p+p+spFac1,-p+p+spFac1-rW,rW+meshLevel*mm, lcWireMesh * mm};                            // left circle

l1b1_1 = newc; Circle(l1b1_1) = {p1b1_1, p1b1_0, p1b1_2};
l1b1_2 = newc; Circle(l1b1_2) = {p1b1_2, p1b1_0, p1b1_3};
l1b1s_3 = newl; Line(l1b1s_3) = {p1b1_1, p1b1_3};

ll1b1_1 = newll; Line Loop(ll1b1_1) = {l1b1_1, l1b1_2, -l1b1s_3};

s1b1_1 = news; Plane Surface(s1b1_1) = {ll1b1_1};

tmpba_1_1[] = {s1b1_1};
sbaf_1_1_1[] = tmpba_1_1[0];


ptsbb_1_1 = newp; Point(ptsbb_1_1) = { -p+p +1.00*spFac1, -p+p +1.00*spFac1, -Rtn+rW +meshLevel*mm };

tmpbb_1_1[] = Extrude {{x1SPWindFac2,0,0}, {0,1,0}, { -p+p +1.00*spFac1, -p+p +1.00*spFac1, -Rtn+rW +meshLevel*mm }, 1.00*alpha } {
  Surface{tmpba_1_1[0]};
};

sbb_1_1[] = {};
sbb_1_1[] += tmpbb_1_1[{2:4}];


//___________________________________________________________________________________________________________________________________________________
// Wire 1b2


ptsbb_1_2 = newp; Point(ptsbb_1_2) = { p+p +1.00*spFac1, -p+p +1.00*spFac1, Rtp-rW +meshLevel*mm };

tmpbb_1_2[] = Extrude {{x2SPWindFac2,0,0}, {0,-1,0}, { p+p +1.00*spFac1, -p+p +1.00*spFac1, Rtp-rW +meshLevel*mm }, 1.00*alpha} {
  Surface{tmpbb_1_1[0]};
}; 

sbb_1_2[] = {};
sbb_1_2[] += tmpbb_1_2[{2:4}];

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


p2a1_0 = newp; Point(p2a1_0) = {p+p+spFac2,p+p+spFac2,-rW+meshLevel*mm, lcWireMesh * mm};                                   // centre circle
p2a1_1 = newp; Point(p2a1_1) = {p+p+spFac2,p+p+spFac2,-2*rW+meshLevel*mm, lcWireMesh * mm};                                 // bottom circle
// p2a1_2 = newp; Point(p2a1_2) = {p+p+spFac2+rW,p+p+spFac2,-rW+meshLevel*mm, lcWireMesh * mm};                             // right circle
p2a1_3 = newp; Point(p2a1_3) = {p+p+spFac2,p+p+spFac2,0+meshLevel*mm, lcWireMesh * mm};                                     // top circle
p2a1_4 = newp; Point(p2a1_4) = {p+p+spFac2-rW,p+p+spFac2,-rW+meshLevel*mm, lcWireMesh * mm};                                // left circle

l2a1s_1 = newl; Line(l2a1s_1) = {p2a1_1, p2a1_3};
l2a1_2 = newl; Circle(l2a1_2) = {p2a1_3, p2a1_0, p2a1_4};
l2a1_3 = newl; Circle(l2a1_3) = {p2a1_4, p2a1_0, p2a1_1};

ll2a1_1 = newll; Line Loop(ll2a1_1) = {l2a1_2, l2a1_3, l2a1s_1};

s2a1_1 = news; Plane Surface(s2a1_1) = {ll2a1_1};

tmpaa_2_1[] = {s2a1_1};
saaf_2_2_1[] = tmpaa_2_1[0];


ptsab_2_1 = newp; Point(ptsab_2_1) = { p+p +1.00*spFac2, p+p +1.00*spFac2, Rtp-rW +meshLevel*mm };

tmpab_2_1[] = Extrude {{y1SPWindFac2,0,0}, {-1,0,0}, { p+p +1.00*spFac2, p+p +1.00*spFac2, Rtp-rW +meshLevel*mm }, 1.00*alpha } {
  Surface{tmpaa_2_1[0]};
};

sab_2_1[] = {};
sab_2_1[] += tmpab_2_1[{2:4}];


//___________________________________________________________________________________________________________________________________________________
// Wire 2a2


ptsab_2_2 = newp; Point(ptsab_2_2) = { p+p +1.00*spFac2, -p+p +1.00*spFac2, -Rtn+rW +meshLevel*mm };

tmpab_2_2[] = Extrude {{y2SPWindFac2,0,0}, {1,0,0}, { p+p +1.00*spFac2, -p+p +1.00*spFac2, -Rtn+rW +meshLevel*mm }, 1.00*alpha} {
  Surface{tmpab_2_1[0]};
};

sab_2_2[] = {};
sab_2_2[] += tmpab_2_2[{2:4}];


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
/////// Face 1a - half wire (y - z) extrude in x direction - Corner 3 to Corner 4
//___________________________________________________________________________________________________________________________________________________
// Wire 1a1


p1a1_0 = newp; Point(p1a1_0) = {p+p+spFac1,p+p+spFac1,rW+meshLevel*mm, lcWireMesh * mm};                                    // centre circle
p1a1_1 = newp; Point(p1a1_1) = {p+p+spFac1,p+p+spFac1,0+meshLevel*mm, lcWireMesh * mm};                                     // bottom circle
// p1a1_2 = newp; Point(p1a1_2) = {p+p+spFac1,p+p+spFac1+rW,rW+meshLevel*mm, lcWireMesh * mm};                              // right circle
p1a1_3 = newp; Point(p1a1_3) = {p+p+spFac1,p+p+spFac1,2*rW+meshLevel*mm, lcWireMesh * mm};                                  // top circle
p1a1_4 = newp; Point(p1a1_4) = {p+p+spFac1,p+p+spFac1-rW,rW+meshLevel*mm, lcWireMesh * mm};                                 // left circle

l1a1s_1 = newl; Line(l1a1s_1) = {p1a1_1, p1a1_3};
l1a1_2 = newl; Circle(l1a1_2) = {p1a1_3, p1a1_0, p1a1_4};
l1a1_3 = newl; Circle(l1a1_3) = {p1a1_4, p1a1_0, p1a1_1};

ll1a1_1 = newll; Line Loop(ll1a1_1) = {l1a1_2, l1a1_3, l1a1s_1};

s1a1_1 = news; Plane Surface(s1a1_1) = {ll1a1_1};

tmpaa_1_1[] = {s1a1_1};
saaf_1_1_1[] = tmpaa_1_1[0];


ptsab_1_1 = newp; Point(ptsab_1_1) = { p+p +1.00*spFac1, p+p +1.00*spFac1, -Rtn+rW +meshLevel*mm };

tmpab_1_1[] = Extrude {{x1SPWindFac2,0,0}, {0,-1,0}, { p+p +1.00*spFac1, p+p +1.00*spFac1, -Rtn+rW +meshLevel*mm }, 1.00*alpha } {
  Surface{tmpaa_1_1[0]};
};

sab_1_1[] = {};
sab_1_1[] += tmpab_1_1[{2:4}];


//___________________________________________________________________________________________________________________________________________________
// Wire 1a2


ptsab_1_2 = newp; Point(ptsab_1_2) = { -p+p +1.00*spFac1, p+p +1.00*spFac1, Rtp-rW +meshLevel*mm };

tmpab_1_2[] = Extrude {{x2SPWindFac2,0,0}, {0,1,0}, { -p+p +1.00*spFac1, p+p +1.00*spFac1, Rtp-rW +meshLevel*mm }, 1.00*alpha} {
  Surface{tmpab_1_1[0]};
};

sab_1_2[] = {};
sab_1_2[] += tmpab_1_2[{2:4}];


sabf_2_1_2[] = {};
sabf_2_1_2[] = {tmpab_1_2[0]};

sl_wire_exterior_surface_1a[] = newreg; Surface Loop(sl_wire_exterior_surface_1a) = { saaf_1_1_1[0], sab_1_1[], sab_1_2[], sabf_2_1_2[0] };
// vol_1a_wire = newv; Compound Volume(vol_1a_wire) = { tmpab_1_1[1], tmpab_1_2[1] };
vol_1a_wire = newreg; Volume(vol_1a_wire) = sl_wire_exterior_surface_1a[];

physvol_1a_wire = newreg; Physical Volume(physvol_1a_wire) = vol_1a_wire;
physsurf_1a_wire = newreg; Physical Surface(physsurf_1a_wire) = { saaf_1_1_1[0], sab_1_1[], sab_1_2[], sabf_2_1_2[0] };


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////// Face 2b - half wire (x - z) extrude in y direction - Corner 1 to Corner 4
//___________________________________________________________________________________________________________________________________________________
// Wire 2b1


p2b1_0 = newp; Point(p2b1_0) = {-p+p+spFac2,-p+p+spFac2,-rW+meshLevel*mm, lcWireMesh * mm};                                 // centre circle
p2b1_1 = newp; Point(p2b1_1) = {-p+p+spFac2,-p+p+spFac2,-2*rW+meshLevel*mm, lcWireMesh * mm};                               // bottom circle
p2b1_2 = newp; Point(p2b1_2) = {-p+p+spFac2+rW,-p+p+spFac2,-rW+meshLevel*mm, lcWireMesh * mm};                              // right circle
p2b1_3 = newp; Point(p2b1_3) = {-p+p+spFac2,-p+p+spFac2,0+meshLevel*mm, lcWireMesh * mm};                                   // top circle
// p2b1_4 = newp; Point(p2b1_4) = {-p+p+spFac2-rW,-p+p+spFac2,-rW+meshLevel*mm, lcWireMesh * mm};                           // left circle

l2b1_1 = newl; Circle(l2b1_1) = {p2b1_1, p2b1_0, p2b1_2};
l2b1_2 = newl; Circle(l2b1_2) = {p2b1_2, p2b1_0, p2b1_3};
l2b1s_3 = newl; Line(l2b1s_3) = {p2b1_1, p2b1_3};

ll2b1_1 = newll; Line Loop(ll2b1_1) = {l2b1_1, l2b1_2, -l2b1s_3};

s2b1_1 = news; Plane Surface(s2b1_1) = {ll2b1_1};

tmpba_2_1[] = {s2b1_1};
sbaf_2_2_1[] = tmpba_2_1[0];


ptsab_2_1 = newp; Point(ptsab_2_1) = { -p+p +1.00*spFac2, -p+p +1.00*spFac2, Rtp-rW +meshLevel*mm };

tmpbb_2_1[] = Extrude {{y1SPWindFac2,0,0}, {1,0,0}, { -p+p +1.00*spFac2, -p+p +1.00*spFac2, Rtp-rW +meshLevel*mm }, 1.00*alpha } {
  Surface{tmpba_2_1[0]};
};

sbb_2_1[] = {};
sbb_2_1[] += tmpbb_2_1[{2:4}];


//___________________________________________________________________________________________________________________________________________________
// Wire 2b2


ptsab_2_2 = newp; Point(ptsab_2_2) = { -p+p +1.00*spFac2, p+p +1.00*spFac2, -Rtn+rW +meshLevel*mm };

tmpbb_2_2[] = Extrude {{y2SPWindFac2,0,0}, {-1,0,0}, { -p+p +1.00*spFac2, p+p +1.00*spFac2, -Rtn+rW +meshLevel*mm }, 1.00*alpha} {
  Surface{tmpbb_2_1[0]};
};

sbb_2_2[] = {};
sbb_2_2[] += tmpbb_2_2[{2:4}];


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


    EndFor
  EndFor
EndFor


For q In {2:2}
  For r In {1:2}
    For s In {2:2}

    // *******************************
    //
    // End faces - sabf_2_1_2, sabf_2_2_2, sbbf_2_1_2, sbbf_2_2_2
    //
    // *******************************

    If(q == 2 && s == 2) 

      ll_facesabb~{q}~{r}~{s}[] = {};

        ll_facesabb~{q}~{r}~{s}[] += Boundary{ Surface{sabf~{q}~{r}~{s}[]}; };

        ptsabfmin~{q}~{r}~{s}() = PointsOf{ Line{Abs(ll_facesabb~{q}~{r}~{s}[0])}; };
        psabfminchz~{q}~{r}~{s}[] = Point{ ptsabfmin~{q}~{r}~{s}[0] };

      For w In {1:1}

        pfacesabbc1chz~{q}~{r}~{s}[] = Point{ ptsabfmin~{q}~{r}~{s}[w] };

        If( pfacesabbc1chz~{q}~{r}~{s}[2] < psabfminchz~{q}~{r}~{s}[2] )

          psabfminz[] = pfacesabbc1chz~{q}~{r}~{s}[2];
          ptsabfminz~{q}~{r}~{s}() = ptsabfmin~{q}~{r}~{s}[w];

        Else

          psabfminz[]  = psabfminchz~{q}~{r}~{s}[2];
          ptsabfminz~{q}~{r}~{s}() = ptsabfmin~{q}~{r}~{s}[0];

        EndIf

      EndFor


      For t In {1:2}

        ptfacesabbc1~{q}~{r}~{s}() = PointsOf{ Line{Abs(ll_facesabb~{q}~{r}~{s}[t])}; };

        For l In {0:1}

          pfacesabbc1chz~{q}~{r}~{s}[] = Point{ ptfacesabbc1~{q}~{r}~{s}[l] };

          If( pfacesabbc1chz~{q}~{r}~{s}[2] < psabfminchz~{q}~{r}~{s}[2] )

            psabfminz[]  = pfacesabbc1chz~{q}~{r}~{s}[2];
            ptsabfminz~{q}~{r}~{s}() = ptfacesabbc1~{q}~{r}~{s}[l];

          EndIf

        EndFor

      EndFor


      ll_facesabb~{q}~{r}~{s}[] = {};

        ll_facesabb~{q}~{r}~{s}[] += Boundary{ Surface{sabf~{q}~{r}~{s}[]}; };

        ptsabfmin~{q}~{r}~{s}() = PointsOf{ Line{Abs(ll_facesabb~{q}~{r}~{s}[0])}; };
        psabfminchz~{q}~{r}~{s}[] = Point{ ptsabfmin~{q}~{r}~{s}[0] };

      For w In {1:1}

        pfacesabbc1chz~{q}~{r}~{s}[] = Point{ ptsabfmin~{q}~{r}~{s}[w] };

        If( pfacesabbc1chz~{q}~{r}~{s}[2] < psabfminchz~{q}~{r}~{s}[2] )

          psabfminz[]  = pfacesabbc1chz~{q}~{r}~{s}[2];
          ptsabfminz~{q}~{r}~{s}() = ptsabfmin~{q}~{r}~{s}[w];

        Else

          psabfminz[]  = psabfminchz~{q}~{r}~{s}[2];
          ptsabfminz~{q}~{r}~{s}() = ptsabfmin~{q}~{r}~{s}[0];

        EndIf

      EndFor

      
      For t In {1:2}

        ptfacesabbc1~{q}~{r}~{s}() = PointsOf{ Line{Abs(ll_facesabb~{q}~{r}~{s}[t])}; };

        For l In {0:1}

          pfacesabbc1chz~{q}~{r}~{s}[] = Point{ ptfacesabbc1~{q}~{r}~{s}[l] };
    
          If( pfacesabbc1chz~{q}~{r}~{s}[2] < psabfminchz~{q}~{r}~{s}[2] )

            psabfminz[]  = pfacesabbc1chz~{q}~{r}~{s}[2];
            ptsabfminz~{q}~{r}~{s}() = ptfacesabbc1~{q}~{r}~{s}[l];

          EndIf

        EndFor

      EndFor


      ll_facesabb~{q}~{r}~{s}[] = {};

        ll_facesabb~{q}~{r}~{s}[] += Boundary{ Surface{sabf~{q}~{r}~{s}[]}; };

        ptsabfmax~{q}~{r}~{s}() = PointsOf{ Line{Abs(ll_facesabb~{q}~{r}~{s}[0])}; };
        psabfmaxchz~{q}~{r}~{s}[] = Point{ ptsabfmax~{q}~{r}~{s}[0] };

      For w In {1:1}

        pfacesabbc1chz~{q}~{r}~{s}[] = Point{ ptsabfmax~{q}~{r}~{s}[w] };

        If( pfacesabbc1chz~{q}~{r}~{s}[2] > psabfmaxchz~{q}~{r}~{s}[2] )

          psabfmaxz~{q}~{r}~{s}[] = pfacesabbc1chz~{q}~{r}~{s}[2];
          ptsabfmaxz~{q}~{r}~{s}() = ptsabfmax~{q}~{r}~{s}[w];

        Else

          psabfmaxz~{q}~{r}~{s}[] = psabfmaxchz~{q}~{r}~{s}[2];
          ptsabfmaxz~{q}~{r}~{s}() = ptsabfmax~{q}~{r}~{s}[0];

        EndIf

      EndFor

      
      For t In {1:2}

        ptfacesabbc1~{q}~{r}~{s}() = PointsOf{ Line{Abs(ll_facesabb~{q}~{r}~{s}[t])}; };

        For l In {0:1}

          pfacesabbc1chz~{q}~{r}~{s}[] = Point{ ptfacesabbc1~{q}~{r}~{s}[l] };
    
          If( pfacesabbc1chz~{q}~{r}~{s}[2] > psabfmaxchz~{q}~{r}~{s}[2] )

            psabfmaxz~{q}~{r}~{s}[] = pfacesabbc1chz~{q}~{r}~{s}[2];
            ptsabfmaxz~{q}~{r}~{s}() = ptfacesabbc1~{q}~{r}~{s}[l];

          EndIf

        EndFor

      EndFor


      ll_facesabb~{q}~{r}~{s}[] = {};

        ll_facesabb~{q}~{r}~{s}[] += Boundary{ Surface{sabf~{q}~{r}~{s}[]}; };

        ptsabfmax~{q}~{r}~{s}() = PointsOf{ Line{Abs(ll_facesabb~{q}~{r}~{s}[0])}; };
        psabfmaxchz~{q}~{r}~{s}[] = Point{ ptsabfmax~{q}~{r}~{s}[0] };

      For w In {1:1}

        pfacesabbc1chz~{q}~{r}~{s}[] = Point{ ptsabfmax~{q}~{r}~{s}[w] };

        If( pfacesabbc1chz~{q}~{r}~{s}[2] > psabfmaxchz~{q}~{r}~{s}[2] )

          psabfmaxz~{q}~{r}~{s}[] = pfacesabbc1chz~{q}~{r}~{s}[2];
          ptsabfmaxz~{q}~{r}~{s}() = ptsabfmax~{q}~{r}~{s}[w];

        Else

          psabfmaxz~{q}~{r}~{s}[] = psabfmaxchz~{q}~{r}~{s}[2];
          ptsabfmaxz~{q}~{r}~{s}() = ptsabfmax~{q}~{r}~{s}[0];

        EndIf

      EndFor

      
      For t In {1:2}

        ptfacesabbc1~{q}~{r}~{s}() = PointsOf{ Line{Abs(ll_facesabb~{q}~{r}~{s}[t])}; };

        For l In {0:1}

          pfacesabbc1chz~{q}~{r}~{s}[] = Point{ ptfacesabbc1~{q}~{r}~{s}[l] };
    
          If( pfacesabbc1chz~{q}~{r}~{s}[2] > psabfmaxchz~{q}~{r}~{s}[2] )

            psabfmaxz~{q}~{r}~{s}[] = pfacesabbc1chz~{q}~{r}~{s}[2];
            ptsabfmaxz~{q}~{r}~{s}() = ptfacesabbc1~{q}~{r}~{s}[l];

          EndIf

        EndFor

      EndFor
    
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


      ll_facesbbb~{q}~{r}~{s}[] = {};

        ll_facesbbb~{q}~{r}~{s}[] += Boundary{ Surface{sbbf~{q}~{r}~{s}[]}; };

        ptsbbfmin~{q}~{r}~{s}() = PointsOf{ Line{Abs(ll_facesbbb~{q}~{r}~{s}[0])}; };
        psbbfminchz~{q}~{r}~{s}[] = Point{ ptsbbfmin~{q}~{r}~{s}[0] };

      For w In {1:1}

        pfacesbbbc1chz~{q}~{r}~{s}[] = Point{ ptsbbfmin~{q}~{r}~{s}[w] };

        If( pfacesbbbc1chz~{q}~{r}~{s}[2] < psbbfminchz~{q}~{r}~{s}[2] )

          psbbfminz~{q}~{r}~{s}[] = pfacesbbbc1chz~{q}~{r}~{s}[2];
          ptsbbfminz~{q}~{r}~{s}() = ptsbbfmin~{q}~{r}~{s}[w];

        Else

          psbbfminz~{q}~{r}~{s}[] = psbbfminchz~{q}~{r}~{s}[2];
          ptsbbfminz~{q}~{r}~{s}() = ptsbbfmin~{q}~{r}~{s}[0];

        EndIf

      EndFor

      
      For t In {1:2}

        ptfacesbbbc1~{q}~{r}~{s}() = PointsOf{ Line{Abs(ll_facesbbb~{q}~{r}~{s}[t])}; };

        For l In {0:1}

          pfacesbbbc1chz~{q}~{r}~{s}[] = Point{ ptfacesbbbc1~{q}~{r}~{s}[l] };
    
          If( pfacesbbbc1chz~{q}~{r}~{s}[2] < psbbfminchz~{q}~{r}~{s}[2] )

            psbbfminz~{q}~{r}~{s}[] = pfacesbbbc1chz~{q}~{r}~{s}[2];
            ptsbbfminz~{q}~{r}~{s}() = ptfacesbbbc1~{q}~{r}~{s}[l];

          EndIf

        EndFor

      EndFor


      ll_facesbbb~{q}~{r}~{s}[] = {};

        ll_facesbbb~{q}~{r}~{s}[] += Boundary{ Surface{sbbf~{q}~{r}~{s}[]}; };

        ptsbbfmin~{q}~{r}~{s}() = PointsOf{ Line{Abs(ll_facesbbb~{q}~{r}~{s}[0])}; };
        psbbfminchz~{q}~{r}~{s}[] = Point{ ptsbbfmin~{q}~{r}~{s}[0] };

      For w In {1:1}

        pfacesbbbc1chz~{q}~{r}~{s}[] = Point{ ptsbbfmin~{q}~{r}~{s}[w] };

        If( pfacesbbbc1chz~{q}~{r}~{s}[2] < psbbfminchz~{q}~{r}~{s}[2] )

          psbbfminz~{q}~{r}~{s}[] = pfacesbbbc1chz~{q}~{r}~{s}[2];
          ptsbbfminz~{q}~{r}~{s}() = ptsbbfmin~{q}~{r}~{s}[w];

        Else

          psbbfminz~{q}~{r}~{s}[] = psbbfminchz~{q}~{r}~{s}[2];
          ptsbbfminz~{q}~{r}~{s}() = ptsbbfmin~{q}~{r}~{s}[0];

        EndIf

      EndFor

    
      For t In {1:2}

        ptfacesbbbc1~{q}~{r}~{s}() = PointsOf{ Line{Abs(ll_facesbbb~{q}~{r}~{s}[t])}; };

        For l In {0:1}

          pfacesbbbc1chz~{q}~{r}~{s}[] = Point{ ptfacesbbbc1~{q}~{r}~{s}[l] };
    
          If( pfacesbbbc1chz~{q}~{r}~{s}[2] < psbbfminchz~{q}~{r}~{s}[2] )

            psbbfminz~{q}~{r}~{s}[] = pfacesbbbc1chz~{q}~{r}~{s}[2];
            ptsbbfminz~{q}~{r}~{s}() = ptfacesbbbc1~{q}~{r}~{s}[l];

          EndIf

        EndFor

      EndFor


      ll_facesbbb~{q}~{r}~{s}[] = {};

        ll_facesbbb~{q}~{r}~{s}[] += Boundary{ Surface{sbbf~{q}~{r}~{s}[]}; };

        ptsbbfmax~{q}~{r}~{s}() = PointsOf{ Line{Abs(ll_facesbbb~{q}~{r}~{s}[0])}; };
        psbbfmaxchz~{q}~{r}~{s}[] = Point{ ptsbbfmax~{q}~{r}~{s}[0] };

      For w In {1:1}

        pfacesbbbc1chz~{q}~{r}~{s}[] = Point{ ptsbbfmax~{q}~{r}~{s}[w] };

        If( pfacesbbbc1chz~{q}~{r}~{s}[2] > psbbfmaxchz~{q}~{r}~{s}[2] )

          psbbfmaxz~{q}~{r}~{s}[] = pfacesbbbc1chz~{q}~{r}~{s}[2];
          ptsbbfmaxz~{q}~{r}~{s}() = ptsbbfmax~{q}~{r}~{s}[w];

        Else

          psbbfmaxz~{q}~{r}~{s}[] = psbbfmaxchz~{q}~{r}~{s}[2];
          ptsbbfmaxz~{q}~{r}~{s}() = ptsbbfmax~{q}~{r}~{s}[0];

        EndIf

      EndFor

      
      For t In {1:2}

        ptfacesbbbc1~{q}~{r}~{s}() = PointsOf{ Line{Abs(ll_facesbbb~{q}~{r}~{s}[t])}; };

        For l In {0:1}

          pfacesbbbc1chz~{q}~{r}~{s}[] = Point{ ptfacesbbbc1~{q}~{r}~{s}[l] };
    
          If( pfacesbbbc1chz~{q}~{r}~{s}[2] > psbbfmaxchz~{q}~{r}~{s}[2] )

            psbbfmaxz~{q}~{r}~{s}[] = pfacesbbbc1chz~{q}~{r}~{s}[2];
            ptsbbfmaxz~{q}~{r}~{s}() = ptfacesbbbc1~{q}~{r}~{s}[l];

          EndIf

        EndFor

      EndFor


      ll_facesbbb~{q}~{r}~{s}[] = {};

        ll_facesbbb~{q}~{r}~{s}[] += Boundary{ Surface{sbbf~{q}~{r}~{s}[]}; };

        ptsbbfmax~{q}~{r}~{s}() = PointsOf{ Line{Abs(ll_facesbbb~{q}~{r}~{s}[0])}; };
        psbbfmaxchz~{q}~{r}~{s}[] = Point{ ptsbbfmax~{q}~{r}~{s}[0] };

      For w In {1:1}

        pfacesbbbc1chz~{q}~{r}~{s}[] = Point{ ptsbbfmax~{q}~{r}~{s}[w] };

        If( pfacesbbbc1chz~{q}~{r}~{s}[2] > psbbfmaxchz~{q}~{r}~{s}[2] )

          psbbfmaxz~{q}~{r}~{s}[] = pfacesbbbc1chz~{q}~{r}~{s}[2];
          ptsbbfmaxz~{q}~{r}~{s}() = ptsbbfmax~{q}~{r}~{s}[w];

        Else

          psbbfmaxz~{q}~{r}~{s}[] = psbbfmaxchz~{q}~{r}~{s}[2];
          ptsbbfmaxz~{q}~{r}~{s}() = ptsbbfmax~{q}~{r}~{s}[0];

        EndIf

      EndFor

      
      For t In {1:2}

        ptfacesbbbc1~{q}~{r}~{s}() = PointsOf{ Line{Abs(ll_facesbbb~{q}~{r}~{s}[t])}; };

        For l In {0:1}

          pfacesbbbc1chz~{q}~{r}~{s}[] = Point{ ptfacesbbbc1~{q}~{r}~{s}[l] };
    
          If( pfacesbbbc1chz~{q}~{r}~{s}[2] > psbbfmaxchz~{q}~{r}~{s}[2] )

            psbbfmaxz~{q}~{r}~{s}[] = pfacesbbbc1chz~{q}~{r}~{s}[2];
            ptsbbfmaxz~{q}~{r}~{s}() = ptfacesbbbc1~{q}~{r}~{s}[l];

          EndIf

        EndFor

      EndFor

    EndIf

    EndFor
  EndFor
EndFor


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////// Window Definition

// WINDOWS / PITCH

// --------------------------------------------------------------------------

// *******************************
// Corner 1
// *******************************
pc1_1 = newp; Point(pc1_1) = { geofx*(+rW-p+p+0*spFac1) + geofx*m*(+rW-p+p+0*spFac1), geofy*(+rW-p+p+0*spFac1) + geofy*n*(+rW-p+p+0*spFac1), 0, lcCopperPlateBdry };


// *******************************
// Corner 2
// *******************************
pc1_2 = newp; Point(pc1_2) = { geofx*(-rW+p+p+0*spFac1) + geofx*m*(-rW+p+p+0*spFac1), geofy*(+rW-p+p+0*spFac1) + geofy*n*(+rW-p+p+0*spFac1), 0, lcCopperPlateBdry };


// *******************************
// Corner 3
// *******************************
pc1_3 = newp; Point(pc1_3) = { geofx*(-rW+p+p+0*spFac1) + geofx*m*(-rW+p+p+0*spFac1), geofy*(-rW+p+p+0*spFac1) + geofy*n*(-rW+p+p+0*spFac1), 0, lcCopperPlateBdry }; 


// *******************************
// Corner 4
// *******************************
pc1_4 = newp; Point(pc1_4) = { geofx*(+rW-p+p+0*spFac1) + geofx*m*(+rW-p+p+0*spFac1), geofy*(-rW+p+p+0*spFac1) + geofy*n*(-rW+p+p+0*spFac1), 0, lcCopperPlateBdry };


// UPPER SQUARE

// *******************************
// Corner 1
// *******************************
ptR1_0 = newp; Point(ptR1_0) = {-p+p+1*spFac1, -p+p+1*spFac1, R-rW, lcCopperPlateBdry};


// *******************************
// Corner 1
// *******************************
ptpR1_1 = newp; Point(ptpR1_1) = {-p+p, -p+p+1*spFac1, R-rW, lcCopperPlateBdry};


// *******************************
// Corner 2
// *******************************
ptpR1_2 = newp; Point(ptpR1_2) = {p+p+1*spFac1, -p+p, R-rW, lcCopperPlateBdry};


// *******************************
// Corner 3
// *******************************
ptpR1_3 = newp; Point(ptpR1_3) = {p+p+1*spFac1, p+p+1*spFac1, R-rW, lcCopperPlateBdry};


// *******************************
// Corner 4
// *******************************
ptpR1_4 = newp; Point(ptpR1_4) = {-p+p+1*spFac1, p+p+1*spFac1, R-rW, lcCopperPlateBdry};


// BOTTOM SQUARE

// *******************************
// Corner 1
// *******************************
ptnR1_1 = newp; Point(ptnR1_1) = {-p+p+1*spFac1, -p+p, -R+rW, lcCopperPlateBdry};


// *******************************
// Corner 2
// *******************************
ptnR1_2 = newp; Point(ptnR1_2) = {p+p+1*spFac1, -p+p+1*spFac1, -R+rW, lcCopperPlateBdry};


// *******************************
// Corner 3
// *******************************
ptnR1_3 = newp; Point(ptnR1_3) = {p+p+1*spFac1, p+p+1*spFac1, -R+rW, lcCopperPlateBdry};


// *******************************
// Corner 4
// *******************************
ptnR1_4 = newp; Point(ptnR1_4) = {-p+p, p+p+1*spFac1, -R+rW, lcCopperPlateBdry};


// End



