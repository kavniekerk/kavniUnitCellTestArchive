Include "mmgas_component_label.pro";
Include "mmgas_meshing_options.pro";

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
/// MMGAS_STR_WIRE GEOMETRY MODULE
//

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
/// GENERAL INFORMATION
//
// mmegas_square_straight_10c.geo
//
// Description
//
// References
//
// See also
//

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
// CONSTANTS

mesh_level = 0.025;                                     // mesh level, in mm
mesh_window = 0.05;                                     // mesh window, in mm

//----------------------------------------------------------
// pillar parameters

r0 = 0.005;                                             // the pillar radius, in mm
pil_f_x = 0;                                            // pillar co-ordinates, multiplication factor in x, 1.25
pil_f_y = 0;                                            // pillar co-ordinates, multiplication factor in y, 1.25
pil_c_x = -0.025*0 - 0/4;                               // pillar co-ordinates, constant factor in x, -0.025
pil_c_y = -0.025*0 - 0/4;                               // pillar co-ordinates, constant factor in y, -0.025

ttD = ( mesh_level - 0.005 );                           // top of dielectric pillar
tteD1 = ( ttD - 0.001 );                                // etched pillar level 1
tteD2 = ( 0.5 * (ttD - tteD1) + tteD1 );                // etched pillar level 2

a = 0.05;                                               // the "pitch", or distance between GEM pillars, in mm

//----------------------------------------------------------
// vertical parameters

r1 = 0.001;                                             // the etching amount (etch radius = r0 + r1), in mm
tlgC = ( 0.0000 ) / 2;                                  // lower electrode ground copper thickness, in mm
tltC = ( 0.0035 ) / 2;                                  // lower electrode top copper thickness, in mm
tM = ( 0.0035 ) / 2;                                    // dielectric thickness, in mm
tubC = ( tltC + tM );                                   // higher electrode bottom copper thickness, in mm
tutC = tubC + ( 0.0035 ) / 2;                           // higher electrode top copper thickness, in mm
lE = 0.5;                                               // distance from GEM plates to upper exterior electrode, in mm
lP = 0.1;                                               // distance from lower LEM plate to pad (readout) plane, in mm

//----------------------------------------------------------
// mesh window and wire parameters

mwf = 1;                                                // mesh_window_factor
mm = 1;                                                 // geometrical scaling
r_w = 0.0025 * mm;                                      // radius of Wiremesh, in microns
p_0 = 0.025;                                            // pitch of the window, in mm
p = 0.025 * mm - 0 * r_w/mwf * mm;                      // pitch of the window, in microns
R = (p * p + r_w * r_w)/(2 * r_w);                      // radius
alpha = Asin((p/R));                                    // angle in radians
Total_Grid_size = (a - 0.01)/2;                         // total grid size, in mm, 0.4

h_f = 0*r_w;                                            // heightfactor

Number_Wires = ( (Total_Grid_size) / (p_0) ) / 2;       // number of wires
Wire_length = ( Total_Grid_size / Number_Wires ) / 2;   // wire length

Number_Units_x = 0;                                     // number of units, 1
Number_Units_y = 0;                                     // number of units, 1

geo_wc_xr = 2*r_w;                                      // y-direction wire in x radial direction
geo_wc_yr = 2*r_w;                                      // x-direction wire in y radial direction

geo_wc_xd = 2*r_w;                                      // x-direction wire in x-direction
geo_wc_yd = 2*r_w;                                      // y-direction wire in y-direction

n_1 = 0;
m_1 = 0;
n_2 = 1;
m_2 = 1;

i_t_x = Number_Wires;
j_t_x = Number_Wires + 1;

i_t_y = Number_Wires;
j_t_y = Number_Wires + 1;

//----------------------------------------------------------
// shell parameters

geo_f_x = 1;                                            // geometric_factor
geo_f_y = 1;                                            // geometric_factor

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
/// GEOMETRY PARAMETERS

//----------------------------------------------------------
// Extrusion Precision

Geometry.ExtrudeSplinePoints = 3;
Geometry.Points = 0;

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
/// MESHING PARAMETERS

//----------------------------------------------------------
// Characteristic lengths - characterization of mesh

// current best dimensions for mesh characteristic lengths

lcDielectricpillar = 0.0025;                            // characterization of dielectric 
lcEtchingpillar = 0.0025;                               // characterization of dielectric etching
lcCopperPlateBdry = 0.0025;                             // characterization of metal surfaces / anode
lcExtElectrodeBdry = 0.025;                             // characterization of external electrode / cathode
lcWireMesh = 0.0025;                                    // characterization of wire electrode

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
/// GEOMETRY MODULE

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
// EXTERIOR SHELL STRUCTURE

//----------------------------------------------------------
// Corner 1

pc1_1 = newp; Point(pc1_1) = {geo_f_x*0+geo_f_x*m_1*a, geo_f_y*0+geo_f_y*n_1*a, tlgC, lcCopperPlateBdry};        // level 1 - bottom lower layer - bottom lower copper electrode
pc2_1 = newp; Point(pc2_1) = {geo_f_x*0+geo_f_x*m_1*a, geo_f_y*0+geo_f_y*n_1*a, tltC, lcCopperPlateBdry};        // level 2 - bottom upper layer
pc3_1 = newp; Point(pc3_1) = {geo_f_x*0+geo_f_x*m_1*a, geo_f_y*0+geo_f_y*n_1*a, tubC, lcCopperPlateBdry};        // level 3 - top lower layer
pc4_1 = newp; Point(pc4_1) = {geo_f_x*0+geo_f_x*m_1*a, geo_f_y*0+geo_f_y*n_1*a, tutC, lcCopperPlateBdry};        // level 4 - top upper layer - top upper copper electrode

//----------------------------------------------------------
// Corner 2

pc1_2 = newp; Point(pc1_2) = {geo_f_x*a+geo_f_x*m_1*a, geo_f_y*0+geo_f_y*n_1*a, tlgC, lcCopperPlateBdry};        // level 1 - bottom lower layer - bottom lower copper electrode
pc2_2 = newp; Point(pc2_2) = {geo_f_x*a+geo_f_x*m_1*a, geo_f_y*0+geo_f_y*n_1*a, tltC, lcCopperPlateBdry};        // level 2 - bottom upper layer
pc3_2 = newp; Point(pc3_2) = {geo_f_x*a+geo_f_x*m_1*a, geo_f_y*0+geo_f_y*n_1*a, tubC, lcCopperPlateBdry};        // level 3 - top lower layer
pc4_2 = newp; Point(pc4_2) = {geo_f_x*a+geo_f_x*m_1*a, geo_f_y*0+geo_f_y*n_1*a, tutC, lcCopperPlateBdry};        // level 4 - top upper layer - top upper copper electrode

//----------------------------------------------------------
// Corner 3

pc1_3 = newp; Point(pc1_3) = {geo_f_x*a+geo_f_x*m_1*a, geo_f_y*a+geo_f_y*n_1*a, tlgC, lcCopperPlateBdry};        // level 1 - bottom lower layer - bottom lower copper electrode
pc2_3 = newp; Point(pc2_3) = {geo_f_x*a+geo_f_x*m_1*a, geo_f_y*a+geo_f_y*n_1*a, tltC, lcCopperPlateBdry};        // level 2 - bottom upper layer
pc3_3 = newp; Point(pc3_3) = {geo_f_x*a+geo_f_x*m_1*a, geo_f_y*a+geo_f_y*n_1*a, tubC, lcCopperPlateBdry};        // level 3 - top lower layer
pc4_3 = newp; Point(pc4_3) = {geo_f_x*a+geo_f_x*m_1*a, geo_f_y*a+geo_f_y*n_1*a, tutC, lcCopperPlateBdry};        // level 4 - top upper layer - top upper copper electrode

//----------------------------------------------------------
// Corner 4

pc1_4 = newp; Point(pc1_4) = {geo_f_x*0+geo_f_x*m_1*a, geo_f_y*a+geo_f_y*n_1*a, tlgC, lcCopperPlateBdry};        // level 1 - bottom lower layer - bottom lower copper electrode
pc2_4 = newp; Point(pc2_4) = {geo_f_x*0+geo_f_x*m_1*a, geo_f_y*a+geo_f_y*n_1*a, tltC, lcCopperPlateBdry};        // level 2 - bottom upper layer
pc3_4 = newp; Point(pc3_4) = {geo_f_x*0+geo_f_x*m_1*a, geo_f_y*a+geo_f_y*n_1*a, tubC, lcCopperPlateBdry};        // level 3 - top lower layer
pc4_4 = newp; Point(pc4_4) = {geo_f_x*0+geo_f_x*m_1*a, geo_f_y*a+geo_f_y*n_1*a, tutC, lcCopperPlateBdry};        // level 4 - top upper layer - top upper copper electrode

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
// COPPER ELECTRODES

//----------------------------------------------------------
// Points between two half pillars on upper LEM

// ptmc_1 = newp; Point(ptmc_1) = {geo_f_x*a/2+geo_f_x*m_1*a, geo_f_y*0+geo_f_y*n_1*a, tubC, lcCopperPlateBdry};
// ptmd_1 = newp; Point(ptmd_1) = {geo_f_x*a/2+geo_f_x*m_1*a, geo_f_y*0+geo_f_y*n_1*a, tlgC, lcCopperPlateBdry};

// ptmc_2 = newp; Point(ptmc_2) = {geo_f_x*a+geo_f_x*m_1*a, geo_f_y*a/2+geo_f_y*n_1*a, tubC, lcCopperPlateBdry};
// ptmd_2 = newp; Point(ptmd_2) = {geo_f_x*a+geo_f_x*m_1*a, geo_f_y*a/2+geo_f_y*n_1*a, tlgC, lcCopperPlateBdry};

// ptmc_3 = newp; Point(ptmc_3) = {geo_f_x*a/2+geo_f_x*m_1*a, geo_f_y*a+geo_f_y*n_1*a, tubC, lcCopperPlateBdry};
// ptmd_3 = newp; Point(ptmd_3) = {geo_f_x*a/2+geo_f_x*m_1*a, geo_f_y*a+geo_f_y*n_1*a, tlgC, lcCopperPlateBdry};

// ptmc_4 = newp; Point(ptmc_4) = {geo_f_x*0+geo_f_x*m_1*a, geo_f_y*a/2+geo_f_y*n_1*a, tubC, lcCopperPlateBdry};
// ptmd_4 = newp; Point(ptmd_4) = {geo_f_x*0+geo_f_x*m_1*a, geo_f_y*a/2+geo_f_y*n_1*a, tlgC, lcCopperPlateBdry};

//----------------------------------------------------------
// Top lower boundary

// pcptl1 = newp; Point(pcptl1) = {geo_f_x*0+geo_f_x*m_1*a, geo_f_y*0+geo_f_y*n_1*a, tltC,lcCopperPlateBdry};
// pcptl2 = newp; Point(pcptl2) = {geo_f_x*a+geo_f_x*m_1*a, geo_f_y*0+geo_f_y*n_1*a, tltC,lcCopperPlateBdry};
// pcptl3 = newp; Point(pcptl3) = {geo_f_x*a+geo_f_x*m_1*a, geo_f_y*a+geo_f_y*n_1*a, tltC,lcCopperPlateBdry};
// pcptl4 = newp; Point(pcptl4) = {geo_f_x*0+geo_f_x*m_1*a, geo_f_y*a+geo_f_y*n_1*a, tltC,lcCopperPlateBdry};

//----------------------------------------------------------
// Top upper boundary

// pcptu1 = newp; Point(pcptu1) = {geo_f_x*0+geo_f_x*m_1*a, geo_f_y*0+geo_f_y*n_1*a, tutC,lcCopperPlateBdry};
// pcptu2 = newp; Point(pcptu2) = {geo_f_x*a+geo_f_x*m_1*a, geo_f_y*0+geo_f_y*n_1*a, tutC,lcCopperPlateBdry};
// pcptu3 = newp; Point(pcptu3) = {geo_f_x*a+geo_f_x*m_1*a, geo_f_y*a+geo_f_y*n_1*a, tutC,lcCopperPlateBdry};
// pcptu4 = newp; Point(pcptu4) = {geo_f_x*0+geo_f_x*m_1*a, geo_f_y*a+geo_f_y*n_1*a, tutC,lcCopperPlateBdry};

//----------------------------------------------------------
// Border lines

//----------------------------------------------------------
// Upper boundary - level 4

// lcptub1a = newl; Line(lcptub1a) = {pc4_1,pc4_2};
// Transfinite Line { lcptub1a } = lcptub1a;
// lcptub2a = newl; Line(lcptub2a) = {pc4_2,pc4_3};
// Transfinite Line { lcptub2a } = lcptub2a;
// lcptub3a = newl; Line(lcptub3a) = {pc4_3,pc4_4};
// Transfinite Line { lcptub3a } = lcptub3a;
// lcptub4a = newl; Line(lcptub4a) = {pc4_4,pc4_1};
// Transfinite Line { lcptub4a } = lcptub4a;

//----------------------------------------------------------
// Lower boundary - level 3

// lcptlb5a = newl; Line(lcptlb5a) = {pc3_1,pc3_2};
// Transfinite Line { lcptlb5a } = lcptlb5a;
// lcptlb6a = newl; Line(lcptlb6a) = {pc3_2,pc3_3};
// Transfinite Line { lcptlb6a } = lcptlb6a;
// lcptlb7a = newl; Line(lcptlb7a) = {pc3_3,pc3_4};
// Transfinite Line { lcptlb7a } = lcptlb7a;
// lcptlb8a = newl; Line(lcptlb8a) = {pc3_4,pc3_1};
// Transfinite Line { lcptlb8a } = lcptlb8a;

//----------------------------------------------------------
// Connect the upper and lower points with lines to form the plate

// lcptib9 = newl; Line(lcptib9) = {pc4_1, pc3_1};
// Transfinite Line { lcptib9 } = lcptib9;
// lcptib10 = newl; Line(lcptib10) = {pc4_2, pc3_2};
// Transfinite Line { lcptib10 } = lcptib10;
// lcptib11 = newl; Line(lcptib11) = {pc4_3, pc3_3};
// Transfinite Line { lcptib11 } = lcptib11;
// lcptib12 = newl; Line(lcptib12) = {pc4_4, pc3_4};
// Transfinite Line { lcptib12 } = lcptib12;

//----------------------------------------------------------
// Points between two half pillars on lower LEM

// pbmd_1 = newp; Point(pbmd_1) = {geo_f_x*a/2+geo_f_x*m_1*a, geo_f_y*0+geo_f_y*n_1*a, tltC, lcCopperPlateBdry};
// pbmc_1 = newp; Point(pbmc_1) = {geo_f_x*a/2+geo_f_x*m_1*a, geo_f_y*0+geo_f_y*n_1*a, tutC, lcCopperPlateBdry};

// pbmd_2 = newp; Point(pbmd_2) = {geo_f_x*a+geo_f_x*m_1*a, geo_f_y*a/2+geo_f_y*n_1*a, tltC, lcCopperPlateBdry};
// pbmc_2 = newp; Point(pbmc_2) = {geo_f_x*a+geo_f_x*m_1*a, geo_f_y*a/2+geo_f_y*n_1*a, tutC, lcCopperPlateBdry};

// pbmd_3 = newp; Point(pbmd_3) = {geo_f_x*a/2+geo_f_x*m_1*a, geo_f_y*a+geo_f_y*n_1*a, tltC, lcCopperPlateBdry};
// pbmc_3 = newp; Point(pbmc_3) = {geo_f_x*a/2+geo_f_x*m_1*a, geo_f_y*a+geo_f_y*n_1*a, tutC, lcCopperPlateBdry};

// pbmd_4 = newp; Point(pbmd_4) = {geo_f_x*0+geo_f_x*m_1*a, geo_f_y*a/2+geo_f_y*n_1*a, tltC, lcCopperPlateBdry};
// pbmc_4 = newp; Point(pbmc_4) = {geo_f_x*0+geo_f_x*m_1*a, geo_f_y*a/2+geo_f_y*n_1*a, tutC, lcCopperPlateBdry};

//----------------------------------------------------------
// Bottom lower boundary

// pcpbl1 = newp; Point(pcpbl1) = {geo_f_x*0+geo_f_x*m_1*a, geo_f_y*0+geo_f_y*n_1*a, tutC,lcCopperPlateBdry};
// pcpbl2 = newp; Point(pcpbl2) = {geo_f_x*a+geo_f_x*m_1*a, geo_f_y*0+geo_f_y*n_1*a, tutC,lcCopperPlateBdry};
// pcpbl3 = newp; Point(pcpbl3) = {geo_f_x*a+geo_f_x*m_1*a, geo_f_y*a+geo_f_y*n_1*a, tutC,lcCopperPlateBdry};
// pcpbl4 = newp; Point(pcpbl4) = {geo_f_x*0+geo_f_x*m_1*a, geo_f_y*a+geo_f_y*n_1*a, tutC,lcCopperPlateBdry};

//----------------------------------------------------------
// Bottom upper boundary

// pcpbu1 = newp; Point(pcpbu1) = {geo_f_x*0+geo_f_x*m_1*a, geo_f_y*0+geo_f_y*n_1*a, tltC,lcCopperPlateBdry};
// pcpbu2 = newp; Point(pcpbu2) = {geo_f_x*a+geo_f_x*m_1*a, geo_f_y*0+geo_f_y*n_1*a, tltC,lcCopperPlateBdry};
// pcpbu3 = newp; Point(pcpbu3) = {geo_f_x*a+geo_f_x*m_1*a, geo_f_y*a+geo_f_y*n_1*a, tltC,lcCopperPlateBdry};
// pcpbu4 = newp; Point(pcpbu4) = {geo_f_x*0+geo_f_x*m_1*a, geo_f_y*a+geo_f_y*n_1*a, tltC,lcCopperPlateBdry};

//----------------------------------------------------------
// Border lines

//----------------------------------------------------------
// Upper boundary - level 2

lcpbub1a = newl; Line(lcpbub1a) = {pc2_1,pc2_2};
// Transfinite Line { lcpbub1a } = lcpbub1a;
lcpbub2a = newl; Line(lcpbub2a) = {pc2_2,pc2_3};
// Transfinite Line { lcpbub2a } = lcpbub2a;
lcpbub3a = newl; Line(lcpbub3a) = {pc2_3,pc2_4};
// Transfinite Line { lcpbub3a } = lcpbub3a;
lcpbub4a = newl; Line(lcpbub4a) = {pc2_4,pc2_1};
// Transfinite Line { lcpbub4a } = lcpbub4a;

//----------------------------------------------------------
// Lower boundary - level 1

lcpblb5a = newl; Line(lcpblb5a) = {pc1_1,pc1_2};
// Transfinite Line { lcpblb5a } = lcpblb5a;
lcpblb6a = newl; Line(lcpblb6a) = {pc1_2,pc1_3};
// Transfinite Line { lcpblb6a } = lcpblb6a;
lcpblb7a = newl; Line(lcpblb7a) = {pc1_3,pc1_4};
// Transfinite Line { lcpblb7a } = lcpblb7a;
lcpblb8a = newl; Line(lcpblb8a) = {pc1_4,pc1_1};
// Transfinite Line { lcpblb8a } = lcpblb8a;

//----------------------------------------------------------
// Connect the upper and lower points with lines to form the plate

lcpbib9 = newl; Line(lcpbib9) = {pc2_1, pc1_1};
// Transfinite Line { lcpbib9 } = lcpbib9;
lcpbib10 = newl; Line(lcpbib10) = {pc2_2, pc1_2};
// Transfinite Line { lcpbib10 } = lcpbib10;
lcpbib11 = newl; Line(lcpbib11) = {pc2_3, pc1_3};
// Transfinite Line { lcpbib11 } = lcpbib11;
lcpbib12 = newl; Line(lcpbib12) = {pc2_4, pc1_4};
// Transfinite Line { lcpbib12 } = lcpbib12;

//----------------------------------------------------------
// Lines connecting the upper and lower level corners

// lcorner1 = newl; Line(lcorner1) = {pc2_1, pc3_1};
// Transfinite Line { lcorner1 } = lcorner1;
// lcorner2 = newl; Line(lcorner2) = {pc2_2, pc3_2};
// Transfinite Line { lcorner2 } = lcorner2;
// lcorner3 = newl; Line(lcorner3) = {pc2_3, pc3_3};
// Transfinite Line { lcorner3 } = lcorner3;
// lcorner4 = newl; Line(lcorner4) = {pc2_4, pc3_4};
// Transfinite Line { lcorner4 } = lcorner4;

//----------------------------------------------------------
// Lines splitting the LEM in half

// lmid1_1 = newl; Line(lmid1_1) = {ptmc_1, ptmd_1};
// lmid1_2 = newl; Line(lmid1_2) = {ptmd_1, pbmd_1};
// lmid1_3 = newl; Line(lmid1_3) = {pbmd_1, pbmc_1};

// lmid2_1 = newl; Line(lmid2_1) = {ptmc_2, ptmd_2};
// lmid2_2 = newl; Line(lmid2_2) = {ptmd_2, pbmd_2};
// lmid2_3 = newl; Line(lmid2_3) = {pbmd_2, pbmc_2};

// lmid3_1 = newl; Line(lmid3_1) = {ptmc_3, ptmd_3};
// lmid3_2 = newl; Line(lmid3_2) = {ptmd_3, pbmd_3};
// lmid3_3 = newl; Line(lmid3_3) = {pbmd_3, pbmc_3};

// lmid4_1 = newl; Line(lmid4_1) = {ptmc_4, ptmd_4};
// lmid4_2 = newl; Line(lmid4_2) = {ptmd_4, pbmd_4};
// lmid4_3 = newl; Line(lmid4_3) = {pbmd_4, pbmc_4};

//----------------------------------------------------------
// Bottom lower copper boundary - level 1

// lcpblb1 = newl; Line(lcpblb1) = {pc1_1,pc1_2};
// Transfinite Line { lcpblb1 } = lcpblb1;
// lcpblb2 = newl; Line(lcpblb2) = {pc1_2,pc1_3};
// Transfinite Line { lcpblb2 } = lcpblb2;
// lcpblb3 = newl; Line(lcpblb3) = {pc1_3,pc1_4};
// Transfinite Line { lcpblb3 } = lcpblb3;
// lcpblb4 = newl; Line(lcpblb4) = {pc1_4,pc1_1};
// Transfinite Line { lcpblb4 } = lcpblb4;

//----------------------------------------------------------
// Top lower copper boundary - level 2

// lcpblt1 = newl; Line(lcpblt1) = {pc2_1,pc2_2};
// Transfinite Line { lcpblt1 } = lcpblt1;
// lcpblt2 = newl; Line(lcpblt2) = {pc2_2,pc2_3};
// Transfinite Line { lcpblt2 } = lcpblt2;
// lcpblt3 = newl; Line(lcpblt3) = {pc2_3,pc2_4};
// Transfinite Line { lcpblt3 } = lcpblt3;
// lcpblt4 = newl; Line(lcpblt4) = {pc2_4,pc2_1};
// Transfinite Line { lcpblt4 } = lcpblt4;

//----------------------------------------------------------
// Bottom upper copper boundary - level 3

// lcpulb1 = newl; Line(lcpulb1) = {pc3_1,pc3_2};
// Transfinite Line { lcpulb1 } = lcpulb1;
// lcpulb2 = newl; Line(lcpulb2) = {pc3_2,pc3_3};
// Transfinite Line { lcpulb2 } = lcpulb2;
// lcpulb3 = newl; Line(lcpulb3) = {pc3_3,pc3_4};
// Transfinite Line { lcpulb3 } = lcpulb3;
// lcpulb4 = newl; Line(lcpulb4) = {pc3_4,pc3_1};
// Transfinite Line { lcpulb4 } = lcpulb4;

//----------------------------------------------------------
// Top upper copper boundary - level 4

// lcpult1 = newl; Line(lcpult1) = {pc4_1,pc4_2};
// Transfinite Line { lcpult1 } = lcpult1;
// lcpult2 = newl; Line(lcpult2) = {pc4_2,pc4_3};
// Transfinite Line { lcpult2 } = lcpult2;
// lcpult3 = newl; Line(lcpult3) = {pc4_3,pc4_4};
// Transfinite Line { lcpult3 } = lcpult3;
// lcpult4 = newl; Line(lcpult4) = {pc4_4,pc4_1};
// Transfinite Line { lcpult4 } = lcpult4;

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
// Wire Mesh

//----------------------------------------------------------
// First set of wires

//----------------------------------------------------------
// x-direction

// Face 1a - half wire (y - z) extrude in x direction - Corner 3 to Corner 4
// Wire 1a1

p0_1a = newp; Point(p0_1a) = {p+p,p+p,-r_w+mesh_level*mm-h_f, lcWireMesh * mm};            // centre circle
p1a_1_1 = newp; Point(p1a_1_1) = {p+p,p+p,-2*r_w+mesh_level*mm-h_f, lcWireMesh * mm};      // bottom circle
// p2_1a = newp; Point(p2_1a) = {p+p,p+p+r_w,-r_w+mesh_level*mm-h_f, lcWireMesh * mm};       // right circle
p1a_3_1 = newp; Point(p1a_3_1) = {p+p,p+p,0+mesh_level*mm-h_f, lcWireMesh * mm};           // top circle
p4_1a = newp; Point(p4_1a) = {p+p,p+p-r_w,-r_w+mesh_level*mm-h_f, lcWireMesh * mm};        // left circle

l2_1as = newl; Line(l2_1as) = {p1a_1_1, p1a_3_1};
l3_1a = newl; Circle(l3_1a) = {p1a_3_1, p0_1a, p4_1a};
l4_1a = newl; Circle(l4_1a) = {p4_1a, p0_1a, p1a_1_1};

ll1_1a = newll; Line Loop(ll1_1a) = {l3_1a, l4_1a, l2_1as};

s_1_1a = news; Plane Surface(s_1_1a) = {ll1_1a};

sa_1_1[] = {};
tmpa_1_1[] = {};
tmpa_1_1[] = {s_1_1a};

tmpa_1_1[] = Extrude {-a/2, 0, 0} { Surface{ s_1_1a }; 
};

sa_1_1[] += tmpa_1_1[{2:4}];

// Wire 1a2

sa_1_2[] = {};
tmpa_1_2[] = {tmpa_1_1[0]};

tmpa_1_2[] = Extrude {-a/2, 0, 0} { Surface{ tmpa_1_2[0] }; 
};

sa_1_2[] += tmpa_1_2[{2:4}];

s_1_1a1[] = s_1_1a;
s_1_1a2[] = tmpa_1_2[0];

/*
sl_wire_exterior_surface_1a[] = newreg; Surface Loop(sl_wire_exterior_surface_1a) = { s_1_1a1[], sa_1_1[], sa_1_2[], s_1_1a2[] };
vol_1a_wire = newreg; Volume(vol_1a_wire) = sl_wire_exterior_surface_1a[];
Physical Volume(physvol_1a_wire) = vol_1a_wire;
Physical Surface(physsurf_1a_wire) = { s_1_1a1[], sa_1_1[], sa_1_2[], s_1_1a2[] };
*/

//----------------------------------------------------------
// Second set of wires

//----------------------------------------------------------
// x-direction

// Face 1b - half wire (y - z) extrude in x direction - Corner 1 to Corner 2
// Wire 1b1

p0_1b = newp; Point(p0_1b) = {-p+p,-p+p,-r_w+mesh_level*mm-h_f, lcWireMesh * mm};
p1b_1_1 = newp; Point(p1b_1_1) = {-p+p,-p+p,-2*r_w+mesh_level*mm-h_f, lcWireMesh * mm};
p2_1b = newp; Point(p2_1b) = {-p+p,r_w - p+p,-r_w+mesh_level*mm-h_f, lcWireMesh * mm};
p1b_3_1 = newp; Point(p1b_3_1) = {-p+p,-p+p,0+mesh_level*mm-h_f, lcWireMesh * mm};
// p4_1b = newp; Point(p4_1b) = {-p+p,-r_w - p+p,-r_w+mesh_level*mm-h_f, lcWireMesh * mm};

l1_1b = newl; Circle(l1_1b) = {p1b_1_1, p0_1b, p2_1b};
l2_1b = newl; Circle(l2_1b) = {p2_1b, p0_1b, p1b_3_1};
l2_1bs = newl; Line(l2_1bs) = {p1b_1_1, p1b_3_1};

ll1_1b = newll; Line Loop(ll1_1b) = {l1_1b, l2_1b, -l2_1bs};

s_1_1b = news; Plane Surface(s_1_1b) = {ll1_1b};

sb_1_1[] = {};
tmpb_1_1[] = {s_1_1b};

tmpb_1_1[] = Extrude {a/2, 0, 0} { Surface{ s_1_1b }; 
};

sb_1_1[] += tmpb_1_1[{2:4}];

// Wire 1b2

sb_1_2[] = {};
tmpb_1_2[] = {tmpb_1_1[0]};

tmpb_1_2[] = Extrude {a/2, 0, 0} { Surface{ tmpb_1_2[0] }; 
};

sb_1_2[] += tmpb_1_2[{2:4}];

s_1_1b1[] = s_1_1b;
s_1_1b2[] = tmpb_1_2[0];

/*
sl_wire_exterior_surface_1b[] = newreg; Surface Loop(sl_wire_exterior_surface_1b) = { s_1_1b1[], sb_1_1[], sb_1_2[], s_1_1b2[] };
vol_1b_wire = newreg; Volume(vol_1b_wire) = sl_wire_exterior_surface_1b[];
Physical Volume(physvol_1b_wire) = vol_1b_wire;
Physical Surface(physsurf_1b_wire) = { s_1_1b1[], sb_1_1[], sb_1_2[], s_1_1b2[] };
*/

//----------------------------------------------------------
// First set of wires

//----------------------------------------------------------
// y-direction

// Face 2a - half wire (x - z) extrude in y direction - Corner 3 to Corner 2
// Wire 2a1

p0_2a = newp; Point(p0_2a) = {p+p,p+p,r_w+mesh_level*mm+h_f, lcWireMesh * mm};
p2a_1_1 = newp; Point(p2a_1_1) = {p+p,p+p,2*r_w+mesh_level*mm+h_f, lcWireMesh * mm};
// p2_2a = newp; Point(p2_2a) = {p+p+r_w,p+p,r_w+mesh_level*mm+h_f, lcWireMesh * mm};
p2a_3_1 = newp; Point(p2a_3_1) = {p+p,p+p,0+mesh_level*mm+h_f, lcWireMesh * mm};
p4_2a = newp; Point(p4_2a) = {p+p-r_w,p+p,r_w+mesh_level*mm+h_f, lcWireMesh * mm};

l2_2as = newl; Line(l2_2as) = {p2a_1_1, p1a_3_1};
l3_2a = newl; Circle(l3_2a) = {p1a_3_1, p0_2a, p4_2a};
l4_2a = newl; Circle(l4_2a) = {p4_2a, p0_2a, p2a_1_1};

ll1_2a = newll; Line Loop(ll1_2a) = {l3_2a, l4_2a, l2_2as};

s_1_2a = news; Plane Surface(s_1_2a) = {ll1_2a};

sa_2_1[] = {};
tmpa_2_1[] = {s_1_2a};

tmpa_2_1[] = Extrude {0, -a/2, 0} { Surface{ s_1_2a }; 
};

sa_2_1[] += tmpa_2_1[{2:4}];

// Wire 2a2

sa_2_2[] = {};
tmpa_2_2[] = {tmpa_2_1[0]};

tmpa_2_2[] = Extrude {0, -a/2, 0} { Surface{ tmpa_2_2[0] }; 
};

sa_2_2[] += tmpa_2_2[{2:4}];

s_1_2a1[] = s_1_2a;
s_1_2a2[] = tmpa_2_2[0];

/*
sl_wire_exterior_surface_2a[] = newreg; Surface Loop(sl_wire_exterior_surface_2a) = { s_1_2a1[], sa_2_1[], sa_2_2[], s_1_2a2[] };
vol_2a_wire = newreg; Volume(vol_2a_wire) = sl_wire_exterior_surface_2a[];
Physical Volume(physvol_2a_wire) = vol_2a_wire;
Physical Surface(physsurf_2a_wire) = { s_1_2a1[], sa_2_1[], sa_2_2[], s_1_2a2[] };
*/

//----------------------------------------------------------
// Second set of wires

//----------------------------------------------------------
// y-direction

// Face 2b - half wire (x - z) extrude in y direction - Corner 1 to Corner 4
// Wire 2b1

p0_2b = newp; Point(p0_2b) = {-p+p,-p+p,r_w+mesh_level*mm+h_f, lcWireMesh * mm};
p2b_1_1 = newp; Point(p2b_1_1) = {-p+p,-p+p,2*r_w+mesh_level*mm+h_f, lcWireMesh * mm};
p2_2b = newp; Point(p2_2b) = {-p+p+r_w,-p+p,r_w+mesh_level*mm+h_f, lcWireMesh * mm};
p2b_3_1 = newp; Point(p2b_3_1) = {-p+p,-p+p,0+mesh_level*mm+h_f, lcWireMesh * mm};
// p4_2b = newp; Point(p4_2b) = {-p+p-r_w,-p+p,r_w+mesh_level*mm+h_f, lcWireMesh * mm};

l1_2b = newl; Circle(l1_2b) = {p2b_1_1, p0_2b, p2_2b};
l2_2b = newl; Circle(l2_2b) = {p2_2b, p0_2b, p1b_3_1};
l2_2bs = newl; Line(l2_2bs) = {p2b_1_1, p1b_3_1};

ll1_2b = newll; Line Loop(ll1_2b) = {l1_2b, l2_2b, -l2_2bs};

s_1_2b = news; Plane Surface(s_1_2b) = {ll1_2b};

sb_2_1[] = {};
tmpb_2_1[] = {s_1_2b};

tmpb_2_1[] = Extrude {0, a/2, 0} { Surface{ s_1_2b }; 
};

sb_2_1[] += tmpb_2_1[{2:4}];

// Wire 2b2

sb_2_2[] = {};
tmpb_2_2[] = {tmpb_2_1[0]};

tmpb_2_2[] = Extrude {0, a/2, 0} { Surface{ tmpb_2_2[] }; 
};

sb_2_2[] += tmpb_2_2[{2:4}];

s_1_2b1[] = s_1_2b;
s_1_2b2[] = tmpb_2_2[0];

/*
sl_wire_exterior_surface_2b[] = newreg; Surface Loop(sl_wire_exterior_surface_2b) = { s_1_2b1[], sb_2_1[], sb_2_2[], s_1_2b2[] };
vol_2b_wire = newreg; Volume(vol_2b_wire) = sl_wire_exterior_surface_2b[];
Physical Volume(physvol_2b_wire) = vol_2b_wire;
Physical Surface(physsurf_2b_wire) = { s_1_2b1[], sb_2_1[], sb_2_2[], s_1_2b2[] };
*/

// Coherence;

//----------------------------------------------------------
// Comparative IF Statement

For q In {1:2}
  For r In {1:2}
    For s In {1:2}

//----------------------------------------------------------
// Face physsurf_bdh_1_1 (Corner 1 - Corner 2)

  If(q == 1&& r == 1) 

  ll_bdhbt~{q}~{r}~{s}[] = {};
  ll_bdhbb~{q}~{r}~{s}[] = {};

  pbdhbt~{q}~{r}~{s}() = {};
  pbdhbb~{q}~{r}~{s}() = {};

  ll_bdhb~{q}~{r}~{s}[] = {};

    ll_bdhb~{q}~{r}~{s}[] += Boundary{ Surface{sb~{r}~{s}[2]}; };

  For t In {0:3}

    If(t == 0 && s == 2)
      ll_bdhbc0~{q}~{r}~{s}[] = {};
      pbdhbc0~{q}~{r}~{s}() = {};
      ll_bdhbc~{q}~{r}~{s}[] = {};

      ll_bdhbc~{q}~{r}~{s}[] += Boundary{ Surface{tmpb~{r}~{s}[t]}; };
      ll_bdhbc0~{q}~{r}~{s}[] = Abs(ll_bdhbc~{q}~{r}~{s}[t]);
      pbdhbc0~{q}~{r}~{s}() = PointsOf{ Line{Abs(ll_bdhbc~{q}~{r}~{s}[t])}; };
    EndIf

    If(t == 1)
      ll_bdhbb~{q}~{r}~{s}[] = Abs(ll_bdhb~{q}~{r}~{s}[t]);
      pbdhbb~{q}~{r}~{s}() = PointsOf{ Line{Abs(ll_bdhb~{q}~{r}~{s}[t])}; };
    EndIf  

    If(t == 1 && s == 2)
      ll_bdhbc1~{q}~{r}~{s}[] = {};
      pbdhbc1~{q}~{r}~{s}() = {};

      ll_bdhbc1~{q}~{r}~{s}[] = Abs(ll_bdhbc~{q}~{r}~{s}[t]);
      pbdhbc1~{q}~{r}~{s}() = PointsOf{ Line{Abs(ll_bdhbc~{q}~{r}~{s}[t])}; };
    EndIf

    If(t == 2 && s == 2)
      ll_bdhb1~{q}~{r}~{s}[] = {};
      pbdhb1~{q}~{r}~{s}() = {};

      ll_bdhb1~{q}~{r}~{s}[] = Abs(ll_bdhb~{q}~{r}~{s}[t]);
      pbdhb1~{q}~{r}~{s}() = PointsOf{ Line{Abs(ll_bdhb~{q}~{r}~{s}[t])}; };
    EndIf

    If(t == 3)  
      ll_bdhbt~{q}~{r}~{s}[] = Abs(ll_bdhb~{q}~{r}~{s}[t]);
      pbdhbt~{q}~{r}~{s}() = PointsOf{ Line{Abs(ll_bdhb~{q}~{r}~{s}[t])}; };
    EndIf

  EndFor

  EndIf

//----------------------------------------------------------
// Face physsurf_bdh_1_2 (Corner 2 - Corner 3)

  If(q == 1&& r == 2)

  ll_bdhat~{q}~{r}~{s}[] = {};
  ll_bdhab~{q}~{r}~{s}[] = {}; 

  pbdhat~{q}~{r}~{s}() = {};
  pbdhab~{q}~{r}~{s}() = {};

  ll_bdha~{q}~{r}~{s}[] = {};

    ll_bdha~{q}~{r}~{s}[] += Boundary{ Surface{sa~{r}~{s}[2]}; };

  For t In {0:3}

    If(t == 0 && s == 2)
      ll_bdhac0~{q}~{r}~{s}[] = {};
      pbdhac0~{q}~{r}~{s}() = {};
      ll_bdhac~{q}~{r}~{s}[] = {};

      ll_bdhac~{q}~{r}~{s}[] += Boundary{ Surface{tmpa~{r}~{s}[t]}; };
      ll_bdhac0~{q}~{r}~{s}[] = Abs(ll_bdhac~{q}~{r}~{s}[t]);
      pbdhac0~{q}~{r}~{s}() = PointsOf{ Line{Abs(ll_bdhac~{q}~{r}~{s}[t])}; };
    EndIf

    If(t == 1)
      ll_bdhab~{q}~{r}~{s}[] = Abs(ll_bdha~{q}~{r}~{s}[t]);
      pbdhab~{q}~{r}~{s}() = PointsOf{ Line{Abs(ll_bdha~{q}~{r}~{s}[t])}; };
    EndIf

    If(t == 1 && s == 2)
      ll_bdhac1~{q}~{r}~{s}[] = {};
      pbdhac1~{q}~{r}~{s}() = {};

      ll_bdhac1~{q}~{r}~{s}[] = Abs(ll_bdhac~{q}~{r}~{s}[t]);
      pbdhac1~{q}~{r}~{s}() = PointsOf{ Line{Abs(ll_bdhac~{q}~{r}~{s}[t])}; };
    EndIf

    If(t == 2 && s == 2)
      ll_bdha1~{q}~{r}~{s}[] = {}; 
      pbdha1~{q}~{r}~{s}() = {};

      ll_bdha1~{q}~{r}~{s}[] = Abs(ll_bdha~{q}~{r}~{s}[t]);
      pbdha1~{q}~{r}~{s}() = PointsOf{ Line{Abs(ll_bdha~{q}~{r}~{s}[t])}; };
    EndIf

    If(t == 3)
      ll_bdhat~{q}~{r}~{s}[] = Abs(ll_bdha~{q}~{r}~{s}[t]);
      pbdhat~{q}~{r}~{s}() = PointsOf{ Line{Abs(ll_bdha~{q}~{r}~{s}[t])}; };
    EndIf

  EndFor

  EndIf

//----------------------------------------------------------
// Face physsurf_bdh_2_1 (Corner 3 - Corner 4)

  If(q == 2&& r == 1)

  ll_bdhat~{q}~{r}~{s}[] = {};
  ll_bdhab~{q}~{r}~{s}[] = {}; 

  pbdhat~{q}~{r}~{s}() = {};
  pbdhab~{q}~{r}~{s}() = {};

  ll_bdha~{q}~{r}~{s}[] = {};

    ll_bdha~{q}~{r}~{s}[] += Boundary{ Surface{sa~{r}~{s}[2]}; };

  For t In {0:3}

    If(t == 0 && s == 2)
      ll_bdhac0~{q}~{r}~{s}[] = {};
      pbdhac0~{q}~{r}~{s}() = {};
      ll_bdhac~{q}~{r}~{s}[] = {};

      ll_bdhac~{q}~{r}~{s}[] += Boundary{ Surface{tmpa~{r}~{s}[t]}; };
      ll_bdhac0~{q}~{r}~{s}[] = Abs(ll_bdhac~{q}~{r}~{s}[t]);
      pbdhac0~{q}~{r}~{s}() = PointsOf{ Line{Abs(ll_bdhac~{q}~{r}~{s}[t])}; };
    EndIf
  
    If(t == 3)
      ll_bdhab~{q}~{r}~{s}[] = Abs(ll_bdha~{q}~{r}~{s}[t]);
      pbdhab~{q}~{r}~{s}() = PointsOf{ Line{Abs(ll_bdha~{q}~{r}~{s}[t])}; };
    EndIf

    If(t == 2 && s == 2)
      ll_bdha1~{q}~{r}~{s}[] = {};
      pbdha1~{q}~{r}~{s}() = {};

      ll_bdha1~{q}~{r}~{s}[] = Abs(ll_bdha~{q}~{r}~{s}[t]);
      pbdha1~{q}~{r}~{s}() = PointsOf{ Line{Abs(ll_bdha~{q}~{r}~{s}[t])}; };
    EndIf

    If(t == 1)
      ll_bdhat~{q}~{r}~{s}[] = Abs(ll_bdha~{q}~{r}~{s}[t]);
      pbdhat~{q}~{r}~{s}() = PointsOf{ Line{Abs(ll_bdha~{q}~{r}~{s}[t])}; };
    EndIf

    If(t == 1 && s == 2)
      ll_bdhac1~{q}~{r}~{s}[] = {};
      pbdhac1~{q}~{r}~{s}() = {};

      ll_bdhac1~{q}~{r}~{s}[] = Abs(ll_bdhac~{q}~{r}~{s}[t]);
      pbdhac1~{q}~{r}~{s}() = PointsOf{ Line{Abs(ll_bdhac~{q}~{r}~{s}[t])}; };
    EndIf

  EndFor

  EndIf

//----------------------------------------------------------
// Face physsurf_bdh_2_2 (Corner 4 - Corner 1)

  If(q == 2&& r == 2) 

  ll_bdhbt~{q}~{r}~{s}[] = {};
  ll_bdhbb~{q}~{r}~{s}[] = {};
  pbdhbt~{q}~{r}~{s}() = {};
  pbdhbb~{q}~{r}~{s}() = {};

  ll_bdhb~{q}~{r}~{s}[] = {};

    ll_bdhb~{q}~{r}~{s}[] += Boundary{ Surface{sb~{r}~{s}[2]}; };

  For t In {0:3}
 
    If(t == 0 && s == 2)
      ll_bdhbc0~{q}~{r}~{s}[] = {};
      pbdhbc0~{q}~{r}~{s}() = {};
      ll_bdhbc~{q}~{r}~{s}[] = {};

      ll_bdhbc~{q}~{r}~{s}[] += Boundary{ Surface{tmpb~{r}~{s}[t]}; };
      ll_bdhbc0~{q}~{r}~{s}[] = Abs(ll_bdhbc~{q}~{r}~{s}[t]);
      pbdhbc0~{q}~{r}~{s}() = PointsOf{ Line{Abs(ll_bdhbc~{q}~{r}~{s}[t])}; };
    EndIf

    If(t == 3)
      ll_bdhbb~{q}~{r}~{s}[] = Abs(ll_bdhb~{q}~{r}~{s}[t]);
      pbdhbb~{q}~{r}~{s}() = PointsOf{ Line{Abs(ll_bdhb~{q}~{r}~{s}[t])}; };
    EndIf

    If(t == 2 && s == 2)
      ll_bdhb1~{q}~{r}~{s}[] = {};
      pbdhb1~{q}~{r}~{s}() = {};

      ll_bdhb1~{q}~{r}~{s}[] = Abs(ll_bdhb~{q}~{r}~{s}[t]);
      pbdhb1~{q}~{r}~{s}() = PointsOf{ Line{Abs(ll_bdhb~{q}~{r}~{s}[t])}; };
    EndIf

    If(t == 1)  
      ll_bdhbt~{q}~{r}~{s}[] = Abs(ll_bdhb~{q}~{r}~{s}[t]);
      pbdhbt~{q}~{r}~{s}() = PointsOf{ Line{Abs(ll_bdhb~{q}~{r}~{s}[t])}; };
    EndIf

    If(t == 1 && s == 2)  
      ll_bdhbc1~{q}~{r}~{s}[] = {};
      pbdhbc1~{q}~{r}~{s}() = {};

      ll_bdhbc1~{q}~{r}~{s}[] = Abs(ll_bdhbc~{q}~{r}~{s}[t]);
      pbdhbc1~{q}~{r}~{s}() = PointsOf{ Line{Abs(ll_bdhbc~{q}~{r}~{s}[t])}; };
    EndIf

  EndFor

  EndIf

  EndFor
 EndFor
EndFor

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
// PILLAR ARCHITECTURE

// For m In {0:Number_Units_x}
//  For n In {0:Number_Units_y}

//----------------------------------------------------------
// pillar 1f (full pillar)

//----------------------------------------------------------
// Center

// pc1_1f~{m}~{n} = newp; Point(pc1_1f~{m}~{n}) = {pil_c_x+geo_f_x*a/2+pil_f_x*geo_f_x*m*a, pil_c_y+geo_f_y*a/2+pil_f_y*geo_f_y*n*a, tlgC,lcDielectricpillar};
// pc2_1f~{m}~{n} = newp; Point(pc2_1f~{m}~{n}) = {pil_c_x+geo_f_x*a/2+pil_f_x*geo_f_x*m*a, pil_c_y+geo_f_y*a/2+pil_f_y*geo_f_y*n*a, tltC,lcDielectricpillar};
// pc3_1f~{m}~{n} = newp; Point(pc3_1f~{m}~{n}) = {pil_c_x+geo_f_x*a/2+pil_f_x*geo_f_x*m*a, pil_c_y+geo_f_y*a/2+pil_f_y*geo_f_y*n*a, tteD1,lcDielectricpillar};
// pc4_1f~{m}~{n} = newp; Point(pc4_1f~{m}~{n}) = {pil_c_x+geo_f_x*a/2+pil_f_x*geo_f_x*m*a, pil_c_y+geo_f_y*a/2+pil_f_y*geo_f_y*n*a, ttD,lcDielectricpillar};

//----------------------------------------------------------
// Dielectric pillar

//----------------------------------------------------------
// Top

// pth1_1f~{m}~{n} = newp; Point(pth1_1f~{m}~{n}) = {pil_c_x+geo_f_x*a/2+pil_f_x*geo_f_x*m*a-1*(r0), pil_c_y+geo_f_y*a/2+pil_f_y*geo_f_y*n*a, ttD,lcDielectricpillar};
// pth2_1f~{m}~{n} = newp; Point(pth2_1f~{m}~{n}) = {pil_c_x+geo_f_x*a/2+pil_f_x*geo_f_x*m*a, pil_c_y+geo_f_y*a/2+pil_f_y*geo_f_y*n*a-1*(r0), ttD,lcDielectricpillar};
// pth3_1f~{m}~{n} = newp; Point(pth3_1f~{m}~{n}) = {pil_c_x+geo_f_x*a/2+pil_f_x*geo_f_x*m*a+1*(r0), pil_c_y+geo_f_y*a/2+pil_f_y*geo_f_y*n*a, ttD,lcDielectricpillar};
// pth4_1f~{m}~{n} = newp; Point(pth4_1f~{m}~{n}) = {pil_c_x+geo_f_x*a/2+pil_f_x*geo_f_x*m*a, pil_c_y+geo_f_y*a/2+pil_f_y*geo_f_y*n*a+1*(r0), ttD,lcDielectricpillar};

// cth1_1f~{m}~{n} = newc; Circle(cth1_1f~{m}~{n}) = {pth1_1f~{m}~{n}, pc4_1f~{m}~{n}, pth2_1f~{m}~{n}};
// cth2_1f~{m}~{n} = newc; Circle(cth2_1f~{m}~{n}) = {pth2_1f~{m}~{n}, pc4_1f~{m}~{n}, pth3_1f~{m}~{n}};
// cth3_1f~{m}~{n} = newc; Circle(cth3_1f~{m}~{n}) = {pth3_1f~{m}~{n}, pc4_1f~{m}~{n}, pth4_1f~{m}~{n}};
// cth4_1f~{m}~{n} = newc; Circle(cth4_1f~{m}~{n}) = {pth4_1f~{m}~{n}, pc4_1f~{m}~{n}, pth1_1f~{m}~{n}};

//----------------------------------------------------------
// Bottom

// pbh1_1f~{m}~{n} = newp; Point(pbh1_1f~{m}~{n}) = {pil_c_x+geo_f_x*a/2+pil_f_x*geo_f_x*m*a-1*(r0), pil_c_y+geo_f_y*a/2+pil_f_y*geo_f_y*n*a, tlgC,lcDielectricpillar};
// pbh2_1f~{m}~{n} = newp; Point(pbh2_1f~{m}~{n}) = {pil_c_x+geo_f_x*a/2+pil_f_x*geo_f_x*m*a, pil_c_y+geo_f_y*a/2+pil_f_y*geo_f_y*n*a-1*(r0), tlgC,lcDielectricpillar};
// pbh3_1f~{m}~{n} = newp; Point(pbh3_1f~{m}~{n}) = {pil_c_x+geo_f_x*a/2+pil_f_x*geo_f_x*m*a+1*(r0), pil_c_y+geo_f_y*a/2+pil_f_y*geo_f_y*n*a, tlgC,lcDielectricpillar};
// pbh4_1f~{m}~{n} = newp; Point(pbh4_1f~{m}~{n}) = {pil_c_x+geo_f_x*a/2+pil_f_x*geo_f_x*m*a, pil_c_y+geo_f_y*a/2+pil_f_y*geo_f_y*n*a+1*(r0), tlgC,lcDielectricpillar};

//----------------------------------------------------------
// Upper Etching

//----------------------------------------------------------
// Top

// ptue1_1f~{m}~{n} = newp; Point(ptue1_1f~{m}~{n}) = {pil_c_x+geo_f_x*a/2+pil_f_x*geo_f_x*m*a-1*(r0+r1), pil_c_y+geo_f_y*a/2+pil_f_y*geo_f_y*n*a, tteD1,lcEtchingpillar};
// ptue2_1f~{m}~{n} = newp; Point(ptue2_1f~{m}~{n}) = {pil_c_x+geo_f_x*a/2+pil_f_x*geo_f_x*m*a, pil_c_y+geo_f_y*a/2+pil_f_y*geo_f_y*n*a-1*(r0+r1), tteD1,lcEtchingpillar};
// ptue3_1f~{m}~{n} = newp; Point(ptue3_1f~{m}~{n}) = {pil_c_x+geo_f_x*a/2+pil_f_x*geo_f_x*m*a+1*(r0+r1), pil_c_y+geo_f_y*a/2+pil_f_y*geo_f_y*n*a, tteD1,lcEtchingpillar};
// ptue4_1f~{m}~{n} = newp; Point(ptue4_1f~{m}~{n}) = {pil_c_x+geo_f_x*a/2+pil_f_x*geo_f_x*m*a, pil_c_y+geo_f_y*a/2+pil_f_y*geo_f_y*n*a+1*(r0+r1), tteD1,lcEtchingpillar};

//----------------------------------------------------------
// Middle Top

// pttue1_1f~{m}~{n} = newp; Point(pttue1_1f~{m}~{n}) = {pil_c_x+geo_f_x*a/2+pil_f_x*geo_f_x*m*a-1*(r0+0.5*r1), pil_c_y+geo_f_y*a/2+pil_f_y*geo_f_y*n*a, tteD2,lcEtchingpillar};
// pttue2_1f~{m}~{n} = newp; Point(pttue2_1f~{m}~{n}) = {pil_c_x+geo_f_x*a/2+pil_f_x*geo_f_x*m*a,pil_c_y+geo_f_y*a/2+pil_f_y*geo_f_y*n*a-1*(r0+0.5*r1), tteD2,lcEtchingpillar};
// pttue3_1f~{m}~{n} = newp; Point(pttue3_1f~{m}~{n}) = {pil_c_x+geo_f_x*a/2+pil_f_x*geo_f_x*m*a+1*(r0+0.5*r1), pil_c_y+geo_f_y*a/2+pil_f_y*geo_f_y*n*a, tteD2,lcEtchingpillar};
// pttue4_1f~{m}~{n} = newp; Point(pttue4_1f~{m}~{n}) = {pil_c_x+geo_f_x*a/2+pil_f_x*geo_f_x*m*a,pil_c_y+geo_f_y*a/2+pil_f_y*geo_f_y*n*a+1*(r0+0.5*r1), tteD2,lcEtchingpillar};

//----------------------------------------------------------
// Circular boundary

// ctue1_1f~{m}~{n} = newc; Circle(ctue1_1f~{m}~{n}) = {ptue1_1f~{m}~{n}, pc3_1f~{m}~{n}, ptue2_1f~{m}~{n}};
// ctue2_1f~{m}~{n} = newc; Circle(ctue2_1f~{m}~{n}) = {ptue2_1f~{m}~{n}, pc3_1f~{m}~{n}, ptue3_1f~{m}~{n}};
// ctue3_1f~{m}~{n} = newc; Circle(ctue3_1f~{m}~{n}) = {ptue3_1f~{m}~{n}, pc3_1f~{m}~{n}, ptue4_1f~{m}~{n}};
// ctue4_1f~{m}~{n} = newc; Circle(ctue4_1f~{m}~{n}) = {ptue4_1f~{m}~{n}, pc3_1f~{m}~{n}, ptue1_1f~{m}~{n}};

// lue1_1f~{m}~{n} = newl; Line(lue1_1f~{m}~{n}) = {ptue1_1f~{m}~{n}, pth1_1f~{m}~{n}};
// Transfinite Line { lue1_1f~{m}~{n} } = lue1_1f~{m}~{n};
// lue2_1f~{m}~{n} = newl; Line(lue2_1f~{m}~{n}) = {ptue2_1f~{m}~{n}, pth2_1f~{m}~{n}};
// Transfinite Line { lue2_1f~{m}~{n} } = lue2_1f~{m}~{n};
// lue3_1f~{m}~{n} = newl; Line(lue3_1f~{m}~{n}) = {ptue3_1f~{m}~{n}, pth3_1f~{m}~{n}};
// Transfinite Line { lue3_1f~{m}~{n} } = lue3_1f~{m}~{n};
// lue4_1f~{m}~{n} = newl; Line(lue4_1f~{m}~{n}) = {ptue4_1f~{m}~{n}, pth4_1f~{m}~{n}};
// Transfinite Line { lue4_1f~{m}~{n} } = lue4_1f~{m}~{n};

//----------------------------------------------------------
// Lower Etching

//----------------------------------------------------------
// Top

// ptle1_1f~{m}~{n} = newp; Point(ptle1_1f~{m}~{n}) = {pil_c_x+geo_f_x*a/2+pil_f_x*geo_f_x*m*a-1*(r0+r1), pil_c_y+geo_f_y*a/2+pil_f_y*geo_f_y*n*a, tltC,lcEtchingpillar}; // +r1
// ptle2_1f~{m}~{n} = newp; Point(ptle2_1f~{m}~{n}) = {pil_c_x+geo_f_x*a/2+pil_f_x*geo_f_x*m*a, pil_c_y+geo_f_y*a/2+pil_f_y*geo_f_y*n*a-1*(r0+r1), tltC,lcEtchingpillar}; // +r1
// ptle3_1f~{m}~{n} = newp; Point(ptle3_1f~{m}~{n}) = {pil_c_x+geo_f_x*a/2+pil_f_x*geo_f_x*m*a+1*(r0+r1), pil_c_y+geo_f_y*a/2+pil_f_y*geo_f_y*n*a, tltC,lcEtchingpillar}; // +r1
// ptle4_1f~{m}~{n} = newp; Point(ptle4_1f~{m}~{n}) = {pil_c_x+geo_f_x*a/2+pil_f_x*geo_f_x*m*a, pil_c_y+geo_f_y*a/2+pil_f_y*geo_f_y*n*a+1*(r0+r1), tltC,lcEtchingpillar}; // +r1

//----------------------------------------------------------
// Bottom

// pble1_1f~{m}~{n} = newp; Point(pble1_1f~{m}~{n}) = {pil_c_x+geo_f_x*a/2+pil_f_x*geo_f_x*m*a-1*(r0+r1), pil_c_y+geo_f_y*a/2+pil_f_y*geo_f_y*n*a, tlgC,lcEtchingpillar}; // +r1
// pble2_1f~{m}~{n} = newp; Point(pble2_1f~{m}~{n}) = {pil_c_x+geo_f_x*a/2+pil_f_x*geo_f_x*m*a, pil_c_y+geo_f_y*a/2+pil_f_y*geo_f_y*n*a-1*(r0+r1), tlgC,lcEtchingpillar}; // +r1
// pble3_1f~{m}~{n} = newp; Point(pble3_1f~{m}~{n}) = {pil_c_x+geo_f_x*a/2+pil_f_x*geo_f_x*m*a+1*(r0+r1), pil_c_y+geo_f_y*a/2+pil_f_y*geo_f_y*n*a, tlgC,lcEtchingpillar}; // +r1
// pble4_1f~{m}~{n} = newp; Point(pble4_1f~{m}~{n}) = {pil_c_x+geo_f_x*a/2+pil_f_x*geo_f_x*m*a, pil_c_y+geo_f_y*a/2+pil_f_y*geo_f_y*n*a+1*(r0+r1), tlgC,lcEtchingpillar}; // +r1

//----------------------------------------------------------
// Circular boundaries

// ctle1_1f~{m}~{n} = newc; Circle(ctle1_1f~{m}~{n}) = {ptle1_1f~{m}~{n}, pc2_1f~{m}~{n}, ptle2_1f~{m}~{n}};
// ctle2_1f~{m}~{n} = newc; Circle(ctle2_1f~{m}~{n}) = {ptle2_1f~{m}~{n}, pc2_1f~{m}~{n}, ptle3_1f~{m}~{n}};
// ctle3_1f~{m}~{n} = newc; Circle(ctle3_1f~{m}~{n}) = {ptle3_1f~{m}~{n}, pc2_1f~{m}~{n}, ptle4_1f~{m}~{n}};
// ctle4_1f~{m}~{n} = newc; Circle(ctle4_1f~{m}~{n}) = {ptle4_1f~{m}~{n}, pc2_1f~{m}~{n}, ptle1_1f~{m}~{n}};

// cble1_1f~{m}~{n} = newc; Circle(cble1_1f~{m}~{n}) = {pble1_1f~{m}~{n}, pc1_1f~{m}~{n}, pble2_1f~{m}~{n}};
// cble2_1f~{m}~{n} = newc; Circle(cble2_1f~{m}~{n}) = {pble2_1f~{m}~{n}, pc1_1f~{m}~{n}, pble3_1f~{m}~{n}};
// cble3_1f~{m}~{n} = newc; Circle(cble3_1f~{m}~{n}) = {pble3_1f~{m}~{n}, pc1_1f~{m}~{n}, pble4_1f~{m}~{n}};
// cble4_1f~{m}~{n} = newc; Circle(cble4_1f~{m}~{n}) = {pble4_1f~{m}~{n}, pc1_1f~{m}~{n}, pble1_1f~{m}~{n}};

// lle1_1f~{m}~{n} = newl; Line(lle1_1f~{m}~{n}) = {ptle1_1f~{m}~{n}, pble1_1f~{m}~{n}};
// Transfinite Line { lle1_1f~{m}~{n} } = lle1_1f~{m}~{n};
// lle2_1f~{m}~{n} = newl; Line(lle2_1f~{m}~{n}) = {ptle2_1f~{m}~{n}, pble2_1f~{m}~{n}};
// Transfinite Line { lle2_1f~{m}~{n} } = lle2_1f~{m}~{n};
// lle3_1f~{m}~{n} = newl; Line(lle3_1f~{m}~{n}) = {ptle3_1f~{m}~{n}, pble3_1f~{m}~{n}};
// Transfinite Line { lle3_1f~{m}~{n} } = lle3_1f~{m}~{n};
// lle4_1f~{m}~{n} = newl; Line(lle4_1f~{m}~{n}) = {ptle4_1f~{m}~{n}, pble4_1f~{m}~{n}};
// Transfinite Line { lle4_1f~{m}~{n} } = lle4_1f~{m}~{n};

//----------------------------------------------------------
// Lines connecting top and bottom

// lconn5_1f~{m}~{n} = newl; Line(lconn5_1f~{m}~{n}) = {ptle1_1f~{m}~{n}, ptue1_1f~{m}~{n}}; // pth1_1f~{m}~{n}
// Transfinite Line { lconn5_1f~{m}~{n} } = lconn5_1f~{m}~{n};
// lconn6_1f~{m}~{n} = newl; Line(lconn6_1f~{m}~{n}) = {ptle2_1f~{m}~{n}, ptue2_1f~{m}~{n}}; // pth2_1f~{m}~{n} 
// Transfinite Line { lconn6_1f~{m}~{n} } = lconn6_1f~{m}~{n};
// lconn7_1f~{m}~{n} = newl; Line(lconn7_1f~{m}~{n}) = {ptle3_1f~{m}~{n}, ptue3_1f~{m}~{n}}; // pth3_1f~{m}~{n}
// Transfinite Line { lconn7_1f~{m}~{n} } = lconn7_1f~{m}~{n};
// lconn8_1f~{m}~{n} = newl; Line(lconn8_1f~{m}~{n}) = {ptle4_1f~{m}~{n}, ptue4_1f~{m}~{n}}; // pth4_1f~{m}~{n}
// Transfinite Line { lconn8_1f~{m}~{n} } = lconn8_1f~{m}~{n};

//----------------------------------------------------------
// Copper planes

//----------------------------------------------------------
// Connect the upper and lower points with lines to form the plate

// lcpbib17~{m}~{n} = newl; Line(lcpbib17~{m}~{n}) = {ptle1_1f~{m}~{n}, pble1_1f~{m}~{n}};
// Transfinite Line { lcpbib17~{m}~{n} } = lcpbib17~{m}~{n};
// lcpbib18~{m}~{n} = newl; Line(lcpbib18~{m}~{n}) = {ptle2_1f~{m}~{n}, pble2_1f~{m}~{n}};
// Transfinite Line { lcpbib18~{m}~{n} } = lcpbib18~{m}~{n};
// lcpbib19~{m}~{n} = newl; Line(lcpbib19~{m}~{n}) = {ptle3_1f~{m}~{n}, pble3_1f~{m}~{n}};
// Transfinite Line { lcpbib19~{m}~{n} } = lcpbib19~{m}~{n};
// lcpbib20~{m}~{n} = newl; Line(lcpbib20~{m}~{n}) = {ptle4_1f~{m}~{n}, pble4_1f~{m}~{n}};
// Transfinite Line { lcpbib20~{m}~{n} } = lcpbib20~{m}~{n};

//----------------------------------------------------------
// Define surfaces

//----------------------------------------------------------
// Surfaces to which voltages will be applied

// llcp_low_rim_1a~{m}~{n} = newreg; Line Loop(llcp_low_rim_1a~{m}~{n}) = { lle1_1f~{m}~{n}, cble1_1f~{m}~{n}, -lle2_1f~{m}~{n}, -ctle1_1f~{m}~{n} };
// llcp_low_rim_1b~{m}~{n} = newreg; Line Loop(llcp_low_rim_1b~{m}~{n}) = { lle2_1f~{m}~{n}, cble2_1f~{m}~{n}, -lle3_1f~{m}~{n}, -ctle2_1f~{m}~{n} };
// llcp_low_rim_1c~{m}~{n} = newreg; Line Loop(llcp_low_rim_1c~{m}~{n}) = { lle3_1f~{m}~{n}, cble3_1f~{m}~{n}, -lle4_1f~{m}~{n}, -ctle3_1f~{m}~{n} };
// llcp_low_rim_1d~{m}~{n} = newreg; Line Loop(llcp_low_rim_1d~{m}~{n}) = { lle4_1f~{m}~{n}, cble4_1f~{m}~{n}, -lle1_1f~{m}~{n}, -ctle4_1f~{m}~{n} };

//----------------------------------------------------------
// Surfaces to which voltages will be applied

// ps_lower_cp1 = newreg; Surface(ps_lower_cp1) = { llcp_low_rim_1a~{m}~{n} };
// surf_lower_cp1[] += ps_lower_cp1;
// ps_lower_cp2 = newreg; Surface(ps_lower_cp2) = { llcp_low_rim_1b~{m}~{n} };
// surf_lower_cp2[] += ps_lower_cp2;
// ps_lower_cp3 = newreg; Surface(ps_lower_cp3) = { llcp_low_rim_1c~{m}~{n} };
// surf_lower_cp3[] += ps_lower_cp3;
// ps_lower_cp4 = newreg; Surface(ps_lower_cp4) = { llcp_low_rim_1d~{m}~{n} };
// surf_lower_cp4[] += ps_lower_cp4;

//----------------------------------------------------------
// Gas & dielectric surfaces

// ll_cyl_dielectric1b~{m}~{n} = newreg; Line Loop(ll_cyl_dielectric1b~{m}~{n}) = { lconn5_1f~{m}~{n}, ctue1_1f~{m}~{n}, -lconn6_1f~{m}~{n}, -ctle1_1f~{m}~{n} }; // cth1_1f~{m}~{n},
// ps_cyl_dielectric1 = newreg; Surface(ps_cyl_dielectric1) = { ll_cyl_dielectric1b~{m}~{n} };
// surf_cyl_dielectric1[] += ps_cyl_dielectric1;
// ll_cyl_dielectric2b~{m}~{n} = newreg; Line Loop(ll_cyl_dielectric2b~{m}~{n}) = { lconn6_1f~{m}~{n}, ctue2_1f~{m}~{n}, -lconn7_1f~{m}~{n}, -ctle2_1f~{m}~{n} }; // cth2_1f~{m}~{n},
// ps_cyl_dielectric2 = newreg; Surface(ps_cyl_dielectric2) = { ll_cyl_dielectric2b~{m}~{n} };
// surf_cyl_dielectric2[] += ps_cyl_dielectric2;
// ll_cyl_dielectric3b~{m}~{n} = newreg; Line Loop(ll_cyl_dielectric3b~{m}~{n}) = { lconn7_1f~{m}~{n}, ctue3_1f~{m}~{n}, -lconn8_1f~{m}~{n}, -ctle3_1f~{m}~{n} }; // cth3_1f~{m}~{n},
// ps_cyl_dielectric3 = newreg; Surface(ps_cyl_dielectric3) = { ll_cyl_dielectric3b~{m}~{n} };
// surf_cyl_dielectric3[] += ps_cyl_dielectric3;
// ll_cyl_dielectric4b~{m}~{n} = newreg; Line Loop(ll_cyl_dielectric4b~{m}~{n}) = { lconn8_1f~{m}~{n}, ctue4_1f~{m}~{n}, -lconn5_1f~{m}~{n}, -ctle4_1f~{m}~{n} }; // cth4_1f~{m}~{n},
// ps_cyl_dielectric4 = newreg; Surface(ps_cyl_dielectric4) = { ll_cyl_dielectric4b~{m}~{n} };
// surf_cyl_dielectric4[] += ps_cyl_dielectric4;

// ll_top_gas2~{m}~{n} = newreg; Line Loop( ll_top_gas2~{m}~{n} ) = { cth1_1f~{m}~{n}, cth2_1f~{m}~{n}, cth3_1f~{m}~{n}, cth4_1f~{m}~{n} };
// ps_top_gas1 = news; Plane Surface(news) = { ll_top_gas2~{m}~{n} };
// surf_top_gas1[] += ps_top_gas1;
// Transfinite Surface { surf_top_gas1[] };
// Recombine Surface { surf_top_gas1[] };

// ll_top_gas3~{m}~{n} = newreg; Line Loop( ll_top_gas3~{m}~{n} ) = { lue1_1f~{m}~{n}, cth1_1f~{m}~{n}, -lue2_1f~{m}~{n}, -ctue1_1f~{m}~{n} };
// ps_top_gas2 = news; Surface(news) = { ll_top_gas3~{m}~{n} };
// surf_top_gas2[] += ps_top_gas2;
// Transfinite Surface { surf_top_gas2[] };
// Recombine Surface { surf_top_gas2[] };

// ll_top_gas4~{m}~{n} = newreg; Line Loop( ll_top_gas4~{m}~{n} ) = { lue2_1f~{m}~{n}, cth2_1f~{m}~{n}, -lue3_1f~{m}~{n}, -ctue2_1f~{m}~{n}  };
// ps_top_gas3 = news; Surface(news) = { ll_top_gas4~{m}~{n} };
// surf_top_gas3[] += ps_top_gas3;
// Transfinite Surface { surf_top_gas3[] };
// Recombine Surface { surf_top_gas3[] };

// ll_top_gas5~{m}~{n} = newreg; Line Loop( ll_top_gas5~{m}~{n} ) = { lue3_1f~{m}~{n}, cth3_1f~{m}~{n}, -lue4_1f~{m}~{n}, -ctue3_1f~{m}~{n} };
// ps_top_gas4 = news; Surface(news) = { ll_top_gas5~{m}~{n} };
// surf_top_gas4[] += ps_top_gas4;
// Transfinite Surface { surf_top_gas4[] };
// Recombine Surface { surf_top_gas4[] };

// ll_top_gas6~{m}~{n} = newreg; Line Loop( ll_top_gas6~{m}~{n} ) = { lue4_1f~{m}~{n}, cth4_1f~{m}~{n}, -lue1_1f~{m}~{n}, -ctue4_1f~{m}~{n} };
// ps_top_gas5 = news; Surface(news) = { ll_top_gas6~{m}~{n} };
// surf_top_gas5[] += ps_top_gas5;
// Transfinite Surface { surf_top_gas5[] };
// Recombine Surface { surf_top_gas5[] };

// ll_top_cp1a1 = newreg; Line Loop( ll_top_cp1a1 ) = { ctle1_1f~{m}~{n}, ctle2_1f~{m}~{n}, ctle3_1f~{m}~{n}, ctle4_1f~{m}~{n} };
// ll_top_cp2a[] += { ll_top_cp1a1 };

// ll_bottom_cp1a1 = newreg; Line Loop( ll_bottom_cp1a1 ) = { cble1_1f~{m}~{n}, cble2_1f~{m}~{n}, cble3_1f~{m}~{n}, cble4_1f~{m}~{n} };
// ll_bottom_cp2a[] += { ll_bottom_cp1a1 };

//   Delete { Point { pc1_1f~{m}~{n} }; }
//   Delete { Point { pc2_1f~{m}~{n} }; }
//   Delete { Point { pc3_1f~{m}~{n} }; }
//   Delete { Point { pc4_1f~{m}~{n} }; }

// EndFor
// EndFor

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
// EXTERNAL ELECTRODES

//----------------------------------------------------------
// Top electrode

pexet1 = newp; Point(pexet1) = {geo_f_x*0+geo_f_x*m_1*a, geo_f_y*0+geo_f_y*n_1*a, tutC+lE,lcExtElectrodeBdry};
// pexet2 = newp; Point(pexet2) = {geo_f_x*a/2+geo_f_x*m_1*a, geo_f_y*0+geo_f_y*n_1*a, tutC+lE,lcExtElectrodeBdry};
pexet3 = newp; Point(pexet3) = {geo_f_x*a+geo_f_x*m_1*a, geo_f_y*0+geo_f_y*n_1*a, tutC+lE,lcExtElectrodeBdry};
pexet4 = newp; Point(pexet4) = {geo_f_x*a+geo_f_x*m_1*a, geo_f_y*a+geo_f_y*n_1*a, tutC+lE,lcExtElectrodeBdry};
// pexet5 = newp; Point(pexet5) = {geo_f_x*a/2+geo_f_x*m_1*a, geo_f_y*a+geo_f_y*n_1*a, tutC+lE,lcExtElectrodeBdry};
pexet6 = newp; Point(pexet6) = {geo_f_x*0+geo_f_x*m_1*a, geo_f_y*a+geo_f_y*n_1*a, tutC+lE,lcExtElectrodeBdry};

//----------------------------------------------------------
// Top electrode lines

lexet1 = newl; Line(lexet1) = {pexet1, pexet3};
// Transfinite Line { lexet1 } = lexet1;
// lexet2 = newl; Line(lexet2) = {pexet2, pexet3};
// Transfinite Line { lexet2 } = lexet2;
lexet3 = newl; Line(lexet3) = {pexet3, pexet4};
// Transfinite Line { lexet3 } = lexet3;
lexet4 = newl; Line(lexet4) = {pexet4, pexet6};
// Transfinite Line { lexet4 } = lexet4;
// lexet5 = newl; Line(lexet5) = {pexet5, pexet6};
// Transfinite Line { lexet5 } = lexet5;
lexet6 = newl; Line(lexet6) = {pexet6, pexet1};
// Transfinite Line { lexet6 } = lexet6;

//----------------------------------------------------------
// Upper electrode boundary - transfinite

// lexett1 = newl; Line(lexett1) = {pexet1, pexet3};
// Transfinite Line { lexett1 } = lexett1; 
// lexett2 = newl; Line(lexett2) = {pexet3, pexet4};
// Transfinite Line { lexett2 } = lexett2;
// lexett3 = newl; Line(lexett3) = {pexet4, pexet6};
// Transfinite Line { lexett3 } = lexett3;
// lexett4 = newl; Line(lexett4) = {pexet6, pexet1};
// Transfinite Line { lexett4 } = lexett4;

//----------------------------------------------------------
// Connect the top electrode to the LEM.

lexetc1 = newl; Line(lexetc1) = {pexet1, pc2_1}; // pc4_1
// Transfinite Line { lexetc1 } = lexetc1;
// lexetc2 = newl; Line(lexetc2) = {pexet2, ptmc_1};
// Transfinite Line { lexetc2 } = lexetc2;
lexetc3 = newl; Line(lexetc3) = {pexet3, pc2_2}; // pc4_2
// Transfinite Line { lexetc3 } = lexetc3;
lexetc4 = newl; Line(lexetc4) = {pexet4, pc2_3}; // pc4_3
// Transfinite Line { lexetc4 } = lexetc4;
// lexetc5 = newl; Line(lexetc5) = {pexet5, ptmc_3};
// Transfinite Line { lexetc5 } = lexetc5;
lexetc6 = newl; Line(lexetc6) = {pexet6, pc2_4}; // pc4_4
// Transfinite Line { lexetc6 } = lexetc6;

//----------------------------------------------------------
// Bottom electrode

// pexeb1 = newp; Point(pexeb1) = {geo_f_x*0+geo_f_x*m_1*a, geo_f_y*0+geo_f_y*n_1*a, mesh_level * mm,lcExtElectrodeBdry};
// pexeb2 = newp; Point(pexeb2) = {geo_f_x*a/2+geo_f_x*m_1*a, geo_f_y*0+geo_f_y*n_1*a, mesh_level * mm,lcExtElectrodeBdry};
// pexeb3 = newp; Point(pexeb3) = {geo_f_x*a+geo_f_x*m_1*a, geo_f_y*0+geo_f_y*n_1*a, mesh_level * mm,lcExtElectrodeBdry};
// pexeb4 = newp; Point(pexeb4) = {geo_f_x*a+geo_f_x*m_1*a, geo_f_y*a+geo_f_y*n_1*a, mesh_level * mm,lcExtElectrodeBdry};
// pexeb5 = newp; Point(pexeb5) = {geo_f_x*a/2+geo_f_x*m_1*a, geo_f_y*a+geo_f_y*n_1*a, mesh_level * mm,lcExtElectrodeBdry};
// pexeb6 = newp; Point(pexeb6) = {geo_f_x*0+geo_f_x*m_1*a, geo_f_y*a+geo_f_y*n_1*a, mesh_level * mm,lcExtElectrodeBdry};

//----------------------------------------------------------
// Bottom wire level boundary - transfinite

// lexebt1 = newl; Line(lexebt1) = {pexeb1, pexeb3};
// Transfinite Line { lexebt1 } = lexebt1;
// lexebt2 = newl; Line(lexebt2) = {pexeb3, pexeb4};
// Transfinite Line { lexebt2 } = lexebt2;
// lexebt3 = newl; Line(lexebt3) = {pexeb4, pexeb6};
// Transfinite Line { lexebt3 } = lexebt3;
// lexebt4 = newl; Line(lexebt4) = {pexeb6, pexeb1};
// Transfinite Line { lexebt4 } = lexebt4;

//----------------------------------------------------------
// Copper plate surfaces

// llcp_up_border1 = newreg; Line Loop(llcp_up_border1) = {lcptlb5a, -lcptib10, -lcptub1a, lcptib9}; // lcptlb5b, -lcptub1b
// pscp_up_border1 = newreg; Plane Surface(pscp_up_border1) = {llcp_up_border1};
// Transfinite Surface { pscp_up_border1 };
// Recombine Surface { pscp_up_border1 };

// llcp_up_border2 = newreg; Line Loop(llcp_up_border2) = {lcptlb6a, -lcptib11, -lcptub2a, lcptib10}; // lcptlb6b, -lcptub2b
// pscp_up_border2 = newreg; Plane Surface(pscp_up_border2) = {llcp_up_border2};
// Transfinite Surface { pscp_up_border2 };
// Recombine Surface { pscp_up_border2 };

// llcp_up_border3 = newreg; Line Loop(llcp_up_border3) = {lcptlb7a, -lcptib12, -lcptub3a, lcptib11}; // lcptlb7b, -lcptub3b
// pscp_up_border3 = newreg; Plane Surface(pscp_up_border3) = {llcp_up_border3};
// Transfinite Surface { pscp_up_border3 };
// Recombine Surface { pscp_up_border3 };

// llcp_up_border4 = newreg; Line Loop(llcp_up_border4) = {lcptlb8a, -lcptib9, -lcptub4a, lcptib12}; // lcptlb8b, -lcptub4b
// pscp_up_border4 = newreg; Plane Surface(pscp_up_border4) = {llcp_up_border4};
// Transfinite Surface { pscp_up_border4 };
// Recombine Surface { pscp_up_border4 };

llcp_low_border1 = newreg; Line Loop(llcp_low_border1) = {lcpblb5a, -lcpbib10, -lcpbub1a, lcpbib9}; // lcpblb5b,  -lcpbub1b
pscp_low_border1 = newreg; Plane Surface(pscp_low_border1) = {llcp_low_border1};
// Transfinite Surface { pscp_low_border1 };
// Recombine Surface { pscp_low_border1 };

llcp_low_border2 = newreg; Line Loop(llcp_low_border2) = {lcpblb6a, -lcpbib11, -lcpbub2a, lcpbib10}; // lcpblb6b, -lcpbub2b  
pscp_low_border2 = newreg; Plane Surface(pscp_low_border2) = {llcp_low_border2};
// Transfinite Surface { pscp_low_border2 };
// Recombine Surface { pscp_low_border2 };

llcp_low_border3 = newreg; Line Loop(llcp_low_border3) = {lcpblb7a, -lcpbib12, -lcpbub3a, lcpbib11}; // lcpblb7b, -lcpbub3b
pscp_low_border3 = newreg; Plane Surface(pscp_low_border3) = {llcp_low_border3};
// Transfinite Surface { pscp_low_border3 };
// Recombine Surface { pscp_low_border3 };

llcp_low_border4 = newreg; Line Loop(llcp_low_border4) = {lcpblb8a, -lcpbib9, -lcpbub4a, lcpbib12}; // lcpblb8b,  -lcpbub4b
pscp_low_border4 = newreg; Plane Surface(pscp_low_border4) = {llcp_low_border4};
// Transfinite Surface { pscp_low_border4 };
// Recombine Surface { pscp_low_border4 };

// llcp_face1 = newreg; Line Loop(llcp_face1) = {lcptub1a, lcptub2a, lcptub3a, lcptub4a}; // lcptub1b, lcptub2b, lcptub3b, lcptub4b 
// llcp_face3 = newreg; Line Loop(llcp_face3) = {lcpbub1a, lcpbub2a, lcpbub3a, lcpbub4a}; // lcpbub1b, lcpbub2b, lcpbub3b, lcpbub4b

//----------------------------------------------------------
// Copper plate surfaces

// ll_side_gas1a = newreg; Line Loop(ll_side_gas1a) = {lcptlb5a, -lcorner2, -lcpbub1a, lcorner1}; // lmid1_2,
// ps_side_gas1a = newreg; Plane Surface(ps_side_gas1a) = {ll_side_gas1a};
// Transfinite Surface { ps_side_gas1a };
// Recombine Surface { ps_side_gas1a };

// ll_side_gas2a = newreg; Line Loop(ll_side_gas2a) = {lcptlb6a, -lcorner3, -lcpbub2a, lcorner2}; // lmid2_2,
// ps_side_gas2a = newreg; Plane Surface(ps_side_gas2a) = {ll_side_gas2a};
// Transfinite Surface { ps_side_gas2a };
// Recombine Surface { ps_side_gas2a };

// ll_side_gas3a = newreg; Line Loop(ll_side_gas3a) = {lcptlb7a, -lcorner4, -lcpbub3a, lcorner3}; // lmid3_2,
// ps_side_gas3a = newreg; Plane Surface(ps_side_gas3a) = {ll_side_gas3a};
// Transfinite Surface { ps_side_gas3a };
// Recombine Surface { ps_side_gas3a };

// ll_side_gas4a = newreg; Line Loop(ll_side_gas4a) = {lcptlb8a, -lcorner1, -lcpbub4a, lcorner4}; // lmid4_2,
// ps_side_gas4a = newreg; Plane Surface(ps_side_gas4a) = {ll_side_gas4a};
// Transfinite Surface { ps_side_gas4a };
// Recombine Surface { ps_side_gas4a };

// ll_side_gas1b = newreg; Line Loop(ll_side_gas1b) = {lcptlb5b, lcorner2, -lcpblb5b, lcorner1}; // -lmid1_2
// ps_side_gas1b = newreg; Plane Surface(ps_side_gas1b) = {ll_side_gas1b};
// ll_side_gas2b = newreg; Line Loop(ll_side_gas2b) = {lcptlb6b, lcorner3, -lcpblb6b, lcorner2}; // -lmid2_2
// ps_side_gas2b = newreg; Plane Surface(ps_side_gas2b) = {ll_side_gas2b};
// ll_side_gas3b = newreg; Line Loop(ll_side_gas3b) = {lcptlb7b, lcorner4, -lcpblb7b, lcorner3}; // -lmid3_2
// ps_side_gas3b = newreg; Plane Surface(ps_side_gas3b) = {ll_side_gas3b};
// ll_side_gas4b = newreg; Line Loop(ll_side_gas4b) = {lcptlb8b, lcorner1, -lcpblb8b, lcorner4}; // -lmid4_2
// ps_side_gas4b = newreg; Plane Surface(ps_side_gas4b) = {ll_side_gas4b};

//----------------------------------------------------------
// Bounding and intersecting surfaces

//----------------------------------------------------------
// Face physsurf_bdh_1_1 (Corner 1 - Corner 2)

l1bdh_1_1_bsurft1 = newl; Line(l1bdh_1_1_bsurft1) = {pexet3, pbdhat_1_2_2[1]};
l2bdh_1_1_bsurft1 = newl; Line(l2bdh_1_1_bsurft1) = {pexet1, p2b_1_1};
l1bdh_1_1_bsurfb1 = newl; Line(l1bdh_1_1_bsurfb1) = {pc2_2, pbdhbb_1_1_2[1]};
l2bdh_1_1_bsurfb1 = newl; Line(l2bdh_1_1_bsurfb1) = {pc2_1, p1b_1_1};

llbdh_1_1_bsurf1t = newreg; Line Loop(llbdh_1_1_bsurf1t) = {lexet1, l1bdh_1_1_bsurft1, -ll_bdhac1_1_2_2[0], -ll_bdhac0_1_2_2[0], -ll_bdhbt_1_1_2[0], -ll_bdhbt_1_1_1[0], -l2_2b, -l1_2b, -l2bdh_1_1_bsurft1};
llbdh_1_1_bsurf1b = newreg; Line Loop(llbdh_1_1_bsurf1b) = {lcpbub1a, l1bdh_1_1_bsurfb1, -ll_bdhbb_1_1_2[0], -ll_bdhbb_1_1_1[0], -l2bdh_1_1_bsurfb1};

psbdh_1_1_bsurf1t = newreg; Plane Surface(psbdh_1_1_bsurf1t) = {llbdh_1_1_bsurf1t};
psbdh_1_1_bsurf1b = newreg; Plane Surface(psbdh_1_1_bsurf1b) = {llbdh_1_1_bsurf1b};

//----------------------------------------------------------
// Face physsurf_bdh_1_2 (Corner 2 - Corner 3)

l1bdh_1_2_bsurft1 = newl; Line(l1bdh_1_2_bsurft1) = {pexet4, p2a_1_1};
l1bdh_1_2_bsurfb1 = newl; Line(l1bdh_1_2_bsurfb1) = {pc2_3, p1a_1_1};

llbdh_1_2_bsurf3t = newreg; Line Loop(llbdh_1_2_bsurf3t) = {lexet3, l1bdh_1_2_bsurft1, ll_bdhat_1_2_1[0], ll_bdhat_1_2_2[0], -l1bdh_1_1_bsurft1};
llbdh_1_2_bsurf3b = newreg; Line Loop(llbdh_1_2_bsurf3b) = {-lcpbub2a, l1bdh_1_1_bsurfb1, ll_bdhbc0_1_1_2[0], ll_bdhbc1_1_1_2[0], -ll_bdhab_1_2_2[0], -ll_bdhab_1_2_1[0], l3_1a, l4_1a, -l1bdh_1_2_bsurfb1};

psbdh_1_2_bsurf3t = newreg; Plane Surface(psbdh_1_2_bsurf3t) = {llbdh_1_2_bsurf3t};
psbdh_1_2_bsurf3b = newreg; Plane Surface(psbdh_1_2_bsurf3b) = {llbdh_1_2_bsurf3b};

//----------------------------------------------------------
// Face physsurf_bdh_2_1 (Corner 3 - Corner 4)

l1bdh_2_1_bsurft4 = newl; Line(l1bdh_2_1_bsurft4) = {pexet6, pbdhbt_2_2_2[1]};
l1bdh_2_1_bsurfb4 = newl; Line(l1bdh_2_1_bsurfb4) = {pc2_4, pbdha1_2_1_2[0]};

llbdh_2_1_bsurf4t = newreg; Line Loop(llbdh_2_1_bsurf4t) = {lexet4, l1bdh_2_1_bsurft4, ll_bdhbc0_2_2_2[0], ll_bdhbc1_2_2_2[0], -ll_bdhat_2_1_2[0], -ll_bdhat_2_1_1[0], l3_2a, l4_2a, -l1bdh_1_2_bsurft1};
llbdh_2_1_bsurf4b = newreg; Line Loop(llbdh_2_1_bsurf4b) = {lcpbub3a, l1bdh_2_1_bsurfb4, -ll_bdhab_2_1_2[0], -ll_bdhab_2_1_1[0], -l1bdh_1_2_bsurfb1};

psbdh_2_1_bsurf4t = newreg; Plane Surface(psbdh_2_1_bsurf4t) = {llbdh_2_1_bsurf4t};
psbdh_2_1_bsurf4b = newreg; Plane Surface(psbdh_2_1_bsurf4b) = {llbdh_2_1_bsurf4b};

//----------------------------------------------------------
// Face physsurf_bdh_2_2 (Corner 4 - Corner 1)

llbdh_2_2_bsurf6t = newreg; Line Loop(llbdh_2_2_bsurf6t) = {lexet6, l2bdh_1_1_bsurft1, ll_bdhbt_2_2_2[0], ll_bdhbt_2_2_1[0], -l1bdh_2_1_bsurft4};
llbdh_2_2_bsurf6b = newreg; Line Loop(llbdh_2_2_bsurf6b) = {lcpbub4a, l2bdh_1_1_bsurfb1, l1_1b, l2_1b, ll_bdhbb_2_2_1[0], ll_bdhbb_2_2_2[0], ll_bdhac0_2_1_2[0], ll_bdhac1_2_1_2[0], -l1bdh_2_1_bsurfb4};

psbdh_2_2_bsurf6t = newreg; Plane Surface(psbdh_2_2_bsurf6t) = {llbdh_2_2_bsurf6t};
psbdh_2_2_bsurf6b = newreg; Plane Surface(psbdh_2_2_bsurf6b) = {llbdh_2_2_bsurf6b};

//----------------------------------------------------------
// Bounding surfaces

ll_bsurf7 = newreg; Line Loop(ll_bsurf7) = {lexet1, lexet3, lexet4, lexet6}; // lexet2, lexet5,
ps_bsurf7 = newreg; Plane Surface( ps_bsurf7 ) = { ll_bsurf7 };
// Transfinite Surface { ps_bsurf7 };
// Recombine Surface { ps_bsurf7 };

ll_top_cp1a2 = newreg; Line Loop(ll_top_cp1a2) = {lcpbub1a, lcpbub2a, lcpbub3a, lcpbub4a}; // lcpbub1b, lcpbub2b, lcpbub3b, lcpbub4b

ps_top_cp2a = news; Plane Surface(news) = {ll_top_cp1a2}; // ll_top_cp2a
surf_top_cp[] += {ps_top_cp2a};
// Transfinite Surface { surf_top_cp[] };
// Recombine Surface { surf_top_cp[] };

// ps_bottom_dielectric1a1 = news; Plane Surface(news) = {ll_bottom_cp2a};
// surf_bottom_dielectric[] = {ps_bottom_dielectric1a1};
// Transfinite Surface { surf_bottom_dielectric[] };
// Recombine Surface { surf_bottom_dielectric[] };

ll_bottom_cp1a2 = newreg; Line Loop(ll_bottom_cp1a2) = {lcpblb5a, lcpblb6a, lcpblb7a, lcpblb8a}; // lcpblb5b, lcpblb6b, lcpblb7b, lcpblb8b
ps_bottom_cp2a = news; Plane Surface(news) = {ll_bottom_cp1a2}; // ll_bottom_cp2a
surf_bottom_cp[] += {ps_bottom_cp2a};
// Transfinite Surface { surf_bottom_cp[] };
// Recombine Surface { surf_bottom_cp[] };

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
//// FINAL DEFINITIONS (SURFACES / VOLUMES)


//------------------------------------------------------------------------------------------
/// SURFACE LOOPS

//----------------------------------------------------------
// Wire Gas Interior Surface Loop - interior wire gas surface loop

sl_wire_gas_total_surface = newreg; Surface Loop(sl_wire_gas_total_surface) = { ps_bsurf7, psbdh_1_1_bsurf1t, psbdh_1_1_bsurf1b, psbdh_1_2_bsurf3t, psbdh_1_2_bsurf3b, psbdh_2_1_bsurf4t, psbdh_2_1_bsurf4b, psbdh_2_2_bsurf6t, psbdh_2_2_bsurf6b, -surf_top_cp[], -sa_1_1[0], -sa_1_1[1], -sa_1_2[0], -sa_1_2[1], -sb_1_1[0], -sb_1_1[1], -sb_1_2[0], -sb_1_2[1], -sa_2_1[0], -sa_2_1[1], -sa_2_2[0], -sa_2_2[1], -sb_2_1[0], -sb_2_1[1], -sb_2_2[0], -sb_2_2[1] };
// -surf_top_gas1[], -surf_top_gas2[], -surf_top_gas3[], -surf_top_gas4[], -surf_top_gas5[], -surf_cyl_dielectric1[], -surf_cyl_dielectric2[], -surf_cyl_dielectric3[], -surf_cyl_dielectric4[],
// pscp_up_border1, pscp_up_border2, pscp_up_border3, pscp_up_border4,
// ps_side_gas1a, ps_side_gas2a, ps_side_gas3a, ps_side_gas4a,
total_sl_wire_gas_total_surface[] += sl_wire_gas_total_surface;

//----------------------------------------------------------
// Gas Exterior Surface Loop - exterior gas surface loop

// total_sl_gas_exterior_surface[0] = newreg; Surface Loop(total_sl_gas_exterior_surface[0]) = { ps_side_gas1b, ps_side_gas2b, ps_side_gas3b, ps_side_gas4b, ps_bsurf2, ps_bsurf5 };

//----------------------------------------------------------
// Dielectric Surface Loop - dielectric surface loop

// sl_dielectric = newreg; Surface Loop(sl_dielectric) = { surf_top_gas1[], surf_top_gas2[], surf_top_gas3[], surf_top_gas4[], surf_top_gas5[], surf_cyl_dielectric1[], surf_cyl_dielectric2[], surf_cyl_dielectric3[], surf_cyl_dielectric4[], -surf_lower_cp1[], 
// -surf_lower_cp2[], -surf_lower_cp3[], -surf_lower_cp4[], surf_bottom_dielectric[] };
// total_sl_dielectric[] += sl_dielectric;

//----------------------------------------------------------
// Wire Volume Surface Loop - wire electrode

// sl_wire = newreg; Surface Loop(sl_wire) = { s_1_2b, sb_1_1[2], sb_1_2[2], tmpa_2_2[0], s_1_1a, sa_2_1[2], sa_2_2[2], tmpb_1_2[0], s_1_2a, sa_1_1[2], sa_1_2[2], tmpb_2_2[0], s_1_1b, sb_2_1[2], sb_2_2[2], tmpa_1_2[0], sa_1_1[0], sa_1_1[1], sa_1_2[0], 
// sa_1_2[1], sb_1_1[0], sb_1_1[1], sb_1_2[0], sb_1_2[1], sa_2_1[0], sa_2_1[1], sa_2_2[0], sa_2_2[1], sb_2_1[0], sb_2_1[1], sb_2_2[0], sb_2_2[1] };
// vol_wire = newreg; Volume(vol_wire) = {sl_wire};

//----------------------------------------------------------
// Lower Electrode Surface Loop - lower electrode surface loop

sl_lower_cp = newreg; Surface Loop(sl_lower_cp) = { surf_top_cp[], pscp_low_border1, pscp_low_border2, pscp_low_border3, pscp_low_border4, surf_bottom_cp[] }; // surf_lower_cp1[], surf_lower_cp2[], surf_lower_cp3[], surf_lower_cp4[],
total_sl_lower_cp[] += sl_lower_cp;

//------------------------------------------------------------------------------------------
/// CONTAINER VOLUME

// vol_container = newreg; Volume(vol_container) = { total_sl_gas_exterior_surface[0] };

//------------------------------------------------------------------------------------------
/// WIRE VOLUME

// total_vol_wire = newreg; Compound Volume(total_vol_wire) = { vol_x1_wire, vol_x2_wire, vol_y1_wire, vol_y2_wire };
// total_vol_wire = newreg; Volume(total_vol_wire) = { total_sl_wire_gas_interior_surface[] };

//------------------------------------------------------------------------------------------
/// GAS VOLUME

vol_gas = newreg; Volume(vol_gas) = { total_sl_wire_gas_total_surface[] }; //  sl_wire_exterior_surface_1a[], sl_wire_exterior_surface_1b[], sl_wire_exterior_surface_2a[], sl_wire_exterior_surface_2b[]
// vol_gas = newreg; Volume(vol_gas) = { total_sl_wire_gas_total_surface[] }; // total_sl_gas_exterior_surface[0]


//------------------------------------------------------------------------------------------
/// COMPONENT VOLUMES

// vol_dielectric = newreg; Volume(vol_dielectric) = total_sl_dielectric[];
vol_lower_cp = newreg; Volume(vol_lower_cp) = total_sl_lower_cp[];


//------------------------------------------------------------------------------------------
/// PHYSICAL SURFACES

//----------------------------------------------------------
// Physical Surfaces - periodic boundary conditions

physsurf_bdh_1_1 = newreg; Physical Surface(physsurf_bd1h1) = { psbdh_1_1_bsurf1t, s_1_2b1[], sb_1_1[2], sb_1_2[2], s_1_2a2[], psbdh_1_1_bsurf1b };             // ps_side_gas1b, ps_bsurf2, pscp_up_border1, ps_side_gas1a,
physsurf_bdh_1_2 = newreg; Physical Surface(physsurf_bd1h2) = { psbdh_1_2_bsurf3t, s_1_1a1[], sa_2_1[2], sa_2_2[2], s_1_1b2[], psbdh_1_2_bsurf3b };             // ps_side_gas2b, pscp_up_border2, ps_side_gas2a,
physsurf_bdh_2_1 = newreg; Physical Surface(physsurf_bd2h1) = { psbdh_2_1_bsurf4t, s_1_2a1[], sa_1_1[2], sa_1_2[2], s_1_2b2[], psbdh_2_1_bsurf4b };             // ps_side_gas3b, ps_bsurf5, pscp_up_border3, ps_side_gas3a,
physsurf_bdh_2_2 = newreg; Physical Surface(physsurf_bd2h2) = { psbdh_2_2_bsurf6t, s_1_1b1[], sb_2_1[2], sb_2_2[2], s_1_1a2[], psbdh_2_2_bsurf6b };             // ps_side_gas4b, pscp_up_border4, ps_side_gas4a,

//----------------------------------------------------------
// Physical surfaces - container surface

// physsurf_container = newreg; Physical Surface(physsurf_container) = { surf_top_cp[], ps_side_gas1a, ps_side_gas2a, ps_side_gas3a, ps_side_gas4a,
// ps_bsurf1, ps_bsurf3, ps_bsurf4, ps_bsurf6, ps_bsurf7, pscp_up_border1, pscp_up_border2, pscp_up_border3, pscp_up_border4, ps_side_gas1b, ps_side_gas2b, ps_side_gas3b, ps_side_gas4b, ps_bsurf2, ps_bsurf5 };

//----------------------------------------------------------
// Physical surfaces - wire surface

// physsurf_wire = newreg; Physical Surface(physsurf_wire) = { sa_1_1[0], sa_1_1[1], sa_1_2[0], sa_1_2[1], sb_1_1[0], sb_1_1[1], sb_1_2[0], sb_1_2[1], sa_2_1[0], sa_2_1[1], sa_2_2[0], sa_2_2[1], sb_2_1[0], sb_2_1[1], sb_2_2[0], sb_2_2[1] };
// physsurf_wire = newreg; Physical Surface(physsurf_wire) = { physsurf_1a_wire, physsurf_1b_wire, physsurf_2a_wire, physsurf_2b_wire };

// physsurf_wire = newreg; Physical Surface(physsurf_wire) = { s_1_2b, sb_1_1[2], sb_1_2[2], tmpa_2_2[0], s_1_1a, sa_2_1[2], sa_2_2[2], tmpb_1_2[0], s_1_2a, sa_1_1[2], sa_1_2[2], tmpb_2_2[0], s_1_1b, sb_2_1[2], sb_2_2[2], tmpa_1_2[0], sa_1_1[0], sa_1_1[1], 
// sa_1_2[0], sa_1_2[1], sb_1_1[0], sb_1_1[1], sb_1_2[0], sb_1_2[1], sa_2_1[0], sa_2_1[1], sa_2_2[0], sa_2_2[1], sb_2_1[0], sb_2_1[1], sb_2_2[0], sb_2_2[1] };

//----------------------------------------------------------
// Physical surfaces - gas exterior / interior surface

physsurf_gas = newreg; Physical Surface(physsurf_gas) = { ps_bsurf7, psbdh_1_1_bsurf1t, psbdh_1_1_bsurf1b, psbdh_1_2_bsurf3t, psbdh_1_2_bsurf3b, psbdh_2_1_bsurf4t, psbdh_2_1_bsurf4b, psbdh_2_2_bsurf6t, psbdh_2_2_bsurf6b, -surf_top_cp[], -sa_1_1[0], -sa_1_1[1], -sa_1_2[0], -sa_1_2[1], -sb_1_1[0], -sb_1_1[1], -sb_1_2[0], -sb_1_2[1], -sa_2_1[0], -sa_2_1[1], -sa_2_2[0], -sa_2_2[1], -sb_2_1[0], -sb_2_1[1], -sb_2_2[0], -sb_2_2[1] };
// Physical Surface(physsurf_gas) = { physsurf_container, -physsurf_x1_wire, -physsurf_x2_wire, -physsurf_y1_wire, -physsurf_y2_wire };

//----------------------------------------------------------
// Physical surfaces - dielectric surface

// physsurf_dielectric = newreg; Physical Surface(physsurf_dielectric) = total_sl_dielectric[];
// physsurf_dielectric = newreg; Physical Surface(physsurf_dielectric) = { surf_top_gas1[], surf_top_gas2[], surf_top_gas3[], surf_top_gas4[], surf_top_gas5[], surf_cyl_dielectric1[], surf_cyl_dielectric2[], surf_cyl_dielectric3[], surf_cyl_dielectric4[], 
// -surf_lower_cp1[], -surf_lower_cp2[], -surf_lower_cp3[], -surf_lower_cp4[], surf_bottom_dielectric[] };

//----------------------------------------------------------
// Physical surfaces - lower electrode surface

physsurf_lower_cp = newreg; Physical Surface(physsurf_lower_cp) = { surf_top_cp[], pscp_low_border1, pscp_low_border2, pscp_low_border3, pscp_low_border4, surf_bottom_cp[] }; // surf_lower_cp1[], surf_lower_cp2[], surf_lower_cp3[], surf_lower_cp4[], 

//----------------------------------------------------------
// Physical surfaces - upper electrode surface

physsurf_upper_el = newreg; Physical Surface(physsurf_upper_el) = { ps_bsurf7 };

//------------------------------------------------------------------------------------------
/// PHYSICAL VOLUMES

//----------------------------------------------------------
// Physical Volumes - Container Volume

// physvol_container = newreg; Physical Volume(physvol_container) = vol_container;

//----------------------------------------------------------
// Physical Volumes - Wire Volume

// physvol_wire = newreg; Physical Volume(physvol_wire) = total_vol_wire;

//----------------------------------------------------------
// Physical Volumes - Gas Volume

physvol_gas = newreg; Physical Volume(physvol_gas) = { vol_gas };
// Physical Volume(physvol_gas) = { vol_gas };

//----------------------------------------------------------
// Physical Volumes - Component Volumes

// physvol_dielectric = newreg; Physical Volume(physvol_dielectric) = { vol_dielectric };
physvol_lower_cp = newreg; Physical Volume(physvol_lower_cp) = { vol_lower_cp };

//----------------------------------------------------------
// FEATURE / ELEMENT REMOVAL

// Delete { Volume { vol_gas }; }

// Coherence;
// Geometry.AutoCoherence = 1;



