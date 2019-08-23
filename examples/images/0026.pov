#include "colors.inc"
#include "finish.inc"

global_settings {assumed_gamma 1 max_trace_level 6}
background {color White}
camera {orthographic
  right -25.25*x up 76.08*y
  direction 1.00*z
  location <0,0,50.00> look_at <0,0,0>}
light_source {<  2.00,   3.00,  40.00> color White
  area_light <0.70, 0, 0>, <0, 0.70, 0>, 3, 3
  adaptive 1 jitter}

#declare simple = finish {phong 0.7}
#declare pale = finish {ambient .5 diffuse .85 roughness .001 specular 0.200 }
#declare intermediate = finish {ambient 0.3 diffuse 0.6 specular 0.10 roughness 0.04 }
#declare vmd = finish {ambient .0 diffuse .65 phong 0.1 phong_size 40. specular 0.500 }
#declare jmol = finish {ambient .2 diffuse .6 specular 1 roughness .001 metallic}
#declare ase2 = finish {ambient 0.05 brilliance 3 diffuse 0.6 metallic specular 0.70 roughness 0.04 reflection 0.15}
#declare ase3 = finish {ambient .15 brilliance 2 diffuse .6 metallic specular 1. roughness .001 reflection .0}
#declare glass = finish {ambient .05 diffuse .3 specular 1. roughness .001}
#declare glass2 = finish {ambient .0 diffuse .3 specular 1. reflection .25 roughness .001}
#declare Rcell = 0.050;
#declare Rbond = 0.100;

#macro atom(LOC, R, COL, TRANS, FIN)
  sphere{LOC, R texture{pigment{color COL transmit TRANS} finish{FIN}}}
#end
#macro constrain(LOC, R, COL, TRANS FIN)
union{torus{R, Rcell rotate 45*z texture{pigment{color COL transmit TRANS} finish{FIN}}}
      torus{R, Rcell rotate -45*z texture{pigment{color COL transmit TRANS} finish{FIN}}}
      translate LOC}
#end

cylinder {<-10.70, -32.52, -33.83>, <  0.13, -36.23, -23.65>, Rcell pigment {Black}}
cylinder {<  0.13, -28.82, -44.00>, < 10.96, -32.52, -33.83>, Rcell pigment {Black}}
cylinder {<  0.13,  36.23, -20.33>, < 10.96,  32.52, -10.15>, Rcell pigment {Black}}
cylinder {<-10.70,  32.52, -10.15>, <  0.13,  28.82,   0.03>, Rcell pigment {Black}}
cylinder {<-10.70, -32.52, -33.83>, <  0.13, -28.82, -44.00>, Rcell pigment {Black}}
cylinder {<  0.13, -36.23, -23.65>, < 10.96, -32.52, -33.83>, Rcell pigment {Black}}
cylinder {<  0.13,  28.82,   0.03>, < 10.96,  32.52, -10.15>, Rcell pigment {Black}}
cylinder {<-10.70,  32.52, -10.15>, <  0.13,  36.23, -20.33>, Rcell pigment {Black}}
cylinder {<-10.70, -32.52, -33.83>, <-10.70,  32.52, -10.15>, Rcell pigment {Black}}
cylinder {<  0.13, -36.23, -23.65>, <  0.13,  28.82,   0.03>, Rcell pigment {Black}}
cylinder {< 10.96, -32.52, -33.83>, < 10.96,  32.52, -10.15>, Rcell pigment {Black}}
cylinder {<  0.13, -28.82, -44.00>, <  0.13,  36.23, -20.33>, Rcell pigment {Black}}
atom(<-10.70, -20.31, -29.38>, 1.32, rgb <0.78, 0.50, 0.20>, 0.0, ase3) // #0 
atom(< -8.90, -20.92, -27.68>, 1.32, rgb <0.78, 0.50, 0.20>, 0.0, ase3) // #1 
atom(< -7.09, -21.54, -25.99>, 1.32, rgb <0.78, 0.50, 0.20>, 0.0, ase3) // #2 
atom(< -8.90, -19.69, -31.08>, 1.32, rgb <0.78, 0.50, 0.20>, 0.0, ase3) // #3 
atom(< -7.09, -20.31, -29.38>, 1.32, rgb <0.78, 0.50, 0.20>, 0.0, ase3) // #4 
atom(< -5.29, -20.92, -27.68>, 1.32, rgb <0.78, 0.50, 0.20>, 0.0, ase3) // #5 
atom(< -7.09, -19.07, -32.77>, 1.32, rgb <0.78, 0.50, 0.20>, 0.0, ase3) // #6 
atom(< -5.29, -19.69, -31.08>, 1.32, rgb <0.78, 0.50, 0.20>, 0.0, ase3) // #7 
atom(< -3.48, -20.31, -29.38>, 1.32, rgb <0.78, 0.50, 0.20>, 0.0, ase3) // #8 
atom(< -8.90, -18.61, -28.76>, 1.32, rgb <0.78, 0.50, 0.20>, 0.0, ase3) // #9 
atom(< -7.09, -19.23, -27.07>, 1.32, rgb <0.78, 0.50, 0.20>, 0.0, ase3) // #10 
atom(< -5.29, -19.85, -25.37>, 1.32, rgb <0.78, 0.50, 0.20>, 0.0, ase3) // #11 
atom(< -7.09, -17.99, -30.46>, 1.32, rgb <0.78, 0.50, 0.20>, 0.0, ase3) // #12 
atom(< -5.29, -18.61, -28.76>, 1.32, rgb <0.78, 0.50, 0.20>, 0.0, ase3) // #13 
atom(< -3.48, -19.23, -27.07>, 1.32, rgb <0.78, 0.50, 0.20>, 0.0, ase3) // #14 
atom(< -5.29, -17.38, -32.15>, 1.32, rgb <0.78, 0.50, 0.20>, 0.0, ase3) // #15 
atom(< -3.48, -17.99, -30.46>, 1.32, rgb <0.78, 0.50, 0.20>, 0.0, ase3) // #16 
atom(< -1.68, -18.61, -28.76>, 1.32, rgb <0.78, 0.50, 0.20>, 0.0, ase3) // #17 
atom(<  9.14, -15.22, -26.12>, 1.32, rgb <0.78, 0.50, 0.20>, 0.0, ase3) // #18 
atom(< -8.76, -17.30, -26.49>, 1.32, rgb <0.78, 0.50, 0.20>, 0.0, ase3) // #19 
atom(<  3.77, -14.54, -34.86>, 1.32, rgb <0.78, 0.50, 0.20>, 0.0, ase3) // #20 
atom(<  1.23, -18.84, -19.65>, 1.32, rgb <0.78, 0.50, 0.20>, 0.0, ase3) // #21 
atom(< -6.82, -16.84, -27.61>, 1.32, rgb <0.78, 0.50, 0.20>, 0.0, ase3) // #22 
atom(< -5.37, -18.33, -26.77>, 1.32, rgb <0.78, 0.50, 0.20>, 0.0, ase3) // #23 
atom(<  5.29, -13.01, -20.96>, 1.32, rgb <0.78, 0.50, 0.20>, 0.0, ase3) // #24 
atom(< -4.86, -16.34, -29.59>, 1.32, rgb <0.78, 0.50, 0.20>, 0.0, ase3) // #25 
atom(< -3.90, -16.75, -27.74>, 1.32, rgb <0.78, 0.50, 0.20>, 0.0, ase3) // #26 
atom(< -0.46,  -9.61, -28.20>, 0.66, rgb <1.00, 0.05, 0.05>, 0.0, ase3) // #27 
atom(< -0.75, -17.27, -25.15>, 0.76, rgb <0.56, 0.56, 0.56>, 0.0, ase3) // #28 
atom(<-10.70,  12.22, -17.54>, 1.32, rgb <0.78, 0.50, 0.20>, 0.0, ase3) // #29 
atom(< -8.90,  11.60, -15.85>, 1.32, rgb <0.78, 0.50, 0.20>, 0.0, ase3) // #30 
atom(< -7.09,  10.98, -14.15>, 1.32, rgb <0.78, 0.50, 0.20>, 0.0, ase3) // #31 
atom(< -8.90,  12.83, -19.24>, 1.32, rgb <0.78, 0.50, 0.20>, 0.0, ase3) // #32 
atom(< -7.09,  12.22, -17.54>, 1.32, rgb <0.78, 0.50, 0.20>, 0.0, ase3) // #33 
atom(< -5.29,  11.60, -15.85>, 1.32, rgb <0.78, 0.50, 0.20>, 0.0, ase3) // #34 
atom(< -7.09,  13.45, -20.93>, 1.32, rgb <0.78, 0.50, 0.20>, 0.0, ase3) // #35 
atom(< -5.29,  12.83, -19.24>, 1.32, rgb <0.78, 0.50, 0.20>, 0.0, ase3) // #36 
atom(< -3.48,  12.22, -17.54>, 1.32, rgb <0.78, 0.50, 0.20>, 0.0, ase3) // #37 
atom(< -8.90,  13.91, -16.92>, 1.32, rgb <0.78, 0.50, 0.20>, 0.0, ase3) // #38 
atom(< -7.09,  13.29, -15.23>, 1.32, rgb <0.78, 0.50, 0.20>, 0.0, ase3) // #39 
atom(< -5.29,  12.68, -13.53>, 1.32, rgb <0.78, 0.50, 0.20>, 0.0, ase3) // #40 
atom(< -7.09,  14.53, -18.62>, 1.32, rgb <0.78, 0.50, 0.20>, 0.0, ase3) // #41 
atom(< -5.29,  13.91, -16.92>, 1.32, rgb <0.78, 0.50, 0.20>, 0.0, ase3) // #42 
atom(< -3.48,  13.29, -15.23>, 1.32, rgb <0.78, 0.50, 0.20>, 0.0, ase3) // #43 
atom(< -5.29,  15.15, -20.32>, 1.32, rgb <0.78, 0.50, 0.20>, 0.0, ase3) // #44 
atom(< -3.48,  14.53, -18.62>, 1.32, rgb <0.78, 0.50, 0.20>, 0.0, ase3) // #45 
atom(< -1.68,  13.91, -16.92>, 1.32, rgb <0.78, 0.50, 0.20>, 0.0, ase3) // #46 
atom(<  9.14,  17.30, -14.29>, 1.32, rgb <0.78, 0.50, 0.20>, 0.0, ase3) // #47 
atom(< -8.76,  15.23, -14.65>, 1.32, rgb <0.78, 0.50, 0.20>, 0.0, ase3) // #48 
atom(<  3.77,  17.99, -23.02>, 1.32, rgb <0.78, 0.50, 0.20>, 0.0, ase3) // #49 
atom(<  1.23,  13.68,  -7.81>, 1.32, rgb <0.78, 0.50, 0.20>, 0.0, ase3) // #50 
atom(< -6.82,  15.69, -15.77>, 1.32, rgb <0.78, 0.50, 0.20>, 0.0, ase3) // #51 
atom(< -5.37,  14.19, -14.93>, 1.32, rgb <0.78, 0.50, 0.20>, 0.0, ase3) // #52 
atom(<  5.29,  19.52,  -9.13>, 1.32, rgb <0.78, 0.50, 0.20>, 0.0, ase3) // #53 
atom(< -4.86,  16.19, -17.76>, 1.32, rgb <0.78, 0.50, 0.20>, 0.0, ase3) // #54 
atom(< -3.90,  15.78, -15.91>, 1.32, rgb <0.78, 0.50, 0.20>, 0.0, ase3) // #55 
atom(< -0.46,  22.92, -16.36>, 0.66, rgb <1.00, 0.05, 0.05>, 0.0, ase3) // #56 
atom(< -0.75,  15.25, -13.31>, 0.76, rgb <0.56, 0.56, 0.56>, 0.0, ase3) // #57 
atom(< -5.29, -18.45, -34.47>, 1.32, rgb <0.78, 0.50, 0.20>, 0.0, ase3) // #58 
atom(< -3.48, -19.07, -32.77>, 1.32, rgb <0.78, 0.50, 0.20>, 0.0, ase3) // #59 
atom(< -1.68, -19.69, -31.08>, 1.32, rgb <0.78, 0.50, 0.20>, 0.0, ase3) // #60 
atom(< -3.48, -17.84, -36.16>, 1.32, rgb <0.78, 0.50, 0.20>, 0.0, ase3) // #61 
atom(< -1.68, -18.45, -34.47>, 1.32, rgb <0.78, 0.50, 0.20>, 0.0, ase3) // #62 
atom(<  0.13, -19.07, -32.77>, 1.32, rgb <0.78, 0.50, 0.20>, 0.0, ase3) // #63 
atom(< -1.68, -17.22, -37.86>, 1.32, rgb <0.78, 0.50, 0.20>, 0.0, ase3) // #64 
atom(<  0.13, -17.84, -36.16>, 1.32, rgb <0.78, 0.50, 0.20>, 0.0, ase3) // #65 
atom(<  1.93, -18.45, -34.47>, 1.32, rgb <0.78, 0.50, 0.20>, 0.0, ase3) // #66 
atom(< -3.48, -16.76, -33.85>, 1.32, rgb <0.78, 0.50, 0.20>, 0.0, ase3) // #67 
atom(< -1.68, -17.38, -32.15>, 1.32, rgb <0.78, 0.50, 0.20>, 0.0, ase3) // #68 
atom(<  0.13, -17.99, -30.46>, 1.32, rgb <0.78, 0.50, 0.20>, 0.0, ase3) // #69 
atom(< -1.68, -16.14, -35.55>, 1.32, rgb <0.78, 0.50, 0.20>, 0.0, ase3) // #70 
atom(<  0.13, -16.76, -33.85>, 1.32, rgb <0.78, 0.50, 0.20>, 0.0, ase3) // #71 
atom(<  1.93, -17.38, -32.15>, 1.32, rgb <0.78, 0.50, 0.20>, 0.0, ase3) // #72 
atom(<  0.13, -15.52, -37.24>, 1.32, rgb <0.78, 0.50, 0.20>, 0.0, ase3) // #73 
atom(<  1.93, -16.14, -35.55>, 1.32, rgb <0.78, 0.50, 0.20>, 0.0, ase3) // #74 
atom(<  3.74, -16.76, -33.85>, 1.32, rgb <0.78, 0.50, 0.20>, 0.0, ase3) // #75 
atom(<  3.73, -17.07, -21.04>, 1.32, rgb <0.78, 0.50, 0.20>, 0.0, ase3) // #76 
atom(< -3.35, -15.45, -31.58>, 1.32, rgb <0.78, 0.50, 0.20>, 0.0, ase3) // #77 
atom(< -1.65, -16.39, -29.77>, 1.32, rgb <0.78, 0.50, 0.20>, 0.0, ase3) // #78 
atom(<  6.65, -16.99, -24.74>, 1.32, rgb <0.78, 0.50, 0.20>, 0.0, ase3) // #79 
atom(< -1.40, -14.99, -32.70>, 1.32, rgb <0.78, 0.50, 0.20>, 0.0, ase3) // #80 
atom(<  0.04, -16.48, -31.86>, 1.32, rgb <0.78, 0.50, 0.20>, 0.0, ase3) // #81 
atom(< 10.70, -11.16, -26.05>, 1.32, rgb <0.78, 0.50, 0.20>, 0.0, ase3) // #82 
atom(<  0.56, -14.48, -34.68>, 1.32, rgb <0.78, 0.50, 0.20>, 0.0, ase3) // #83 
atom(<  1.51, -14.89, -32.83>, 1.32, rgb <0.78, 0.50, 0.20>, 0.0, ase3) // #84 
atom(< -5.87, -11.46, -23.11>, 0.66, rgb <1.00, 0.05, 0.05>, 0.0, ase3) // #85 
atom(<  4.66, -15.42, -30.24>, 0.76, rgb <0.56, 0.56, 0.56>, 0.0, ase3) // #86 
atom(< -5.29,  14.07, -22.63>, 1.32, rgb <0.78, 0.50, 0.20>, 0.0, ase3) // #87 
atom(< -3.48,  13.45, -20.93>, 1.32, rgb <0.78, 0.50, 0.20>, 0.0, ase3) // #88 
atom(< -1.68,  12.83, -19.24>, 1.32, rgb <0.78, 0.50, 0.20>, 0.0, ase3) // #89 
atom(< -3.48,  14.69, -24.33>, 1.32, rgb <0.78, 0.50, 0.20>, 0.0, ase3) // #90 
atom(< -1.68,  14.07, -22.63>, 1.32, rgb <0.78, 0.50, 0.20>, 0.0, ase3) // #91 
atom(<  0.13,  13.45, -20.93>, 1.32, rgb <0.78, 0.50, 0.20>, 0.0, ase3) // #92 
atom(< -1.68,  15.30, -26.02>, 1.32, rgb <0.78, 0.50, 0.20>, 0.0, ase3) // #93 
atom(<  0.13,  14.69, -24.33>, 1.32, rgb <0.78, 0.50, 0.20>, 0.0, ase3) // #94 
atom(<  1.93,  14.07, -22.63>, 1.32, rgb <0.78, 0.50, 0.20>, 0.0, ase3) // #95 
atom(< -3.48,  15.76, -22.01>, 1.32, rgb <0.78, 0.50, 0.20>, 0.0, ase3) // #96 
atom(< -1.68,  15.15, -20.32>, 1.32, rgb <0.78, 0.50, 0.20>, 0.0, ase3) // #97 
atom(<  0.13,  14.53, -18.62>, 1.32, rgb <0.78, 0.50, 0.20>, 0.0, ase3) // #98 
atom(< -1.68,  16.38, -23.71>, 1.32, rgb <0.78, 0.50, 0.20>, 0.0, ase3) // #99 
atom(<  0.13,  15.76, -22.01>, 1.32, rgb <0.78, 0.50, 0.20>, 0.0, ase3) // #100 
atom(<  1.93,  15.15, -20.32>, 1.32, rgb <0.78, 0.50, 0.20>, 0.0, ase3) // #101 
atom(<  0.13,  17.00, -25.41>, 1.32, rgb <0.78, 0.50, 0.20>, 0.0, ase3) // #102 
atom(<  1.93,  16.38, -23.71>, 1.32, rgb <0.78, 0.50, 0.20>, 0.0, ase3) // #103 
atom(<  3.74,  15.76, -22.01>, 1.32, rgb <0.78, 0.50, 0.20>, 0.0, ase3) // #104 
atom(<  3.73,  15.45,  -9.20>, 1.32, rgb <0.78, 0.50, 0.20>, 0.0, ase3) // #105 
atom(< -3.35,  17.08, -19.74>, 1.32, rgb <0.78, 0.50, 0.20>, 0.0, ase3) // #106 
atom(< -1.65,  16.13, -17.93>, 1.32, rgb <0.78, 0.50, 0.20>, 0.0, ase3) // #107 
atom(<  6.65,  15.53, -12.90>, 1.32, rgb <0.78, 0.50, 0.20>, 0.0, ase3) // #108 
atom(< -1.40,  17.54, -20.86>, 1.32, rgb <0.78, 0.50, 0.20>, 0.0, ase3) // #109 
atom(<  0.04,  16.04, -20.02>, 1.32, rgb <0.78, 0.50, 0.20>, 0.0, ase3) // #110 
atom(< 10.70,  21.37, -14.21>, 1.32, rgb <0.78, 0.50, 0.20>, 0.0, ase3) // #111 
atom(<  0.56,  18.04, -22.85>, 1.32, rgb <0.78, 0.50, 0.20>, 0.0, ase3) // #112 
atom(<  1.51,  17.63, -21.00>, 1.32, rgb <0.78, 0.50, 0.20>, 0.0, ase3) // #113 
atom(< -5.87,  21.06, -11.28>, 0.66, rgb <1.00, 0.05, 0.05>, 0.0, ase3) // #114 
atom(<  4.66,  17.10, -18.40>, 0.76, rgb <0.56, 0.56, 0.56>, 0.0, ase3) // #115 
atom(< -5.29, -22.16, -24.29>, 1.32, rgb <0.78, 0.50, 0.20>, 0.0, ase3) // #116 
atom(< -3.48, -22.78, -22.59>, 1.32, rgb <0.78, 0.50, 0.20>, 0.0, ase3) // #117 
atom(< -1.68, -23.39, -20.90>, 1.32, rgb <0.78, 0.50, 0.20>, 0.0, ase3) // #118 
atom(< -3.48, -21.54, -25.99>, 1.32, rgb <0.78, 0.50, 0.20>, 0.0, ase3) // #119 
atom(< -1.68, -22.16, -24.29>, 1.32, rgb <0.78, 0.50, 0.20>, 0.0, ase3) // #120 
atom(<  0.13, -22.78, -22.59>, 1.32, rgb <0.78, 0.50, 0.20>, 0.0, ase3) // #121 
atom(< -1.68, -20.92, -27.68>, 1.32, rgb <0.78, 0.50, 0.20>, 0.0, ase3) // #122 
atom(<  0.13, -21.54, -25.99>, 1.32, rgb <0.78, 0.50, 0.20>, 0.0, ase3) // #123 
atom(<  1.93, -22.16, -24.29>, 1.32, rgb <0.78, 0.50, 0.20>, 0.0, ase3) // #124 
atom(< -3.48, -20.46, -23.67>, 1.32, rgb <0.78, 0.50, 0.20>, 0.0, ase3) // #125 
atom(< -1.68, -21.08, -21.98>, 1.32, rgb <0.78, 0.50, 0.20>, 0.0, ase3) // #126 
atom(<  0.13, -21.70, -20.28>, 1.32, rgb <0.78, 0.50, 0.20>, 0.0, ase3) // #127 
atom(< -1.68, -19.85, -25.37>, 1.32, rgb <0.78, 0.50, 0.20>, 0.0, ase3) // #128 
atom(<  0.13, -20.46, -23.67>, 1.32, rgb <0.78, 0.50, 0.20>, 0.0, ase3) // #129 
atom(<  1.93, -21.08, -21.98>, 1.32, rgb <0.78, 0.50, 0.20>, 0.0, ase3) // #130 
atom(<  0.13, -19.23, -27.07>, 1.32, rgb <0.78, 0.50, 0.20>, 0.0, ase3) // #131 
atom(<  1.93, -19.85, -25.37>, 1.32, rgb <0.78, 0.50, 0.20>, 0.0, ase3) // #132 
atom(<  3.74, -20.46, -23.67>, 1.32, rgb <0.78, 0.50, 0.20>, 0.0, ase3) // #133 
atom(<  3.73, -13.37, -31.21>, 1.32, rgb <0.78, 0.50, 0.20>, 0.0, ase3) // #134 
atom(< -3.35, -19.15, -21.40>, 1.32, rgb <0.78, 0.50, 0.20>, 0.0, ase3) // #135 
atom(<  9.18, -16.39, -29.77>, 1.32, rgb <0.78, 0.50, 0.20>, 0.0, ase3) // #136 
atom(< -4.18, -16.99, -24.74>, 1.32, rgb <0.78, 0.50, 0.20>, 0.0, ase3) // #137 
atom(< -1.40, -18.69, -22.52>, 1.32, rgb <0.78, 0.50, 0.20>, 0.0, ase3) // #138 
atom(<  0.04, -20.18, -21.68>, 1.32, rgb <0.78, 0.50, 0.20>, 0.0, ase3) // #139 
atom(< -0.13, -11.16, -26.05>, 1.32, rgb <0.78, 0.50, 0.20>, 0.0, ase3) // #140 
atom(<  0.56, -18.19, -24.51>, 1.32, rgb <0.78, 0.50, 0.20>, 0.0, ase3) // #141 
atom(<  1.51, -18.60, -22.66>, 1.32, rgb <0.78, 0.50, 0.20>, 0.0, ase3) // #142 
atom(<  4.96, -11.46, -23.11>, 0.66, rgb <1.00, 0.05, 0.05>, 0.0, ase3) // #143 
atom(< -6.17, -15.42, -30.24>, 0.76, rgb <0.56, 0.56, 0.56>, 0.0, ase3) // #144 
atom(< -5.29,  10.36, -12.45>, 1.32, rgb <0.78, 0.50, 0.20>, 0.0, ase3) // #145 
atom(< -3.48,   9.75, -10.76>, 1.32, rgb <0.78, 0.50, 0.20>, 0.0, ase3) // #146 
atom(< -1.68,   9.13,  -9.06>, 1.32, rgb <0.78, 0.50, 0.20>, 0.0, ase3) // #147 
atom(< -3.48,  10.98, -14.15>, 1.32, rgb <0.78, 0.50, 0.20>, 0.0, ase3) // #148 
atom(< -1.68,  10.36, -12.45>, 1.32, rgb <0.78, 0.50, 0.20>, 0.0, ase3) // #149 
atom(<  0.13,   9.75, -10.76>, 1.32, rgb <0.78, 0.50, 0.20>, 0.0, ase3) // #150 
atom(< -1.68,  11.60, -15.85>, 1.32, rgb <0.78, 0.50, 0.20>, 0.0, ase3) // #151 
atom(<  0.13,  10.98, -14.15>, 1.32, rgb <0.78, 0.50, 0.20>, 0.0, ase3) // #152 
atom(<  1.93,  10.36, -12.45>, 1.32, rgb <0.78, 0.50, 0.20>, 0.0, ase3) // #153 
atom(< -3.48,  12.06, -11.84>, 1.32, rgb <0.78, 0.50, 0.20>, 0.0, ase3) // #154 
atom(< -1.68,  11.44, -10.14>, 1.32, rgb <0.78, 0.50, 0.20>, 0.0, ase3) // #155 
atom(<  0.13,  10.83,  -8.44>, 1.32, rgb <0.78, 0.50, 0.20>, 0.0, ase3) // #156 
atom(< -1.68,  12.68, -13.53>, 1.32, rgb <0.78, 0.50, 0.20>, 0.0, ase3) // #157 
atom(<  0.13,  12.06, -11.84>, 1.32, rgb <0.78, 0.50, 0.20>, 0.0, ase3) // #158 
atom(<  1.93,  11.44, -10.14>, 1.32, rgb <0.78, 0.50, 0.20>, 0.0, ase3) // #159 
atom(<  0.13,  13.29, -15.23>, 1.32, rgb <0.78, 0.50, 0.20>, 0.0, ase3) // #160 
atom(<  1.93,  12.68, -13.53>, 1.32, rgb <0.78, 0.50, 0.20>, 0.0, ase3) // #161 
atom(<  3.74,  12.06, -11.84>, 1.32, rgb <0.78, 0.50, 0.20>, 0.0, ase3) // #162 
atom(<  3.73,  19.16, -19.37>, 1.32, rgb <0.78, 0.50, 0.20>, 0.0, ase3) // #163 
atom(< -3.35,  13.37,  -9.56>, 1.32, rgb <0.78, 0.50, 0.20>, 0.0, ase3) // #164 
atom(<  9.18,  16.13, -17.93>, 1.32, rgb <0.78, 0.50, 0.20>, 0.0, ase3) // #165 
atom(< -4.18,  15.53, -12.90>, 1.32, rgb <0.78, 0.50, 0.20>, 0.0, ase3) // #166 
atom(< -1.40,  13.83, -10.68>, 1.32, rgb <0.78, 0.50, 0.20>, 0.0, ase3) // #167 
atom(<  0.04,  12.34,  -9.84>, 1.32, rgb <0.78, 0.50, 0.20>, 0.0, ase3) // #168 
atom(< -0.13,  21.37, -14.21>, 1.32, rgb <0.78, 0.50, 0.20>, 0.0, ase3) // #169 
atom(<  0.56,  14.33, -12.67>, 1.32, rgb <0.78, 0.50, 0.20>, 0.0, ase3) // #170 
atom(<  1.51,  13.92, -10.82>, 1.32, rgb <0.78, 0.50, 0.20>, 0.0, ase3) // #171 
atom(<  4.96,  21.06, -11.28>, 0.66, rgb <1.00, 0.05, 0.05>, 0.0, ase3) // #172 
atom(< -6.17,  17.10, -18.40>, 0.76, rgb <0.56, 0.56, 0.56>, 0.0, ase3) // #173 
atom(<  0.13, -20.31, -29.38>, 1.32, rgb <0.78, 0.50, 0.20>, 0.0, ase3) // #174 
atom(<  1.93, -20.92, -27.68>, 1.32, rgb <0.78, 0.50, 0.20>, 0.0, ase3) // #175 
atom(<  3.74, -21.54, -25.99>, 1.32, rgb <0.78, 0.50, 0.20>, 0.0, ase3) // #176 
atom(<  1.93, -19.69, -31.08>, 1.32, rgb <0.78, 0.50, 0.20>, 0.0, ase3) // #177 
atom(<  3.74, -20.31, -29.38>, 1.32, rgb <0.78, 0.50, 0.20>, 0.0, ase3) // #178 
atom(<  5.54, -20.92, -27.68>, 1.32, rgb <0.78, 0.50, 0.20>, 0.0, ase3) // #179 
atom(<  3.74, -19.07, -32.77>, 1.32, rgb <0.78, 0.50, 0.20>, 0.0, ase3) // #180 
atom(<  5.54, -19.69, -31.08>, 1.32, rgb <0.78, 0.50, 0.20>, 0.0, ase3) // #181 
atom(<  7.35, -20.31, -29.38>, 1.32, rgb <0.78, 0.50, 0.20>, 0.0, ase3) // #182 
atom(<  1.93, -18.61, -28.76>, 1.32, rgb <0.78, 0.50, 0.20>, 0.0, ase3) // #183 
atom(<  3.74, -19.23, -27.07>, 1.32, rgb <0.78, 0.50, 0.20>, 0.0, ase3) // #184 
atom(<  5.54, -19.85, -25.37>, 1.32, rgb <0.78, 0.50, 0.20>, 0.0, ase3) // #185 
atom(<  3.74, -17.99, -30.46>, 1.32, rgb <0.78, 0.50, 0.20>, 0.0, ase3) // #186 
atom(<  5.54, -18.61, -28.76>, 1.32, rgb <0.78, 0.50, 0.20>, 0.0, ase3) // #187 
atom(<  7.35, -19.23, -27.07>, 1.32, rgb <0.78, 0.50, 0.20>, 0.0, ase3) // #188 
atom(<  5.54, -17.38, -32.15>, 1.32, rgb <0.78, 0.50, 0.20>, 0.0, ase3) // #189 
atom(<  7.35, -17.99, -30.46>, 1.32, rgb <0.78, 0.50, 0.20>, 0.0, ase3) // #190 
atom(<  9.15, -18.61, -28.76>, 1.32, rgb <0.78, 0.50, 0.20>, 0.0, ase3) // #191 
atom(< -1.69, -15.22, -26.12>, 1.32, rgb <0.78, 0.50, 0.20>, 0.0, ase3) // #192 
atom(<  2.07, -17.30, -26.49>, 1.32, rgb <0.78, 0.50, 0.20>, 0.0, ase3) // #193 
atom(<  3.77, -18.24, -24.68>, 1.32, rgb <0.78, 0.50, 0.20>, 0.0, ase3) // #194 
atom(<  1.23, -15.14, -29.82>, 1.32, rgb <0.78, 0.50, 0.20>, 0.0, ase3) // #195 
atom(<  4.01, -16.84, -27.61>, 1.32, rgb <0.78, 0.50, 0.20>, 0.0, ase3) // #196 
atom(<  5.46, -18.33, -26.77>, 1.32, rgb <0.78, 0.50, 0.20>, 0.0, ase3) // #197 
atom(<  5.29,  -9.30, -31.14>, 1.32, rgb <0.78, 0.50, 0.20>, 0.0, ase3) // #198 
atom(<  5.97, -16.34, -29.59>, 1.32, rgb <0.78, 0.50, 0.20>, 0.0, ase3) // #199 
atom(<  6.93, -16.75, -27.74>, 1.32, rgb <0.78, 0.50, 0.20>, 0.0, ase3) // #200 
atom(< -0.46, -13.31, -18.02>, 0.66, rgb <1.00, 0.05, 0.05>, 0.0, ase3) // #201 
atom(< -0.75, -13.57, -35.32>, 0.76, rgb <0.56, 0.56, 0.56>, 0.0, ase3) // #202 
atom(<  0.13,  12.22, -17.54>, 1.32, rgb <0.78, 0.50, 0.20>, 0.0, ase3) // #203 
atom(<  1.93,  11.60, -15.85>, 1.32, rgb <0.78, 0.50, 0.20>, 0.0, ase3) // #204 
atom(<  3.74,  10.98, -14.15>, 1.32, rgb <0.78, 0.50, 0.20>, 0.0, ase3) // #205 
atom(<  1.93,  12.83, -19.24>, 1.32, rgb <0.78, 0.50, 0.20>, 0.0, ase3) // #206 
atom(<  3.74,  12.22, -17.54>, 1.32, rgb <0.78, 0.50, 0.20>, 0.0, ase3) // #207 
atom(<  5.54,  11.60, -15.85>, 1.32, rgb <0.78, 0.50, 0.20>, 0.0, ase3) // #208 
atom(<  3.74,  13.45, -20.93>, 1.32, rgb <0.78, 0.50, 0.20>, 0.0, ase3) // #209 
atom(<  5.54,  12.83, -19.24>, 1.32, rgb <0.78, 0.50, 0.20>, 0.0, ase3) // #210 
atom(<  7.35,  12.22, -17.54>, 1.32, rgb <0.78, 0.50, 0.20>, 0.0, ase3) // #211 
atom(<  1.93,  13.91, -16.92>, 1.32, rgb <0.78, 0.50, 0.20>, 0.0, ase3) // #212 
atom(<  3.74,  13.29, -15.23>, 1.32, rgb <0.78, 0.50, 0.20>, 0.0, ase3) // #213 
atom(<  5.54,  12.68, -13.53>, 1.32, rgb <0.78, 0.50, 0.20>, 0.0, ase3) // #214 
atom(<  3.74,  14.53, -18.62>, 1.32, rgb <0.78, 0.50, 0.20>, 0.0, ase3) // #215 
atom(<  5.54,  13.91, -16.92>, 1.32, rgb <0.78, 0.50, 0.20>, 0.0, ase3) // #216 
atom(<  7.35,  13.29, -15.23>, 1.32, rgb <0.78, 0.50, 0.20>, 0.0, ase3) // #217 
atom(<  5.54,  15.15, -20.32>, 1.32, rgb <0.78, 0.50, 0.20>, 0.0, ase3) // #218 
atom(<  7.35,  14.53, -18.62>, 1.32, rgb <0.78, 0.50, 0.20>, 0.0, ase3) // #219 
atom(<  9.15,  13.91, -16.92>, 1.32, rgb <0.78, 0.50, 0.20>, 0.0, ase3) // #220 
atom(< -1.69,  17.30, -14.29>, 1.32, rgb <0.78, 0.50, 0.20>, 0.0, ase3) // #221 
atom(<  2.07,  15.23, -14.65>, 1.32, rgb <0.78, 0.50, 0.20>, 0.0, ase3) // #222 
atom(<  3.77,  14.28, -12.85>, 1.32, rgb <0.78, 0.50, 0.20>, 0.0, ase3) // #223 
atom(<  1.23,  17.38, -17.99>, 1.32, rgb <0.78, 0.50, 0.20>, 0.0, ase3) // #224 
atom(<  4.01,  15.69, -15.77>, 1.32, rgb <0.78, 0.50, 0.20>, 0.0, ase3) // #225 
atom(<  5.46,  14.19, -14.93>, 1.32, rgb <0.78, 0.50, 0.20>, 0.0, ase3) // #226 
atom(<  5.29,  23.22, -19.30>, 1.32, rgb <0.78, 0.50, 0.20>, 0.0, ase3) // #227 
atom(<  5.97,  16.19, -17.76>, 1.32, rgb <0.78, 0.50, 0.20>, 0.0, ase3) // #228 
atom(<  6.93,  15.78, -15.91>, 1.32, rgb <0.78, 0.50, 0.20>, 0.0, ase3) // #229 
atom(< -0.46,  19.21,  -6.19>, 0.66, rgb <1.00, 0.05, 0.05>, 0.0, ase3) // #230 
atom(< -0.75,  18.95, -23.49>, 0.76, rgb <0.56, 0.56, 0.56>, 0.0, ase3) // #231 
