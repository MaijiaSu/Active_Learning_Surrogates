wipe						                                                      
model basic -ndm 1 -ndf 1				                                                  
if { [file exists output] == 0 } {                                                         
  file mkdir output
}
  
#-----define geometry
node 1 0      				
node 2 1   
node 3 2
node 4 3
#node 5 4
#node 6 5
#node 7 6

#-----define boundary condition
fix 1 1 

#-----define mass
mass 2 1.0e6
mass 3 1.0e6
mass 4 1.0e6
#mass 5 30000
#mass 6 30000
#mass 7 30000 

#-----define material
#uniaxialMaterial Elastic      1        10;
#uniaxialMaterial Hardening      1        10 0.1 0.05 0;

#uniaxialMaterial KikuchiAikenHDR 1 X0.6 1.0 1.0;

uniaxialMaterial BoucWen 1 .1 3.0e6 5 [expr 1/(2.0*pow(0.04,5))] [expr 1/(2.0*pow(0.04,5))] 1 0 0 0;
uniaxialMaterial BoucWen 2 .1 2.8e6 5 [expr 1/(2.0*pow(0.04,5))] [expr 1/(2.0*pow(0.04,5))] 1 0 0 0;
uniaxialMaterial BoucWen 3 .1 1.5e6 5 [expr 1/(2.0*pow(0.04,5))] [expr 1/(2.0*pow(0.04,5))] 1 0 0 0;
#uniaxialMaterial BoucWen 4 .1 5.5e7 5 [expr 1/(2.0*pow(0.01,5))] [expr 1/(2.0*pow(0.01,5))] 1 0 0 0;
#uniaxialMaterial BoucWen 5 .1 4.0e7 5 [expr 1/(2.0*pow(0.01,5))] [expr 1/(2.0*pow(0.01,5))] 1 0 0 0;
#uniaxialMaterial BoucWen 6 .1 2.0e7 5 [expr 1/(2.0*pow(0.01,5))] [expr 1/(2.0*pow(0.01,5))] 1 0 0 0;

#uniaxialMaterial Steel01 1 12.0e4 3.0e6 0.1;
#uniaxialMaterial Steel01 2 11.2e4 2.8e6 0.1;
#uniaxialMaterial Steel01 3 6.0e4 1.5e6 0.1;
#uniaxialMaterial Steel01 4 5.5e5 5.5e7 0.1;
#uniaxialMaterial Steel01 5 4.0e5 4.0e7 0.1;
#uniaxialMaterial Steel01 6 2.0e5 2.0e7 0.1;


#-----define element  
element truss 1 1 2 1 1 -doRayleigh 1
element truss 2 2 3 1 2 -doRayleigh 1
element truss 3 3 4 1 3 -doRayleigh 1
#element truss 4 4 5 1 4 -doRayleigh 1
#element truss 5 5 6 1 5 -doRayleigh 1
#element truss 6 6 7 1 6 -doRayleigh 1


#-----define recorder
recorder Node    -file output/disp_1.out     -precision 15          -node  2       -dof 1      disp;
recorder Drift    -file output/disp_2.out     -precision 15        -iNode  2  -jNode  3     -dof 1      -perpDirn 1;	
recorder Drift    -file output/disp_3.out     -precision 15        -iNode  3  -jNode  4     -dof 1      -perpDirn 1;	
#recorder Drift    -file output/disp_4.out     -precision 15        -iNode  4  -jNode  5     -dof 1      -perpDirn 1;	
#recorder Drift    -file output/disp_5.out     -precision 15        -iNode  5  -jNode  6     -dof 1      -perpDirn 1;	
#recorder Drift    -file output/disp_6.out     -precision 15        -iNode  6  -jNode  7     -dof 1      -perpDirn 1;
recorder Element  -file output/force_1.out     -precision 15        -ele   1       axialForce;
recorder Element  -file output/force_3.out     -precision 15        -ele   3       axialForce;

#recorder Node     -file output/disp_2.out     -precision 15        -node    3     -dof 1      disp;	
#recorder Node     -file output/disp_3.out     -precision 15        -node    4     -dof 1      disp;	
#recorder Node     -file output/disp_4.out     -precision 15        -node    5     -dof 1      disp;
#recorder Node     -file output/disp_5.out     -precision 15        -node    6     -dof 1      disp;	
#recorder Node     -file output/disp_6.out     -precision 15        -node    7     -dof 1      disp;
	
#recorder Node    -file output/accel_1.out    -precision 15        -node  2       -dof 1      accel;
#recorder Node    -file output/accel_2.out    -precision 15        -node  3       -dof 1      accel;
#recorder Node    -file output/accel_3.out    -precision 15        -node  4       -dof 1      accel;
#recorder Node    -file output/accel_4.out    -precision 15        -node  5       -dof 1      accel;
#recorder Node    -file output/accel_5.out    -precision 15        -node  6       -dof 1      accel;
#recorder Node    -file output/accel_6.out    -precision 15        -node  7       -dof 1      accel;

#recorder Node    -file output/vel_1.out    -precision 15        -node  2       -dof 1      vel;
#recorder Node    -file output/vel_2.out    -precision 15        -node  3       -dof 1      vel;
#recorder Node    -file output/vel_3.out    -precision 15        -node  4       -dof 1      vel;
#recorder Node    -file output/vel_4.out    -precision 15        -node  5       -dof 1      vel;
#recorder Node    -file output/vel_5.out    -precision 15        -node  6       -dof 1      vel;
#recorder Node    -file output/vel_6.out    -precision 15        -node  7       -dof 1      vel;

#recorder Node    -file output/accel_1.out    -precision 15        -node  2       -dof 1      accel;
#recorder Node    -file output/disp_1.out     -precision 15       -time    -node  2       -dof 1      disp;	
#recorder Node    -file output/accel_1.out    -precision 15       -time    -node  2       -dof 1      accel;

#pattern Plain 1 Linear {
#load 2  1;
#load 3  1;
#load 4  1;
#load 5  1;
#load 6  1;
#load 7  1;

#load 2  142596;
#load 3  142596;
#load 4  142596;
#load 5  142596;
#load 6  142596;
#load 7  142596;
#}

#constraints Plain;     				
#numberer Plain;					
#system BandGeneral;				
#test NormDispIncr 1.0e-8 1000; 				
#algorithm Newton;					
#integrator LoadControl 0.01;
#integrator DisplacementControl 2 1 0.0001;					
#analysis Static;				
#analyze 181;	

source in.tcl







