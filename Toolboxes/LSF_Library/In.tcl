# define DAMPING--------------------------------------------------------------------------------------
# apply Rayleigh DAMPING from $xDamp
# D=$alphaM*M + $betaKcurr*Kcurrent + $betaKcomm*KlastCommit + $beatKinit*$Kinitial
set xDamp 0.05;				# 5% damping ratio
set lambdaN [eigen 2]; 			# eigenvalue 
set lambdaI [lindex $lambdaN [expr 0]];
set lambdaJ [lindex $lambdaN [expr 1]];
set w1 [expr pow($lambdaI,0.5)];					# w1 (1st mode circular frequency)
set w2 [expr pow($lambdaJ,0.5)];					# w2 (2nd mode circular frequency)
set alphaM    [expr $xDamp*(2.*$w1*$w2)/($w1+$w2)];       # M-prop. damping; D = alphaM*M
set betaKcurr 0.;         			# K-proportional damping;      +beatKcurr*KCurrent
set betaKcomm [expr $xDamp*2./($w1+$w2)];   	# K-prop. damping parameter;   +betaKcomm*KlastCommitt
set betaKinit 0.;         			# initial-stiffness proportional damping      +beatKinit*Kini
rayleigh $alphaM $betaKcurr $betaKinit $betaKcomm; 		# RAYLEIGH damping

# variables for ground motion
set dt 0.01;
set Tol 1.0e-8;
set DtAnalysis	$dt;	# time-step Dt for lateral analysis
set TmaxAnalysis 10;	

# define acceleration loading 
set accelSeries "Series -dt $dt -filePath accel.txt -factor 1";

pattern UniformExcitation 1 1 -accel $accelSeries

#pattern MultipleSupport 1  {	
#		             
#		groundMotion 1 Plain -disp $accelSeries;   
#	     	imposedMotion 1  1 1
#	
#};	

source LibAnalysisDynamicParameters.tcl;

set Nsteps [expr int($TmaxAnalysis/$DtAnalysis)];
set ok [analyze $Nsteps $DtAnalysis];			# actually perform analysis; returns ok=0 if analysis was successful	 
                                                  

if {$ok != 0} {      ;					# analysis was not successful.
	# --------------------------------------------------------------------------------------------------
	# change some analysis parameters to achieve convergence
	# performance is slower inside this loop
	#    Time-controlled analysis
	set ok 0;
	set controlTime [getTime];
	while {$controlTime < $TmaxAnalysis && $ok == 0} {
		set controlTime [getTime]
		set ok [analyze 1 $DtAnalysis]
		if {$ok != 0} {
			puts "Trying Newton with Initial Tangent .."
			test NormDispIncr   $Tol 1000  0
			algorithm Newton -initial
			set ok [analyze 1 $DtAnalysis]
			test $testTypeDynamic $TolDynamic $maxNumIterDynamic  0
			algorithm $algorithmTypeDynamic
		}
		if {$ok != 0} {
			puts "Trying Broyden .."
			algorithm Broyden 8
			set ok [analyze 1 $DtAnalysis]
			algorithm $algorithmTypeDynamic
		}
		if {$ok != 0} {
			puts "Trying NewtonWithLineSearch .."
			algorithm NewtonLineSearch .8
			set ok [analyze 1 $DtAnalysis]
			algorithm $algorithmTypeDynamic
		}
	}
};      # end if ok !0
