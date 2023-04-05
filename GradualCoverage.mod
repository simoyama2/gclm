#param t0;
#let t0 := time();

param clients integer;
param locations integer; 
param radii integer; 

set I := 1..clients by 1;
set J := 1..locations by 1;
set R := 1..radii by 1;

param rho{R}; 
param distance{I,J}; 
param alpha{I};
param weights{I}; 
param populations{I};
param b; 
param dmin; 
param dmax; 
param prob{i in I, j in J, r in R} = 
	if distance[i,j] <= dmin then 1
	else if distance[i,j] > dmax then 0
	else 1 / (1+10**((distance[i,j] - dmin) / (rho[r] - dmin)-1/b) );

param budget integer; #NEW



var x{j in J,r in R} binary;
var P{I}; 
var devplus{I} >= 0;
var devminus{I} >=0;

minimize z: sum{i in I} weights[i]*devminus[i];

subject to probability{i in I}:
	 (1 - (prod{j in J,r in R} (1 - prob[i,j,r] * x[j,r])))  = P[i];

subject to sums{j in J}:
	sum{r in R} x[j,r] <= 1;

subject to deviations{i in I}:
	P[i] + devminus[i] - devplus[i] = alpha[i];

subject to nequip: #NEW
	sum{j in J, r in R} x[j,r] = budget;

#subject to trava: #NEW
#	sum{j in J} y[j] >= x[j,r];



#subject to fermi {i in I, j in J, r in R: distance[i,j] <= dmin}:
#	prob[i,j,r]  = 1;

#subject to fermi2 {i in I, j in J, r in R: dmin < distance[i,j] <= dmax}:
#	prob[i,j,r] = 1 / (1+10**((distance[i,j] - dmin) / (rho[r] - dmin)-1/b) );
	
#subject to fermi3 {i in I, j in J, r in R: distance[i,j] > dmax}:
#	prob[i,j,r] = 0;
	

#subject to covering{i in I}:
#	sum{j in I} a[i,j]*x[j] >= y[i];
	
#subject to nequip:
#	sum{j in I} x[j]=equipments;


#solve;
#display _total_solve_time;
#display x;
#display y;
#display z;





