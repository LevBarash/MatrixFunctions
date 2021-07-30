//
// This program calculates matrix elements of the exponential
// of the transverse-field Ising model Hamiltonian using divided differences
//
// This program is introduced in the paper:
// Lev Barash, Stefan Guttel, Itay Hen, Calculating Elements of Matrix Functions using Divided Differences
//
// This program is licensed under a Creative Commons Attribution 4.0 International License:
// http://creativecommons.org/licenses/by/4.0/
//
// ExExFloat datatype and calculation of divided differences are described in the paper:
// L. Gupta, L. Barash, I. Hen, Calculating the divided differences of the exponential function by addition and removal of inputs, Computer Physics Communications 254, 107385 (2020)
//

#include<stdio.h>
#include<cstdlib>
#include"divdiff.h"

#define qmax 8        // maximal addition to the initial expansion order, i.e., Q = Ndifferentbits + qmax
#define L 8           // linear size of lattice
#define n (L*L)       // n = L^2

#define tolerance 1e-15    // error tolerance
// note that this is absolute error rather than relative error
// note that actual error is determined by the remainder and can be much smaller than the error estimate by the smallest term

double beta=1;
double Gamma=0.01;

unsigned long long zi = 16210525687446977967ULL; // index of canonical basis vector z_i
unsigned long long zf = 16210525687446977967ULL; // index of canonical basis vector z_f

int q; divdiff* d;

double ddsum,ddkahanc,coefficient=1; unsigned long long npaths; int currEnergy; int currs=1;

int k[n]; int differentbits[n]; int Ndifferentbits=0;
int currk[n]; int currseq[qmax+n+1]; int currq=0;

int lattice[n];

int spin(int x, int y){ return lattice[x*L+y];}

void flip(int k){ lattice[k] *= -1;} // perform the spin flip

int flipDeltaE(int k){ // perform the spin flip and calculate the corresponding energy change
	int x, y, x1,x2,x3,x4,y1,y2,y3,y4;
	x = k/L; y = k-x*L;
	x1=(x+L-1)%L; x3=(x+1)%L; x2=x4=x;
	y4=(y+L-1)%L; y2=(y+1)%L; y1=y3=y;
	lattice[k] *= -1;
	return 2*lattice[k]*(spin(x1,y1)+spin(x2,y2)+spin(x3,y3)+spin(x4,y4));
}

int IsingEnergy(){ // calculate energy of a given configuration of spins
	int x,y,x1,y1,x2,y2,x3,y3,x4,y4,energy=0;
	for(x=0;x<L;x++) for(y=0;y<L;y++){
		x1=(x+L-1)%L; x3=(x+1)%L; x2=x4=x;
		y4=(y+L-1)%L; y2=(y+1)%L; y1=y3=y;
		energy+=spin(x,y)*(spin(x1,y1)+spin(x2,y2)+spin(x3,y3)+spin(x4,y4));
	}
	return energy/2;
}

void output(){ // process the path found
	double y = d->divdiffs[q].get_double() - ddkahanc; // employing Kahan summation to accurately perform
	volatile double t = ddsum + y;                     // the operation ddsum += d->divdiffs[q].get_double();
	volatile double z = t - ddsum; // volatile keywords are used to inhibit compiler optimizations that are dangerous for Kahan summation
	ddkahanc = z - y;
	ddsum = t;
	npaths++;
}

void sequences(){ // for a given k_0,...,k_{n-1} consider all paths such that number of flips of spin #i is k_i.
	int i, prevEnergy;
	if(currq<q){ for(i=0;i<n;i++) if(currk[i]>0){
		currk[i]--; currseq[currq]=i; currq++;
		prevEnergy=currEnergy; currEnergy+=flipDeltaE(i); d->AddElement(-beta*currEnergy,currs);
		sequences();
		currq--; flip(i); currEnergy=prevEnergy; d->RemoveElement(); currk[i]++;
	}} else	output();
}

void allsequences(){
	for(int i=0;i<n;i++) currk[i]=k[i]; currq=0; sequences();
}

void processSubset(int walls[]){  // process subsets of size (n-1) of {0,1,2,...,(q-Ndifferentbits)/2+n-2}.
        int curr;       // (k_i+1) / 2 or k_i/2+1 = walls[i] - walls[i-1], where walls[n-1]=(q-Ndifferentbits)/2+n-1 and walls[-1]=-1. Then, k_0+...+k_{n-1}=q.
	for(int i=0;i<n;i++){
                curr =  (i<n-1) ? walls[i] : (q-Ndifferentbits)/2+n-1; curr -= i>0 ? walls[i-1] : -1;
                k[i] =  differentbits[i] ? curr*2-1 : (curr-1)*2;
	}
	allsequences();
}

void generateSubsets(int m){	// generate all possible subsets of size n-1 of {0,1,2,...,m-1}
	int i,j,r=n-1,index,data[n-1];
	for(i=0;i<r;i++) data[i] = i;
	while(1){
		processSubset(data); 
		index = r-1; while(index>=0 && data[index]==m-r+index) index--;
		if(index<0) break;
		for(j=r-1;j>=index;j--) data[j]=data[index]+j-index+1;
	}
}

void allk(){	// generate all arrays k_0,k_1,...,k_{n-1} such that k_0+...+k_{n-1}=q, and k_i has the same evenness as differentbits[i]
        generateSubsets((q-Ndifferentbits)/2+n-1);
}

int main(int argc, char** argv){
	if(argc>=2) zi = strtoull(argv[1],NULL,10);
	if(argc>=3) zf = strtoull(argv[2],NULL,10);
	int i; double sum,totalsum; divdiff_init();
        unsigned long long k=1; for(i=0;i<n;i++) { lattice[i]= (zi&k)==0?-1:1; k*=2;}
        for(i=0;i<n;i++) if(((zi>>i)&1) == ((zf>>i)&1)) differentbits[i]=0; else { differentbits[i]=1; Ndifferentbits++;}
        printf("Number of different bits = %d\n",Ndifferentbits);
	divdiff dd(qmax+Ndifferentbits+4,500); d = &dd; currEnergy=IsingEnergy();
        for(q=Ndifferentbits;q<=Ndifferentbits+qmax;q+=2){
	        if(q>0) { coefficient*=(beta*beta*Gamma*Gamma); coefficient/=(q*(q-1)); currs=(int)ceil(4.0*q/3.5);}
		dd.CurrentLength=0; dd.AddElement(-beta*currEnergy,currs);
		npaths=0; ddsum=0.0; ddkahanc=0.0; allk();
		sum = coefficient*ddsum;
		if(q==Ndifferentbits) totalsum=sum; else totalsum+=sum;
		printf("q =%3d, number of paths = %15llu, sum = %.17g\n",q,npaths,sum); fflush(stdout);
		if(tolerance>=fabs(sum)) break;
	}
	printf("total sum = %.17g\n",totalsum);
	divdiff_clear_up();
}
