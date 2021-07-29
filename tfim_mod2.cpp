//
// This program calculates matrix elements of the exponential
// of the transverse-field Ising model Hamiltonian modulo two using divided differences
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

#define tolerance 1e-15      // error tolerance
// note that this is absolute error rather than relative error
// note that actual error is determined by the remainder and can be much smaller than the error estimate by the smallest term

double beta=1;
double Gamma=0.05;

unsigned long long zi = 16210525687446977967ULL; // index of canonical basis vector z_i
unsigned long long zf = 16210525687446977967ULL; // index of canonical basis vector z_f

int q;

int k[n]; int differentbits[n]; int Ndifferentbits=0;
int currk[n]; int currseq[qmax+n]; int currE[qmax+n+1]; int currq;

unsigned long long E1occurences[qmax+n+2];
ExExFloat divdiffs[qmax+n+2]; ExExFloat coefficient, Tolerance;

int lattice[n];

int spin(int x, int y){ return lattice[x*L+y];}

void flip(int k){ lattice[k] *= -1;} // perform the spin flip

int IsingEnergy(){ // calculate energy of a given configuration of spins
	int x,y,x1,y1,x2,y2,x3,y3,x4,y4,energy=0;
	for(x=0;x<L;x++) for(y=0;y<L;y++){
		x1=(x+L-1)%L; x3=(x+1)%L; x2=x4=x;
		y4=(y+L-1)%L; y2=(y+1)%L; y1=y3=y;
		energy+=spin(x,y)*(spin(x1,y1)+spin(x2,y2)+spin(x3,y3)+spin(x4,y4));
	}
	return (std::abs(energy)/8)%2;
}

void output(){ // process the path found
	int k=0;
	for(int i=0;i<q+1;i++) if(currE[i]==1) k++;
	E1occurences[k]++;
}

void sequences(){ // for a given k_0,...,k_{n-1} consider all paths such that number of flips of spin #i is k_i.
	int i;
	if(currq<q){ for(i=0;i<n;i++) if(currk[i]>0){
		currk[i]--; currseq[currq]=i; currq++;
		sequences();
		currq--; currk[i]++;
	}} else{ 
		for(i=0;i<q+1;i++){
			currE[i]=IsingEnergy();
			if(i<q) flip(currseq[i]);
		}
		for(i=0;i<q;i++) if(differentbits[currseq[i]]) flip(currseq[i]);
		output();
	}
}

void allsequences(){
	for(int i=0;i<n;i++) currk[i]=k[i]; currq=0; sequences();
}

void processSubset(int walls[]){  // process subsets of size (n-1) of {0,1,2,...,(q-Ndifferentbits)/2+n-2}.
	int curr;	// (k_i+1) / 2 or k_i/2+1 = walls[i] - walls[i-1], where walls[n-1]=(q-Ndifferentbits)/2+n-1 and walls[-1]=-1. Then, k_0+...+k_{n-1}=q.
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

void fill_divdiffs(){
	int i, E1=3;
	divdiff d(qmax+Ndifferentbits+2,500); ExExFloat coeff; // auto initialized as coeff=1
	for(E1=0;E1<=q+1;E1++){
		for(i=0;i<E1;i++) d.z[i]=-beta;
		for(i=E1;i<=q;i++) d.z[i]=0;
		d.AddAll(q+1); divdiffs[E1] = d.divdiffs[q];
		E1occurences[E1] = 0;
	}
	for(i=1;i<=q;i++){ coeff*=(beta*Gamma);	coeff/=i;}
	coefficient = coeff;
}

int main(int argc, char** argv){
	if(argc>=2) zi = strtoull(argv[1],NULL,10);
	if(argc>=3) zf = strtoull(argv[2],NULL,10);
	unsigned long long npaths; int i;
	ExExFloat sum,totalsum,summand; divdiff_init(); Tolerance = tolerance;
	unsigned long long k=1; for(i=0;i<n;i++) { lattice[i]= (zi&k)==0?-1:1; k*=2;}
	for(i=0;i<n;i++) if(((zi>>i)&1) == ((zf>>i)&1)) differentbits[i]=0; else { differentbits[i]=1; Ndifferentbits++;}
	printf("Number of different bits = %d\n",Ndifferentbits);
	for(q=Ndifferentbits;q<=Ndifferentbits+qmax;q+=2){
		fill_divdiffs(); npaths=0;
		allk();
		for(i=0;i<=q+1;i++){
			npaths+=E1occurences[i];
			summand = coefficient*divdiffs[i]*(double)E1occurences[i];
			if(i==0) sum=summand; else sum+=summand;
		}
		if(q==Ndifferentbits) totalsum=sum; else totalsum+=sum;
		printf("q =%3d, number of paths = %15llu, sum = ",q,npaths); sum.print(); printf("\n"); fflush(stdout);
		if(Tolerance>=sum.abs()) break;
	}
	printf("total sum = "); totalsum.print(); printf("\n");
	divdiff_clear_up();
}
