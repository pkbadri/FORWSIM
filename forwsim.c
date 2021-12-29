#include <math.h>
#include <stdio.h>
#include <fstream>
#include <iostream>
#include <time.h>
#include <ctype.h>
#include <assert.h>
#include "mtrand.h"
using namespace std;

#define     ZERO  0         // Constant. Do not change.
int       SSIZE;            // Number of samples to output from final population. 
int     POPSIZE;            // Total number of chromsomes in the diploid population.
int   SEQLENGTH;            // Length of the DNA sequence. Choose larger than DELETE*mu*POPSIZE.
double       re;            // Per-generation per-sequence rate of recombination.
double       mu;            // Per-generation per-sequence rate of mutation.
double        s;            // Probability of self-fertilization.
int         GEN;            // Number of generations.
int      DELETE;            // Intervals after which fixed mutations are removed.  


unsigned long init[4] = {0x123, 0x234, 0x345, 0x456},length = 4;
unsigned long duplicate[4];
MTRand drand;//double in [0, 1) generator, already init  
MTRand mt;

static int *indicator,*multhit;
static int *n2count,*n2r,*n3count,*n3r,*n4count,*n4r,*n5count,*n5r;
static int *n6count,*n6r,*n7count,*n7r,*n8count,*n8r,*n9count,*n9r;
static int *current, *size, *siz, *temp1;
static int *indic, *temp2, *recount, flag1,flag2;

static inline double  minimum(double i,double j){return(i>j) ? j:i;}
static inline double  maximum(double i,double j){return(i>j) ? i:j;}
static inline int poisson(float mean){int count=0;double q,bound;bound = exp(-mean);for(q=1;q>=bound;q*=drand()){++count;}return(count-1);}

//Ancestry information for a generation
class next{
public:
double *np;
int *n;
};



//Integer arrays that represent chromosomes
class member{
public:
int *q;//Choose a suitable value depending on the expected number of mutations per individual
};static member *temp;



//Pointers to chromosome arrays.
struct node{
member *chr;
};static node *odd,*even;




main(int argc, char *argv[])
{ 
SSIZE = atoi(argv[2]);POPSIZE = atoi(argv[4]);SEQLENGTH = atoi(argv[6]);re = atof(argv[8]);
mu = atof(argv[10]);s = atof(argv[12]);GEN = atoi(argv[14]);DELETE = atoi(argv[16]);
cout<<SSIZE<<" ";cout<<POPSIZE<<" ";cout<<SEQLENGTH<<" ";cout<<re<<" ";cout<<mu<<" ";cout<<s<<" ";cout<<GEN<<" ";cout<<DELETE<<"\n";

if(POPSIZE % 2 ==1){cout<<"Error. pop must be even number.\n"; return(0);}

int i,j,k,l,m,n,o,r,trp,trk,start,end,ind;double p,rec,x,y;

n2count = new int[POPSIZE];n2r = new int[POPSIZE];
n3count = new int[POPSIZE];n3r = new int[POPSIZE];
n4count = new int[POPSIZE];n4r = new int[POPSIZE];
n5count = new int[POPSIZE];n5r = new int[POPSIZE];
n6count = new int[POPSIZE];n6r = new int[POPSIZE];
n7count = new int[POPSIZE];n7r = new int[POPSIZE];
n8count = new int[POPSIZE];n8r = new int[POPSIZE];
n9count = new int[POPSIZE];n9r = new int[POPSIZE];
current = new int[POPSIZE];size = new int[POPSIZE];siz = new int[POPSIZE];
temp1 = new int[POPSIZE];temp2 = new int[POPSIZE];
recount = new int[POPSIZE];indic = new int [POPSIZE];


int finished, array_size, flag;
next *gen1,*gen2,*gen3,*gen4,*gen5,*gen6,*gen7,*gen8,*gen9,*gtemp;
double a,b,c;
gen1 = new next;gen1->np = new double[POPSIZE];gen1->n = new int[POPSIZE];
gen2 = new next;gen2->np = new double[POPSIZE];gen2->n = new int[POPSIZE];
gen3 = new next;gen3->np = new double[POPSIZE];gen3->n = new int[POPSIZE];
gen4 = new next;gen4->np = new double[POPSIZE];gen4->n = new int[POPSIZE];
gen5 = new next;gen5->np = new double[POPSIZE];gen5->n = new int[POPSIZE];
gen6 = new next;gen6->np = new double[POPSIZE];gen6->n = new int[POPSIZE];
gen7 = new next;gen7->np = new double[POPSIZE];gen7->n = new int[POPSIZE];
gen8 = new next;gen8->np = new double[POPSIZE];gen8->n = new int[POPSIZE];
gen9 = new next;gen9->np = new double[POPSIZE];gen9->n = new int[POPSIZE];
gtemp = new next;gtemp->np = new double[POPSIZE];gtemp->n = new int[POPSIZE];
temp = new member;odd = new node[POPSIZE];even = new node[POPSIZE];
indicator = new int[SEQLENGTH];multhit = new int[SEQLENGTH];


rec = 1 - exp(-re);finished=0;
for(i=0;i<length;i++){init[i]=time(NULL);duplicate[i] = init[i];}mt.seed(init,length);

cout<<init[0]<<"\t"<<init[1]<<"\t"<<init[2]<<"\t"<<init[3]<<"\n";

array_size = 2500 + int(0.25*POPSIZE*mu);

while (finished == 0){

memset(&size[0],'\0',4*POPSIZE);memset(&siz[0],'\0',4*POPSIZE);memset(&multhit[0],'\0',4*SEQLENGTH);
a=0;b=0;c=0;flag=0;array_size += 100;
for(i=0;i<POPSIZE;i++){odd[i].chr = new member;odd[i].chr->q = new int[array_size];even[i].chr = new member;even[i].chr->q = new int[array_size];}
 for(i=0;i<length;i++){init[i]=duplicate[i];}mt.seed(init,length);

m=0; 
while(m<POPSIZE){
gen1->np[m] = drand();gen1->np[m+1] = drand();x=drand();
if(x>=s){gen1->n[m]=int(POPSIZE*drand());gen1->n[m+1]=int(POPSIZE*drand());}
else{y=drand();if(y<=0.50){gen1->n[m]=int(POPSIZE*drand());gen1->n[m+1]= gen1->n[m];}
else{gen1->n[m]=int(POPSIZE*drand());if((gen1->n[m])%2==0){gen1->n[m+1]=(gen1->n[m])+1;}if((gen1->n[m])%2==1){gen1->n[m+1]=(gen1->n[m])-1;}}
}

gen2->np[m] = drand();gen2->np[m+1] = drand();x=drand();
if(x>=s){gen2->n[m]=int(POPSIZE*drand());gen2->n[m+1]=int(POPSIZE*drand());}
else{y=drand();if(y<=0.50){gen2->n[m]=int(POPSIZE*drand());gen2->n[m+1]= gen2->n[m];}
else{gen2->n[m]=int(POPSIZE*drand());if((gen2->n[m])%2==0){gen2->n[m+1]=(gen2->n[m])+1;}if((gen2->n[m])%2==1){gen2->n[m+1]=(gen2->n[m])-1;}}
}

gen3->np[m] = drand();gen3->np[m+1] = drand();x=drand();
if(x>=s){gen3->n[m]=int(POPSIZE*drand());gen3->n[m+1]=int(POPSIZE*drand());}
else{y=drand();if(y<=0.50){gen3->n[m]=int(POPSIZE*drand());gen3->n[m+1]= gen3->n[m];}
else{gen3->n[m]=int(POPSIZE*drand());if((gen3->n[m])%2==0){gen3->n[m+1]=(gen3->n[m])+1;}if((gen3->n[m])%2==1){gen3->n[m+1]=(gen3->n[m])-1;}}
}

gen4->np[m] = drand();gen4->np[m+1] = drand();x=drand();
if(x>=s){gen4->n[m]=int(POPSIZE*drand());gen4->n[m+1]=int(POPSIZE*drand());}
else{y=drand();if(y<=0.50){gen4->n[m]=int(POPSIZE*drand());gen4->n[m+1]= gen4->n[m];}
else{gen4->n[m]=int(POPSIZE*drand());if((gen4->n[m])%2==0){gen4->n[m+1]=(gen4->n[m])+1;}if((gen4->n[m])%2==1){gen4->n[m+1]=(gen4->n[m])-1;}}
}

gen5->np[m] = drand();gen5->np[m+1] = drand();x=drand();
if(x>=s){gen5->n[m]=int(POPSIZE*drand());gen5->n[m+1]=int(POPSIZE*drand());}
else{y=drand();if(y<=0.50){gen5->n[m]=int(POPSIZE*drand());gen5->n[m+1]= gen5->n[m];}
else{gen5->n[m]=int(POPSIZE*drand());if((gen5->n[m])%2==0){gen5->n[m+1]=(gen5->n[m])+1;}if((gen5->n[m])%2==1){gen5->n[m+1]=(gen5->n[m])-1;}}
}

gen6->np[m] = drand();gen6->np[m+1] = drand();x=drand();
if(x>=s){gen6->n[m]=int(POPSIZE*drand());gen6->n[m+1]=int(POPSIZE*drand());}
else{y=drand();if(y<=0.50){gen6->n[m]=int(POPSIZE*drand());gen6->n[m+1]= gen6->n[m];}
else{gen6->n[m]=int(POPSIZE*drand());if((gen6->n[m])%2==0){gen6->n[m+1]=(gen6->n[m])+1;}if((gen6->n[m])%2==1){gen6->n[m+1]=(gen6->n[m])-1;}}
}

gen7->np[m] = drand();gen7->np[m+1] = drand();x=drand();
if(x>=s){gen7->n[m]=int(POPSIZE*drand());gen7->n[m+1]=int(POPSIZE*drand());}
else{y=drand();if(y<=0.50){gen7->n[m]=int(POPSIZE*drand());gen7->n[m+1]= gen7->n[m];}
else{gen7->n[m]=int(POPSIZE*drand());if((gen7->n[m])%2==0){gen7->n[m+1]=(gen7->n[m])+1;}if((gen7->n[m])%2==1){gen7->n[m+1]=(gen7->n[m])-1;}}
}

gen8->np[m] = drand();gen8->np[m+1] = drand();x=drand();
if(x>=s){gen8->n[m]=int(POPSIZE*drand());gen8->n[m+1]=int(POPSIZE*drand());}
else{y=drand();if(y<=0.50){gen8->n[m]=int(POPSIZE*drand());gen8->n[m+1]= gen8->n[m];}
else{gen8->n[m]=int(POPSIZE*drand());if((gen8->n[m])%2==0){gen8->n[m+1]=(gen8->n[m])+1;}if((gen8->n[m])%2==1){gen8->n[m+1]=(gen8->n[m])-1;}}
}

gen9->np[m] = drand();gen9->np[m+1] = drand();x=drand();
if(x>=s){gen9->n[m]=int(POPSIZE*drand());gen9->n[m+1]=int(POPSIZE*drand());}
else{y=drand();if(y<=0.50){gen9->n[m]=int(POPSIZE*drand());gen9->n[m+1]= gen9->n[m];}
else{gen9->n[m]=int(POPSIZE*drand());if((gen9->n[m])%2==0){gen9->n[m+1]=(gen9->n[m])+1;}if((gen9->n[m])%2==1){gen9->n[m+1]=(gen9->n[m])-1;}}
}
m=m+2;}

memset(&n2count[0],'\0',4*POPSIZE);memset(&n2r[0],'\0',4*POPSIZE);
memset(&n3count[0],'\0',4*POPSIZE);memset(&n3r[0],'\0',4*POPSIZE);
memset(&n4count[0],'\0',4*POPSIZE);memset(&n4r[0],'\0',4*POPSIZE);
memset(&n5count[0],'\0',4*POPSIZE);memset(&n5r[0],'\0',4*POPSIZE);
memset(&n6count[0],'\0',4*POPSIZE);memset(&n6r[0],'\0',4*POPSIZE);
memset(&n7count[0],'\0',4*POPSIZE);memset(&n7r[0],'\0',4*POPSIZE);
memset(&n8count[0],'\0',4*POPSIZE);memset(&n8r[0],'\0',4*POPSIZE);
memset(&n9count[0],'\0',4*POPSIZE);memset(&n9r[0],'\0',4*POPSIZE);



for(i=1;i<=GEN;i++){
//Create odd generation
if(i%2 == 1){ 
memset(&size[0],'\0',4*POPSIZE);memset(&indic[0],'\0',4*POPSIZE);memset(&recount[0],'\0',4*POPSIZE);
for(m=0;m<POPSIZE;m++){
temp1[m] = gen2->n[m];n2count[temp1[m]] = 1; 
if(gen2->np[m]<rec){
if(gen2->n[m]%2==0){n2r[gen2->n[m]+1] = 1;}
if(gen2->n[m]%2==1){n2r[gen2->n[m]-1] = 1;}
}current[m] = m;}
  
for(m=0;m<POPSIZE;m++){
temp2[m] = temp1[gen3->n[m]];n3count[temp2[m]] = 1; 
if(gen3->np[m]<rec){
if(gen3->n[m]%2==0){n3r[temp1[gen3->n[m]+1]] = 1;}
if(gen3->n[m]%2==1){n3r[temp1[gen3->n[m]-1]] = 1;}
}} 

for(m=0;m<POPSIZE;m++){
temp1[m] = temp2[gen4->n[m]];n4count[temp1[m]] = 1; 
if(gen4->np[m]<rec){
if(gen4->n[m]%2==0){n4r[temp2[gen4->n[m]+1]] = 1;}
if(gen4->n[m]%2==1){n4r[temp2[gen4->n[m]-1]] = 1;}
}} 

for(m=0;m<POPSIZE;m++){
temp2[m] = temp1[gen5->n[m]];n5count[temp2[m]] = 1; 
if(gen5->np[m]<rec){
if(gen5->n[m]%2==0){n5r[temp1[gen5->n[m]+1]] = 1;}
if(gen5->n[m]%2==1){n5r[temp1[gen5->n[m]-1]] = 1;}
}} 

for(m=0;m<POPSIZE;m++){
temp1[m] = temp2[gen6->n[m]];n6count[temp1[m]] = 1;
if(gen6->np[m]<rec){
if(gen6->n[m]%2==0){n6r[temp2[gen6->n[m]+1]] = 1;}
if(gen6->n[m]%2==1){n6r[temp2[gen6->n[m]-1]] = 1;}
}}

for(m=0;m<POPSIZE;m++){
temp2[m] = temp1[gen7->n[m]];n7count[temp2[m]] = 1;
if(gen7->np[m]<rec){
if(gen7->n[m]%2==0){n7r[temp1[gen7->n[m]+1]] = 1;}
if(gen7->n[m]%2==1){n7r[temp1[gen7->n[m]-1]] = 1;}
}}

for(m=0;m<POPSIZE;m++){
temp1[m] = temp2[gen8->n[m]];n8count[temp1[m]] = 1;
if(gen8->np[m]<rec){
if(gen8->n[m]%2==0){n8r[temp2[gen8->n[m]+1]] = 1;}
if(gen8->n[m]%2==1){n8r[temp2[gen8->n[m]-1]] = 1;}
}}

for(m=0;m<POPSIZE;m++){
temp2[m] = temp1[gen9->n[m]];n9count[temp2[m]] = 1;
if(gen9->np[m]<rec){
if(gen9->n[m]%2==0){n9r[temp1[gen9->n[m]+1]] = 1;}
if(gen9->n[m]%2==1){n9r[temp1[gen9->n[m]-1]] = 1;}
}}


m=0;
while(m<POPSIZE){
flag1=0;
while(ZERO==0){
if(n2count[m]==0){flag1=1;break;}
if(n3r[m]==0){
if(n3count[m]==0){flag1=1;break;}
if(n4r[m]==0){
if(n4count[m]==0){flag1=1;break;}
if(n5r[m]==0){
if(n5count[m]==0){flag1=1;break;}
if(n6r[m]==0){
if(n6count[m]==0){flag1=1;break;}
if(n7r[m]==0){
if(n7count[m]==0){flag1=1;break;}
if(n8r[m]==0){
if(n8count[m]==0){flag1=1;break;}
if(n9r[m]==0){
if(n9count[m]==0){flag1=1;break;}
}}}}}}}
break;}

flag2=0;
while(ZERO==0){
if(n2count[m+1]==0){flag2=1;break;}
if(n3r[m+1]==0){
if(n3count[m+1]==0){flag2=1;break;}
if(n4r[m+1]==0){
if(n4count[m+1]==0){flag2=1;break;}
if(n5r[m+1]==0){
if(n5count[m+1]==0){flag2=1;break;}
if(n6r[m+1]==0){
if(n6count[m+1]==0){flag2=1;break;}
if(n7r[m+1]==0){
if(n7count[m+1]==0){flag2=1;break;}
if(n8r[m+1]==0){
if(n8count[m+1]==0){flag2=1;break;}
if(n9r[m+1]==0){
if(n9count[m+1]==0){flag2=1;break;}
}}}}}}}
break;}

if(flag1==1 && flag2==1){indic[m]=1;indic[m+1]=1;}
else{
if(flag1==1 && n2r[m]==0){indic[m]=1;}
if(flag2==1 && n2r[m+1]==0){indic[m+1]=1;}
}

if(indic[m]==0 || (i%DELETE>=DELETE-10 || i%DELETE==0 || i>=GEN-10)){
p = gen1->np[m];j = gen1->n[m];
if(p < rec){if(j%2==0){k = j+1;}if(j%2==1){k = j-1;}++recount[(j+k)/2];}
else{
if(current[j]!=j){memcpy(&odd[m].chr->q[0],&odd[current[j]-POPSIZE].chr->q[0],4*siz[j]);size[m]=siz[j];++c;}
if(current[j]==j){temp = odd[m].chr;odd[m].chr = even[j].chr;even[j].chr = temp;current[j]=m + POPSIZE;size[m] = siz[j];}
}}
  
if(indic[m+1]==0 || (i%DELETE>=DELETE-10 || i%DELETE==0 || i>=GEN-10)){
p = gen1->np[m+1];j = gen1->n[m+1];
if(p < rec){if(j%2==0){k = j+1;}if(j%2==1){k = j-1;}++recount[(j+k)/2];}
else{
if(current[j]!=j){memcpy(&odd[m+1].chr->q[0],&odd[current[j]-POPSIZE].chr->q[0],4*siz[j]);size[m+1]=siz[j];++c;}
if(current[j]==j){temp = odd[m+1].chr;odd[m+1].chr = even[j].chr;even[j].chr = temp;current[j]=m+1+POPSIZE;size[m+1] = siz[j];}
}}
m=m+2;}


for(m=0;m<POPSIZE;m++){
p = gen1->np[m];j = gen1->n[m];
if(i%DELETE > 0  && i%DELETE < DELETE-10 && i < GEN-10){if(indic[m]==1){continue;}}

if(p<rec){
r = int((SEQLENGTH-1)*drand());if(j%2==0){k = j+1;}if(j%2==1){k = j-1;}--recount[(j+k)/2];
if(recount[(j+k)/2]==0){if(current[j]!=j && current[k]==k){o=j;j=k;k=o;}}
if(current[j]==j && siz[j] > 0){
if(r>=even[j].chr->q[siz[j]-1]){
if(recount[(j+k)/2]==0){temp=odd[m].chr;odd[m].chr=even[j].chr;even[j].chr=temp;size[m]=siz[j];}
else{memcpy(&odd[m].chr->q[0],&even[j].chr->q[0],4*siz[j]);size[m]=siz[j];a +=1.0;}
}
  
else{
start=0;end=siz[j]-1;trk=end;ind=0;
while(trk > 0 && trk <=siz[j]-1){
if(r >=even[j].chr->q[trk-1] && r <even[j].chr->q[trk]){ind=1;break;}
if(r >=even[j].chr->q[trk]){start=trk;trk=(start+end)/2;if(trk==start){++trk;}continue;}
if(r <even[j].chr->q[trk]){end=trk;trk=(start+end)/2;}}
if(ind==1){
if(recount[(j+k)/2]==0){temp=odd[m].chr;odd[m].chr=even[j].chr;even[j].chr=temp;size[m]=trk;}
else{memcpy(&odd[m].chr->q[0],&even[j].chr->q[0],4*trk);size[m]=trk;a += 0.5;}
}}
}
 
if(current[j]!=j && siz[j] > 0){
trp=current[j]-POPSIZE;  
if(r>=odd[trp].chr->q[siz[j]-1]){memcpy(&odd[m].chr->q[0],&odd[trp].chr->q[0],4*siz[j]);size[m]=siz[j]; a += 1.0;}
else{
start=0;end=siz[j]-1;trk=end;ind=0;
while(trk > 0 && trk<=siz[j]-1){
if(r >=odd[trp].chr->q[trk-1] && r <odd[trp].chr->q[trk]){ind=1;break;}
if(r >=odd[trp].chr->q[trk]){start=trk;trk=(start+end)/2;if(trk==start){++trk;}continue;}
if(r <odd[trp].chr->q[trk]){end=trk;trk=(start+end)/2;}}
if(ind==1){memcpy(&odd[m].chr->q[0],&odd[trp].chr->q[0],4*trk);size[m]=trk;a +=0.5;}}}

if(current[k]==k && siz[k]>0){
if(r<even[k].chr->q[0]){memcpy(&odd[m].chr->q[size[m]],&even[k].chr->q[0],4*siz[k]);size[m]=size[m]+siz[k];a += 1.0;}
else{
start=0;end=siz[k]-1;trk=end;ind=0;
while(trk > 0 && trk <=siz[k]-1){
if(r >=even[k].chr->q[trk-1] && r <even[k].chr->q[trk]){ind=1;break;}
if(r >=even[k].chr->q[trk]){start=trk;trk=(start+end)/2;if(trk==start){++trk;}continue;}
if(r <even[k].chr->q[trk]){end=trk;trk=(start+end)/2;}}
if(ind==1){memcpy(&odd[m].chr->q[size[m]],&even[k].chr->q[trk],4*(siz[k]-trk));size[m]=size[m]+siz[k]-trk;a += 0.5;}}}

if(current[k]!=k && siz[k]>0){trp=current[k]-POPSIZE;
if(r<odd[trp].chr->q[0]){memcpy(&odd[m].chr->q[size[m]],&odd[trp].chr->q[0],4*siz[k]);size[m]=size[m]+siz[k];a += 1.0;}
else{
start=0;end=siz[k]-1;trk=end;ind=0;
while(trk > 0 && trk<=siz[k]-1){
if(r >=odd[trp].chr->q[trk-1] && r <odd[trp].chr->q[trk]){ind=1;break;}
if(r >=odd[trp].chr->q[trk]){start=trk;trk=(start+end)/2;if(trk==start){++trk;}continue;}
if(r <odd[trp].chr->q[trk]){end=trk;trk=(start+end)/2;}}
if(ind==1){memcpy(&odd[m].chr->q[size[m]],&odd[trp].chr->q[trk],4*(siz[k]-trk));size[m]=size[m]+siz[k]-trk;a += 0.5;}}}
}}  

for(m=0;m<POPSIZE;m++){ 

if(i%DELETE < DELETE-10 && i%DELETE > 0 && i < GEN-10){if(indic[m]==1){continue;}}

//INSERT NEW MUTATIONS
k = poisson(mu);
if(k + size[m] > array_size){cout<<"Out of Memory. Increasing array size and rerunning.\n";break;}

for(l=0;l<k;l++){
trp = int(SEQLENGTH*drand());while(multhit[trp]==1){trp = int(SEQLENGTH*drand());}multhit[trp]=1;
if((size[m] > 0 && trp>=odd[m].chr->q[size[m]-1]) || size[m]==0){odd[m].chr->q[size[m]]=trp;}
else{
if(trp<=odd[m].chr->q[0]){memmove(&odd[m].chr->q[1],&odd[m].chr->q[0],4*size[m]);odd[m].chr->q[0]=trp;++b;}
else{
start=0;end=size[m]-1;trk=end;ind=0;
while(trk > 0 && trk<=size[m]-1){
if(trp > odd[m].chr->q[trk-1] && trp <=odd[m].chr->q[trk]){ind=1;break;} 
if(trp > odd[m].chr->q[trk]){start=trk;trk=(start+end)/2;if(trk==start){++trk;}continue;} 
if(trp <= odd[m].chr->q[trk]){end=trk;trk=(start+end)/2;continue;}} 
if(ind==1){memmove(&odd[m].chr->q[trk+1],&odd[m].chr->q[trk],4*(size[m]-trk));odd[m].chr->q[trk]=trp;b +=0.5;}
if(ind==0){cout<<"Error\n";return(0);}}
}++size[m];}
}

if(m < POPSIZE){for(m=0;m<POPSIZE;m++){delete odd[m].chr->q;delete odd[m].chr;delete even[m].chr->q;delete even[m].chr;}flag=1;break;}

 
gtemp = gen1;gen1 = gen2;gen2 = gen3;gen3 = gen4;gen4 = gen5;gen5 = gen6;gen6 = gen7;gen7 = gen8;gen8=gen9;gen9 = gtemp;
m=0;
while(m<POPSIZE){
gen9->np[m] = drand();gen9->np[m+1] = drand();x=drand();
if(x>=s){gen9->n[m]=int(POPSIZE*drand());gen9->n[m+1]=int(POPSIZE*drand());}
else{y=drand();if(y<=0.50){gen9->n[m]=int(POPSIZE*drand());gen9->n[m+1]= gen9->n[m];}
else{gen9->n[m]=int(POPSIZE*drand());if((gen9->n[m])%2==0){gen9->n[m+1]=(gen9->n[m])+1;}
else{gen9->n[m+1]=(gen9->n[m])-1;}}
}
m=m+2;}

memset(&n2count[0],'\0',4*POPSIZE);memset(&n2r[0],'\0',4*POPSIZE);
memset(&n3count[0],'\0',4*POPSIZE);memset(&n3r[0],'\0',4*POPSIZE);
memset(&n4count[0],'\0',4*POPSIZE);memset(&n4r[0],'\0',4*POPSIZE);
memset(&n5count[0],'\0',4*POPSIZE);memset(&n5r[0],'\0',4*POPSIZE);
memset(&n6count[0],'\0',4*POPSIZE);memset(&n6r[0],'\0',4*POPSIZE);
memset(&n7count[0],'\0',4*POPSIZE);memset(&n7r[0],'\0',4*POPSIZE);
memset(&n8count[0],'\0',4*POPSIZE);memset(&n8r[0],'\0',4*POPSIZE);
memset(&n9count[0],'\0',4*POPSIZE);memset(&n9r[0],'\0',4*POPSIZE);
}

  
if(i%2 == 0){ 
memset(&siz[0],'\0',4*POPSIZE);memset(&indic[0],'\0',4*POPSIZE);memset(&recount[0],'\0',4*POPSIZE);
for(m=0;m<POPSIZE;m++){
temp1[m] = gen2->n[m];n2count[temp1[m]] = 1; 
if(gen2->np[m]<rec){
if(gen2->n[m]%2==0){n2r[gen2->n[m]+1] = 1;}
if(gen2->n[m]%2==1){n2r[gen2->n[m]-1] = 1;}
}current[m] = m;} 

for(m=0;m<POPSIZE;m++){
temp2[m] = temp1[gen3->n[m]];n3count[temp2[m]] = 1; 
if(gen3->np[m]<rec){
if(gen3->n[m]%2==0){n3r[temp1[gen3->n[m]+1]] = 1;}
if(gen3->n[m]%2==1){n3r[temp1[gen3->n[m]-1]] = 1;}
}} 
  
for(m=0;m<POPSIZE;m++){
temp1[m] = temp2[gen4->n[m]];n4count[temp1[m]] = 1; 
if(gen4->np[m]<rec){
if(gen4->n[m]%2==0){n4r[temp2[gen4->n[m]+1]] = 1;}
if(gen4->n[m]%2==1){n4r[temp2[gen4->n[m]-1]] = 1;}
}} 

for(m=0;m<POPSIZE;m++){
temp2[m] = temp1[gen5->n[m]];n5count[temp2[m]] = 1; 
if(gen5->np[m]<rec){
if(gen5->n[m]%2==0){n5r[temp1[gen5->n[m]+1]] = 1;}
if(gen5->n[m]%2==1){n5r[temp1[gen5->n[m]-1]] = 1;}
}} 

for(m=0;m<POPSIZE;m++){
temp1[m] = temp2[gen6->n[m]];n6count[temp1[m]] = 1;
if(gen6->np[m]<rec){
if(gen6->n[m]%2==0){n6r[temp2[gen6->n[m]+1]] = 1;}
if(gen6->n[m]%2==1){n6r[temp2[gen6->n[m]-1]] = 1;}
}}

for(m=0;m<POPSIZE;m++){
temp2[m] = temp1[gen7->n[m]];n7count[temp2[m]] = 1;
if(gen7->np[m]<rec){
if(gen7->n[m]%2==0){n7r[temp1[gen7->n[m]+1]] = 1;}
if(gen7->n[m]%2==1){n7r[temp1[gen7->n[m]-1]] = 1;}
}}

for(m=0;m<POPSIZE;m++){
temp1[m] = temp2[gen8->n[m]];n8count[temp1[m]] = 1;
if(gen8->np[m]<rec){
if(gen8->n[m]%2==0){n8r[temp2[gen8->n[m]+1]] = 1;}
if(gen8->n[m]%2==1){n8r[temp2[gen8->n[m]-1]] = 1;}
}}

for(m=0;m<POPSIZE;m++){
temp2[m] = temp1[gen9->n[m]];n9count[temp2[m]] = 1;
if(gen9->np[m]<rec){
if(gen9->n[m]%2==0){n9r[temp1[gen9->n[m]+1]] = 1;}
if(gen9->n[m]%2==1){n9r[temp1[gen9->n[m]-1]] = 1;}
}}


m=0;
while(m<POPSIZE){
flag1=0;
while(ZERO==0){
if(n2count[m]==0){flag1=1;break;}
if(n3r[m]==0){
if(n3count[m]==0){flag1=1;break;}
if(n4r[m]==0){
if(n4count[m]==0){flag1=1;break;}
if(n5r[m]==0){
if(n5count[m]==0){flag1=1;break;}
if(n6r[m]==0){
if(n6count[m]==0){flag1=1;break;}
if(n7r[m]==0){
if(n7count[m]==0){flag1=1;break;}
if(n8r[m]==0){
if(n8count[m]==0){flag1=1;break;}
if(n9r[m]==0){
if(n9count[m]==0){flag1=1;break;}
}}}}}}}
break;}


flag2=0;
while(ZERO==0){
if(n2count[m+1]==0){flag2=1;break;}
if(n3r[m+1]==0){
if(n3count[m+1]==0){flag2=1;break;}
if(n4r[m+1]==0){
if(n4count[m+1]==0){flag2=1;break;}
if(n5r[m+1]==0){
if(n5count[m+1]==0){flag2=1;break;}
if(n6r[m+1]==0){
if(n6count[m+1]==0){flag2=1;break;}
if(n7r[m+1]==0){
if(n7count[m+1]==0){flag2=1;break;}
if(n8r[m+1]==0){
if(n8count[m+1]==0){flag2=1;break;}
if(n9r[m+1]==0){
if(n9count[m+1]==0){flag2=1;break;}
}}}}}}}
break;}

if(flag1==1 && flag2==1){indic[m]=1;indic[m+1]=1;}
else{
if(flag1==1 && n2r[m]==0){indic[m]=1;}
if(flag2==1 && n2r[m+1]==0){indic[m+1]=1;}}
  

if(indic[m]==0 || (i%DELETE>=DELETE-10 || i%DELETE==0 || i>=GEN-10)){
p = gen1->np[m];j = gen1->n[m];if(p < rec){if(j%2==0){k = j+1;}if(j%2==1){k = j-1;}++recount[(j+k)/2];}
else{
if(current[j]!=j){memcpy(&even[m].chr->q[0],&even[current[j]-POPSIZE].chr->q[0],4*size[j]);siz[m]=size[j];}
if(current[j]==j){temp = even[m].chr;even[m].chr = odd[j].chr;odd[j].chr = temp;current[j]=m+POPSIZE;siz[m] = size[j];}
}}

if(indic[m+1]==0 || (i%DELETE>=DELETE-10 || i%DELETE==0 || i>=GEN-10)){
p = gen1->np[m+1];j = gen1->n[m+1];if(p < rec){if(j%2==0){k = j+1;}if(j%2==1){k = j-1;}++recount[(j+k)/2];}
else{
if(current[j]!=j){memcpy(&even[m+1].chr->q[0],&even[current[j]-POPSIZE].chr->q[0],4*size[j]);siz[m+1]=size[j];}
if(current[j]==j){temp = even[m+1].chr;even[m+1].chr = odd[j].chr;odd[j].chr = temp;current[j]=m+1+POPSIZE;siz[m+1] = size[j];}
}}
m=m+2;}


for(m=0;m<POPSIZE;m++){
p = gen1->np[m];j = gen1->n[m];
if(i%DELETE < DELETE-10 && i%DELETE > 0 && i < GEN-10){if(indic[m]==1){continue;}}
  
if(p<rec){
r = int((SEQLENGTH-1)*drand());if(j%2==0){k = j+1;}if(j%2==1){k = j-1;}--recount[(j+k)/2];  
if(recount[(j+k)/2]==0){if(current[j]!=j && current[k]==k){o=j;j=k;k=o;}}
if(current[j]==j && size[j] > 0){
if(r>=odd[j].chr->q[size[j]-1]){
if(recount[(j+k)/2]==0){temp = even[m].chr;even[m].chr = odd[j].chr;odd[j].chr = temp;current[j] = m+POPSIZE;siz[m] = size[j];}
else{memcpy(&even[m].chr->q[0],&odd[j].chr->q[0],4*size[j]);siz[m]=size[j];}
}
  
else{
start=0;end=size[j]-1;trk=end;ind=0;
while(trk > 0 && trk <=size[j]-1){
if(r >=odd[j].chr->q[trk-1] && r <odd[j].chr->q[trk]){ind=1;break;}
if(r >=odd[j].chr->q[trk]){start=trk;trk=(start+end)/2;if(trk==start){++trk;}continue;}
if(r <odd[j].chr->q[trk]){end=trk;trk=(start+end)/2;}
}

if(ind==1){
if(recount[(j+k)/2]==0){temp = even[m].chr;even[m].chr = odd[j].chr;odd[j].chr = temp;current[j]=m+POPSIZE;siz[m] = trk;}  
else{memcpy(&even[m].chr->q[0],&odd[j].chr->q[0],4*trk);siz[m]=trk;}
}}
}

if(current[j]!=j && size[j]>0){
trp=current[j]-POPSIZE;
if(r>=even[trp].chr->q[size[j]-1]){memcpy(&even[m].chr->q[0],&even[trp].chr->q[0],4*size[j]);siz[m]=size[j];}
else{
start=0;end=size[j]-1;trk=end;ind=0;
while(trk > 0 && trk<=size[j]-1){
if(r >=even[trp].chr->q[trk-1] && r <even[trp].chr->q[trk]){ind=1;break;}
if(r >=even[trp].chr->q[trk]){start=trk;trk=(start+end)/2;if(trk==start){++trk;}continue;}
if(r <even[trp].chr->q[trk]){end=trk;trk=(start+end)/2;}}
if(ind==1){memcpy(&even[m].chr->q[0],&even[trp].chr->q[0],4*trk);siz[m]=trk;}}}

if(current[k]==k && size[k]>0){
if(r<odd[k].chr->q[0]){memcpy(&even[m].chr->q[siz[m]],&odd[k].chr->q[0],4*size[k]);siz[m]=siz[m]+size[k];}
else{
start=0;end=size[k]-1;trk=end;ind=0;
while(trk > 0 && trk <=size[k]-1){
if(r >=odd[k].chr->q[trk-1] && r <odd[k].chr->q[trk]){ind=1;break;}
if(r >=odd[k].chr->q[trk]){start=trk;trk=(start+end)/2;if(trk==start){++trk;}continue;}
if(r <odd[k].chr->q[trk]){end=trk;trk=(start+end)/2;}}
if(ind==1){memcpy(&even[m].chr->q[siz[m]],&odd[k].chr->q[trk],4*(size[k]-trk));siz[m]=siz[m]+size[k]-trk;}}}

if(current[k]!=k && size[k]>0){trp=current[k]-POPSIZE;
if(r<even[trp].chr->q[0]){memcpy(&even[m].chr->q[siz[m]],&even[trp].chr->q[0],4*size[k]);siz[m]=siz[m]+size[k];}
else{
start=0;end=size[k]-1;trk=end;ind=0;
while(trk > 0 && trk<=size[k]-1){
if(r >=even[trp].chr->q[trk-1] && r <even[trp].chr->q[trk]){ind=1;break;}
if(r >=even[trp].chr->q[trk]){start=trk;trk=(start+end)/2;if(trk==start){++trk;}continue;}
if(r <even[trp].chr->q[trk]){end=trk;trk=(start+end)/2;}}
if(ind==1){memcpy(&even[m].chr->q[siz[m]],&even[trp].chr->q[trk],4*(size[k]-trk));siz[m]=siz[m]+size[k]-trk;}}}
}}
  

for(m=0;m<POPSIZE;m++){

if(i%DELETE < DELETE-10 && i%DELETE > 0 && i < GEN-10){if(indic[m]==1){continue;}}

//INSERT NEW MUTATIONS
k = poisson(mu);
if(k + siz[m] > array_size){cout<<"Out of Memory. Increasing array size and rerunning.\n";break;}

for(l=0;l<k;l++){
trp = int(SEQLENGTH*drand());while(multhit[trp]==1){trp = int(SEQLENGTH*drand());}multhit[trp]=1;
if((siz[m]>0 && trp>=even[m].chr->q[siz[m]-1]) || siz[m]==0){even[m].chr->q[siz[m]]=trp;}
else{
if(trp<=even[m].chr->q[0]){memmove(&even[m].chr->q[1],&even[m].chr->q[0],4*siz[m]);even[m].chr->q[0]=trp;}
else{
start=0;end=siz[m]-1;trk=end;ind=0;
while(trk > 0 && trk<=siz[m]-1){
if(trp > even[m].chr->q[trk-1] && trp <=even[m].chr->q[trk]){ind=1;break;}
if(trp > even[m].chr->q[trk]){start=trk;trk=(start+end)/2;if(trk==start){++trk;}continue;}
if(trp <=even[m].chr->q[trk]){end=trk;trk=(start+end)/2;continue;}
}
if(ind==1){memmove(&even[m].chr->q[trk+1],&even[m].chr->q[trk],4*(siz[m]-trk));even[m].chr->q[trk]=trp;}
if(ind==0){cout<<"Error\n";return(0);}}
}++siz[m];}
}

if(m < POPSIZE){for(m=0;m<POPSIZE;m++){delete odd[m].chr->q;delete odd[m].chr;delete even[m].chr->q;delete even[m].chr;}flag=1;break;}

gtemp = gen1;gen1 = gen2;gen2 = gen3;gen3 = gen4;gen4 = gen5;gen5 = gen6;gen6 = gen7;gen7 = gen8;gen8=gen9;gen9 = gtemp;
m=0;
while(m<POPSIZE){
gen9->np[m] = drand();gen9->np[m+1] = drand();x=drand();
if(x>=s){gen9->n[m]=int(POPSIZE*drand());gen9->n[m+1]=int(POPSIZE*drand());}
else{y=drand();if(y<=0.50){gen9->n[m]=int(POPSIZE*drand());gen9->n[m+1]= gen9->n[m];}
else{gen9->n[m]=int(POPSIZE*drand());if((gen9->n[m])%2==0){gen9->n[m+1]=(gen9->n[m])+1;}
else{gen9->n[m+1]=(gen9->n[m])-1;}}
}
m=m+2;}
 
memset(&n2count[0],'\0',4*POPSIZE);memset(&n2r[0],'\0',4*POPSIZE);
memset(&n3count[0],'\0',4*POPSIZE);memset(&n3r[0],'\0',4*POPSIZE);
memset(&n4count[0],'\0',4*POPSIZE);memset(&n4r[0],'\0',4*POPSIZE);
memset(&n5count[0],'\0',4*POPSIZE);memset(&n5r[0],'\0',4*POPSIZE);
memset(&n6count[0],'\0',4*POPSIZE);memset(&n6r[0],'\0',4*POPSIZE);
memset(&n7count[0],'\0',4*POPSIZE);memset(&n7r[0],'\0',4*POPSIZE);
memset(&n8count[0],'\0',4*POPSIZE);memset(&n8r[0],'\0',4*POPSIZE);
memset(&n9count[0],'\0',4*POPSIZE);memset(&n9r[0],'\0',4*POPSIZE);
}

//REMOVE FIXED MUTATIONS
if(i%DELETE==0 || i==GEN){

if((i&1)==0){
memset(&indicator[0],'\0',4*SEQLENGTH);
for(m=0;m<POPSIZE;m++){
for(l=0;l<siz[m];l++){++indicator[even[m].chr->q[l]];}
}

for(l=0;l<SEQLENGTH;l++){
if(indicator[l]==0 || indicator[l]==POPSIZE){multhit[l]=0;}
if(indicator[l] > POPSIZE || indicator[l]<0){cout<<i<<" "<<"Error2222\n";return(0);}
}

int count=0;
for(m=0;m<POPSIZE;m++){
k=0;trp=siz[m];l=0;
while(l< trp){
trk = even[m].chr->q[l];
if(indicator[trk]==POPSIZE){++l;continue;}
else{even[m].chr->q[k] = trk;++k;++l;continue;}
}siz[m]=k;count=count+siz[m];}
cout<<"Generation "<<i<<" Total number of mutations: "<<count<<"\n";
}

else{
memset(&indicator[0],'\0',4*SEQLENGTH);
for(m=0;m<POPSIZE;m++){
for(l=0;l<size[m];l++){++indicator[odd[m].chr->q[l]];}
}


for(l=0;l<SEQLENGTH;l++){
if(indicator[l]==0 || indicator[l]==POPSIZE){multhit[l]=0;}
if(indicator[l] > POPSIZE || indicator[l]<0){cout<<i<<" "<<"Error2222\n";return(0);}
}

int count=0;
for(m=0;m<POPSIZE;m++){
k=0;trp=size[m];l=0;
while(l< trp){
trk = odd[m].chr->q[l];
if(indicator[trk]==POPSIZE){++l;continue;}
else{odd[m].chr->q[k] = trk;++k;++l;continue;}
}size[m]=k;count=count+size[m];}
cout<<"Generation "<<i<<" Total number of mutations: "<<count<<"\n";
}
}

}


if(flag == 1){continue;}


ofstream output1("finalpopulation.txt",ios::out);
for(m=0;m<SSIZE;m++){
for(l=0;l<siz[m];l++){
if(GEN%2==0){output1<<even[m].chr->q[l]<<" ";}
else{output1<<odd[m].chr->q[l]<<" ";}
}output1<<"\n";}
output1.close();


finished=1;}







}

  




  
  






