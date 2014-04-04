#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h> 

 
//#define square(a) (a*a) 
#define cube(a) (a*a*a) 

#define frac (0.95)

#define OVD (200.)
#define rho_crit (27.755e10)

#define _R200
#define _ONLYMUTUAL
#define _VERB
#define _DEBUG
//#define _ULTRADEBUG
#define _PERIODIC
#define _CRIT
//#define _ENOUGH
#define _FRACT1
//#define _SQUARE
#define _CUBE
//#define _RNDCELL
#define MAXTRIALS 10
//#define _NOCHECK
#define _ONLYBIG

float rho_ref;

long select_cell(); 
long select_cell_rnd(long *, long *, long *, float *, int, float); 
long select_heaviest_cell(long *, long *, long *, float *, int, float); 
long select_part(long );
int exclude_cell(long,float , float *, float *, float *, long ,long, long);
void exclude(long,float,float *,float *,float *,long ,long,long);
int check_HaloR(long ,float *,float *,float *,float *);
int check_HaloR_in_mesh(long,float *, float *, float * , float *,long,long,long);
int check_HaloR_in_cell(long ,float *, float *, float * , float *,long ,long,long);



long **ListOfPart, **ListOfHalos, *Nexcluded, *excluded, *NPartPerCell, *NHalosPerCell;
float *MassLeft,lcell,Lbox;
long NCells,NTotCells;
float tol,mpart;




//provisional
float pi = 3.1416; 
float square(float a){
	return a*a;
}

/*long cube(long a){
	return a*a*a;
}
*/

float R200_from_mass(float Mass) {
	return  (float) pow((3./(4.*OVD*rho_ref*pi)*Mass),(1./3.));
}



void place_halos_Mglobal(long NHalosTot, float *HaloMass, long Nlin, long NTotPart, float *PartX, float *PartY, float *PartZ, float L, float mp, float fract, float Nmin, ,float *HaloX, float *HaloY, float *HaloZ){

fprintf(stderr,"Hi! This is place_halos.c\n");
fprintf(stdout,"Hi! This is place_halos.c\n");


//Initiallising -------------------------------------------------
	long i,j,k,lin_ijk,check;
	long *count, tol_nhalos=3;
	long ihalo,ilong, ipart, Nhalos;
	double invL = 1./L,diff;
	float R, Mmin;	
	time_t t0,t1,t2,t3,t4,t5;
	NCells = Nlin;
	lcell=L/NCells;
	Lbox = L;
	
	t0=time(NULL);
	NTotCells = cube(NCells);

	//temp
	float *HaloR;	
	
	//Allocate memory for the arrays 
	excluded = (long *) calloc(NTotPart, sizeof(long));
	NPartPerCell = (long *) calloc(NCells*NCells*NCells,sizeof(long));
	NHalosPerCell = (long *) calloc(NCells*NCells*NCells,sizeof(long));
	MassLeft = (float *) calloc(NCells*NCells*NCells,sizeof(float));
	count = (long *) calloc(NCells*NCells*NCells,sizeof(long));
	Nexcluded = (long *) calloc(cube(NCells),sizeof(long)); 

	HaloR = (float *) calloc(NHalosTot,sizeof(float));

	//Initiallise random numbers
	srand (time(NULL));

	#ifdef _CRIT
	rho_ref = rho_crit;
	#else
	rho_ref = mp/(L*L*L)*NTotPart;
	#endif	



#ifdef _CRIT
	fprintf(stderr,"#def _CRIT\n");
#endif
#ifdef _EXHALOS
	fprintf(stderr,"#def _EXHALOS\n");
#endif
#ifdef _VERB
	fprintf(stderr,"#def _VERB\n");
#endif 
#ifdef _DEBUG
	fprintf(stderr,"#def _DEBUG \n");
#endif
#ifdef _ULTRADEBUG
	fprintf(stderr,"#def _ULTRADEBUG \n");
#endif
#ifdef _PERIODIC
	fprintf(stderr,"#def _PERIODIC\n");
#endif
#ifdef _FRACT1
	fprintf(stderr,"#def _FRACT1\n");
#endif
#ifdef _SQUARE
	fprintf(stderr,"#def _SQUARE\n");
#endif
#ifdef _CUBE
	fprintf(stderr,"#def _CUBE\n");
#endif
#ifdef _RNDCELL
	fprintf(stderr,"#def _RNDCELL\n");
#endif
#ifdef _NOCHECK
	fprintf(stderr,"#def _NOCHECK\n");
#endif
#ifdef _ONLYBIG
	fprintf(stderr,"#def _ONLYBIG\n");
#endif

	
	#ifdef _VERB
	fprintf(stderr,"\nMassFunction computed globally. Particles and Halos placed in %ld^3 cells\n",NCells);
	fprintf(stderr,"Exclusion done only with haloes.\n",NCells);
	fprintf(stderr,"BOX = %f  lcell =%f   rho_ref = %e  invL %f\n",L,lcell,rho_ref,invL);
	fprintf(stderr,"Nhalos = %ld NPart = %ld\n",NHalosTot, NTotPart);
	fprintf(stderr,"\nRAND_MAX=%d\n",RAND_MAX);
	fprintf(stderr,"X[0] = %f Y[0] = %f Z[0] = %f\n",PartX[0],PartY[0],PartZ[0]);
	fprintf(stderr,"X[1] = %f Y[1] = %f Z[1] = %f\n",PartX[1],PartY[1],PartZ[1]);
	fprintf(stderr,"M[0] = %f \n",HaloMass[0]);
	fprintf(stderr,"M[1] = %f \n",HaloMass[1]);
	fprintf(stderr,"\nExclusion done only with haloes. Minimmum mass= %e. Minimum part per halo = %f. Effective mp (not the real one) %e\n",HaloMass[NHalosTot-1],Nmin,mp);
	#endif	
	tol=tol_nhalos*HaloMass[NHalosTot-1];
	fprintf(stderr," Usiing mass tolerance %e, corresponging to %ld halos\n",tol,tol_nhalos);


	//Mmin = HaloMass[NHalosTot-1];
	if (Nmin<1){
		fprintf(stderr,"ERROR: Halo mass can not me lower than particle mass\n");
		exit(0);
	}
	t1=time(NULL);
 	diff = difftime(t1,t0);
	fprintf(stderr,"time of initialisation %f\n",diff);
// ------------------------------------------------- Initiallised



#ifdef _VERB
	fprintf(stderr,"Assigning particles to grid ...\n");
#endif


//Assign particles to grid ------------------------------------
	//count particles per cell
	for (ilong=0;ilong<NTotPart;ilong++) {
		i = (long) (invL * PartX[ilong]*NCells);
		j = (long) (invL * PartY[ilong]*NCells);
		k = (long) (invL * PartZ[ilong]*NCells);
		if (i<0 || i>=NCells || j<0 || j>=NCells || k<0 || k>=NCells){
			fprintf(stderr,"ERROR: Particle %ld at [%f,%f,%f] seems to be out of the right box interval [0.,%f), placed at cell [%ld,%ld,%ld]\n",ilong,PartX[ilong],PartY[ilong],PartZ[ilong],L,i,j,k);	
			exit(0);
		}
		lin_ijk = k+j*NCells+i*NCells*NCells;
		NPartPerCell[lin_ijk]++;
#ifdef _DEBUG
		if(ilong<10 || ilong > NTotPart -10 || ilong==243666)
			fprintf(stderr,"ipart=%ld  cell: %ld=[%ld,%ld,%ld] N=%ld, Pos= [%f,%f,%f]\n",ilong,lin_ijk,i,j,k,NPartPerCell[lin_ijk],PartX[ilong],PartY[ilong],PartZ[ilong]);
#endif
	}
#ifdef _DEBUG
	fprintf(stderr,"... particles counted ...\n");
	t2=time(NULL);
 	diff = difftime(t2,t1);
	fprintf(stderr,"time counting %f\n",diff);
#endif
	//Alloc Enough Memory
	ListOfPart = (long **) calloc(NCells*NCells*NCells,sizeof(long *));
	ListOfHalos = (long **) calloc(NCells*NCells*NCells,sizeof(long *));
	for (i=0;i<NCells;i++){
	for (j=0;j<NCells;j++){
	for (k=0;k<NCells;k++){
		lin_ijk = k+j*NCells+i*NCells*NCells;
		ListOfPart[lin_ijk] = (long *) calloc(NPartPerCell[lin_ijk],sizeof(long));
		Nhalos = (long) (NPartPerCell[lin_ijk]/Nmin);
		ListOfHalos[lin_ijk] = (long *) calloc(Nhalos,sizeof(long));
		MassLeft[lin_ijk] = NPartPerCell[lin_ijk]*mp; 
#ifdef _ULTRADEBUG
		if (lin_ijk<10 || lin_ijk > (NCells*NCells*NCells) - 10){
			fprintf(stderr,"Allocated %ld (longs) in ListOfPart(%ld=[%ld,%ld,%ld])\n",NPartPerCell[lin_ijk],lin_ijk,i,j,k);
			fprintf(stderr,"Allocated %ld (longs) in ListOfHalos(%ld=[%ld,%ld,%ld])\n",Nhalos,lin_ijk,i,j,k);
		}
#endif		
	}	
	}
	}
#ifdef _DEBUG
	fprintf(stderr,"... memory allocated ...\n");
	t3=time(NULL);
 	diff = difftime(t3,t2);
	fprintf(stderr,"time allocating %f\n",diff);
#endif

	for (ilong=0;ilong<NTotPart;ilong++) {
		i = (long) (invL * PartX[ilong]*NCells);
		j = (long) (invL * PartY[ilong]*NCells);
		k = (long) (invL * PartZ[ilong]*NCells);
		lin_ijk = k+j*NCells+i*NCells*NCells;
		ListOfPart[lin_ijk][count[lin_ijk]] = ilong;
		count[lin_ijk]++;
	}
//----------------------------------- Particles assigned to grid

#ifdef _VERB
	fprintf(stderr," ...done\n\n");
	t4=time(NULL);
 	diff = difftime(t4,t3);
	fprintf(stderr,"time of the actual assignment %f\n",diff);
#endif
#ifdef _DEBUG
	fprintf(stderr," Mass Function\n");
	fprintf(stderr,"%e \n",HaloMass[0]);
	fprintf(stdout,"%e \n",HaloMass[1]);
	for (ihalo=0;ihalo<15;ihalo++){
		fprintf(stderr,"halo %ld: ",ihalo);
		fprintf(stderr,"M=%e\n",HaloMass[ihalo]);
	}
#endif

#ifdef _VERB
	fprintf(stderr,"\nPlacing Halos...\n\n");
#endif

	for (ihalo=0;ihalo<NHalosTot;ihalo++){

		#ifdef _DEBUG
		fprintf(stderr,"\n- Halo %ld ",ihalo);
		#endif
		#ifdef _RNDCELL
		lin_ijk = select_cell_rnd(&i,&j,&k,MassLeft,NCells,HaloMass[ihalo]);
		#else
		lin_ijk = select_cell();
		#endif
		if(lin_ijk<0) {
			fprintf(stderr,"Maximum Mass cell was %e\n",MassLeft[select_heaviest_cell(&i,&j,&k,MassLeft,NCells,HaloMass[ihalo])]);
			break;
		}
		k=lin_ijk%(NCells);
		j=((lin_ijk-k)/NCells)%NCells;
		i=(lin_ijk-k-j*NCells)/(NCells*NCells);
		#ifdef _DEBUG
		fprintf(stderr,"\n assigned to cell %ld=[%ld,%ld,%ld]\n Mass in Cell: Before M=%e    ",lin_ijk,i,j,k,(MassLeft)[lin_ijk]);
		#endif
		(MassLeft)[lin_ijk]=(MassLeft)[lin_ijk]-HaloMass[ihalo];

		if ((MassLeft[lin_ijk])<0.){
			(MassLeft)[lin_ijk]=0.;
			fprintf(stderr,"WARNING: getting negative masses\n");
		}
		#ifdef _DEBUG
		fprintf(stderr,"After M=%e    ",(MassLeft)[lin_ijk]);
		#endif

		do {
			ipart = select_part(lin_ijk);		
               		HaloX[ihalo] = PartX[ipart];
               		HaloY[ihalo] = PartY[ipart];
               		HaloZ[ihalo] = PartZ[ipart];
			R=R200_from_mass(HaloMass[ihalo]);
			HaloR[ihalo]=R;
			#ifndef _NOCHECK
			check = check_HaloR_in_mesh(ihalo,HaloX,HaloY,HaloZ,HaloR,i,j,k);
			#else
			check = 1;
			#endif
			if (check==0){
				#ifdef _DEBUG
				fprintf(stderr,"Refused part : %ld\n",ipart);
				#endif
			}

		} while (check==0);
		#ifdef _DEBUG
		fprintf(stderr,"halo %ld assigned to particle %ld at [%f,%f,%f]. R= %f, M= %e\n",ihalo,ipart,HaloX[ihalo],HaloY[ihalo],HaloZ[ihalo],R,HaloMass[ihalo]);
		#endif
		ListOfHalos[lin_ijk][NHalosPerCell[lin_ijk]]=ihalo;
		NHalosPerCell[lin_ijk]++;
	}

#ifdef _VERB
	t5=time(NULL);
 	diff = difftime(t5,t4);
	fprintf(stderr,"time placing %f\n",diff);
 	diff = difftime(t5,t0);
	fprintf(stderr,"total time in .c %f\n",diff);
	fprintf(stderr,"\nEverything done!!!\n");
#endif
	free(count); free(NPartPerCell); free(ListOfPart); free(excluded);
}









void place_halos_ranked(long NHalosTot, float *HaloMass, long Nlin, long NTotPart, float *PartX, float *PartY, float *PartZ, float L, float mp, float fract, float Nmin, float *HaloX, float *HaloY, float *HaloZ){
	
		

//Initiallising -------------------------------------------------
	long i,j,k,lin_ijk,check,trial;
	long *count;
	long ihalo,ilong, ipart;
	double invL = 1./L, diff;
	float R;	
	long Nhalos;
	time_t t0, t1, t2,t3,t4,t5;
	t0=time(NULL);
#ifdef _CRIT
	rho_ref = rho_crit;
#else
	rho_ref = (HaloMass[NHalosTot]/Nmin)/(L*L*L)*NTotPart;
#endif
	
#ifdef _FRACT1
	mpart=mp;
#else
	mpart=fract*mp;
#endif
	if (Nlin!=0){
		NCells=Nlin;
		lcell = L/NCells;
		fprintf(stderr,"Using RANKED approach, with lcell approx %f\n",lcell);
	}
	else {	
		lcell = 2*R200_from_mass(HaloMass[0]);
		NCells = (long) (L/lcell);	
		lcell = L/NCells;
		fprintf(stderr,"Using RANKED approach, with lcell computed automatically. M[0]=%e, R[0]=%f, lcell aprox %f\n",HaloMass[0],R200_from_mass(HaloMass[0]),lcell);
	}


	fprintf(stderr,"Box size: %f. Number of Cells: %ld. Actual length of cell lcell=%f\n",L,NCells,lcell);

	Lbox = L;
	NTotCells = cube(NCells);

	//temp
	float *HaloR;	
	
	//Allocate memory for the arrays 
	excluded = (long *) calloc(NTotPart, sizeof(long));
	NPartPerCell = (long *) calloc(NCells*NCells*NCells,sizeof(long));
	MassLeft = (float *) calloc(NCells*NCells*NCells,sizeof(float));
	count = (long *) calloc(NCells*NCells*NCells,sizeof(long));
	Nexcluded = (long *) calloc(cube(NCells),sizeof(long)); 
	
	HaloR = (float *) calloc(NHalosTot,sizeof(float));


	//Initiallise random numbers
	srand (time(NULL));


#ifdef _R200
	fprintf(stderr,"#def _R200\n");
#endif
#ifdef _CRIT
	fprintf(stderr,"#def _CRIT\n");
#endif
#ifdef _MUTUAL
	fprintf(stderr,"#def _MUTUAL\n");
#endif
#ifdef _ONLYMUTUAL
	fprintf(stderr,"#def _ONLYMUTUAL\n");
#endif
#ifdef _VERB
	fprintf(stderr,"#def _VERB\n");
#endif 
#ifdef _DEBUG
	fprintf(stderr,"#def _DEBUG \n");
#endif
#ifdef _ULTRADEBUG
	fprintf(stderr,"#def _ULTRADEBUG \n");
#endif
#ifdef _PERIODIC
	fprintf(stderr,"#def _PERIODIC\n");
#endif
#ifdef _EXHALOS
	fprintf(stderr,"#def _EXHALOS");
#endif
#ifdef _FRACT1
	fprintf(stderr,"#def _FRACT1");
#else
	fprintf(stderr,"using f=%f\n",fract);
#endif
#ifdef _SQUARE
	fprintf(stderr,"#def _SQUARE\n");
#endif
#ifdef _CUBE
	fprintf(stderr,"#def _CUBE\n");
#endif
#ifdef _ONLYBIG
	fprintf(stderr,"#def _ONLYBIG\n");
#endif
#ifdef _NOCHECK
	fprintf(stderr,"#def _NOCHECK\n");
#endif

	
	#ifdef _VERB
	fprintf(stderr,"\nMassFunction computed globally, Particles placed in %ld^3 cells\n",NCells);
	fprintf(stderr,"BOX = %f  lcell =%f   rho_ref = %e  invL %f\n",L,lcell,rho_ref,invL);
	fprintf(stderr,"mp = %f  fract =%f  mp used =\n",mp,fract,mpart);
	fprintf(stderr,"RAND_MAX=%d\n",RAND_MAX);
	fprintf(stderr,"Nhalos = %ld NPart = %ld\n",NHalosTot, NTotPart);
	fprintf(stderr,"X[0] = %f Y[0] = %f Z[0] = %f\n",PartX[0],PartY[0],PartZ[0]);
	fprintf(stderr,"X[1] = %f Y[1] = %f Z[1] = %f\n",PartX[1],PartY[1],PartZ[1]);
	fprintf(stderr,"M[0] = %f \n",HaloMass[0]);
	fprintf(stderr,"M[1] = %f \n",HaloMass[1]);
	#endif
	t1=time(NULL);
 	diff = difftime(t1,t0);
	fprintf(stderr,"time of initialisation %f\n",diff);
// ------------------------------------------------- Initiallised



//Assign particles to grid ------------------------------------
	//count particles per cell
	for (ilong=0;ilong<NTotPart;ilong++) {
		i = (long) (invL * PartX[ilong]*NCells);
		j = (long) (invL * PartY[ilong]*NCells);
		k = (long) (invL * PartZ[ilong]*NCells);
		if (i<0 || i>=NCells || j<0 || j>=NCells || k<0 || k>=NCells){
			fprintf(stderr,"ERROR: Particle %ld at [%f,%f,%f] seems to be out of the right box interval [0.,%f), placed at cell [%ld,%ld,%ld]\n",ilong,PartX[ilong],PartY[ilong],PartZ[ilong],L,i,j,k);	
		}
		lin_ijk = k+j*NCells+i*NCells*NCells;
		NPartPerCell[lin_ijk]++;
#ifdef _DEBUG
		if(ilong<10 || ilong > NTotPart -10 || ilong==243666)
			fprintf(stderr,"ipart=%ld  cell: %ld=[%ld,%ld,%ld] N=%ld, Pos= [%f,%f,%f]\n",ilong,lin_ijk,i,j,k,NPartPerCell[lin_ijk],PartX[ilong],PartY[ilong],PartZ[ilong]);
#endif
	}
#ifdef _DEBUG
	fprintf(stderr,"... particles counted ...\n");
	t2=time(NULL);
 	diff = difftime(t2,t1);
	fprintf(stderr,"time counting %f\n",diff);
#endif
	//Alloc Enough Memory
	ListOfPart = (long **) calloc(NCells*NCells*NCells,sizeof(long *));
	ListOfHalos = (long **) calloc(NCells*NCells*NCells,sizeof(long *));
	NHalosPerCell = (long *) calloc(NCells*NCells*NCells,sizeof(long *));
	for (i=0;i<NCells;i++){
	for (j=0;j<NCells;j++){
	for (k=0;k<NCells;k++){
		lin_ijk = k+j*NCells+i*NCells*NCells;
		ListOfPart[lin_ijk] = (long *) calloc(NPartPerCell[lin_ijk],sizeof(long));
		Nhalos = (long) (NPartPerCell[lin_ijk]/Nmin);
		ListOfHalos[lin_ijk] = (long *) calloc(Nhalos,sizeof(long));
		MassLeft[lin_ijk] = NPartPerCell[lin_ijk]*mpart; 

#ifdef _ULTRADEBUG
		if (lin_ijk<10 || lin_ijk > (NCells*NCells*NCells) - 10 || lin_ijk ==16943336){
			fprintf(stderr,"Allocated %ld (longs) in ListOfPart(%ld=[%ld,%ld,%ld])\n",NPartPerCell[lin_ijk],lin_ijk,i,j,k);
			fprintf(stderr,"Allocated %ld (longs) in ListOfHalos(%ld=[%ld,%ld,%ld])\n",Nhalos,lin_ijk,i,j,k);
		}
#endif
	}	
	}
	}
#ifdef _DEBUG
	fprintf(stderr,"... memory allocated ...\n");
	t3=time(NULL);
 	diff = difftime(t3,t2);
	fprintf(stderr,"time allocating %f\n",diff);
#endif

	for (ilong=0;ilong<NTotPart;ilong++) {
		i = (long) (invL * PartX[ilong]*NCells);
		j = (long) (invL * PartY[ilong]*NCells);
		k = (long) (invL * PartZ[ilong]*NCells);
		lin_ijk = k+j*NCells+i*NCells*NCells;
		ListOfPart[lin_ijk][count[lin_ijk]] = ilong;
		count[lin_ijk]++;
	}
//----------------------------------- Particles assigned to grid

	t4=time(NULL);
 	diff = difftime(t4,t3);
	fprintf(stderr,"time of the actual assignment %f\n",diff);
	fprintf(stderr,"Placing halos down...\n");
//Placing Halos-------------------------------------------------
	for (ihalo=0;ihalo<NHalosTot;ihalo++){
		trial = 0;
#ifdef _DEBUG
		fprintf(stderr,"\n- Halo %ld with mass %e",ihalo,HaloMass[ihalo]);
#endif
		do {
			lin_ijk = select_heaviest_cell(&i,&j,&k,MassLeft,NCells,HaloMass[ihalo]);	
#ifdef _DEBUG
			fprintf(stderr," assigned to cell %ld = [%ld,%ld,%ld]",lin_ijk,i,j,k);
#endif
			ipart = select_part(lin_ijk);		
#ifdef _DEBUG
			fprintf(stderr,", to particle %ld",ipart);
#endif
               		HaloX[ihalo] = PartX[ipart];
              		HaloY[ihalo] = PartY[ipart];
               		HaloZ[ihalo] = PartZ[ipart];
#ifdef _DEBUG
			fprintf(stderr," at [%f,%f,%f]. ",HaloX[ihalo],HaloY[ihalo],HaloZ[ihalo]);
#endif
			#ifndef _NOCHECK
			check = check_HaloR_in_mesh(ihalo,HaloX,HaloY,HaloZ,HaloR,i,j,k);
			#else 
			check=1;
			#endif
			if (check==0){
				#ifdef _DEBUG
				fprintf(stderr,"Refused part : %ld, check = %ld\n",ipart,check);
				#endif
				trial++;
			}
			if (trial==MAXTRIALS){
				fprintf(stderr,"WARNING: after %ld trials, couldnt place a halo here. Decreasing mass in cell\n",trial);
				trial =0;
				MassLeft[lin_ijk]-=HaloMass[ihalo];
			}
		} while (check==0);
		ListOfHalos[lin_ijk][NHalosPerCell[lin_ijk]]=ihalo;
		NHalosPerCell[lin_ijk]++;
		R=R200_from_mass(HaloMass[ihalo]);
	#ifdef _DEBUG
		fprintf(stderr,"Mass in cell before %e, ",MassLeft[lin_ijk]);
	#endif
	#ifdef _ONLYMUTUAL
		MassLeft[lin_ijk]-=HaloMass[ihalo];
		if (MassLeft[lin_ijk]<0.)
			MassLeft[lin_ijk]=0.;
	#else
		exclude(ipart,R,PartX,PartY,PartZ,i,j,k);
	#endif
	#ifdef _DEBUG
		fprintf(stderr,"Mass in cell after %e, ",MassLeft[lin_ijk]);
	#endif
	
	#ifdef _DEBUG
		fprintf(stderr,"after %e. \n",MassLeft[lin_ijk]);
	#endif
	}


	t5=time(NULL);
 	diff = difftime(t5,t4);
	fprintf(stderr,"time placing %f\n",diff);
 	diff = difftime(t5,t0);
	fprintf(stderr,"total time in .c %f\n",diff);
	fprintf(stderr,"Everything done!\n");

}




void exclude(long ipart,float R,float *PartX,float *PartY,float *PartZ,long i,long j,long k) {

	float X = PartX[ipart];
	float Y = PartY[ipart];
        float Z = PartZ[ipart];
	float R2 = R*R;  
	long faces=0, edges=0,vertices=0;

		if (2.*R>lcell){
			fprintf(stderr,"WARNING: there are haloes with diameter D= 2R = %f greater than the size of the cell l = %f.\nPlease change the size of the cell\n",2.*R,lcell);
		}

		exclude_cell(ipart,R2,PartX,PartY,PartZ,i,j,k);
			

		//Boundaries: faces
		if (X+R>lcell*(i+1)) {
			faces++;
			#ifdef _DEBUG
			fprintf(stderr,"Radius extends beyond face X+\n");
			#endif
			exclude_cell(ipart,R2,PartX,PartY,PartZ,i+1,j,k);
		}
		if (X-R<lcell*i) {
			faces++;
			#ifdef _DEBUG
			fprintf(stderr,"Radius extends beyond face X-\n");
			#endif
			exclude_cell(ipart,R2,PartX,PartY,PartZ,i-1,j,k);
		}
		if (Y+R>lcell*(j+1)) {
			faces++;
			#ifdef _DEBUG
			fprintf(stderr,"Radius extends beyond face Y+\n");
			#endif
			exclude_cell(ipart,R2,PartX,PartY,PartZ,i,j+1,k);
		}
		if (Y-R<lcell*j) {
			faces++;
			#ifdef _DEBUG
			fprintf(stderr,"Radius extends beyond face Y-\n");
			#endif
			exclude_cell(ipart,R2,PartX,PartY,PartZ,i,j-1,k);
		}
		if (Z+R>lcell*(k+1)) {
			faces++;
			#ifdef _DEBUG
			fprintf(stderr,"Radius extends beyond face Z+\n");
			#endif
			exclude_cell(ipart,R2,PartX,PartY,PartZ,i,j,k+1);
		}
		if (Z-R<lcell*k) {
			faces++;
			#ifdef _DEBUG
			fprintf(stderr,"Radius extends beyond face Z-\n");
			#endif
			exclude_cell(ipart,R2,PartX,PartY,PartZ,i,j,k-1);
		}
		

		//Boundaries: edges		
		  //XY
		  if (square(X-lcell*(i+1))+square(Y-lcell*(j+1))<R2){
			edges++;
			#ifdef _DEBUG
			fprintf(stderr,"%f + %f < %f\n",square(X-lcell*(i+1)),square(Y-lcell*(j+1)),R2);
			fprintf(stderr,"Radius extends beyond edge X+Y+\n");
			#endif
			exclude_cell(ipart,R2,PartX,PartY,PartZ,i+1,j+1,k);
		  }
		  if (square(X-lcell*(i+1))+square(Y-lcell*(j))<R2){
			edges++;
			#ifdef _DEBUG
			fprintf(stderr,"Radius extends beyond edge X+Y-\n");
			#endif
			exclude_cell(ipart,R2,PartX,PartY,PartZ,i+1,j-1,k);
		  }
		  if (square(X-lcell*(i))+square(Y-lcell*(j+1))<R2){
			edges++;
			#ifdef _DEBUG
			fprintf(stderr,"Radius extends beyond edge X-Y+\n");
			#endif
			exclude_cell(ipart,R2,PartX,PartY,PartZ,i-1,j+1,k);
		  }
		  if (square(X-lcell*(i))+square(Y-lcell*(j))<R2){
			edges++;
			#ifdef _DEBUG
			fprintf(stderr,"Radius extends beyond edge X-Y-\n");
			#endif
			exclude_cell(ipart,R2,PartX,PartY,PartZ,i-1,j-1,k);
		  }
		  //YZ
		  if (square(Z-lcell*(k+1))+square(Y-lcell*(j+1))<R2){
			edges++;
			#ifdef _DEBUG
			fprintf(stderr,"Radius extends beyond edge Y+Z+\n");
			#endif
			exclude_cell(ipart,R2,PartX,PartY,PartZ,i,j+1,k+1);
		  }
		  if (square(Z-lcell*(k+1))+square(Y-lcell*(j))<R2){
			edges++;
			#ifdef _DEBUG
			fprintf(stderr,"Radius extends beyond edge Z+Y-\n");
			#endif
			exclude_cell(ipart,R2,PartX,PartY,PartZ,i,j-1,k+1);
		  }
		  if (square(Z-lcell*(k))+square(Y-lcell*(j+1))<R2){
			edges++;
			#ifdef _DEBUG
			fprintf(stderr,"Radius extends beyond edge Z-Y+\n");
			#endif
			exclude_cell(ipart,R2,PartX,PartY,PartZ,i,j+1,k-1);
		  }
		  if (square(Z-lcell*(k))+square(Y-lcell*(j))<R2){
			edges++;
			#ifdef _DEBUG
			fprintf(stderr,"Radius extends beyond edge Z-Y-\n");
			#endif
			exclude_cell(ipart,R2,PartX,PartY,PartZ,i,j-1,k-1);
		  }
		  //XZ
		  if (square(X-lcell*(i+1))+square(Z-lcell*(k+1))<R2){
			edges++;
			#ifdef _DEBUG
			fprintf(stderr,"Radius extends beyond edge X+Z+\n");
			#endif
			exclude_cell(ipart,R2,PartX,PartY,PartZ,i+1,j,k+1);
		  }
		  if (square(X-lcell*(i+1))+square(Z-lcell*(k))<R2){
			edges++;
			#ifdef _DEBUG
			fprintf(stderr,"Radius extends beyond edge X+Z-\n");
			#endif
			exclude_cell(ipart,R2,PartX,PartY,PartZ,i+1,j,k-1);
		  }
		  if (square(X-lcell*(i))+square(Z-lcell*(k+1))<R2){
			edges++;
			#ifdef _DEBUG
			fprintf(stderr,"Radius extends beyond edge X-Z+\n");
			#endif
			exclude_cell(ipart,R2,PartX,PartY,PartZ,i-1,j,k+1);
		  }
		  if (square(X-lcell*(i))+square(Z-lcell*(k))<R2){
			edges++;
			#ifdef _DEBUG
			fprintf(stderr,"Radius extends beyond edge X-Z-\n");
			#endif
			exclude_cell(ipart,R2,PartX,PartY,PartZ,i-1,j,k-1);
		  }	


		//Boundaries: vertices
		if (square(X-lcell*(i+1))+square(Y-lcell*(j+1))+square(Z-lcell*(k+1))<R2){
			vertices++;
			#ifdef _DEBUG
			fprintf(stderr,"Radius extends beyond vertex X+Y+Z+\n");
			#endif
			exclude_cell(ipart,R2,PartX,PartY,PartZ,i+1,j+1,k+1);
		}
		if (square(X-lcell*(i))+square(Y-lcell*(j+1))+square(Z-lcell*(k+1))<R2){
			vertices++;
			#ifdef _DEBUG
			fprintf(stderr,"Radius extends beyond vertex X-Y+Z+\n");
			#endif
			exclude_cell(ipart,R2,PartX,PartY,PartZ,i-1,j+1,k+1);
		}
		if (square(X-lcell*(i+1))+square(Y-lcell*(j))+square(Z-lcell*(k+1))<R2){
			vertices++;
			#ifdef _DEBUG
			fprintf(stderr,"Radius extends beyond vertex X+Y-Z+\n");
			#endif
			exclude_cell(ipart,R2,PartX,PartY,PartZ,i+1,j-1,k+1);
		}
		if (square(X-lcell*(i+1))+square(Y-lcell*(j+1))+square(Z-lcell*(k))<R2){
			vertices++;
			#ifdef _DEBUG
			fprintf(stderr,"Radius extends beyond vertex X+Y+Z-\n");
			#endif
			exclude_cell(ipart,R2,PartX,PartY,PartZ,i+1,j+1,k-1);
		}
		if (square(X-lcell*(i))+square(Y-lcell*(j))+square(Z-lcell*(k+1))<R2){
			vertices++;
			#ifdef _DEBUG
			fprintf(stderr,"Radius extends beyond vertex X-Y-Z+\n");
			#endif
			exclude_cell(ipart,R2,PartX,PartY,PartZ,i-1,j-1,k+1);
		}
		if (square(X-lcell*(i))+square(Y-lcell*(j+1))+square(Z-lcell*(k))<R2){
			vertices++;
			#ifdef _DEBUG
			fprintf(stderr,"Radius extends beyond vertex X-Y+Z-\n");
			#endif
			exclude_cell(ipart,R2,PartX,PartY,PartZ,i-1,j+1,k-1);
		}
		if (square(X-lcell*(i+1))+square(Y-lcell*(j))+square(Z-lcell*(k))<R2){
			vertices++;
			#ifdef _DEBUG
			fprintf(stderr,"Radius extends beyond vertex X+Y-Z-\n");
			#endif
			exclude_cell(ipart,R2,PartX,PartY,PartZ,i+1,j-1,k-1);
		}
		if (square(X-lcell*(i))+square(Y-lcell*(j))+square(Z-lcell*(k))<R2){
			vertices++;
			#ifdef _DEBUG
			fprintf(stderr,"Radius extends beyond vertex X-Y-Z-\n");
			#endif
			exclude_cell(ipart,R2,PartX,PartY,PartZ,i-1,j-1,k-1);
		}

		if (edges==1 && faces<2) fprintf(stderr,"ERROR1: Something wrong with the boundaries of the cells\n");
		if (edges>=2 && faces!=3) fprintf(stderr,"ERROR2: Something wrong with the boundaries of the cells: %ld edges   %ld faces \n",edges,faces);
		if (vertices==1 && (faces!=3 || edges!=3) ) fprintf(stderr,"ERROR3: Something wrong with the boundaries of the cells:   %ld vertices   %ld edges   %ld faces \n",vertices,edges,faces);
		if (faces>3)  fprintf(stderr,"ERROR4: Something wrong with the boundaries of the cells\n");
		if (edges>3)  fprintf(stderr,"ERROR5: Something wrong with the boundaries of the cells %ld vertices   %ld edges   %ld faces \n",vertices,edges,faces);
		if (vertices>1)  fprintf(stderr,"ERROR6: Something wrong with the boundaries of the cells\n");	

}

int exclude_cell(long ipart,float R2, float *PartX, float *PartY, float *PartZ, long i,long j, long k){
	long ilong, jpart;
	float X=PartX[ipart],Y=PartY[ipart],Z=PartZ[ipart];

#ifdef _PERIODIC
		if (i==-1){
			i = NCells -1;
			X = X + Lbox;
		}
		else if (i==NCells){
			i = 0;
			X = X - Lbox;
		}

		if (j==-1){
			j = NCells -1;
			Y = Y + Lbox;
		}
		else if (j==NCells){
			j = 0;
			Y = Y - Lbox;
		}
		if (k==-1){
			k = NCells -1;
			Z = Z + Lbox;
		}
		else if (k==NCells){
			k = 0;
			Z = Z - Lbox;
		}
#endif
	long lin_ijk = k+j*NCells+i*NCells*NCells,Nbef;
	if (MassLeft[lin_ijk]==0.){
		#ifdef _DEBUG
		fprintf(stderr,"No mass left in cell [%ld,%ld,%ld]  M=%f \n",i,j,k,MassLeft[lin_ijk]);
		#endif
		return 0;
	}
	Nbef = Nexcluded[lin_ijk];
	#ifdef _DEBUG
	fprintf(stderr,"Before: Nexcluded = %ld  . Excluding cell [%ld,%ld,%ld]...",Nexcluded[lin_ijk],i,j,k);
	#endif
	if (i>=0 && i<NCells && j>=0 && j<NCells && k>=0 && k<NCells){
		for (ilong=0; ilong<NPartPerCell[lin_ijk];ilong++){
                                jpart = ListOfPart[lin_ijk][ilong];
//				if (lin_ijk==0)
//						fprintf(stderr,"Particle i=%ld/%ld    ID=%ld at [%f,%f,%f]",ilong,NPartPerCell[lin_ijk],jpart,PartX[jpart],PartY[jpart],PartZ[jpart]);
                                if (excluded[jpart]==0){
                                        if (R2>(square(PartX[jpart]-X)+square(PartY[jpart]-Y)+square(PartZ[jpart]-Z))){
                                                excluded[jpart]=1;
                                                Nexcluded[lin_ijk]++;
//						if (lin_ijk==0)
//							fprintf(stderr," just excluded: %ld \n",excluded[jpart]);
                                        }
//					else 
//					if (lin_ijk==0)
//						fprintf(stderr,"  not excluded: %ld \n",excluded[jpart]);
					
					
                                }
				
//				else if (lin_ijk==0)
//						fprintf(stderr,"  was already excluded: %ld \n",excluded[jpart]);
                }
		
	}
	#ifdef _DEBUG
	fprintf(stderr,"After: Nexcluded = %ld  .  Out of %ld particles\n",Nexcluded[lin_ijk], NPartPerCell[lin_ijk]);
	#endif
	#ifndef _ONLYMUTUAL
	MassLeft[lin_ijk]=MassLeft[lin_ijk]-mpart*(Nexcluded[lin_ijk]-Nbef);
	#endif
	if (Nexcluded[lin_ijk]>frac*NPartPerCell[lin_ijk]){
		#ifdef _VERB
		fprintf(stderr,"WARNING: Cell %ld was completely excluded when there was still %e M_sun to be assigned\n",lin_ijk,MassLeft[lin_ijk]);
		fprintf(stderr,"WARNING: NPartPerCell %ld         Nexcluded %ld\n",Nexcluded[lin_ijk],NPartPerCell[lin_ijk]);
		#endif
		MassLeft[lin_ijk]=0.;
	}			
	return 0;
}



long select_part(long ijk){
	long i_rnd,ipart;
         do {
		i_rnd = (long) (NPartPerCell[ijk] * ((double)rand()/(RAND_MAX+1.0)));
                ipart = ListOfPart[ijk][i_rnd];
         } while (excluded[ipart]==1);	
	return ipart;
}

long select_heaviest_cell(long *x, long *y, long *z, float *MassArray, int N, float M) {
	long i,j,k,lin_ijk, out_ijk;
	float max=0.0;	
	for (i=0;i<N;i++){
        for (j=0;j<N;j++){
        for (k=0;k<N;k++){
		lin_ijk = k+j*N+i*N*N;
		if (max<(MassArray)[lin_ijk]){
			max=(MassArray)[lin_ijk];
			*x=i;
			*y=j;
			*z=k;
			out_ijk = lin_ijk;	
		}
	}
	}
	}
	return out_ijk;
}
long select_cell_rnd(long *x, long *y, long *z, float *MassArray, int N, float M) {
	long i,j,k;
	i = rand()%N;
	j = rand()%N;
	k = rand()%N;
	return k+j*N+i*N*N;
}
long select_cell() {
	long i;
	double TotMass = 0., d_rand;

	#ifdef _ENOUGH
		float threshold; 
		threshold=M-tol;
	#else
	float threshold = 6.4; //test
	#endif

	for (i=0;i<NTotCells;i++){
		#ifdef _ENOUGH
		if (MassLeft[i]>threshold)
		#endif


		#ifndef _SQUARE
		 #ifndef _CUBE
		TotMass += (double) MassLeft[i];
		 #else
		TotMass += (double) (MassLeft[i]*MassLeft[i]*1.0e-30*MassLeft[i]);
		 #endif
		#else
		TotMass += (double) (MassLeft[i]*MassLeft[i]);
		#endif	
	}
	d_rand = TotMass * ((double)rand()/(RAND_MAX+1.0));
#ifdef _ULTRADEBUG
	fprintf(stderr,"Rand: %e TotMass: %e  HaloMass: %e.    ",d_rand,TotMass,M);
#endif
	TotMass = 0.;

	for (i=0;i<NTotCells;i++){
	#ifdef _ENOUGH
	if ((MassLeft)[i]>threshold)
	#endif
		#ifndef _SQUARE
		 #ifndef _CUBE
		TotMass += (double) MassLeft[i];
		 #else
		TotMass += (double) (MassLeft[i]*MassLeft[i]*1.0e-30*MassLeft[i]);
		 #endif
		#else
		TotMass += (double) (MassLeft[i]*MassLeft[i]);
		#endif	
		#ifdef _ULTRADEBUG
		fprintf(stderr,"%e, ",TotMass);
		#endif
		if (TotMass > d_rand){
			return i;
		}
	}
	fprintf(stderr,"ERROR: no cell selected: TotMass = %e, rnd = %e. th= , tol =%e\n",TotMass,d_rand,tol);
	return -1;		
}

int check_HaloR_in_cell(long ipart,float *PartX, float *PartY, float *PartZ, float *PartR, long i,long j, long k){
        long jpart,jj;
        float X=PartX[ipart],Y=PartY[ipart],Z=PartZ[ipart],R=PartR[ipart];

#ifdef _PERIODIC
                if (i==-1){
                        i = NCells -1;
                        X = X + Lbox;
                }
                else if (i==NCells){
                        i = 0;
                        X = X - Lbox;
                }

                if (j==-1){
                        j = NCells -1;
                        Y = Y + Lbox;
                }
                else if (j==NCells){
                        j = 0;
                        Y = Y - Lbox;
                }
                if (k==-1){
                        k = NCells -1;
                        Z = Z + Lbox;
                }
                else if (k==NCells){
                        k = 0;
                        Z = Z - Lbox;
                }
#endif
  if (i>=0 && i<NCells && j>=0 && j<NCells && k>=0 && k<NCells){
#ifdef _ULTRADEBUG
       fprintf(stderr,"Checking cell [%ld,%ld,%ld] for halo %ld",i,j,k,ipart);

#endif
        long lin_ijk = k+j*NCells+i*NCells*NCells;
#ifdef _ULTRADEBUG
	fprintf(stderr,"= %ld.",lin_ijk);
#endif
	for (jj=0; jj<NHalosPerCell[lin_ijk]; jj++){
#ifdef _ULTRADEBUG
		fprintf(stderr,"jj=%ld/%ld ",jj,NHalosPerCell[lin_ijk]);
#endif
		jpart=ListOfHalos[lin_ijk][jj];
#ifdef _ULTRADEBUG
		fprintf(stderr,"jpart=%ld ",jpart);
#endif
		#ifdef _ONLYBIG
		if ((square(X-PartX[jpart])+square(Y-PartY[jpart])+square(Z-PartZ[jpart]))<square(PartR[jpart])) {
		#else
		if ((square(X-PartX[jpart])+square(Y-PartY[jpart])+square(Z-PartZ[jpart]))<square(R+PartR[jpart])) {
		#endif
#ifdef _DEBUG
                        fprintf(stderr,"\nChecking cell [%ld,%ld,%ld] for halo %ld",i,j,k,ipart);
                        fprintf(stderr," lin_ijk= %ld.",lin_ijk);
                        fprintf(stderr,"jj=%ld/%ld ",jj,NHalosPerCell[lin_ijk]);
                        fprintf(stderr,"jpart=%ld ",jpart);
                        fprintf(stderr,"refused!\n");
#endif
			return 0;
		}
	}
	return 1;
  }
  else {
	fprintf(stderr,"WARNING: Computing distances outside the box\n");
	return 1;
  }
}




int check_HaloR_in_mesh(long ihalo,float *X, float *Y, float *Z , float *R,long i,long j,long k){
	int l,m,n;
	for (l=i-1;l<=i+1;l++){
	for (m=j-1;m<=j+1;m++){
	for (n=k-1;n<=k+1;n++){
		if (check_HaloR_in_cell(ihalo,X,Y,Z,R,l,m,n)==0)
			return 0;

	}
	}
	}

	return 1;
}
int check_HaloR(long i,float *X, float *Y, float *Z , float *R){
	long j;
	for (j=0; j<i; j++){
		if ((square(X[i]-X[j])+square(Y[i]-Y[j])+square(Z[i]-Z[j]))<square(R[i]+R[j])){
			#ifdef _DEBUG
			fprintf(stderr,"\nOverlap between previous halo %ld placed at [%f,%f,%f] with R=%f and the selected particle at [%f,%f,%f] with R=%f. ",j,X[j],Y[j],Z[j],R[j],X[i],Y[i],Z[i],R[i]);
			#endif
			return 0;	
		}
	}
	return 1;
}
