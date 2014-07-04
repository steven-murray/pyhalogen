#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h> 
 
//#define square(a) (a*a) 
//#define cube(a) (a*a*a) 
#define frac (0.95)

#define OVD (200.)
#define rho_crit (27.755e10)
//#define mp (5.29384e+09/20)
//#define M20 (5.29384e+09)
//#define M20 (5.58e11*20.) // change anything related to this




#define _R200
#define _MUTUAL
#define _VERB
#define _DEBUG
#define _PERIODIC
#define _CRIT
//#define _ENOUGH

float rho_ref;

long select_cell(long *, long *, long *, float *, int, float); 
long select_part(long );
int exclude_cell(long,float , float *, float *, float *, long ,long, long);
void exclude(long,float,float *,float *,float *,long ,long,long);
int check_HaloR(long ,float *,float *,float *,float *);

long **ListOfPart, *Nexcluded, *excluded, *NPartPerCell;
float *MassLeft,lcell,Lbox;
long NCells;





float pi = 3.1416; 
float square(float a){
	return a*a;
}
float cube(float a){
	return a*a*a;
}



float R200_from_mass(float Mass) {
	return  (float) pow((3./(4.*OVD*rho_ref*pi)*Mass),(1./3.));
}





 
void place_halos_samecells(long NHalosTot, float *HaloMass, long NtotCells, int *FirstHaloInCell, long NTotPart, float *PartX, float *PartY, float *PartZ, float L, float mp, float *HaloX, float *HaloY, float *HaloZ){


fprintf(stderr,"Hi! This is place_halos.c  (no glob)\n");
fprintf(stdout,"Hi! This is place_halos.c (no glob)\n");
 
//Initiallising -------------------------------------------------
	long i,j,k,lin_ijk,rnd, starthalo, endhalo;
	long *NPartPerCell, *count, *excluded;
	long ihalo,ilong, ipart,jpart, Nexcluded;
	long **ListOfPart;
	float  R, R2, *MassLeft;
	float invL = 1./L;
	long NCells=(long) pow(NtotCells+1.0,1./3); //This should be changed!!

	
	//Allocate memory for the arrays 
	excluded = (long *) calloc(NTotPart, sizeof(long));
	NPartPerCell = (long *) calloc(NCells*NCells*NCells,sizeof(long));
	MassLeft = (float *) calloc(NCells*NCells*NCells,sizeof(float));
	count = (long *) calloc(NCells*NCells*NCells,sizeof(long));
	
	//Initiallise random numbers
	srand (time(NULL));

#ifdef _CRIT
	rho_ref = rho_crit;
#else
	rho_ref = mp/(L*L*L)*NTotPart;
#endif	

	#ifdef _VERB
	 fprintf(stderr,"BOX = %f    rho_mean = %e\n",L,rho_ref);
	 fprintf(stderr,"RAND_MAX=%d\n",RAND_MAX);
	 fprintf(stderr,"Nlincells = %ld Ntotcells = %ld Nhalos = %ld\n",NCells,NtotCells,NHalosTot);
	#endif
// ------------------------------------------------- Initiallised



//Assign particles to grid ------------------------------------

	//count particles per node
	for (ilong=0;ilong<NTotPart;ilong++) {
		i = (long) (invL * PartX[ilong]*NCells);
		j = (long) (invL * PartY[ilong]*NCells);
		k = (long) (invL * PartZ[ilong]*NCells);

		lin_ijk = k+j*NCells+i*NCells*NCells;
		NPartPerCell[lin_ijk]++;
	}
	//Alloc Enough Memory
	ListOfPart = (long **) calloc(NCells*NCells*NCells,sizeof(long *));
	for (i=0;i<NCells;i++){
	for (j=0;j<NCells;j++){
	for (k=0;k<NCells;k++){
		lin_ijk = k+j*NCells+i*NCells*NCells;
		ListOfPart[lin_ijk] = (long *) calloc(NPartPerCell[lin_ijk],sizeof(long *));
		MassLeft[lin_ijk] = NPartPerCell[lin_ijk]*mp; 
	}	
	}
	}

	for (ilong=0;ilong<NTotPart;ilong++) {
		i = (long) (invL * PartX[ilong]*NCells);
		j = (long) (invL * PartY[ilong]*NCells);
		k = (long) (invL * PartZ[ilong]*NCells);
	
		lin_ijk = k+j*NCells+i*NCells*NCells;
		ListOfPart[lin_ijk][count[lin_ijk]] = ilong;
		count[lin_ijk]++;
	}
//----------------------------------- Particles assigned to grid




//PLACING HALOES-------------------------------------------------
	//Walk the grid
	for (lin_ijk=0;lin_ijk<NtotCells;lin_ijk++){
		starthalo = FirstHaloInCell[lin_ijk];
		Nexcluded=0;
		if (lin_ijk==NtotCells-1)
			endhalo=NHalosTot;
		else
			endhalo = FirstHaloInCell[lin_ijk+1];

		#ifdef _VERB
		fprintf(stderr,"Node %ld, halo start %ld, halo end %ld\nPartciles in Cell: %ld\n",lin_ijk,starthalo,endhalo,NPartPerCell[lin_ijk]);
		#endif
		for(ihalo = starthalo; ihalo<endhalo; ihalo++){


			//Now we actually place them!
			do {	
			  rnd = (long) NPartPerCell[lin_ijk] * ((float)rand()/(RAND_MAX+1.0));	
			  ipart = ListOfPart[lin_ijk][rnd];

			} while (excluded[ipart]==1);


			HaloX[ihalo] = PartX[ipart];
			HaloY[ihalo] = PartY[ipart];
			HaloZ[ihalo] = PartZ[ipart];

			#ifdef _R200
			  R = R200_from_mass(HaloMass[ihalo]);
			  R2 = R*R;
			  for (j=0; j<NPartPerCell[lin_ijk];j++){
				jpart = ListOfPart[lin_ijk][j];
				if (excluded[jpart]==0){
					if (R2>((PartX[jpart]-HaloX[ihalo])*(PartX[jpart]-HaloX[ihalo])+(PartY[jpart]-HaloY[ihalo])*(PartY[jpart]-HaloY[ihalo])+(PartZ[jpart]-HaloZ[ihalo])*(PartZ[jpart]-HaloZ[ihalo]))){
						excluded[jpart]=1;
						Nexcluded++;
					} 
				}
			  }
			  if ((( (float) Nexcluded * 1.0/NPartPerCell[lin_ijk])>frac)){
				fprintf(stderr,"WARNING: Reached the  %f %% of the particles excluded for the mass %e for the node %ld\n",100*frac,HaloMass[ihalo],lin_ijk);
				break;
			  }


			#else
			excluded[ipart]=1;
			Nexcluded++;
			#endif

		}

		#ifdef _VERB
		fprintf(stderr,"Done with %ld (excluded %ld particles)\n",lin_ijk,Nexcluded);
		#endif

	}
//------------------------------------------------- HALOES PLACED



	#ifdef _VERB
		fprintf(stderr,"All haloes placed\n");
	#endif

	free(count); free(NPartPerCell); free(ListOfPart); free(excluded);
}





void place_halos_Mglobal(long NHalosTot, float *HaloMass, long Nlin, long NTotPart, float *PartX, float *PartY, float *PartZ, float L, float mp, float *HaloX, float *HaloY, float *HaloZ){

fprintf(stderr,"Hi! This is place_halos.c\n");
fprintf(stdout,"Hi! This is place_halos.c\n");


//Initiallising -------------------------------------------------
	long i,j,k,lin_ijk;
	long *count;
	long ihalo,ilong, ipart;
	double invL = 1./L;
	float R;	
	NCells = Nlin;
	lcell=L/NCells;
	Lbox = L;
	//long NTotCells = cube(NCells);

	//temp
	float *HaloR;	
	


	//Allocate memory for the arrays 
	excluded = (long *) calloc(NTotPart, sizeof(long));
	NPartPerCell = (long *) calloc(NCells*NCells*NCells,sizeof(long));
	MassLeft = (float *) calloc(NCells*NCells*NCells,sizeof(float));
	count = (long *) calloc(NCells*NCells*NCells,sizeof(long));
	Nexcluded = (long *) calloc(cube(NCells),sizeof(long)); 
	//do we need to allocate the output?
	HaloR = (float *) calloc(NHalosTot,sizeof(float));


	//Initiallise random numbers
	srand (time(NULL));

#ifdef _CRIT
	rho_ref = rho_crit;
#else
	rho_ref = mp/(L*L*L)*NTotPart;
#endif	



#ifdef _R200
	fprintf(stderr,"#def _R200\n");
#endif
#ifdef _CRIT
	fprintf(stderr,"#def _CRIT\n");
#endif
#ifdef _MUTUAL
	fprintf(stderr,"#def _MUTUAL\n");
#endif
#ifdef _VERB
	fprintf(stderr,"#def _VERB\n");
#endif 
#ifdef _DEBUG
	fprintf(stderr,"#def _DEBUG \n");
#endif
#ifdef _PERIODIC
	fprintf(stderr,"#def _PERIODIC\n");
#endif

	
	#ifdef _VERB
	fprintf(stderr,"\nMassFunction computed globally, Particles placed in %ld^3 cells\n",NCells);
	fprintf(stderr,"BOX = %f  lcell =%f   rho_ref = %e  invL %f\n",L,lcell,rho_ref,invL);
	fprintf(stderr,"RAND_MAX=%d\n",RAND_MAX);
	fprintf(stderr,"Nhalos = %ld NPart = %ld\n",NHalosTot, NTotPart);
	fprintf(stderr,"X[0] = %f Y[0] = %f Z[0] = %f\n",PartX[0],PartY[0],PartZ[0]);
	fprintf(stderr,"X[1] = %f Y[1] = %f Z[1] = %f\n",PartX[1],PartY[1],PartZ[1]);
	fprintf(stderr,"M[0] = %f \n",HaloMass[0]);
	fprintf(stderr,"M[1] = %f \n",HaloMass[1]);
	#endif
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
#endif
	//Alloc Enough Memory
	ListOfPart = (long **) calloc(NCells*NCells*NCells,sizeof(long *));
	for (i=0;i<NCells;i++){
	for (j=0;j<NCells;j++){
	for (k=0;k<NCells;k++){
		lin_ijk = k+j*NCells+i*NCells*NCells;
		ListOfPart[lin_ijk] = (long *) calloc(NPartPerCell[lin_ijk],sizeof(long));
		MassLeft[lin_ijk] = NPartPerCell[lin_ijk]*mp; 
#ifdef _DEBUG
		if (lin_ijk<10 || lin_ijk > (NCells*NCells*NCells) - 10)
			fprintf(stderr,"Allocated %ld in lin_ijk = %ld [%ld,%ld,%ld]\n",NPartPerCell[lin_ijk],lin_ijk,i,j,k);
#endif
	}	
	}
	}
#ifdef _DEBUG
	fprintf(stderr,"... memory allocated ...\n");
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
	fprintf(stderr," Mass Function\n");
	fprintf(stderr,"%e \n",HaloMass[0]);
	fprintf(stdout,"%e \n",HaloMass[1]);
	for (ihalo=0;ihalo<15;ihalo++){
		fprintf(stderr,"halo %ld: ",ihalo);
		fprintf(stderr,"M=%e\n",HaloMass[ihalo]);
	}
	fprintf(stderr," Placing Halos...\n\n");
#endif


	for (ihalo=0;ihalo<NHalosTot;ihalo++){
	

#ifdef _DEBUG
		fprintf(stderr,"\n- Halo %ld ",ihalo);
#endif
	
		//Select a Cell to place the halo (weighted with the mass)
		lin_ijk = select_cell(&i,&j,&k,MassLeft,NCells,HaloMass[ihalo]);
#ifdef _DEBUG
		fprintf(stderr,"ijk=%ld  ",lin_ijk);
		fprintf(stderr,", after: M=%e   \n",MassLeft[lin_ijk]);
#endif

		//Place the halo in one particle of that cell
		ipart = select_part(lin_ijk);		
               	HaloX[ihalo] = PartX[ipart];
               	HaloY[ihalo] = PartY[ipart];
               	HaloZ[ihalo] = PartZ[ipart];
		R=R200_from_mass(HaloMass[ihalo]);
		HaloR[ihalo]=R;
		#ifdef _MUTUAL
		while(check_HaloR(ihalo,HaloX,HaloY,HaloZ,HaloR)==0) {
			#ifdef _DEBUG
			fprintf(stderr,"Refused part : %ld. Excluded = %ld\n",ipart,excluded[ipart]);
			#endif
			ipart = select_part(lin_ijk);		
                	HaloX[ihalo] = PartX[ipart];
                	HaloY[ihalo] = PartY[ipart];
                	HaloZ[ihalo] = PartZ[ipart];
			R=R200_from_mass(HaloMass[ihalo]);
			HaloR[ihalo]=R;
		}
		#endif
		//Exclude particles in a R200 radius for future selections
#ifdef _DEBUG
		fprintf(stderr,"halo %ld assigned to particle %ld at [%f,%f,%f]. R= %f, M= %e\n",ihalo,ipart,HaloX[ihalo],HaloY[ihalo],HaloZ[ihalo],R,HaloMass[ihalo]);
#endif
		exclude(ipart,R,PartX,PartY,PartZ,i,j,k);            


	}
#ifdef _VERB
	fprintf(stderr,"\nEverything done!!!\n");
#endif
	free(count); free(NPartPerCell); free(ListOfPart); free(excluded);
}


void exclude(long ipart,float R,float *PartX,float *PartY,float *PartZ,long i,long j,long k) {

	float X = PartX[ipart];
	float Y = PartY[ipart];
        float Z = PartZ[ipart];
	float R2 = R*R;  
	long faces=0, edges=0,vertices=0;

		if (2.*R>lcell){
			fprintf(stderr,"ERROR: there are haloes with diameter D= 2R = %f greater than the size of the cell l = %f.\nPlease change the size of the cell\n",2.*R,lcell);
			exit(0);
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
	long lin_ijk = k+j*NCells+i*NCells*NCells;
	if (MassLeft[lin_ijk]==0.){
		#ifdef _DEBUG
		fprintf(stderr,"No mass left in cell [%ld,%ld,%ld]  M=%f \n",i,j,k,MassLeft[lin_ijk]);
		#endif
		return 0;
	}
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

	if (Nexcluded[lin_ijk]>frac*NPartPerCell[lin_ijk]){
		#ifdef _VERB
		fprintf(stderr,"WARNING: Cell %ld was completely excluded when there was still %e M_sun to be assigned\n",lin_ijk,MassLeft[lin_ijk]);
		fprintf(stderr,"WARNING: NPartPerCell %ld         Nexcluded %ld\n",Nexcluded[lin_ijk],NPartPerCell[lin_ijk]);
		#endif
		MassLeft[lin_ijk]=0.;
	}			
	return 0;
}

/*
long select_part(long *Part,long N,long *excl){
	long i_rnd,ipart;
         do {
		i_rnd = (long) (N * ((double)rand()/(RAND_MAX+1.0)));
                ipart = Part[i_rnd];
         } while (excl[ipart]==1);	
	return ipart;
}*/

long select_part(long ijk){
	long i_rnd,ipart;
         do {
		i_rnd = (long) (NPartPerCell[ijk] * ((double)rand()/(RAND_MAX+1.0)));
                ipart = ListOfPart[ijk][i_rnd];
         } while (excluded[ipart]==1);	
	return ipart;
}



long select_cell(long *x, long *y, long *z, float *MassArray, int N, float M) {
	long i,j,k,lin_ijk;
	float TotMass = 0., d_rand;
	float threshold; 

	#ifdef _ENOUGH
		threshold=M;
	#else
		threshold=0;
	#endif

	for (i=0;i<N;i++){
        for (j=0;j<N;j++){
        for (k=0;k<N;k++){
		lin_ijk = k+j*N+i*N*N;
		if ((MassArray)[lin_ijk]>threshold)
			TotMass += (MassArray)[lin_ijk];
	}
	}
	}
	d_rand = TotMass * ((double)rand()/(RAND_MAX+1.0));
#ifdef _DEBUG
			fprintf(stderr,"Rand: %e TotMass: %e  HaloMass: %e.    ",d_rand,TotMass,M);
#endif
	TotMass = 0;

	for (i=0;i<N;i++){
        for (j=0;j<N;j++){
        for (k=0;k<N;k++){
	lin_ijk = k+j*N+i*N*N;
	if ((MassArray)[lin_ijk]>threshold)
			TotMass += (MassArray)[lin_ijk];
#ifdef _DEBUG
//			fprintf(stderr,"%e, ",TotMass);
#endif
		if (TotMass > d_rand){
			*x=i;
			*y=j;
			*z=k;
#ifdef _DEBUG
			fprintf(stderr,"\nassigned to cell %ld=[%ld,%ld,%ld]\n Mass in Cell: Before M=%e    ",lin_ijk,i,j,k,(MassArray)[lin_ijk]);
#endif

			(MassArray)[lin_ijk]=(MassArray)[lin_ijk]-M;

#ifdef _DEBUG
			fprintf(stderr,"\n After M=%e    ",(MassArray)[lin_ijk]);
#endif
			if ((MassArray[lin_ijk])<0.){
				(MassArray)[lin_ijk]=0.;
				#ifdef _ENOUGH
				fprintf(stderr,"ERROR: getting negative masses\n");
				#endif
			}
			return lin_ijk;
		}
	}
	}
	}
	fprintf(stderr,"ERROR: no cell selected\n");
	return -1;		
}


int check_HaloR(long i,float *X, float *Y, float *Z , float *R){
	long j;
	for (j=0; j<i; j++){
		if ((square(X[i]-X[j])+square(Y[i]-Y[j])+square(Z[i]-Z[j]))<square(R[i]+R[j])){
//			#ifdef _DEBUG
			fprintf(stderr,"\nOverlap between previous halo %ld placed at [%f,%f,%f] with R=%f and the selected particle at [%f,%f,%f] with R=%f. ",j,X[j],Y[j],Z[j],R[j],X[i],Y[i],Z[i],R[i]);
//			#endif
			return 0;	
		}
	}
	return 1;
}
