#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <string>

#include "triforce-0.39/triforce.h"



void readPDB(string name){
	
}


main()
{ 
	//printf("it has started....\n");
	DataFile df0("/home/nils/triforce/dev/dataConvex.dat",Binary);
	DataFile df1("/home/nils/triforce/dev/dataConcave.dat",Binary);
	
	DataFile df0forces0("/home/nils/triforce/dev/dataConvex0.dat",Binary);
	DataFile df0forces1("/home/nils/triforce/dev/dataConvex1.dat",Binary);
	DataFile df0forces2("/home/nils/triforce/dev/dataConvex2.dat",Binary);
	
	DataFile df1forces0("/home/nils/triforce/dev/dataConcave0.dat",Binary);
	DataFile df1forces1("/home/nils/triforce/dev/dataConcave1.dat",Binary);
	DataFile df1forces2("/home/nils/triforce/dev/dataConcave2.dat",Binary);

	
	Data3D *dat0;
	Data3D *dat1;
	
	Data3D *dat0forces0;
	Data3D *dat0forces1;
	Data3D *dat0forces2;

	Data3D *dat1forces0;
	Data3D *dat1forces1;
	Data3D *dat1forces2;
	
	dat0 = df0.digest3DBinaryTable();
	dat1 = df1.digest3DBinaryTable();
	
	dat0forces0 = df0forces0.digest3DBinaryTable();
	dat0forces1 = df0forces1.digest3DBinaryTable();
	dat0forces2 = df0forces2.digest3DBinaryTable();
	
	dat1forces0 = df1forces0.digest3DBinaryTable();
	dat1forces1 = df1forces1.digest3DBinaryTable();
	dat1forces2 = df1forces2.digest3DBinaryTable();
	
	Surface3D* surf0 = new Surface3D(dat0);
	Surface3D* surf1 = new Surface3D(dat1);
	
	Surface3D* surf0forces0 = new Surface3D(dat0forces0);
	Surface3D* surf0forces1 = new Surface3D(dat0forces1);
	Surface3D* surf0forces2 = new Surface3D(dat0forces2);

	Surface3D* surf1forces0 = new Surface3D(dat1forces0);
	Surface3D* surf1forces1 = new Surface3D(dat1forces1);
	Surface3D* surf1forces2 = new Surface3D(dat1forces2);
	
	dat2 = df4.digest3DBinaryTable();
	
	//dat1->print(); 
	
	DataFile df2("/home/nils/seawater/nonpolar/fluorene.top",MapCSV);
	Topology *top;
	top = df2.digestTOP();
	
	//b->print();
	
	DataFile df3("/home/nils/seawater/nonpolar/fluorene.gro",MapCSV);
	Molecule *mol;
	mol = df3.digestGRO(*top);
	//mol = new Molecule(*top);

	/*
	//testcase 0 (single occlusion)
	mol->addRealAtom(0,-1.5,0,string("C1"));
	mol->addRealAtom(-1,1.2,1.15,string("C1"));
	mol->addRealAtom(-1,1.2,-1.15,string("C1"));
	mol->addRealAtom(0.8,0,0,string("C1"));
	*/
	
	/*
	//testcase 1 (outside / inside)
	mol->addRealAtom(0,-1.5,0,string("C1"));
	mol->addRealAtom(0,1.3,0.5,string("C1"));
	mol->addRealAtom(0,1.3,-0.5,string("C1"));
	mol->addRealAtom(0,1.3,-1.5,string("C1"));
	*/
	
	/*
	//testcase 2 (both)
	mol->addRealAtom(0,-1.5,0,string("C1"));
	mol->addRealAtom(0,1.3,0.75,string("C1"));
	mol->addRealAtom(0,1.3,-0.75,string("C1"));
	mol->addRealAtom(0,1.2,0,string("C1"));
	*/

	/*
	//testcase 3 (totally inside)
	mol->addRealAtom(0,-1.5,0,string("C1"));
	mol->addRealAtom(0,1.3,1,string("C1"));
	mol->addRealAtom(0,1.3,-1,string("C1"));
	mol->addRealAtom(0,1.85,0,string("C1"));
	*/
	
	
	//testcase 4 (singles)
	//mol->addRealAtom(0,0,0,string("C1"));
	//mol->addRealAtom(1.5,1.3,1,string("C1"));
	//mol->addRealAtom(1.5,1.3,-1,string("C1"));
	//mol->addRealAtom(0,2,0,string("C1"));
	//mol->addRealAtom(2,0,2,string("C1"));
	
	
	/*
	//testcase 4 (triforce)
	mol->addRealAtom(0,-1.5,0,string("C1"));
	mol->addRealAtom(1.2,1.3,1,string("C1"));
	mol->addRealAtom(1.2,1.3,-1,string("C1"));
	mol->addRealAtom(-1.2,1.3,0,string("C1"));
	*/
	
	mol->print();
	
	//c->print();
	
	Interpolation interpolator0(dat0);
	Interpolation interpolator1(dat1);
	
	Interpolation interpolator0forces0(dat0forces0);
	Interpolation interpolator0forces1(dat0forces1);
	Interpolation interpolator0forces2(dat0forces2);
	
	Interpolation interpolator1forces0(dat1forces0);
	Interpolation interpolator1forces1(dat1forces1);
	Interpolation interpolator1forces2(dat1forces2);
	
/*	
	for(int i=0; i<10; ++i){
		Vector v(3);
		double r;
		v(0) = i*1.7/1000.0;
		v(1) = 1.1;
		v(2) = 0.75;
		
		r = interpolator0.interpolate(v);
		
		printf("%f %f\n",v(2),r);
		
	}
	exit(0);         
*/	
	
	//Tessellation tessellation(*mol, &interpolator2);
	Tessellation tessellation(*mol);
	tessellation.build(true);
	
	IntegratorTriforce integrator(&interpolator0, &interpolator0forces0, &interpolator0forces1, &interpolator0forces2, &interpolator1, &interpolator0forces0, &interpolator1forces1, &interpolator2forces2);
	double area = integrator.integrate(mol, &tessellation);
	
	//IntegratorGaussBonnet integrator3;
	//double area3 = integrator3.integrate(mol, &tessellation2);
	
	IntegratorNumerical integrator4(100);
	double area4 = integrator4.integrate(mol, &tessellation);
	
	
	
	printf("total area of molecule: %f %f\n",area,area4);
	
	
	/*
	
	Tessellation tessellation2(*mol);
	tessellation2.build(false);
	
	IntegratorTriforce integrator5(&interpolator0, &interpolator1);
	double area5 = integrator5.integrate(mol, &tessellation2);
	
	
	printf("total area of molecule: %f %f\n",area,area5);
	
	
	IntegratorStatistical integrator2(1000);
	double area2 = integrator2.integrate(mol, &tessellation2);

	IntegratorGaussBonnet integrator3;
	double area3 = integrator3.integrate(mol, &tessellation2);
	
	IntegratorNumerical integrator4(1000);
	double area4 = integrator4.integrate(mol, &tessellation2);
	
	printf("total area of molecule: %f, %f, %f, %f\n",area,area2,area3,area4);
	*/
	
}
