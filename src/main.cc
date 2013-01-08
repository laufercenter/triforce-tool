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
	DataFile df0("/home/nils/triforce/dev/dataConcave.dat",Binary);
	DataFile df1("/home/nils/triforce/dev/dataConcave.dat",Binary);
	DataFile df2("/home/nils/triforce/dev/dataConcave.dat",Binary);
	DataFile df3("/home/nils/triforce/dev/dataConcave.dat",Binary);
	Data3D* dat0 = df0.digest3DBinaryTable();
	Data3D* dat1 = df1.digest3DBinaryTable();
	Data3D* dat2 = df2.digest3DBinaryTable();
	Data3D* dat3 = df3.digest3DBinaryTable();
	
	Surface3D* surf0 = new Surface3D(dat0);
	Surface3D* surf1 = new Surface3D(dat1);
	Surface3D* surf2 = new Surface3D(dat2);
	Surface3D* surf3 = new Surface3D(dat3);
	
	
	
	//dat1->print(); 
	
	DataFile dftop("/home/nils/seawater/nonpolar/fluorene.top",MapCSV);
	Topology *top;
	top = dftop.digestTOP();
	
	//b->print();
	
	DataFile dfgro("/home/nils/seawater/nonpolar/fluorene.gro",MapCSV);
	Molecule *mol, *mol2;
	mol = dfgro.digestGRO(*top);
	mol2 = dfgro.digestGRO(*top);
	//mol = new Molecule(*top);
	//mol2=mol;
	
	/*
	//testcase 0 (single occlusion)
	mol->addRealAtom(0,0,0,string("C1"));
	mol->addRealAtom(0,4.2,2.15,string("C1"));
	mol->addRealAtom(0,4.2,-2.15,string("C1"));
	//mol->addRealAtom(0,4.2,0,string("C1"));
	mol->addRealAtom(-0.6,2,0,6.8,1);
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
	
	Interpolation interpolator0(surf0);
	Interpolation interpolator1(surf1);
	Interpolation interpolator2(surf2);
	Interpolation interpolator3(surf3);
	
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
	
	mol->print();
	
	IntegratorTriforce integrator(&interpolator0, &interpolator1, &interpolator2, &interpolator3);
	double area = integrator.integrate(mol, &tessellation);
	
	//IntegratorGaussBonnet integrator3;
	//double area3 = integrator3.integrate(mol, &tessellation2);
	
	IntegratorNumerical integrator4(30);
	double area4 = integrator4.integrate(mol2);
	
	mol->print();
	mol2->print();
	
	
	
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
