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
	Data3D *dat0;
	Data3D *dat1;
	dat0 = df0.digest3DBinaryTable();
	dat1 = df1.digest3DBinaryTable();
	
	dat1->print();
	
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
	mol->addRealAtom(0,1.2,1.15,string("C1"));
	mol->addRealAtom(0,1.2,-1.15,string("C1"));
	mol->addRealAtom(1,1.3,-0.5,string("C1"));
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
	
	/*
	//testcase 4 (singles)
	mol->addRealAtom(0,-1.5,0,string("C1"));
	mol->addRealAtom(1.5,1.3,1,string("C1"));
	mol->addRealAtom(1.5,1.3,-1,string("C1"));
	mol->addRealAtom(-1.5,1.3,0,string("C1"));
	*/
	
	/*
	//testcase 4 (triforce)
	mol->addRealAtom(0,-1.5,0,string("C1"));
	mol->addRealAtom(1.2,1.3,1,string("C1"));
	mol->addRealAtom(1.2,1.3,-1,string("C1"));
	mol->addRealAtom(-1.2,1.3,0,string("C1"));
	*/
	
	//mol->print();
	
	//c->print();
	
	Interpolation interpolator0(dat0);
	Interpolation interpolator1(dat1);
	
	Tessellation tessellation(*mol);
	tessellation.build(true);
	
	IntegratorTriforce integrator(&interpolator0, &interpolator1);
	double area = integrator.integrate(mol, &tessellation);
	
	
	printf("total area of molecule: %f\n",area);
	
	/*
	Tessellation tessellation2(*mol);
	tessellation2.build(false);
	
	IntegratorStatistical integrator2(10000);
	double area2 = integrator2.integrate(mol, &tessellation2);

	IntegratorGaussBonnet integrator3;
	double area3 = integrator3.integrate(mol, &tessellation2);
	
	IntegratorNumerical integrator4(1000);
	double area4 = integrator4.integrate(mol, &tessellation2);
	
	printf("total area of molecule: %f, %f, %f, %f\n",area,area2,area3,area4);
	*/
	
}
