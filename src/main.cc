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
	DataFile df("/home/nils/triforce/dev/dataConvex.dat",Binary);
	Data3D *dat;
	dat = df.digest3DBinaryTable();
	
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
	mol->addRealAtom(0,1.3,0.5,string("C1"));
	mol->addRealAtom(0,1.3,-0.5,string("C1"));
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
	
	Interpolation interpolator(dat);
	
	Tessellation tessellation(*mol);
	tessellation.build();
	string s("gbonnet.csv");
	tessellation.outputGaussBonnetData(s);
	
	Integrator integrator(&tessellation, &interpolator, mol);
	double area = integrator.integrate();
	
	printf("total area of molecule: %f\n",area);
	
	
}
