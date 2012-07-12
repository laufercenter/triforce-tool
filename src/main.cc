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
	Data3D *a;
	a = df.digestBinaryTable();
	
	DataFile df2("/home/nils/seawater/nonpolar/fluorene.top",MapCSV);
	Topology *b;
	b = df2.digestTOP();
	
	//b->print();
	
	DataFile df3("/home/nils/seawater/nonpolar/fluorene.gro",MapCSV);
	Molecule *c;
	c = df3.digestGRO(*b);
	
	//c->print();
	
	Interpolation *interpolator = new Interpolation(a);
	
	int res=100;
	double x[res];
	double y[res];
	Vector vec(3);
	for(int i=0;i<res;++i){
		double v = M_PI * (double)i/(double)res;
		vec(0)=0.6;
		vec(1)=v;
		vec(2)=0.5;
		y[i]=interpolator->multiPointTaylor(vec);
		x[i]=v;
		printf("%f %f\n",x[i],y[i]);
	}
	
	
}
