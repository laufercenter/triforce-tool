#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <string>

#include "triforce-0.39/triforce.h"


#include "bfgs.h"




DataFile *df0;
DataFile *df1;
DataFile *df2;
DataFile *df3;
Data3D *dat0;
Data3D *dat1;
Data3D *dat2;
Data3D *dat3;

Surface3D* surf0;
Surface3D* surf1;
Surface3D* surf2;
Surface3D* surf3;

Interpolation *interpolator0;
Interpolation *interpolator1;
Interpolation *interpolator2;
Interpolation *interpolator3;



Molecule *mol;
IntegratorTriforce *integrator;
IntegratorNumerical *integratorNumerical;
vector<Vector> coordinates;
vector<double> radii;
vector<vector<double*> > forces;




void triforce_report(int step, double funct_value, Vector x){
	vector<Vector>  coordinates = mol->fetchCoordinates();
	printf("%d %f",step, funct_value);
	//outputting coordinates
	for(int i=0; i<coordinates.size(); ++i){
		printf(" %f",radii[i]);
		for(int j=0; j<3; ++j){
			printf(" %f",x(i*3+j));
		}
	}
	printf("\n");
}
	
	


double triforce_func(Vector m){
	Vector v(3);
	double area;
	vector<Vector>  coordinates = mol->fetchCoordinates();
	for(int i=0; i<coordinates.size(); ++i){
		for(int j=0; j<3; ++j){
			v(j) = m(i*3+j);
		}
		mol->setInternallyStoredAtomCoordinates(i,v);
		
	}
	
	
	Tessellation *t = new Tessellation(*mol);
	t->build(true);
	area = integrator->integrate(mol, t);
	delete(t);
	

	return area;
}

Vector triforce_derivative (Vector unused){
	Vector m;
	forces = mol->fetchForcePointers();
	m = Vector(forces.size()*3);
	
	for(int i=0; i<forces.size(); ++i){
		for(int j=0; j<3; ++j){
			m(i*3+j) = *(forces[i][j]);
		}
	}
	
	return m;
	
	
}


main()
{ 
	
	
	df0 = new DataFile("/home/nils/triforce/dev/dataConcave.dat");
	df1 = new DataFile("/home/nils/triforce/dev/dataConcave0.dat");
	df2 = new DataFile("/home/nils/triforce/dev/dataConcave1.dat");
	df3 = new DataFile("/home/nils/triforce/dev/dataConcave2.dat");
	dat0 = df0->digest3DBinaryTable();
	dat1 = df1->digest3DBinaryTable();
	dat2 = df2->digest3DBinaryTable();
	dat3 = df3->digest3DBinaryTable();
	
	surf0 = new Surface3D(dat0);
	surf1 = new Surface3D(dat1);
	surf2 = new Surface3D(dat2);
	surf3 = new Surface3D(dat3);
	
	interpolator0 = new Interpolation(surf0,TAYLOR_QUADRATIC);
	interpolator1 = new Interpolation(surf1,TAYLOR_QUADRATIC);
	interpolator2 = new Interpolation(surf2,TAYLOR_QUADRATIC);
	interpolator3 = new Interpolation(surf3,TAYLOR_QUADRATIC);
	
	
	
	
	
	
	
	
	
	/*
	//printf("it has started....\n");
	DataFile df0("/home/nils/triforce/dev/dataConcave.dat",Binary);
	DataFile df1("/home/nils/triforce/dev/dataConcave0.dat",Binary);
	DataFile df2("/home/nils/triforce/dev/dataConcave1.dat",Binary);
	DataFile df3("/home/nils/triforce/dev/dataConcave2.dat",Binary);
	Data3D* dat0 = df0.digest3DBinaryTable();
	Data3D* dat1 = df1.digest3DBinaryTable();
	Data3D* dat2 = df2.digest3DBinaryTable();
	Data3D* dat3 = df3.digest3DBinaryTable();
	
	Surface3D* surf0 = new Surface3D(dat0);
	Surface3D* surf1 = new Surface3D(dat1);
	Surface3D* surf2 = new Surface3D(dat2);
	Surface3D* surf3 = new Surface3D(dat3);
	*/
	
/*	
	Interpolation interpolator0(surf0,TAYLOR_LINEAR);
	Interpolation interpolator1(surf1,TAYLOR_LINEAR);
	Interpolation interpolator2(surf2,TAYLOR_LINEAR);
	Interpolation interpolator3(surf3,TAYLOR_LINEAR);
*/	

/*

	Interpolation interpolator0(surf0,TAYLOR_QUADRATIC);
	Interpolation interpolator1(surf1,TAYLOR_QUADRATIC);
	Interpolation interpolator2(surf2,TAYLOR_QUADRATIC);
	Interpolation interpolator3(surf3,TAYLOR_QUADRATIC);
	
	
	
	
	DataFile dftop("/home/nils/seawater/nonpolar/fluorene.top",MapCSV);
	Topology *top;
	top = dftop.digestTOP();
	
	
	DataFile dfgro("/home/nils/seawater/nonpolar/fluorene.gro",MapCSV);
	Molecule *mol_triforce, *mol_numerical, *mol_fd;
	*/
	/*
	mol_triforce = dfgro.digestGRO(*top);
	mol_numerical = dfgro.digestGRO(*top);
	mol_fd = dfgro.digestGRO(*top);
	
	Vector ppt(3);
	ppt.zeros();
	ppt(2) = 0.1;
	
	//mol_triforce->perturbInternallyStoredAtomCoordinates(21, ppt);
	//mol_numerical->perturbInternallyStoredAtomCoordinates(21, ppt);
	*/
	
	/*
	
	mol_triforce = new Molecule(*top);
	mol_numerical = new Molecule(*top);
	mol_fd = new Molecule(*top);
	
	
	//testcase 0 (single occlusion)
	mol_triforce->addInternallyStoredAtom(0,0,0,string("C1"));
	mol_triforce->addInternallyStoredAtom(1.5,4.2,2.15,string("C1"));
	//mol_triforce->addInternallyStoredAtom(-0.5,4.2,-2.15,string("C1"));
	//mol_triforce->addInternallyStoredAtom(0,4.2,0,string("C1"));
	
	mol_numerical->addInternallyStoredAtom(0,0,0,string("C1"));
	mol_numerical->addInternallyStoredAtom(1.5,4.2,2.15,string("C1"));
	//mol_numerical->addInternallyStoredAtom(-0.5,4.2,-2.15,string("C1"));
	//mol_numerical->addInternallyStoredAtom(0,4.2,0,string("C1"));

	mol_fd->addInternallyStoredAtom(0,0,0,string("C1"));
	mol_fd->addInternallyStoredAtom(1.5,4.2,2.15,string("C1"));
	//mol_fd->addInternallyStoredAtom(-0.5,4.2,-2.15,string("C1"));
	//mol_fd->addInternallyStoredAtom(0,4.2,0,string("C1"));
	
	//mol->addInternallyStoredAtom(-0.6,2,0,6.8,1);
	
	*/
	
	/*
	//testcase 0 (single occlusion)
	mol_triforce->addInternallyStoredAtom(0,0,0,string("C1"));
	mol_triforce->addInternallyStoredAtom(2.15,4.2,0,string("C1"));
	mol_triforce->addInternallyStoredAtom(-2.15,4.2,0,string("C1"));
	//mol_triforce->addInternallyStoredAtom(0,4.2,0,string("C1"));
	
	mol_numerical->addInternallyStoredAtom(0,0,0,string("C1"));
	mol_numerical->addInternallyStoredAtom(2.15,4.2,0,string("C1"));
	mol_numerical->addInternallyStoredAtom(-2.15,4.2,0,string("C1"));
	//mol_numerical->addInternallyStoredAtom(0,4.2,0,string("C1"));

	mol_fd->addInternallyStoredAtom(0,0,0,string("C1"));
	mol_fd->addInternallyStoredAtom(2.15,4.2,0,string("C1"));
	mol_fd->addInternallyStoredAtom(-2.15,4.2,0,string("C1"));
	//mol_fd->addInternallyStoredAtom(0,4.2,0,string("C1"));
	*/
	
	
	/*
	//testcase 1 (outside / inside)
	mol->addInternallyStoredAtom(0,-1.5,0,string("C1"));
	mol->addInternallyStoredAtom(0,1.3,0.5,string("C1"));
	mol->addInternallyStoredAtom(0,1.3,-0.5,string("C1"));
	mol->addInternallyStoredAtom(0,1.3,-1.5,string("C1"));
	*/
	
	/*
	//testcase 2 (both)
	mol->addInternallyStoredAtom(0,-1.5,0,string("C1"));
	mol->addInternallyStoredAtom(0,1.3,0.75,string("C1"));
	mol->addInternallyStoredAtom(0,1.3,-0.75,string("C1"));
	mol->addInternallyStoredAtom(0,1.2,0,string("C1"));
	*/

	/*
	//testcase 3 (totally inside)
	mol->addInternallyStoredAtom(0,-1.5,0,string("C1"));
	mol->addInternallyStoredAtom(0,1.3,1,string("C1"));
	mol->addInternallyStoredAtom(0,1.3,-1,string("C1"));
	mol->addInternallyStoredAtom(0,1.85,0,string("C1"));
	*/
	
	
	//testcase 4 (singles)
	//mol->addInternallyStoredAtom(0,0,0,string("C1"));
	//mol->addInternallyStoredAtom(1.5,1.3,1,string("C1"));
	//mol->addInternallyStoredAtom(1.5,1.3,-1,string("C1"));
	//mol->addInternallyStoredAtom(0,2,0,string("C1"));
	//mol->addInternallyStoredAtom(2,0,2,string("C1"));
	
	
	/*
	//testcase 4 (triforce)
	mol->addInternallyStoredAtom(0,-1.5,0,string("C1"));
	mol->addInternallyStoredAtom(1.2,1.3,1,string("C1"));
	mol->addInternallyStoredAtom(1.2,1.3,-1,string("C1"));
	mol->addInternallyStoredAtom(-1.2,1.3,0,string("C1"));
	*/
	
	
	
	
	
	
	/*
	//DISCONTINUITY ******************************************************************************
	vector<double*> areas;
	vector<vector<double*> > forces;
	double l;
	int detail;
	Vector coordinate0(3),coordinate1(3);
	Tessellation *t;
	double area;
	double cl;
	IntegratorTriforce integrator(&interpolator0, &interpolator1, &interpolator2, &interpolator3);

	mol_triforce = new Molecule(*top);
	
	l = 4.2;
	detail = 1000;
	
	mol_triforce->addInternallyStoredAtom(-l/2.0,0,0,1.2,0);
	mol_triforce->addInternallyStoredAtom(l/2.0,0,0,1.2,0);
	
	areas = mol_triforce->fetchAreaPointers();
	forces = mol_triforce->fetchForcePointers();

	
	coordinate0 = Vector(3).zeros();
	coordinate1 = Vector(3).zeros();

	
	for(int i=0; i < detail; ++i){
		cl = l-(l*i/(detail-1));
		
		coordinate0(0) = -cl/2.0;
		coordinate1(0) = cl/2.0;
		
		
		mol_triforce->setInternallyStoredAtomCoordinates(0, coordinate0);
		mol_triforce->setInternallyStoredAtomCoordinates(1, coordinate1);
		t = new Tessellation(*mol_triforce);
		t->build(true);
		area = integrator.integrate(mol_triforce, t);
		delete(t);
		
		printf("%d %f %f %f %f %f %f %f %f %f %f\n",i,cl,area,*(areas[0]),*(areas[1]), *(forces[0][0]), *(forces[0][1]), *(forces[0][2]), *(forces[1][0]), *(forces[1][1]), *(forces[1][2]));
	}
	
	
	
	
	exit(0);
	//DISCONTINUITY ******************************************************************************
	*/
	
	
	
	/*
	//MINIMISATION 0 *****************************************************************************
	vector<double*> areas;
	vector<vector<double*> > forces;
	vector<Vector> coordinates;
	int steps;
	Vector p(3);
	Tessellation *t;
	double area, area2;
	double steplength;
	IntegratorTriforce integrator(&interpolator0, &interpolator1, &interpolator2, &interpolator3);
	IntegratorNumerical integratorNumerical(100);

	mol_triforce = new Molecule(*top);
	
	steps = 50;
	steplength = 0.01;
	
	mol_triforce->addInternallyStoredAtom(0,0,0,2.8,0);
	mol_triforce->addInternallyStoredAtom(3.9,0,0,1.2,0);
	mol_triforce->addInternallyStoredAtom(0,3.9,0,1.2,0);
	
	areas = mol_triforce->fetchAreaPointers();
	forces = mol_triforce->fetchForcePointers();
	coordinates = mol_triforce->fetchCoordinates();
	
	
	p = Vector(3).zeros();

	
	for(int i=0; i < steps; ++i){
		fprintf(stderr,"STEP: %d\n",i);
		
		t = new Tessellation(*mol_triforce);
		t->build(true);
		area = integrator.integrate(mol_triforce, t);
		delete(t);
		
		fprintf(stderr,"triforce\n");
		mol_triforce->print();
		
		
		area2 = integratorNumerical.integrate(mol_triforce);
		
		fprintf(stderr,"numerical\n");
		mol_triforce->print();
		
		
		coordinates = mol_triforce->fetchCoordinates();
		//printf("%d %f %f %f %f %f %f %f %f %f %f %f %f %f\n",i,area,*(areas[0]),*(areas[1]),*(areas[2]),coordinates[0](0),coordinates[0](1),coordinates[0](2),coordinates[1](0),coordinates[1](1),coordinates[1](2),coordinates[2](0),coordinates[2](1),coordinates[2](2));
		//printf("%d %f %f\n",i,area,area2);
		printf("%d %f %f 0 0 %f %f %f %f %f %f %f %f %f\n",i,area,area2,coordinates[0](0),coordinates[0](1),coordinates[0](2),coordinates[1](0),coordinates[1](1),coordinates[1](2),coordinates[2](0),coordinates[2](1),coordinates[2](2));
		
		
		
		
		
		
		for(int j=0; j < forces.size(); ++j){
			for(int k=0; k<3; ++k)
				p(k) = -*(forces[j][k]) * steplength;
			
			mol_triforce->perturbInternallyStoredAtomCoordinates(j, p);
		}
		
		
	}
	//MINIMISATION 0 *****************************************************************************
	*/
	
	
	//MINIMISATION 1 *****************************************************************************
	
	TriforceInterface trii(string("/home/nils/triforce/dev"));
	
	Molecule *mol_triforce, *mol_numerical, *mol_fd;
	
	
	vector<double*> areas;
	vector<vector<double*> > forces;
	vector<Vector> coordinates;
	vector<double> radii;
	int steps;
	Vector p(3);
	Tessellation *t;
	double area, area2;
	double steplength;
	
	DataFile dftop(string("/home/nils/seawater/nonpolar/fluorene.top"));
	Topology *top = dftop.digestTOP();
	
	
	DataFile dfgro(string("/home/nils/seawater/nonpolar/fluorene.gro"));

	//integrator = new IntegratorTriforce(interpolator0, interpolator1, interpolator2, interpolator3);
	integratorNumerical = new IntegratorNumerical(20);
	
	mol_triforce = dfgro.digestGRO(*top);

	
	steps = 300;
	steplength = 0.01;
	
	areas = mol_triforce->fetchAreaPointers();
	forces = mol_triforce->fetchForcePointers();
	coordinates = mol_triforce->fetchCoordinates();
	radii = mol_triforce->fetchRadii();
	
	p = Vector(3).zeros();

	
	for(int i=0; i < steps; ++i){
		fprintf(stderr,"STEP: %d\n",i);

		area = trii.calculateSurfaceArea(*mol_triforce);
		
		coordinates = mol_triforce->fetchCoordinates();
		


		
		fprintf(stderr,"triforce\n");
		
		mol_triforce->print();
		
		area2 = integratorNumerical->integrate(mol_triforce);
		
		fprintf(stderr,"numerical\n");
		
		mol_triforce->print();
		
		printf("%d %f %f",i,area,area2);
		
		for(int j=0; j < forces.size(); ++j){
			printf(" %f %f %f %f",radii[j], coordinates[j](0),coordinates[j](1),coordinates[j](2));
		}
		printf("\n");
		//if(abs(area-area2) > 10) exit(-1);
		
		
		for(int j=0; j < forces.size(); ++j){
			for(int k=0; k<3; ++k)
				p(k) = -*(forces[j][k]) * steplength;
			
			mol_triforce->perturbInternallyStoredAtomCoordinates(j, p);
		}
		
		
	}
	//MINIMISATION 1 *****************************************************************************
	
	
	/*
	//MINIMISATION 2 *****************************************************************************
	DataFile dftop("/home/nils/seawater/nonpolar/fluorene.top",MapCSV);
	Topology *top = dftop.digestTOP();
	
	
	DataFile dfgro("/home/nils/seawater/nonpolar/fluorene.gro",MapCSV);

	mol = dfgro.digestGRO(*top);
	integrator = new IntegratorTriforce(interpolator0, interpolator1, interpolator2, interpolator3);
	
	
	
	coordinates = mol->fetchCoordinates();
	radii = mol->fetchRadii();
	Vector m(coordinates.size()*3);
	
	for(int i=0; i<coordinates.size(); ++i){
		for(int j=0; j<3; ++j){
			m(i*3+j) = coordinates[i](j);
		}
	}
	
	BFGS bfgs;
	
	bfgs.optimise(triforce_func, triforce_derivative, triforce_report, m, 0.000001, 100);
	
	
	
	
	
	
	
	//MINIMISATION 2 *****************************************************************************
	*/
	
	
	
	
	
/*
	//Semi-analytical area and derivatives
	Tessellation tessellation(*mol_triforce);
	tessellation.build(true);
	
	IntegratorTriforce integrator(&interpolator0, &interpolator1, &interpolator2, &interpolator3);
	double area = integrator.integrate(mol_triforce, &tessellation);

	//Numerical area
	IntegratorNumerical integratorNumerical(100);
	double area4 = integratorNumerical.integrate(mol_numerical);
	
	printf("total area of molecule: %f %f\n",area,area4);
	
	mol_triforce->print();
	mol_numerical->print();
	
	
	
	exit(1);
	
	//finite difference derivatives
	vector<Vector> fd_forces;
	vector<Vector> atoms;
	//Tessellation *t;
	atoms = mol_fd->fetchCoordinates();
	Vector coordinate(3);
	Vector force(3);
	Vector pp(3),pn(3);
	double areap,arean;

	for(int i=0; i<atoms.size(); ++i){
		fd_forces.push_back(Vector(3).zeros());
	}
	
#define FD3 0.1
	
	
	for(int i=0; i<atoms.size(); ++i){
		fprintf(stdout,"STEP: %d\n",i);
		coordinate = mol_fd->getInternallyStoredAtomCoordinates(i);
		for(int j=0; j<3; ++j){
			pp = Vector(3).zeros();
			pn = Vector(3).zeros();
			pp(j) = FD3;
			pn(j) = -FD3;
			
			
			mol_fd->perturbInternallyStoredAtomCoordinates(i, pp);
			t = new Tessellation(*mol_fd);
			t->build(true);
			areap = integrator.integrate(mol_fd, t);
			mol_fd->setInternallyStoredAtomCoordinates(i, coordinate);
			delete(t);
			
			mol_fd->perturbInternallyStoredAtomCoordinates(i, pn);
			t = new Tessellation(*mol_fd);
			t->build(true);
			arean = integrator.integrate(mol_fd, t);
			mol_fd->setInternallyStoredAtomCoordinates(i, coordinate);
			delete(t);
			
			force(j) = (areap - arean)/(2*FD3);
			fprintf(stdout,"AREAS: %f %f\n",areap,arean);
		}
		fd_forces[i] = force;
	}
	
	*/
	
/*	
	for(int i=0; i<atoms.size(); ++i){
		fprintf(stdout,"STEP: %d\n",i);
		coordinate = mol_fd->getInternallyStoredAtomCoordinates(i);
		for(int j=0; j<3; ++j){
			pp = Vector(3).zeros();
			pn = Vector(3).zeros();
			pp(j) = FD3;
			pn(j) = -FD3;
			
			
			mol_fd->perturbInternallyStoredAtomCoordinates(i, pp);
			areap = integratorNumerical.integrate(mol_fd);
			mol_fd->setInternallyStoredAtomCoordinates(i, coordinate);
			
			mol_fd->perturbInternallyStoredAtomCoordinates(i, pn);
			arean = integratorNumerical.integrate(mol_fd);
			mol_fd->setInternallyStoredAtomCoordinates(i, coordinate);
			
			force(j) = (areap - arean)/(2*FD3);
			fprintf(stdout,"AREAS: %f %f\n",areap,arean);
		}
		fd_forces[i] = force;
	}
*/	
	/*
	for(int i=0; i<atoms.size(); ++i){
		fprintf(stdout,"FD FORCES [%d]: %f %f %f\n",i,fd_forces[i](0),fd_forces[i](1),fd_forces[i](2));
	}
	
	mol_triforce->print();
	mol_numerical->print();
	*/
	
	
}
