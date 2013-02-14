#ifdef HAVE_CONFIG_H
#include <config.h>
#endif


#include <string>
#include <boost/program_options.hpp>
#include "triforce-0.39/triforce.h"


namespace po = boost::program_options;



Molecule *mol, *mol_numerical;
string libraryPath, top, gro;
bool error;
double area, area_numerical;
int numericalDetail;
double numericalDifference;
int compare;



void progressbar(double percentage){
	int size=80;
	
	cout << "\r" "[";
	for(int i=0; i<size; ++i){
		if((double)i/(double)(size-1) <= percentage)
			cout << "=";
		else
			cout << " ";
	}
	cout << "] ";
	cout << ((int)ceil(percentage*100));
	cout << "%";
	cout << std::flush;
	
	if(percentage>=1) cout << endl;
	
}


int main(int ac, char* av[])
{ 
	try{
		po::options_description desc("Allowed options");
		desc.add_options()
			("help,h", "help")
			("library_path,l", po::value<string>(&libraryPath), "set library path")
			("topology,t", po::value<string>(&top), "set topology file")
			("structure,s", po::value<string>(&gro), "set structure file")
			("numerical_detail,n", po::value<int>(&numericalDetail)->default_value(0), "numerical detail (default off)")
			("numerical_difference,d", po::value<double>(&numericalDifference)->default_value(0), "numerical difference (default off)")
			("compare,c", po::value<int>(&compare)->default_value(0), "show only semi-analytical/numerical difference (default off)");
			
		po::variables_map vm;
		po::store(po::parse_command_line(ac,av,desc), vm);
		po::notify(vm);

		
		
		if(vm.count("help")){
			exit(0);
		}
		
		
		error=false;
		if(!vm.count("library_path")){
			fprintf(stderr, "library path not set\n");
			error=true;
		}
		if(!vm.count("topology")){
			fprintf(stderr, "topology file not set\n");
			error=true;
		}
		if(!vm.count("structure")){
			fprintf(stderr, "structure file not set\n");
			error=true;
		}
		
		if(error){
			cout << desc << "\n";
			
			exit(-1);
		}
	}
	catch(exception& e) {
		cerr << "error: " << e.what() << "\n";
		exit(-1);
	}
	
	
	
	
	DataFile dftop(top);
	Topology *top = dftop.digestTOP();
	DataFile dfgro(gro);
	Molecule *mol;
	mol = dfgro.digestGRO(*top);
	

	TriforceInterface trii(libraryPath);
	double area = trii.calculateSurfaceArea(*mol);
	
	if(compare==0){
		mol->print();
		printf("TOTAL AREA: %f\n",area);
	}
	
	if(numericalDetail>0){
		Molecule *mol_numerical;
		mol_numerical = dfgro.digestGRO(*top);
		
		IntegratorNumerical integratorNumerical(numericalDetail,numericalDifference);
		area_numerical = integratorNumerical.integrate(mol_numerical, -1, &progressbar);
		
		if(compare==0){
			mol_numerical->print();
			printf("TOTAL AREA: %f\n",area_numerical);
		}
		else{
			mol->printDifference(mol_numerical);
			printf("TOTAL AREA: %f\n",area - area_numerical);
		}
	}

	
	
	
}
