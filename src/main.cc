#ifdef HAVE_CONFIG_H
#include <config.h>
#endif


#include <string>
#include <boost/program_options.hpp>
#include "triforce-0.39/triforce.h"


namespace po = boost::program_options;



Molecule *mol;
string libraryPath, top, gro;
bool error;
double area;


int main(int ac, char* av[])
{ 
	try{
		po::options_description desc("Allowed options");
		desc.add_options()
			("help,h", "help")
			("library_path,l", po::value<string>(&libraryPath), "set library path")
			("topology,t", po::value<string>(&top), "set topology file")
			("structure,s", po::value<string>(&gro), "set structure file");
			
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
			fprintf(stderr, "structure not set\n");
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

	mol->print();
	
	printf("TOTAL AREA: %f\n",area);
	
	
}
