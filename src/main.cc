#ifdef HAVE_CONFIG_H
#include <config.h>
#endif


#include <string>
#include <boost/program_options.hpp>
#include <boost/filesystem.hpp>
#include <triforce.h>


namespace po = boost::program_options;

using namespace boost::filesystem;



Molecule *mol, *mol_numerical;
string libraryPath, top, struc;
bool error;
double area, area_numerical;
int numericalDetail;
double numericalDifference;
int compare;
string topMode;
TopologyMode topm;
unsigned int buffer;
unsigned int slack;



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
			("structure,s", po::value<string>(&struc), "set structure file")
			("top_mode,m", po::value<string>(&topMode)->default_value(string("gromacs")), "set topology mode")
			("numerical_detail,n", po::value<int>(&numericalDetail)->default_value(0), "numerical detail (default off)")
			("numerical_difference,d", po::value<double>(&numericalDifference)->default_value(0), "numerical difference (default off)")
			("compare,c", po::value<int>(&compare)->default_value(0), "show only semi-analytical/numerical difference (default off)")
			("mld-buffer,b", po::value<unsigned int>(&buffer)->default_value(0), "use b multi layered depth buffers")
			("slack,k", po::value<unsigned int>(&slack)->default_value(1), "use depth-buffer with slack k");
			
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
			top=libraryPath+"/generic.top";
			topMode=string("generic");
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
	
	
	topm=GROMACS;
	if(topMode.compare("generic")==0){
		topm=GROMACS_GENERIC;
	}

	
		
		
	
	DataFile dftop(top);
	Topology *top = dftop.digestTOP(topm);
	//top->print();
	
	
	DataFile dfgro(struc);
	Molecule *mol;
	
	if(path(struc).extension().compare(string(".gro"))==0){
		mol = dfgro.digestGRO(*top);

	}
	else if(path(struc).extension().compare(string(".pdb"))==0){
		mol = dfgro.digestPDB(*top);
	}
	else if(path(struc).extension().compare(string(".xyzr"))==0){
		mol = dfgro.digestXYZR();
	}
	else{
		printf("unknown structure file extension\n");
		exit(-1);
	}
	
	mol->generateNeighbourList();
	

	TriforceInterface trii(libraryPath, buffer, slack);
	double area = trii.calculateSurfaceArea(*mol);
	
	if(compare==0){
		//mol->printxyz();
		mol->print();
		printf("TOTAL AREA: %f\n",area);
	}
	
	if(numericalDetail>0){
		printf("\nNumerical evaluation\n");
		DataFile dfgro2(struc);
		Molecule *mol_numerical;
		if(path(struc).extension().compare(string(".gro"))==0){
			mol_numerical = dfgro2.digestGRO(*top);

		}
		if(path(struc).extension().compare(string(".pdb"))==0){
			mol_numerical = dfgro2.digestPDB(*top);
		}
		
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
