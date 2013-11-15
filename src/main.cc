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
float area, area_numerical;
int numericalDetail;
float numericalDifference;
int compare;
string topMode;
TopologyMode topm;
unsigned int buffer;
unsigned int slack;
unsigned int hydrogens;
bool useHydrogens;
string output;
FILE* outputfile;
bool usingOutputFile;



void progressbar(float percentage){
	int size=80;
	
	cout << "\r" "[";
	for(int i=0; i<size; ++i){
		if((float)i/(float)(size-1) <= percentage)
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
			("output,o", po::value<string>(&output), "set output file")
			("library_path,l", po::value<string>(&libraryPath), "set library path")
			("topology,t", po::value<string>(&top), "set topology file")
			("structure,s", po::value<string>(&struc), "set structure file")
			("top_mode,m", po::value<string>(&topMode)->default_value(string("gromacs")), "set topology mode")
			("numerical_detail,n", po::value<int>(&numericalDetail)->default_value(0), "numerical detail (default off)")
			("numerical_difference,d", po::value<float>(&numericalDifference)->default_value(0), "numerical difference (default off)")
			("mld-buffer,b", po::value<unsigned int>(&buffer)->default_value(0), "use b multi layered depth buffers (default 0)")
			("slack,k", po::value<unsigned int>(&slack)->default_value(6), "use depth-buffer with slack k (default 6)")
			("hydrogens,y", po::value<unsigned int>(&hydrogens)->default_value(1), "use hydrogens (default on)");
			
		po::variables_map vm;
		po::store(po::parse_command_line(ac,av,desc), vm);
		po::notify(vm);

		
		
		if(vm.count("help")){
			exit(0);
		}
		
		
		
		topm=GROMACS;
		if(topMode.compare("generic")==0){
			topm=GROMACS_GENERIC;
		}
		else if(topMode.compare("elemental")==0){
			topm=GROMACS_ELEMENTAL;
		}
		
		
		
		error=false;
		if(!vm.count("library_path")){
			fprintf(stderr, "library path not set\n");
			error=true;
		}
		if(!vm.count("topology")){
			if(topm==GROMACS_ELEMENTAL){
				top=libraryPath+"/elemental.top";
			}
			else{
				top=libraryPath+"/generic.top";
				topMode=string("generic");
				topm=GROMACS_GENERIC;
			}
		}
		if(!vm.count("structure")){
			fprintf(stderr, "structure file not set\n");
			error=true;
		}
		usingOutputFile=false;
		if(vm.count("output")){
			outputfile = fopen(output.c_str(),"w");
			usingOutputFile=true;
		}
		else outputfile = stdout;
		
		if(error){
			cout << desc << "\n";
			
			exit(-1);
		}
	}
	catch(exception& e) {
		cerr << "error: " << e.what() << "\n";
		exit(-1);
	}
	
	if(hydrogens==0) useHydrogens=false;
	else useHydrogens=true;



	
		
		
	
	DataFile dftop(top);
	Topology *top = dftop.digestTOP(topm);
	//top->print();
	
	
	DataFile dfgro(struc);
	Molecule *mol;
	
	if(path(struc).extension().compare(string(".gro"))==0){
		mol = dfgro.digestGRO(*top,useHydrogens);

	}
	else if(path(struc).extension().compare(string(".pdb"))==0){
		mol = dfgro.digestPDB(*top,useHydrogens);
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
	float area = trii.calculateSurfaceArea(*mol);
	
	//mol->printxyz();
	mol->print(outputfile);
	fprintf(outputfile,"TOTAL AREA: %f\n",area);
	
	if(numericalDetail>0){
		fprintf(outputfile,"\nNumerical evaluation\n");
		DataFile dfgro2(struc);
		Molecule *mol_numerical;
		if(path(struc).extension().compare(string(".gro"))==0){
			mol_numerical = dfgro2.digestGRO(*top,useHydrogens);

		}
		else if(path(struc).extension().compare(string(".pdb"))==0){
			mol_numerical = dfgro2.digestPDB(*top,useHydrogens);
		}
		else{
			printf("unknown structure file extension\n");
			exit(-1);
		}
		
		IntegratorNumerical integratorNumerical(numericalDetail,numericalDifference);
		area_numerical = integratorNumerical.integrate(mol_numerical, -1, &progressbar);
		
		mol_numerical->print(outputfile);
		fprintf(outputfile,"TOTAL AREA: %f\n",area_numerical);
	}
	
	trii.printBenchmark(outputfile);
	
	
	if(usingOutputFile) fclose(outputfile);

	
	
	
}
