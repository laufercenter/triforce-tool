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
string libraryPath, topname, topname1, struc;
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
string output, output2;
FILE *outputfile0, *outputfile1;
bool usingOutputFile0;
bool usingOutputFile1;
Topology *top, *top1;


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
			("output2,g", po::value<string>(&output2), "set output2 file")
			("library_path,l", po::value<string>(&libraryPath), "set library path")
			("topology,t", po::value<string>(&topname), "set topology file")
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
		else if(topMode.compare("generic2")==0){
			topm=GROMACS_GENERIC2;
		}
		
		
		
		error=false;
		if(!vm.count("library_path")){
			fprintf(stderr, "library path not set\n");
			error=true;
		}
		if(!vm.count("topology")){
			if(topm==GROMACS_ELEMENTAL){
				topname=libraryPath+"/elemental.top";
			}
			else if(topm==GROMACS_GENERIC){
				topname=libraryPath+"/generic.top";
				topMode=string("generic");
				topm=GROMACS_GENERIC;
			}
			else if(topm==GROMACS_GENERIC2){
				topname=libraryPath+"/generic2.top";
				topname1=libraryPath+"/elemental2.top";				
				topMode=string("generic2");
			}
		}
		if(!vm.count("structure")){
			fprintf(stderr, "structure file not set\n");
			error=true;
		}
		usingOutputFile0=false;
		if(vm.count("output")){
			outputfile0 = fopen(output.c_str(),"w");
			usingOutputFile0=true;
		}
		else outputfile0 = stdout;
		
		usingOutputFile1=false;
		if(vm.count("output2")){
			outputfile1 = fopen(output2.c_str(),"w");
			usingOutputFile1=true;
		}
		
		if(error){
			cout << desc << "\n";
			
			exit(-1);
		}
	}
	catch(std::exception& e) {
		cerr << "error: " << e.what() << "\n";
		exit(-1);
	}
	
	if(hydrogens==0) useHydrogens=false;
	else useHydrogens=true;



	
		
		
	
	DataFile dftop(topname);
	top = dftop.digestTOP(topm);

	
	if(topm==GROMACS_GENERIC2){
		DataFile dftop(topname1);
		//add elemental topology on top of other one
		top = dftop.digestTOP(GROMACS_ELEMENTAL,top);
	}
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
	
	

	TriforceInterface trii(libraryPath, buffer, slack);
	area = trii.calculateSurfaceArea(*mol);
	
	//mol->printxyz();
	mol->print(outputfile0);
	
	if(usingOutputFile1){
		trii.printSurfaces(*mol,outputfile1);
	}
	
	
	fprintf(stdout,"TOTAL AREA: %f\n",area);
	
	if(numericalDetail>0){
		fprintf(outputfile0,"\nNumerical evaluation\n");
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
		
		mol_numerical->print(outputfile0);
		fprintf(outputfile0,"TOTAL AREA: %f\n",area_numerical);
	}
	
	trii.printBenchmark(stdout);
	
	
	if(usingOutputFile0){
		fclose(outputfile0);
	}
	if(usingOutputFile1){
		fclose(outputfile1);
	}
	

	
	
	
}
