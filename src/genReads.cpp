// ***************************************************************************
// genReads.cpp (c) 2018 Zhenhua Yu <qasim0208@163.com>
// Health Informatics Lab, Ningxia University
// All rights reserved.

#include <iostream>
#include <string>
#include <unistd.h>
#include <getopt.h>
#include <ctime>

#include "MyDefine.h"

void parseArgs(int argc, char *argv[]);
void usage(const char* app);

int main(int argc, char *argv[]) {
	/*** record elapsed time ***/
	time_t start_t, end_t;
	long time_used;
	start_t = time(NULL);
	
	parseArgs(argc, argv);
	
	/*** create thread pool ***/
	threadPool = new ThreadPool(config.getIntPara("threads"));
	threadPool->pool_init();
	
	/*** load sequence data ***/
	genome.loadRefSeq();
	
	/*** create fragments ***/
	malbac.createFrags();
	
	/*** amplify fragments ***/
	malbac.amplify();
	
	/*** load model ***/
	profile.train(config.getStringPara("profile"));
	
	/*** generating reads ***/
	malbac.yieldReads();
	cerr << "\nReads generation done!" << endl;
	
	end_t = time(NULL);
	time_used = end_t-start_t;
	int minutes = time_used/60;
	int seconds = time_used - minutes*60;
	cerr << "\nElapsed time: " << minutes << " minutes and " << seconds << " seconds!\n" << endl;
	
	return 0;
}

void parseArgs(int argc, char *argv[]) {
	string modelFile = "", inputFile = "";
	string outputPrefix = "", layout = "PE";
	
	int primersPerKB = 5;
	double ADOR = 0.01, FPR = 0.05;
	
	int threads = 1, readLength = 0, coverage = 0;
	int isize = 300;
	
	/*** record elapsed time ***/
	time_t start_t, end_t;
	long time_used;
	start_t = time(NULL);
	
	struct option long_options[] = {
		{"help", no_argument, 0, 'h'},
		{"input", required_argument, 0, 'i'},
		{"primers", required_argument, 0, 'p'},
		{"ador", required_argument, 0, 'a'},
		{"fpr", required_argument, 0, 'f'},
		{"model", required_argument, 0, 'm'},
		{"layout", required_argument, 0, 'l'},
		//{"readLength", required_argument, 0, 'L'},
		{"coverage", required_argument, 0, 'c'},
		{"isize", required_argument, 0, 's'},
		{"threads", required_argument, 0, 't'},
		{"output", required_argument, 0, 'o'},
		{0, 0, 0, 0}
	};

	int c;
	//Parse command line parameters
	while((c = getopt_long(argc, argv, "hi:p:a:f:m:l:c:s:t:o:", long_options, NULL)) != -1){
		switch(c){
			case 'h':
				usage(argv[0]);
				exit(0);
			case 'i':
				inputFile = optarg;
				break;
			case 'p':
				primersPerKB = atoi(optarg);
				break;
			case 'a':
				ADOR = atof(optarg);
				break;
			case 'f':
				FPR = atof(optarg);
				break;
			case 'm':
				modelFile = optarg;
				break;
			case 'l':
				layout = optarg;
				break;
			case 'L':
				readLength = atoi(optarg);
				break;
			case 'c':
				coverage = atoi(optarg);
				break;
			case 's':
				isize = atoi(optarg);
				break;
			case 't':
				threads = atoi(optarg);
				break;
			case 'o':
				outputPrefix = optarg;
				break;
			default :
				usage(argv[0]);
				exit(1);
		}
	}
	
	if(inputFile.empty()) {
		cerr << "Error: reference file (.fasta) not specified!" << endl;
		usage(argv[0]);
		exit(1);
	}
	
	/*
	if(primersPerKB < 1 || primersPerKB > 15) {
		cerr << "Error: the value of parameter \"primers\" should be 1~15!" << endl;
		exit(1);
	}
	*/
	
	if(ADOR < 0 || ADOR > 0.5) {
		cerr << "Error: the value of parameter \"ador\" should be 0~0.5!" << endl;
		exit(1);
	}
	
	if(FPR < 0 || FPR > 0.5) {
		cerr << "Error: the value of parameter \"fpr\" should be 0~0.5!" << endl;
		exit(1);
	}
	
	if(modelFile.empty()) {
		cerr << "Error: sequencing profile must be specified!" << endl;
		usage(argv[0]);
		exit(1);
	}
	
	if(outputPrefix.empty()) {
		cerr << "Error: the prefix of output file not specified!" << endl;
		usage(argv[0]);
		exit(1);
	}
	if(layout.empty()) {
		cerr << "Warning: sequence layout not specified!" << endl;
		cerr << "use the default value: \"PE for paired-end\"" << endl;
		layout = "PE";
	}
	else if(layout.compare("SE") != 0 && layout.compare("PE") != 0) {
		cerr << "Error: sequence layout incorrectly specified!" << endl;
		cerr << "should be SE (single end) or PE (paired-end)" << endl;
		exit(1);
	}
	/*
	if(readLength < 50 || readLength > 500) {
		cerr << "Error: read length should be 50~500!" << endl;
		exit(1);
	}
	if(layout.compare("PE") == 0 && isize < readLength) {
		cerr << "Error: insert size should be not smaller than read length!" << endl;
		exit(1);
	}
	*/
	if(coverage < 1) {
		cerr << "Error: sequence coverage should be a positive integer!" << endl;
		exit(1);
	}
	if(threads < 1) {
		cerr << "Error: number of threads should be a positive integer!" << endl;
		exit(1);
	}
	
	config.setStringPara("ref", inputFile);
	config.setIntPara("primers", primersPerKB);
	config.setRealPara("ador", ADOR);
	config.setRealPara("fpr", FPR);
	config.setStringPara("profile", modelFile);
	config.setStringPara("output", outputPrefix);
	config.setStringPara("layout", layout);
	//config.setIntPara("readLength", readLength);
	config.setIntPara("coverage", coverage);
	config.setIntPara("isize", isize);
	config.setIntPara("threads", threads);
}

void usage(const char* app) {
	cerr << "\nVersion: " << current_version << endl << endl;
	cerr << "Usage: " << app << " [options]" << endl
		<< endl
		<< "Options:" << endl
		<< "    -h, --help                      give this information" << endl
		<< "    -i, --input <string>            sequence file (.fasta) generated by simuVars program" << endl
		<< "  MALBAC options:" << endl
		<< "    -p, --primers <int>             average number of primers attached to 1KB-length genomic fragment during amplification [Default:2]" << endl
		<< "    -a, --ador <float>              allele dropout rate [Default:0.01]" << endl
		<< "    -f, --fpr <float>               false positive rate [Default:0.05]" << endl
		<< "  Read simulation options:" << endl
		<< "    -m, --model <string>            profile inferred from real sequencing data" << endl
		<< "    -l, --layout <string>           read layout (SE for single end, PE for paired-end) [Default:PE]" << endl
		//<< "    -L, --readLength <int>          the length of reads to be generated" << endl
		<< "    -c, --coverage <int>            sequencing coverage [Default:5]" << endl
		<< "    -s, --isize <int>               mean insert size for paired-end sequencing [Default:300]" << endl
		<< "    -t, --threads <int>             number of threads to use [Default:1]" << endl
		<< "    -o, --output <string>           the prefix of output file" << endl
		<< endl
		<< "Example:" << endl
		<< "    " << app << " -i /path/to/ref.fa -m /path/to/hiseq2500.profile -t 5 -o /path/to/reads" << endl
		<< endl
		<< "Author: Zhenhua Yu <qasim0208@163.com>\n" << endl;
}

