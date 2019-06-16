// ***************************************************************************
// simuVars.cpp (c) 2018 Zhenhua Yu <qasim0208@163.com>
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
	
	/*** load data ***/
	genome.loadData();
	/*** create and save sequences ***/
	genome.saveSequence();
	
	end_t = time(NULL);
	time_used = end_t-start_t;
	int minutes = time_used/60;
	int seconds = time_used - minutes*60;
	cerr << "\nElapsed time: " << minutes << " minutes and " << seconds << " seconds!\n" << endl;
	
	return 0;
}

void parseArgs(int argc, char *argv[]) {
	string refFile = "", snpFile = "";
	string varFile = "", outFile = "";

	struct option long_options[] = {
		{"help", no_argument, 0, 'h'},
		{"ref", required_argument, 0, 'r'},
		{"snp", required_argument, 0, 's'},
		{"var", required_argument, 0, 'v'},
		{"output", required_argument, 0, 'o'},
		{0, 0, 0, 0}
	};

	int c;
	//Parse command line parameters
	while((c = getopt_long(argc, argv, "hr:s:v:o:", long_options, NULL)) != -1){
		switch(c){
			case 'h':
				usage(argv[0]);
				exit(0);
			case 'r':
				refFile = optarg;
				break;
			case 's':
				snpFile = optarg;
				break;
			case 'v':
				varFile = optarg;
				break;
			case 'o':
				outFile = optarg;
				break;
			default :
				usage(argv[0]);
				exit(1);
		}
	}

	if(refFile.empty()){
		cerr << "Use --ref to specify the reference file (fasta)." << endl;
		usage(argv[0]);
		exit(1);
	}
	
	if(snpFile.empty()){
		cerr << "Warning: SNP file not specified!" << endl;
		cerr << "No SNPs will be inserted into the genome." << endl;
	}
	
	if(varFile.empty()){
		cerr << "Warning: variation file not specified!" << endl;
		cerr << "No variations will be inserted into the genome." << endl;
	}
	
	if(outFile.empty()){
		cerr << "Use --output to specify the output file." << endl;
		usage(argv[0]);
		exit(1);
	}
	
	config.setStringPara("ref", refFile);
	config.setStringPara("snp", snpFile);
	config.setStringPara("var", varFile);
	config.setStringPara("output", outFile);
}

void usage(const char* app) {
	cerr << "\nVersion: " << current_version << endl << endl;
	cerr << "Usage: " << app << " [options]" << endl
		<< endl
		<< "Options:" << endl
		<< "    -h, --help                      give this information" << endl
		<< "    -r, --ref <string>              reference file (.fasta)" << endl
		<< "    -s, --snp <string>              SNP file containing the SNPs to be simulated [Default:null]" << endl
		<< "    -v, --var <string>              variation file containing the genomic variations to be simulated [Default:null]" << endl
		<< "    -o, --output <string>           output file (.fasta) to save generated sequences" << endl
		<< endl
		<< "Example:" << endl
		<< "    " << app << " -r /path/to/hg19.fa -s /path/to/hg19.snp138.1based.txt -v /path/to/variation.txt -o /path/to/results.fa" << endl
		<< endl
		<< "    " << app << " -r /path/to/hg19.fa -v /path/to/variation.txt -o /path/to/results.fa" << endl
		<< endl
		<< "    " << app << " -r /path/to/hg19.fa -s /path/to/hg19.snp138.1based.txt -o /path/to/results.fa" << endl
		<< endl
		<< "Author: Zhenhua Yu <qasim0208@163.com>\n" << endl;
}

