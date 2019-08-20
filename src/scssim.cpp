// ***************************************************************************
// scssim.cpp (c) 2019 Zhenhua Yu <qasim0208@163.com>
// Health Informatics Lab, Ningxia University
// All rights reserved.

#include <iostream>
#include <string>
#include <unistd.h>
#include <getopt.h>
#include <ctime>

#include "MyDefine.h"

void parseArgs_simuVars(int argc, char *argv[]);
void parseArgs_learnProfile(int argc, char *argv[]);
void parseArgs_genReads(int argc, char *argv[]);
void parseArgs(int argc, char *argv[]);
void usage(const char* app);
void usage_simuVars(const char* app);
void usage_learnProfile(const char* app);
void usage_genReads(const char* app);

int main(int argc, char *argv[]) {
	/*** record elapsed time ***/
	time_t start_t, end_t;
	long time_used;
	start_t = time(NULL);
	
	parseArgs(argc, argv);
	
	string subcmd = argv[1];
	
	if(subcmd.compare("simuvars") == 0) {
		/*** load data ***/
		genome.loadData();
		/*** create and save sequences ***/
		genome.saveSequence();
	}
	else if(subcmd.compare("learn") == 0) {
		/*** load data ***/
		genome.loadTrainData();
		/*** profile learning ***/
		profile.init();
		profile.train();
	}
	else { // genreads
		srand(start_t);
		/*** create thread pool ***/
		threadPool = new ThreadPool(config.getIntPara("threads"));
		threadPool->pool_init();
		
		/*** load sequence data ***/
		genome.loadData();
		
		/*** load model ***/
		profile.train(config.getStringPara("profile"));
		
		/*** create fragments ***/
		malbac.createFrags();
		
		/*** amplify fragments ***/
		malbac.amplify();
		
		/*** generating reads ***/
		malbac.yieldReads();
		cerr << "\nReads generation done!" << endl;
	}
	
	end_t = time(NULL);
	time_used = end_t-start_t;
	int minutes = time_used/60;
	int seconds = time_used - minutes*60;
	cerr << "\nElapsed time: " << minutes << " minutes and " << seconds << " seconds!\n" << endl;
	
	return 0;
}

void parseArgs(int argc, char *argv[]) {
	if(argc == 1) {
		usage(argv[0]);
		exit(0);
	}
	string subcmd = argv[1];
	if(subcmd.compare("-h") == 0 || subcmd.compare("--help") == 0) {
		usage(argv[0]);
		exit(0);
	}
	else if(subcmd.compare("-v") == 0 || subcmd.compare("--version") == 0) {
		cerr << "SCSsim version " << current_version << endl;
		exit(0);
	}
	else if(subcmd.compare("simuvars") == 0) {
		parseArgs_simuVars(argc-1, &argv[1]);
	}
	else if(subcmd.compare("learn") == 0) {
		parseArgs_learnProfile(argc-1, &argv[1]);
	}
	else if(subcmd.compare("genreads") == 0) {
		parseArgs_genReads(argc-1, &argv[1]);
	}
	else {
		cerr << "Error: unrecognized subcommand \"" << subcmd << "\"." << endl;
		usage(argv[0]);
		exit(0);
	}
}

void parseArgs_simuVars(int argc, char *argv[]) {
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
				usage_simuVars(argv[0]);
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
				usage_simuVars(argv[0]);
				exit(1);
		}
	}

	if(refFile.empty()){
		cerr << "Use --ref to specify the reference file (fasta)." << endl;
		usage_simuVars(argv[0]);
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
		usage_simuVars(argv[0]);
		exit(1);
	}
	
	config.setStringPara("ref", refFile);
	config.setStringPara("snp", snpFile);
	config.setStringPara("var", varFile);
	config.setStringPara("output", outFile);
}

void parseArgs_learnProfile(int argc, char *argv[]) {
	string bamFile = "", targetFile = "";
	string vcfFile = "", refFile = "";
	string outFile = "", samtools = "";
	int wsize = 1000;
	int kmer = 3;

	struct option long_options[] = {
		{"help", no_argument, 0, 'h'},
		{"bam", required_argument, 0, 'b'},
		{"target", required_argument, 0, 't'},
		{"vcf", required_argument, 0, 'v'},
		{"ref", required_argument, 0, 'r'},
		{"wsize", required_argument, 0, 'w'},
		{"kmer", required_argument, 0, 'k'},
		{"output", required_argument, 0, 'o'},
		{"samtools", required_argument, 0, 's'},
		{0, 0, 0, 0}
	};

	int c;
	//Parse command line parameters
	while((c = getopt_long(argc, argv, "hb:t:v:r:w:k:o:s:", long_options, NULL)) != -1) {
		switch(c){
			case 'h':
				usage_learnProfile(argv[0]);
				exit(0);
			case 'b':
				bamFile = optarg;
				break;
			case 't':
				targetFile = optarg;
				break;
			case 'v':
				vcfFile = optarg;
				break;
			case 'r':
				refFile = optarg;
				break;
			case 'w':
				wsize = atoi(optarg);
				break;
			case 'k':
				kmer = atoi(optarg);
				break;
			case 'o':
				outFile = optarg;
				break;
			case 's':
				samtools = optarg;
				break;
			default :
				usage_learnProfile(argv[0]);
				exit(1);
		}
	}

	if(bamFile.empty()){
		cerr << "Use --bam to specify a normal BAM file." << endl;
		usage_learnProfile(argv[0]);
		exit(1);
	}

	if(!targetFile.empty()){
		cerr << "Target regions specified!" << endl;
	}
	
	if(vcfFile.empty()){
		cerr << "Use --vcf to specify the VCF file generated from the normal BAM." << endl;
		usage_learnProfile(argv[0]);
		exit(1);
	}

	if(refFile.empty()){
		cerr << "Use --ref to specify the reference file(.fasta) to which the reads are aligned." << endl;
		usage_learnProfile(argv[0]);
		exit(1);
	}
	
	if(outFile.empty()){
		cerr << "Use --output to specify the output file." << endl;
		usage_learnProfile(argv[0]);
		exit(1);
	}

	if(samtools.empty()) {
		cerr << "Warning: the path of samtools not specified!" << endl;
		cerr << "Assume the tool has been installed and included in the system PATH!" << endl;
	}
	
	if(wsize < 100) {
		cerr << "Error: the value of parameter \"wsize\" should be at least 100!" << endl;
		usage_learnProfile(argv[0]);
		exit(1);
	}
	
	if(kmer < 1 || kmer > 5) {
		cerr << "Error: parameter \"kmer\" should be a positive integer with maximum value of 5!" << endl;
		exit(1);
	}
	
	config.setStringPara("bam", bamFile);
	config.setStringPara("ref", refFile);
	config.setStringPara("target", targetFile);
	config.setStringPara("vcf", vcfFile);
	config.setStringPara("samtools", samtools);
	config.setStringPara("output", outFile);
	config.setIntPara("fragSize", wsize);
	config.setIntPara("kmer", kmer);
}

void parseArgs_genReads(int argc, char *argv[]) {
	string modelFile = "", inputFile = "";
	string outputPrefix = "", layout = "PE";
	
	long primers = 100000;
	double gamma = 1e-9;
	
	int threads = 1, isize = 260;
	double coverage = 5;
	
	struct option long_options[] = {
		{"help", no_argument, 0, 'h'},
		{"input", required_argument, 0, 'i'},
		{"primers", required_argument, 0, 'p'},
		{"gamma", required_argument, 0, 'r'},
		{"model", required_argument, 0, 'm'},
		{"layout", required_argument, 0, 'l'},
		{"coverage", required_argument, 0, 'c'},
		{"isize", required_argument, 0, 's'},
		{"threads", required_argument, 0, 't'},
		{"output", required_argument, 0, 'o'},
		{0, 0, 0, 0}
	};

	int c;
	//Parse command line parameters
	while((c = getopt_long(argc, argv, "hi:p:r:m:l:c:s:t:o:", long_options, NULL)) != -1){
		switch(c){
			case 'h':
				usage_genReads(argv[0]);
				exit(0);
			case 'i':
				inputFile = optarg;
				break;
			case 'p':
				primers = atol(optarg);
				break;
			case 'r':
				gamma = atof(optarg);
				break;
			case 'm':
				modelFile = optarg;
				break;
			case 'l':
				layout = optarg;
				break;
			case 'c':
				coverage = atof(optarg);
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
				usage_genReads(argv[0]);
				exit(1);
		}
	}
	
	if(inputFile.empty()) {
		cerr << "Error: reference file (.fasta) not specified!" << endl;
		usage_genReads(argv[0]);
		exit(1);
	}
	
	if(primers < 1000) {
		cerr << "Error: the value of parameter \"primers\" should be at least 1000!" << endl;
		exit(1);
	}
	
	if(gamma <= 0 || gamma > 1e-8) {
		cerr << "Error: the value of parameter \"gamma\" should be in 0~1e-8!" << endl;
		exit(1);
	}
	
	if(modelFile.empty()) {
		cerr << "Error: sequencing profile must be specified!" << endl;
		usage_genReads(argv[0]);
		exit(1);
	}
	
	if(outputPrefix.empty()) {
		cerr << "Error: the prefix of output file not specified!" << endl;
		usage_genReads(argv[0]);
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
	if(coverage <= 0) {
		cerr << "Error: sequencing coverage not properly specified!" << endl;
		exit(1);
	}
	if(threads < 1) {
		cerr << "Error: number of threads should be a positive integer!" << endl;
		exit(1);
	}
	
	config.setStringPara("ref", inputFile);
	config.setStringPara("profile", modelFile);
	config.setStringPara("output", outputPrefix);
	config.setStringPara("layout", layout);
	config.setIntPara("primers", primers);
	config.setRealPara("gamma", gamma);
	config.setRealPara("coverage", coverage);
	config.setIntPara("isize", isize);
	config.setIntPara("threads", threads);
}

void usage(const char* app) {
	cerr << "\nSCSsim version: " << current_version << endl;
	cerr << "Usage: " << app << " [subcommand] [options]" << endl
		<< endl
		<< "Optional arguments:" << endl
		<< "    -h, --help                      give this information" << endl
		<< "    -v, --version <string>          print software version" << endl
		<< endl
		<< "Available subcmds:" << endl
		<< "    simuvars          simulate the genome sequence of single cells" << endl
		<< "    learn             learn sequencing profiles from real sequencing data" << endl
		<< "    genreads          simulate sequencing reads of single cell" << endl
		<< endl
		<< "Author: Zhenhua Yu <qasim0208@163.com>\n" << endl;
}

void usage_simuVars(const char* app) {
	cerr << "Usage: scssim " << app << " [options]" << endl
		<< endl
		<< "Options:" << endl
		<< "    -h, --help                      give this information" << endl
		<< "    -r, --ref <string>              reference file (.fasta)" << endl
		<< "    -s, --snp <string>              SNP file containing the SNPs to be simulated [Default:null]" << endl
		<< "    -v, --var <string>              variation file containing the genomic variations to be simulated [Default:null]" << endl
		<< "    -o, --output <string>           output file (.fasta) to save generated sequences" << endl
		<< endl
		<< "Example:" << endl
		<< "    scssim " << app << " -r /path/to/hg19.fa -s /path/to/hg19.snp138.1based.txt -v /path/to/variation.txt -o /path/to/results.fa" << endl
		<< endl
		<< "    scssim " << app << " -r /path/to/hg19.fa -v /path/to/variation.txt -o /path/to/results.fa" << endl
		<< endl
		<< "    scssim " << app << " -r /path/to/hg19.fa -s /path/to/hg19.snp138.1based.txt -o /path/to/results.fa" << endl
		<< endl
		<< "Author: Zhenhua Yu <qasim0208@163.com>\n" << endl;
}

void usage_learnProfile(const char* app) {
	cerr << "Usage: scssim " << app << " [options]" << endl
		<< endl
		<< "Options:" << endl
		<< "    -h, --help                      give this information" << endl
		<< "    -b, --bam <string>              normal BAM file" << endl
		<< "    -t, --target <string>           exome target file (.bed) for whole-exome sequencing[default:null]" << endl
		<< "    -v, --vcf <string>              the VCF file generated from the normal BAM" << endl
		<< "    -r, --ref <string>              genome reference file (.fasta) to which the reads were aligned" << endl
		<< "    -w, --wsize <int>               the length of windows used to infer GC-content bias[default:1000]" << endl
		<< "    -k, --kmer <int>                the length of kmer sequence [default:3]" << endl
		<< "    -o, --output <string>           output file" << endl
		<< "    -s, --samtools <string>         the path of samtools [default:samtools]" << endl
		<< endl
		<< "Example:" << endl
		<< "    scssim " << app << " -b /path/to/normal.bam -t /path/to/normal.bed -v /path/to/normal.vcf -r /path/to/ref.fa > /path/to/results.profile" << endl
		<< endl
		<< "    scssim " << app << " -b /path/to/normal.bam -v /path/to/normal.vcf -r /path/to/ref.fa -o /path/to/results.profile -s /path/to/samtools" << endl
		<< endl
		<< "Author: Zhenhua Yu <qasim0208@163.com>\n" << endl;
}

void usage_genReads(const char* app) {
	cerr << "Usage: scssim " << app << " [options]" << endl
		<< endl
		<< "Options:" << endl
		<< "    -h, --help                      give this information" << endl
		<< "    -i, --input <string>            sequence file (.fasta) generated by simuVars program" << endl
		<< "  MALBAC options:" << endl
		<< "    -p, --primers <int>             the number of primers [Default:100000]" << endl
		<< "    -r, --gamma <float>             a parameter controlling the number of primers used in each cycle [Default:1e-9]" << endl
		<< "  Read simulation options:" << endl
		<< "    -m, --model <string>            profile inferred from real sequencing data" << endl
		<< "    -l, --layout <string>           read layout (SE for single end, PE for paired-end) [Default:PE]" << endl
		<< "    -c, --coverage <float>          sequencing coverage [Default:5]" << endl
		<< "    -s, --isize <int>               mean insert size for paired-end sequencing [Default:260]" << endl
		<< "    -t, --threads <int>             number of threads to use [Default:1]" << endl
		<< "    -o, --output <string>           the prefix of output file" << endl
		<< endl
		<< "Example:" << endl
		<< "    scssim " << app << " -i /path/to/ref.fa -m /path/to/hiseq2500.profile -t 5 -o /path/to/reads" << endl
		<< endl
		<< "Author: Zhenhua Yu <qasim0208@163.com>\n" << endl;
}

