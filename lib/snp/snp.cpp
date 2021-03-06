// ***************************************************************************
// snp.cpp (c) 2013 zhenhua yu <yzh163@mail.ustc.edu.cn>
// HI Lab, University of Science and Technology of China
// All rights reserved.

#include <cstring>

#include "snp.h"
#include "split.h"


SNP::SNP(string name, long long position, string observed, char strand, char ref) {
	this->name = name;
	this->position = position;
	this->observed = observed;
	this->strand = strand;
	this->ref = ref;
	vector<string> strs = split(observed, '/');
	if(strand == '-'){
		ref = getComplement(ref);
	}
	if(strs[0].at(0) == ref) {
	nucleotide = strs[1].at(0);
	} 
	else{
	nucleotide = strs[0].at(0);
	}

	if(strand == '-'){
	nucleotide = getComplement(nucleotide);
	}

	frequency = 0;
	type = unknown;
}

SNP::SNP(string name, long long position, char ref, char nucleotide) {
	this->name = name;
	this->position = position;
	this->ref = ref;
	this->nucleotide = nucleotide;
	
	frequency = 0;
	type = unknown;
}

SNP::SNP(void) {
}

SNP::~SNP(void) {
}

string SNP::getName() {
	return name;
}

long long SNP::getPosition() {
	return position;
}

char SNP::getStrand() {
	return strand;
}

string SNP::getObserved() {
	return observed;
}

char SNP::getRef() {
	return ref;
}

char SNP::getNucleotide() {
	return nucleotide;
}

float SNP::getFrequency() {
	return frequency;
}

varType SNP::getType() {
	return type;
}

char SNP::getComplement(char nucleotide) {
	switch(nucleotide){
		case 'A': return 'T';
		case 'T': return 'A';
		case 'C': return 'G';
		case 'G': return 'C';
		case 'a': return 't';
		case 't': return 'a';
		case 'c': return 'g';
		case 'g': return 'c';
		case 'N': return 'N';
		default : return 'N';
	}
}

void SNP::setFrequency(float value) {
	this->frequency = value;
}

void SNP::setType(varType type) {
	this->type = type;
}


SNPOnChr::SNPOnChr(void) {
	snp_num = 0;
	SNPFile = NULL;
}

SNPOnChr::~SNPOnChr(void) {
	if(SNPFile != NULL) {
		fclose(SNPFile);
	}
}

vector<string> SNPOnChr::getChroms() {
	vector<string> chroms;
	map<string, vector<SNP> >::iterator it;
	for(it = this->begin(); it != this->end(); it++){
		chroms.push_back(it->first);
	}
	return chroms;
}

long SNPOnChr::SNPNumber() {
	return snp_num;
}

string SNPOnChr::aberOfChr(string chromosome) {
	size_t indx = chromosome.find("chrom");
	if(indx == string::npos){
		indx = chromosome.find("chr");
		if(indx != string::npos){
			chromosome = chromosome.substr(indx+3,chromosome.size()-3);
		}
	}
	else{
		chromosome = chromosome.substr(indx+5,chromosome.size()-5);
	}
	return chromosome;
}

void SNPOnChr::readSNPs(string fname) {
	SNPFile = fopen(fname.c_str(),"r");
	if(!(SNPFile = fopen(fname.c_str(),"r"))){
		cerr << "can not open SNP file " << fname << endl;
		exit(-1);
	}
	vector<string> chromosomes;
	vector<string>::iterator it;
	vector<vector<SNP> > AllSNP;
	vector<vector<SNP> >::iterator it2;
	snp_num = 0;
	char buf[1000];
	char * elems[10];
	int elemnum;
	long line_num = 0;
    while(fgets(buf, 1000, SNPFile)){
		line_num++;
		// here assume the snp file is tab-delimited, every line being:
		// SNP id, chromosome, position, letters, strand, c_ref
		elemnum = split(buf, '\t', elems);
		if(elemnum == 6){
			snp_num++;
			//SNP snp(elems[0], atoll(elems[2]), elems[3], *elems[4], *elems[5]);
			SNP snp("", atoll(elems[2]), elems[3], *elems[4], *elems[5]);
			string chromosome = elems[1];
			chromosome = aberOfChr(chromosome);
			for(it = chromosomes.begin(); it < chromosomes.end(); it++){
				if(!chromosome.compare(*it)){
					break;
				}
			}
			if(it < chromosomes.end()){
				it2 = AllSNP.begin();
				it2 += (it-chromosomes.begin());
				(*it2).push_back(snp);
			}
			else{
				vector<SNP> v;
				v.push_back(snp);
				AllSNP.push_back(v);
				chromosomes.push_back(chromosome);
				//cerr << chromosome << endl;
			}
			
		}
		else{
			cerr << "Warning: malformed snp file " << fname <<
                    		", there should be 6 fields @line " << line_num << endl;
		    	cerr << buf << endl;
		    	//exit(1);
        	}
    }
	for(size_t i = 0; i < chromosomes.size(); i++){
		this->insert(make_pair(chromosomes[i],AllSNP[i]));
	}
	
}
