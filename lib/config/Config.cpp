// ***************************************************************************
// Config.cpp (c) 2018 Zhenhua Yu <qasim0208@163.com>
// Health Informatics Lab, Ningxia University
// All rights reserved.

#include <iostream>
#include <fstream>
#include <cassert>

#include "MyDefine.h"
#include "Config.h"

Config::Config() {
	string strParaNames[] = {"bam", "profile", "ref", "target", "var", "snp", 
							"vcf", "samtools", "bases", "output", "layout"};
	
	/*---start default configuration---*/
	
	int n = sizeof(strParaNames)/sizeof(string);
	for(int i = 0; i < n; i++) {
		if(strParaNames[i].compare("layout") == 0) {
			stringParas.insert(make_pair(strParaNames[i], "PE"));
		}
		else if(strParaNames[i].compare("bases") == 0) {
			stringParas.insert(make_pair(strParaNames[i], "ACGT"));
		}
		else {
			stringParas.insert(make_pair(strParaNames[i], ""));
		}
	}
	intParas.insert(make_pair("primers", 100000));
	intParas.insert(make_pair("threads", 1));
	intParas.insert(make_pair("verbose", 1));
	intParas.insert(make_pair("readLength", 0));
	intParas.insert(make_pair("ploidy", 2));
	intParas.insert(make_pair("isize", 260));
	intParas.insert(make_pair("kmer", 3));
	intParas.insert(make_pair("bins", 0));
	intParas.insert(make_pair("ampliconMaxLen", 2000));
	intParas.insert(make_pair("ampliconMinLen", 1000));
	intParas.insert(make_pair("fragSize", 1000));
	
	realParas.insert(make_pair("coverage", 0));
	realParas.insert(make_pair("indelRate", 0.00025));
	realParas.insert(make_pair("gamma", -1));
	realParas.insert(make_pair("ber", 3.4e-4));
	// refer to https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0105585
	/*---end default configuration---*/
}

string Config::getStringPara(string paraName) {
	if(stringParas.find(paraName) != stringParas.end()) {
		return stringParas[paraName];
	}
	else {
		cerr << "Error: unrecognized parameter name \"" << paraName << "\"" << endl;
		exit(1);
	}
	return "";
}

void Config::setStringPara(string paraName, string value) {
	if(stringParas.find(paraName) != stringParas.end()) {
		stringParas[paraName] = value;
	}
	else {
		cerr << "Warning: unrecognized parameter name \"" << paraName << "\"" << endl;
	}
}

long Config::getIntPara(string paraName) {
	if(intParas.find(paraName) != intParas.end()) {
		return intParas[paraName];
	}
	else {
		cerr << "Error: unrecognized parameter name \"" << paraName << "\"" << endl;
		exit(1);
	}
}

void Config::setIntPara(string paraName, long value) {
	if(intParas.find(paraName) != intParas.end()) {
		intParas[paraName] = value;
	}
	else {
		cerr << "Warning: unrecognized parameter name \"" << paraName << "\"" << endl;
	}
}

double Config::getRealPara(string paraName) {
	if(realParas.find(paraName) != realParas.end()) {
		return realParas[paraName];
	}
	else {
		cerr << "Error: unrecognized parameter name \"" << paraName << "\"" << endl;
		exit(1);
	}
}

void Config::setRealPara(string paraName, double value) {
	if(realParas.find(paraName) != realParas.end()) {
		realParas[paraName] = value;
	}
	else {
		cerr << "Warning: unrecognized parameter name \"" << paraName << "\"" << endl;
	}
}
