// ***************************************************************************
// Malbac.cpp (c) 2018 Zhenhua Yu <qasim0208@163.com>
// Health Informatics Lab, Ningxia University
// All rights reserved.

#include <iostream>
#include <fstream>
#include <cmath>
#include <unistd.h>
#include <algorithm>

#include "MyDefine.h"
#include "Malbac.h"

Malbac::Malbac() {
	pthread_mutex_init(&pm, NULL);
}

void Malbac::extendSemiAmplicons(vector<Amplicon>& amplicons) {
	pthread_mutex_lock(&pm);
	semiAmplicons.insert(semiAmplicons.end(),amplicons.begin(),amplicons.end());
	pthread_mutex_unlock(&pm);
}

void Malbac::extendFullAmplicons(vector<Amplicon>& amplicons) {
	pthread_mutex_lock(&pm);
	fullAmplicons.insert(fullAmplicons.end(),amplicons.begin(),amplicons.end());
	pthread_mutex_unlock(&pm);
}

void Malbac::createFrags() {
	genome.splitToFrags(fragments);
}

void Malbac::amplify() {
	cerr << "\nMALBAC amplification..." << endl;
	setPrimers(true);
	amplifyFrags();
	for(int i = 0; i < 5; i++) {
		cerr << "cycle number: " << i+1 << endl;
		setPrimers(false);
		amplifySemiAmplicons();
		cerr << fullAmplicons.size() << endl;
		cerr << "semi amplicon amplification done!" << endl;
		if(i < 4) {
			amplifyFrags();
			cerr << semiAmplicons.size() << endl;
			cerr << "fragment amplification done!" << endl;
		}
	}
}

void Malbac::amplifyAndSaveProducts() {
	string outFile = config.getStringPara("output");
	ofstream ofs;
	ofs.open(outFile.c_str());
	if(!ofs.is_open()) {
		cerr << "Error: cannot open output file " << outFile << "!" << endl;
		exit(-1);
	}
	
	cerr << "\nMALBAC amplification..." << endl;
	setPrimers(true);
	amplifyFrags();
	for(int i = 0; i < 5; i++) {
		cerr << "cycle number: " << i+1 << endl;
		setPrimers(false);
		amplifySemiAmplicons();
		cerr << fullAmplicons.size() << endl;
		cerr << "semi amplicon amplification done!" << endl;
		saveFullAmplicons(ofs);
		if(i < 4) {
			amplifyFrags();
			cerr << semiAmplicons.size() << endl;
			cerr << "fragment amplification done!" << endl;
		}
	}
	
	ofs.close();
}

void Malbac::setPrimers(bool onlyFrags) {
	unsigned long i, k;
	unsigned long fragNum = fragments.size();
	unsigned int length, gcContent;
	double gcSum = 0, totalLen = 0;
	unsigned int primersPerKB = config.getIntPara("primers");
	
	for(i = 0; i < fragNum; i++) {
		gcContent = fragments[i].getGCcontent();
		length = fragments[i].getLength();
		gcSum += 1.0*gcContent/Fragment::maxSize;
		totalLen += length;
	}
	
	if(!onlyFrags) {
		unsigned long semiAmpNum = semiAmplicons.size();
		for(i = 0; i < semiAmpNum; i++) {
			gcContent = semiAmplicons[i].getGCcontent();
			length = semiAmplicons[i].getLength();
			gcSum += 1.0*gcContent/Fragment::maxSize;
			totalLen += length;
		}
	}
	
	unsigned int totalPrimers = totalLen/1000*primersPerKB;
	unsigned int count = 0;
	for(i = 0; i < fragNum; i++) {
		gcContent = fragments[i].getGCcontent();
		k = 1.0*gcContent/Fragment::maxSize*totalPrimers/gcSum;
		fragments[i].setPrimers(k);
		count += k;
	}
	if(!onlyFrags) {
		unsigned long semiAmpNum = semiAmplicons.size();
		for(i = 0; i < semiAmpNum; i++) {
			gcContent = semiAmplicons[i].getGCcontent();
			k = 1.0*gcContent/Fragment::maxSize*totalPrimers/gcSum;
			semiAmplicons[i].setPrimers(k);
		}
	}
}

void Malbac::saveFullAmplicons(ofstream& ofs) {
	unsigned long i, j, k;
	unsigned long ampliconNum = fullAmplicons.size();
	char c;
	char* seq;
	unsigned int sindx, length;
	int width = 100;
	/*
	for(i = 0; i < ampliconNum; i++) {
		ofs << ">" << "amp_" << i+1 << endl;
		seq = fullAmplicons[i].getSequence();
		sindx = 0;
		length = strlen(seq);
		while(sindx < length) {
			if(sindx+width > length) {
				ofs << &(seq[sindx]) << endl;
			}
			else {
				c = seq[sindx+width];
				seq[sindx+width] = '\0';
				ofs << &(seq[sindx]) << endl;
				seq[sindx+width] = c;
			}
			sindx += width;
		}
		delete[] seq;
	}
	*/
	//fullAmplicons.clear();
}

void Malbac::amplifyFrags() {
	unsigned long sindx = 0;
	unsigned long fragNum = fragments.size();
	vector<unsigned long*> threadParas;
	int i = 0;
	unsigned long loadPerThread = max((unsigned long) 10, fragNum/threadPool->getThreadNumber());
	//int loadPerThread = 1000;
	while(sindx < fragNum) {
		unsigned long* indxs = new unsigned long[2];
		if(sindx+loadPerThread > fragNum) {
			indxs[0] = sindx;
			indxs[1] = fragNum-1;
		}
		else {
			indxs[0] = sindx;
			indxs[1] = sindx+loadPerThread-1;
		}
		threadPool->pool_add_work(&Fragment::batchAmplify, indxs, i++);
		threadParas.push_back(indxs);
		sindx += loadPerThread;
	}
	threadPool->wait();
	for(i = 0; i < threadParas.size(); i++) {
		delete[] threadParas[i];
	}
}

void Malbac::amplifySemiAmplicons() {
	unsigned long sindx = 0;
	unsigned long semiAmpliconNum = semiAmplicons.size();
	vector<unsigned long*> threadParas;
	int i = 0;
	unsigned long loadPerThread = max((unsigned long) 10, semiAmpliconNum/threadPool->getThreadNumber());
	while(sindx < semiAmpliconNum) {
		unsigned long* indxs = new unsigned long[2];
		if(sindx+loadPerThread > semiAmpliconNum) {
			indxs[0] = sindx;
			indxs[1] = semiAmpliconNum-1;
		}
		else {
			indxs[0] = sindx;
			indxs[1] = sindx+loadPerThread-1;
		}
		threadPool->pool_add_work(&Amplicon::batchAmplify, indxs, i++);
		threadParas.push_back(indxs);
		sindx += loadPerThread;
	}
	threadPool->wait();
	for(i = 0; i < threadParas.size(); i++) {
		delete[] threadParas[i];
	}
}

void Malbac::setReadCounts(unsigned long reads) {
	unsigned long i;
	unsigned long ampliconNum = fullAmplicons.size();
	double WL = 0;
	for(i = 0; i < ampliconNum; i++) {
		WL += fullAmplicons[i].getWeightedLength();
	}
	
	unsigned long sum = 0;
	for(i = 0; i < ampliconNum; i++) {
		long readCount = (fullAmplicons[i].getWeightedLength()/WL)*reads;
		readNumbers.push_back(readCount);
		sum += readCount;
	}
	unsigned long n = reads-sum;
	while(n > 0) {
		i = randomInteger(0, ampliconNum);
		readNumbers[i]++;
		n--;
	}
}

void Malbac::yieldReads() {	
	int i;
	int refLen = genome.getGenomeLength();
	
	unsigned long reads = refLen*config.getIntPara("coverage")/config.getIntPara("readLength");
	if(config.isVerbose()) {	
		cerr << "\nNumber of reads to sample: " << reads << endl;
	}
	
	srand(time(0));
	setReadCounts(reads);
	
	string fqFilePrefix = config.getStringPara("output");
	if(config.isPairedEnd()) {
		string outFile1 = fqFilePrefix+"_1.fq";
		string outFile2 = fqFilePrefix+"_2.fq";
		swp = new SeqWriter(outFile1, outFile2);
	}
	else {
		string outFile = fqFilePrefix+".fq";
		swp = new SeqWriter(outFile);
	}
	
	
	cerr << "\n*****Producing reads*****" << endl;
	unsigned long ampliconNum = fullAmplicons.size();
	unsigned long sindx = 0;
	vector<unsigned long*> threadParas;
	unsigned long loadPerThread = max((unsigned long) 10, ampliconNum/threadPool->getThreadNumber());
	while(sindx < ampliconNum) {
		unsigned long* indxs = new unsigned long[2];
		if(sindx+loadPerThread > ampliconNum) {
			indxs[0] = sindx;
			indxs[1] = ampliconNum-1;
		}
		else {
			indxs[0] = sindx;
			indxs[1] = sindx+loadPerThread-1;
		}
		threadPool->pool_add_work(&Amplicon::yieldReads, indxs, i++);
		threadParas.push_back(indxs);
		sindx += loadPerThread;
	}
	threadPool->wait();
	for(i = 0; i < threadParas.size(); i++) {
		delete[] threadParas[i];
	}
	
	delete swp;
}
