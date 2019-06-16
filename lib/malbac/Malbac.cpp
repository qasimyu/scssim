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
	pthread_mutex_init(&pm_amp, NULL);
	pthread_mutex_init(&pm_primer, NULL);
	semiAmplicons = createNullLink();
	fullAmplicons = createNullLink();
	primers = NULL;
	readNumbers = NULL;
}

Malbac::~Malbac() {
	if(primers != NULL) {
		for(int i = 0;i < ptypeCount; i++) {
			delete[] primers[i];
		}
		delete[] primers;
	}
	if(readNumbers != NULL) {
		delete[] readNumbers;
	}
}

void Malbac::createPrimers() {
	int i, j, k, n;
	int tmp[8];
	string bases = config.getStringPara("bases");
	int N = bases.length();
	long primerNum = config.getIntPara("primers");
	
	ptypeCount = pow(N, 8);
	primers = new char*[ptypeCount];
	for(i = 0; i < ptypeCount; i++) {
		primers[i] = new char[9];
		primers[i][8] = '\0';
	}
	
	for(i = 0; i < 8; i++) {
		tmp[i] = 0;
	}
	k = 0;
	while(tmp[0] < N) {
		PrimerIndex *sIndx = &rIndex;
		for(i = 0; i < 8; i++) {
			char c = bases[tmp[i]];
			primers[k][i] = c;
			sIndx = &(sIndx->nextIndexs[c]);
		}
		sIndx->index = k;
		sIndx->count = primerNum;
		k++;
		n = 1;
		for(i = 7; i > 0; i--) {
			if(n == 0) {
				break;
			}
			tmp[i] += n;
			if(tmp[i] == N) {
				tmp[i] = 0;
				n = 1;
			}
			else {
				n = 0;
			}
		}
		tmp[0] += n;
	}
	totalPrimers = k*primerNum;
}

long Malbac::getPrimerCount(const char *s) {
	PrimerIndex *sIndx = &rIndex;
	for(int i = 0; s[i] != '\0'; i++) {
		sIndx = &(sIndx->nextIndexs[s[i]]);
	}
	return sIndx->count;
}

int Malbac::updatePrimerCount(const char *s, int n) {
	PrimerIndex *sIndx = &rIndex;
	for(int i = 0; s[i] != '\0'; i++) {
		sIndx = &(sIndx->nextIndexs[s[i]]);
	}
	if(sIndx->count+n < 0) {
		return 0;
	}
	pthread_mutex_lock(&pm_primer);
	sIndx->count += n;
	pthread_mutex_unlock(&pm_primer);
	return 1;
}

void Malbac::extendSemiAmplicons(AmpliconLink amplicons) {
	if(amplicons == NULL) {
		return;
	}
	if(amplicons->link == amplicons) {
		free(amplicons);
		return;
	}
	pthread_mutex_lock(&pm_amp);
	AmpliconLink p = amplicons->link;
	amplicons->link = semiAmplicons->link;
	semiAmplicons->link = p->link;
	semiAmplicons = amplicons;
	Amplicon& amplicon = semiAmplicons->link->amplicon;
	amplicon.setData(amplicon.getData()+p->amplicon.getData()); // update semi amplicon count
	free(p);
	pthread_mutex_unlock(&pm_amp);
}

void Malbac::extendFullAmplicons(AmpliconLink amplicons) {
	if(amplicons == NULL) {
		return;
	}
	if(amplicons->link == amplicons) {
		free(amplicons);
		return;
	}
	pthread_mutex_lock(&pm_amp);
	AmpliconLink p = amplicons->link;
	amplicons->link = fullAmplicons->link;
	fullAmplicons->link = p->link;
	fullAmplicons = amplicons;
	Amplicon& amplicon = fullAmplicons->link->amplicon;
	amplicon.setData(amplicon.getData()+p->amplicon.getData()); // update full amplicon count
	free(p);
	pthread_mutex_unlock(&pm_amp);
}

void Malbac::createFrags() {
	genome.splitToFrags(fragments);
}

unsigned long Malbac::getAmpliconCount(AmpliconLink linkList) {
	return linkList->link->amplicon.getData();
}

void Malbac::amplify() {
	cerr << "\nMALBAC amplification..." << endl;
	
	double gamma = config.getRealPara("gamma");
	if(gamma <= 0 || gamma > 0.001) {
		unsigned long refLen = genome.getGenomeLength()/2;
		gamma = 3.0e-7*refLen/1000000;
		config.setRealPara("gamma", gamma);
		cerr << "the value of parameter \"gamma\" was set to " << gamma << endl;
	}
	
	createPrimers();
	setPrimers(true);
	amplifyFrags();
	//cerr << getAmpliconCount(semiAmplicons) << endl;
	for(int i = 0; i < 5; i++) {
		if(totalPrimers == 0) {
			break;
		}
		cerr << "cycle number: " << i+1 << endl;
		setPrimers(false);
		amplifySemiAmplicons();
		//cerr << getAmpliconCount(fullAmplicons) << endl;
		cerr << "semi amplicon amplification done!" << endl;
		if(i < 4) {
			amplifyFrags();
			//cerr << getAmpliconCount(semiAmplicons) << endl;
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
	createPrimers();
	setPrimers(true);
	cerr << "amplifyFrags" << endl;
	amplifyFrags();
	for(int i = 0; i < 5; i++) {
		if(totalPrimers == 0) {
			break;
		}
		cerr << "cycle number: " << i+1 << endl;
		setPrimers(false);
		amplifySemiAmplicons();
		cerr << getAmpliconCount(fullAmplicons) << endl;
		cerr << "semi amplicon amplification done!" << endl;
		saveFullAmplicons(ofs);
		if(i < 4) {
			amplifyFrags();
			cerr << getAmpliconCount(semiAmplicons) << endl;
			cerr << "fragment amplification done!" << endl;
		}
	}
	
	ofs.close();
}

void Malbac::setPrimers(bool onlyFrags) {
	unsigned long i, j, k;
	unsigned long templateNum = 0, fragNum = fragments.size();
	unsigned int length, gcContent;
	double totalLen = 0;
	
	for(i = 0; i < fragNum; i++) {
		length = fragments[i].getLength();
		totalLen += length;
	}
	templateNum += fragNum;
	
	if(!onlyFrags) {
		AmpliconLink p = semiAmplicons->link->link;
		templateNum += getAmpliconCount(semiAmplicons);
		while(p != semiAmplicons->link) {
			Amplicon& amplicon = p->amplicon; 
			length = amplicon.getLength();
			totalLen += length;
			p = p->link;
		}
	}
	
	double gamma = config.getRealPara("gamma");
	unsigned long usedPrimers = totalPrimers*gamma;
	unsigned long count = 0;
	for(i = 0; i < fragNum; i++) {
		length = fragments[i].getLength();
		k = 1.0*length/totalLen*usedPrimers;
		fragments[i].setPrimers(k);
		count += k;
	}
	if(!onlyFrags) {
		AmpliconLink p = semiAmplicons->link->link;
		while(p != semiAmplicons->link) {
			Amplicon& amplicon = p->amplicon;
			length = amplicon.getLength();
			k = 1.0*length/totalLen*usedPrimers;
			amplicon.setPrimers(k);
			count += k;
			p = p->link;
		}
	}
	long n = usedPrimers-count;
	while(n > 0) {
		j = threadPool->randomInteger(0, 50);
		if(j > n) j = n;
		i = threadPool->randomInteger(0, templateNum);
		if(i < fragNum) {
			fragments[i].setPrimers(fragments[i].getPrimers()+j);
		}
		else {
			i -= fragNum;
			AmpliconLink p = semiAmplicons->link;
			k = 0;
			while(k <= i) {
				p = p->link;
				k++;
			}
			Amplicon& amplicon = p->amplicon;
			amplicon.setPrimers(amplicon.getPrimers()+j);
		}
		n -= j;
	}
	totalPrimers -= usedPrimers;
}

void Malbac::saveFullAmplicons(ofstream& ofs) {
	unsigned long i = 1, j, k;
	unsigned long ampliconNum = getAmpliconCount(fullAmplicons);
	char c;
	char* seq;
	unsigned int sindx, length;
	int width = 100;
	AmpliconLink q, p = fullAmplicons->link->link;
	while(p != fullAmplicons->link) {
		ofs << ">" << "amp_" << i++ << endl;
		seq = p->amplicon.getSequence();
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
		q = p;
		p = p->link;
		free(q);
	}
	p->amplicon.setData(0);
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
	unsigned long semiAmpliconNum = getAmpliconCount(semiAmplicons);
	vector<AmpliconLink*> threadParas;
	int i = 0, j =0;
	unsigned long loadPerThread = max((unsigned long) 10, semiAmpliconNum/threadPool->getThreadNumber());
	AmpliconLink p = semiAmplicons->link->link;
	while(p != semiAmplicons->link) {
		AmpliconLink* pNodeRange = new AmpliconLink[2];
		pNodeRange[0] = p;
		i = 0;
		while(i < loadPerThread && p != semiAmplicons->link) {
			p = p->link;
			i++;
		}
		pNodeRange[1] = p;
		threadPool->pool_add_work(&Amplicon::batchAmplify, pNodeRange, j++);
		threadParas.push_back(pNodeRange);
	}
	threadPool->wait();
	for(i = 0; i < threadParas.size(); i++) {
		delete[] threadParas[i];
	}
}

void Malbac::setReadCounts(long reads) {
	int k = config.isPairedEnd()? 2:1;
	unsigned long i;
	unsigned long ampliconNum = getAmpliconCount(fullAmplicons);
	double WL = 0;
	Matrix<double> wls(1, ampliconNum);
	AmpliconLink p = fullAmplicons->link->link;
	i = 0;
	while(p != fullAmplicons->link) {
		double wl = p->amplicon.getWeightedLength();
		WL += wl;
		wls.set(0, i++, wl);
		p = p->link;
	}
	
	readNumbers = new int[ampliconNum];
	memset(readNumbers, 0, sizeof(int)*ampliconNum);
	long sum = 0;
	for(i = 0; i < ampliconNum; i++) {
		int readCount = (wls.get(0,i)/WL)*reads;
		readCount -= readCount%k;
		readNumbers[i] = readCount;
		sum += readCount;
	}
	
	//wls.normalize(0);
	//Matrix<double> aprob = wls.cumsum();
	wls.clear();
	
	reads -= sum;
	while(reads > 0) {
		i = threadPool->randomInteger(0, ampliconNum);
		//i = randIndx(aprob, true);
		readNumbers[i] += k;
		reads -= k;
	}
}

void Malbac::yieldReads() {	
	int i, j = 0;	
	unsigned long refLen = genome.getGenomeLength()/2;
	unsigned long reads = refLen*config.getIntPara("coverage")/config.getIntPara("readLength");
	if(config.isVerbose()) {	
		cerr << "\nNumber of reads to generate: " << reads << endl;
	}
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
	unsigned long ampliconNum = getAmpliconCount(fullAmplicons);
	vector<AmpliconLink*> threadParas;
	unsigned long loadPerThread = max((unsigned long) 10, ampliconNum/threadPool->getThreadNumber());
	AmpliconLink p = fullAmplicons->link->link;
	while(p != fullAmplicons->link) {
		AmpliconLink* pNodeRange = new AmpliconLink[2];
		pNodeRange[0] = p;
		i = 0;
		while(i < loadPerThread && p != fullAmplicons->link) {
			p = p->link;
			i++;
		}
		pNodeRange[1] = p;
		threadPool->pool_add_work(&Amplicon::yieldReads, pNodeRange, j++);
		threadParas.push_back(pNodeRange);
	}
	threadPool->wait();
	for(i = 0; i < threadParas.size(); i++) {
		delete[] threadParas[i];
	}
	
	delete swp;
}
