// ***************************************************************************
// Fragment.cpp (c) 2018 Zhenhua Yu <qasim0208@163.com>
// Health Informatics Lab, Ningxia University
// All rights reserved.

#include <cstdio>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <unistd.h>

#include "MyDefine.h"
#include "Fragment.h"

int Fragment::maxSize = 100000;
int Fragment::minSize = 10000;

Fragment::Fragment(string chr, long startPos, int length, int strand) {
	assert(!chr.empty());
	assert(startPos >= 0);
	//assert(length >= minSize && length <= maxSize);
	assert(strand == 1 || strand == -1);
	this->chr = chr;
	this->startPos = startPos;
	this->length = length;
	this->strand = strand;
	
	this->sequence = NULL;
}

Fragment::~Fragment() {
	delete[] sequence;
}

void Fragment::describe() {
	cerr << "chr=" << chr << ", spos=" << startPos 
		<< ", length=" << length << ", strand=" << strand << endl;
}

void Fragment::createSequence() {
	if(strand == 1) {
		sequence = genome.getSubSequence(chr, startPos-1, length);
	}
	else {
		long chrLen = genome.getChromLen(chr);
		sequence = genome.getSubSequence(chr, chrLen-startPos-length+1, length);
		sequence = getComplementSeq(sequence);
		reverse(sequence, sequence+strlen(sequence));
	}
	gcContent = countGC(sequence);
}

void Fragment::amplify(vector<Amplicon>& semiAmplicons, unsigned long tmplIndx) {
	unsigned int i, j, k, n;
	unsigned int spos, ampliconLen, curSize;
	double ador = config.getRealPara("ador");
	double fpr = config.getRealPara("fpr");
	char* fragSeq = getSequence();
	string bases = config.getStringPara("bases");
	int maxLen = config.getIntPara("ampliconMaxLen");
	int minLen = config.getIntPara("ampliconMinLen");
	short int* posAttached = new short int[length];
	memset(posAttached, 0, length*sizeof(short int));
	
	for(i = 0; i < primerNum; i++) {
		int tryTimes = 0;
		do {
			spos = threadPool->randomInteger(0, length-34);
			ampliconLen = threadPool->randomInteger(minLen, maxLen+1);
			tryTimes++;
			if(tryTimes > 50) {
				break;
			}
		} while(spos+ampliconLen > length || posAttached[spos+ampliconLen-1] == 1);
		if(tryTimes > 50) {
			break;
		}
		
		posAttached[spos+ampliconLen-1] = 1;

		char c = fragSeq[spos+ampliconLen];
		fragSeq[spos+ampliconLen] = '\0';
		int gcNum = countGC(fragSeq+spos);
		fragSeq[spos+ampliconLen] = c;
		
		vector<AmpError> errs;
		k = 0;
		for(j = 0; j < ampliconLen-8; j++) {
			double p = threadPool->randomDouble(0, 1);
			if(p <= ador) {
				AmpError ampErr(0, j, 0);
				errs.push_back(ampErr);
				k++;
				if(fragSeq[spos+j] == 'C' || fragSeq[spos+j] == 'G') {
					gcNum--;
				}
				continue;
			}
			p = threadPool->randomDouble(0, 1);
			if(p <= fpr) {
				char base = getComplementBase(fragSeq[spos+j]);
				do {
					n = threadPool->randomInteger(0, bases.size());
				} while(bases[n] == base);
				
				if(bases[n] == 'C' || bases[n] == 'G') {
					gcNum++;
				}
				if(base == 'C' || base == 'G') {
					gcNum--;
				}
				
				AmpError ampErr(1, j, n);
				errs.push_back(ampErr);
			}
		}
		
		AmpError *ampErrs = NULL;
		if(!errs.empty()) {
			ampErrs = new AmpError[errs.size()+1];
			for(j = 0; j < errs.size(); j++) {
				memcpy(ampErrs[j].getData(), errs[j].getData(), 4);
			}
			ampErrs[j].setAlt(bases.size());
		}
		//Amplicon tmp(true, this, ampErrs, spos, ampliconLen-k, gcNum);
		Amplicon tmp(true, tmplIndx, ampErrs, spos, ampliconLen-k, gcNum);
		semiAmplicons.push_back(tmp);
	}
	//delete[] fragSeq;
	delete[] posAttached;
}

void* Fragment::batchAmplify(const void* args) {
	unsigned long* indexs = (unsigned long*) args;
	unsigned long sindx = indexs[0];
	unsigned long eindx = indexs[1];
	unsigned long i;
	vector<Fragment>& frags = malbac.getFrags();
	vector<Amplicon> results;
	
	for(i = sindx; i <= eindx; i++) {
		//frags[i].amplify(results);
		frags[i].amplify(results, i);
	}
	
	malbac.extendSemiAmplicons(results);
}

