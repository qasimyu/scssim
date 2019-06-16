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
	sequence = genome.getSubSequence(chr, startPos-1, length);
	if(strand == 1) {
		reverse(sequence, sequence+strlen(sequence));
	}
	else {
		sequence = getComplementSeq(sequence);
	}
	gcContent = countGC(sequence);
	
}

void Fragment::amplify(AmpliconLink& results) {
	unsigned int i, j, k, n;
	unsigned int spos, ampliconLen, curSize;
	double ber = config.getRealPara("ber");
	
	string bases = config.getStringPara("bases");
	int maxLen = config.getIntPara("ampliconMaxLen");
	int minLen = config.getIntPara("ampliconMinLen");
	
	if(length < minLen+27) {
		return;
	}
	
	char* fragSeq = getSequence();
	char* fragSeq_c = new char[strlen(fragSeq)+1];
	strcpy(fragSeq_c, fragSeq);
	fragSeq_c = getComplementSeq(fragSeq_c);
	
	short int* posAttached = new short int[length];
	memset(posAttached, 0, length*sizeof(short int));
	
	for(i = 0; i < primerNum; i++) {
		int tryTimes = 0;
		do {
			spos = threadPool->randomInteger(27, length);
			ampliconLen = threadPool->randomDouble(minLen, maxLen+1);
			tryTimes++;
			if(tryTimes > 50) {
				break;
			}
			if(spos+ampliconLen > length || posAttached[spos] == 1) {
				continue;
			}
			char c = fragSeq_c[spos+8];
			fragSeq_c[spos+8] = '\0';
			j = malbac.updatePrimerCount(&fragSeq_c[spos], -1);
			fragSeq_c[spos+8] = c;
			if(j == 1) {
				break;
			}
		} while(1);
		if(tryTimes > 50) {
			break;
		}
		
		posAttached[spos] = 1;
		
		char c = fragSeq_c[spos+ampliconLen];
		fragSeq_c[spos+ampliconLen] = '\0';
		int gcNum = countGC(fragSeq_c+spos);
		fragSeq_c[spos+ampliconLen] = c;	
		
		vector<AmpError> errs;
		for(j = 8; j < ampliconLen; j++) {
			double p = threadPool->randomDouble(0, 1);
			if(p < ber) {
				char base = fragSeq_c[spos+j];
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
		Amplicon tmp(true, this, ampErrs, spos, ampliconLen, max(0, gcNum));
		insertLinkList(results, tmp);
	}
	delete[] fragSeq_c;
	delete[] posAttached;
}

void* Fragment::batchAmplify(const void* args) {
	unsigned long* indexs = (unsigned long*) args;
	unsigned long sindx = indexs[0];
	unsigned long eindx = indexs[1];
	unsigned long i;
	vector<Fragment>& frags = malbac.getFrags();
	AmpliconLink results = createNullLink();
	
	for(i = sindx; i <= eindx; i++) {
		frags[i].amplify(results);
	}
	
	malbac.extendSemiAmplicons(results);
}

