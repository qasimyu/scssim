// ***************************************************************************
// Amplicon.cpp (c) 2018 Zhenhua Yu <qasim0208@163.com>
// Health Informatics Lab, Ningxia University
// All rights reserved.

#include <iostream>
#include <algorithm>

#include "MyDefine.h"
#include "Fragment.h"
#include "Amplicon.h"

AmpError::AmpError(unsigned char etype, unsigned int pos, unsigned char alt) {
	data[0] = (etype << 6) | ((pos >> 21) & 0x0000003F);
	data[1] = (pos >> 13) & 0x000000FF;
	data[2] = (pos >> 5) & 0x000000FF;
	data[3] = ((pos & 0x0000001F) << 3) | alt;
}

unsigned char AmpError::getErrType() {
	return data[0] >> 6;
}

unsigned int AmpError::getPos() {
	unsigned int pos = data[0];
	pos = (pos & 0x0000003F) << 21;
	unsigned int tmp = data[1];
	tmp = tmp << 13;
	pos |= tmp;
	tmp = data[2];
	tmp = tmp << 5;
	pos |= tmp;
	tmp = data[3];
	tmp = tmp >> 3;
	pos |= tmp;
	return pos;
}

unsigned char AmpError::getAlt() {
	return data[3] & 0x07;
}

void AmpError::setAlt(unsigned char alt) {
	data[3] = (data[3] & 0xF8) | (alt & 0x07);
}

Amplicon::Amplicon(bool isSemi, unsigned long tmplIndx, AmpError *ampErrs, unsigned int startPos, unsigned int length, unsigned int gcContent) {
//Amplicon::Amplicon(bool isSemi, void *tmpl, AmpError *ampErrs, unsigned int startPos, unsigned int length, unsigned int gcContent) {
	//assert(tmpl != NULL);
	
	data[0] = isSemi? 0x80:0;
	data[0] |= (startPos >> 10) & 0x0000007F;
	data[1] = (startPos >> 2) & 0x000000FF;
	data[2] = ((startPos & 0x00000003) << 6) | ((length >> 11) & 0x0000003F);
	data[3] = (length >> 3) & 0x000000FF;
	data[4] = ((length & 0x00000007) << 5) | ((gcContent >> 12) & 0x0000001F);
	data[5] = (gcContent >> 4) & 0x000000FF;
	data[6] = (gcContent & 0x0000000F) << 4;
	
	//this->tmpl = tmpl;
	this->tmplIndx = tmplIndx;
	this->ampErrs = ampErrs;
}

void Amplicon::clear() {
	delete[] ampErrs;
	ampErrs = NULL;
}

bool Amplicon::isSemi() {
	return data[0] & 0x80;
}

void Amplicon::setPrimers(unsigned short int primers) {
	data[6] = (data[6] & 0xE0) | ((primers >> 8) & 0x000F);
	data[7] = primers & 0x00FF;
}

unsigned short int Amplicon::getPrimers() {
	unsigned short int primers = data[6];
	primers = (primers & 0x000F) << 8;
	unsigned int tmp = data[7];
	primers |= tmp;
	return primers;
}

unsigned int Amplicon::getStartPos() {
	unsigned int startPos = data[0];
	startPos = (startPos & 0x0000007F) << 10;
	unsigned int tmp = data[1];
	tmp = tmp << 2;
	startPos |= tmp;
	tmp = data[2];
	tmp = tmp >> 6;
	startPos |= tmp;
	return startPos;
}

unsigned int Amplicon::getLength() {
	unsigned int length = data[2];
	length = (length & 0x0000003F) << 11;
	unsigned int tmp = data[3];
	tmp = tmp << 3;
	length |= tmp;
	tmp = data[4];
	tmp = tmp >> 5;
	length |= tmp;
	return length;
}

void Amplicon::setLength(unsigned int length) {
	data[2] = (data[2] & 0xC0) | ((length >> 11) & 0x0000003F);
	data[3] = (length >> 3) & 0x000000FF;
	data[4] = (data[4] & 0x1F) | ((length & 0x00000007) << 5);
}

unsigned int Amplicon::getGCcontent() {
	unsigned int gcContent = data[4];
	gcContent = (gcContent & 0x0000001F) << 12;
	unsigned int tmp = data[5];
	tmp = tmp << 4;
	gcContent |= tmp;
	tmp = data[6];
	tmp = tmp >> 4;
	gcContent |= tmp;
	return gcContent;
}

void Amplicon::setGCcontent(unsigned int gcContent) {
	data[4] = (data[4] & 0xE0) | ((gcContent >> 12) & 0x0000001F);
	data[5] = (gcContent >> 4) & 0x000000FF;
	data[6] = (data[6] & 0x0F) | ((gcContent & 0x0000000F) << 4);
}

void Amplicon::amplify(vector<Amplicon>& fullAmplicons, unsigned long tmplIndx) {
	unsigned int i, j, k, n;
	unsigned int spos, ampliconLen;
	double ador = config.getRealPara("ador");
	double fpr = config.getRealPara("fpr");
	char* semiAmpSeq = getSequence();
	string bases = config.getStringPara("bases");
	unsigned int length = getLength();
	unsigned short int primerNum = getPrimers();
	int maxLen = config.getIntPara("ampliconMaxLen");
	int minLen = config.getIntPara("ampliconMinLen");
	short int* posAttached = new short int[length];
	memset(posAttached, 0, length*sizeof(unsigned short int));
	
	for(i = 0; i < primerNum; i++) {
		int tryTimes = 0;
		spos = 0;
		do {
			//spos = randomInteger(0, length-34);
			ampliconLen = threadPool->randomDouble(minLen, maxLen+1);
			tryTimes++;
			if(tryTimes > 50) {
				break;
			}
		} while(spos+ampliconLen > length || posAttached[spos+ampliconLen-1] == 1);
		if(tryTimes > 50) {
			break;
		}
		
		posAttached[spos+ampliconLen-1] = 1;
		
		char c = semiAmpSeq[spos+ampliconLen];
		semiAmpSeq[spos+ampliconLen] = '\0';
		unsigned int gcNum = countGC(semiAmpSeq+spos);
		semiAmpSeq[spos+ampliconLen] = c;
		
		vector<AmpError> errs;
		k = 0;
		for(j = 0; j < ampliconLen-8; j++) {
			double p = threadPool->randomDouble(0, 1);
			if(p <= ador) {
				//AmpError ampErr(ADO, j, 0);
				AmpError ampErr(0, j, 0);
				errs.push_back(ampErr);
				k++;
				if(semiAmpSeq[spos+j] == 'C' || semiAmpSeq[spos+j] == 'G') {
					gcNum--;
				}
				continue;
			}
			p = threadPool->randomDouble(0, 1);
			if(p <= fpr) {
				char base = getComplementBase(semiAmpSeq[spos+j]);
				do {
					n = threadPool->randomDouble(0, bases.size());
				} while(bases[n] == base);
				
				if(bases[n] == 'C' || bases[n] == 'G') {
					gcNum++;
				}
				if(base == 'C' || base == 'G') {
					gcNum--;
				}
				//AmpError ampErr(FP, j, bases[n]);
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
		//Amplicon tmp(false, this, ampErrs, spos, ampliconLen-k, gcNum);
		Amplicon tmp(false, tmplIndx, ampErrs, spos, ampliconLen-k, gcNum);
		fullAmplicons.push_back(tmp);
	}
	
	delete[] semiAmpSeq;
	delete[] posAttached;
}

void* Amplicon::batchAmplify(const void* args) {
	unsigned long* indexs = (unsigned long*) args;
	unsigned long sindx = indexs[0];
	unsigned long eindx = indexs[1];
	unsigned long i;
	vector<Amplicon>& semiAmplicons = malbac.getSemiAmplicons();
	vector<Amplicon> results;
	
	for(i = sindx; i <= eindx; i++) {
		//semiAmplicons[i].amplify(results);
		semiAmplicons[i].amplify(results, i);
	}
	
	malbac.extendFullAmplicons(results);
}

char* Amplicon::getSequence() {
	unsigned int i, j, k;
	unsigned int n, m, sindx;
	unsigned int semiLength, length;
	unsigned int spos;
	string bases = config.getStringPara("bases");
	if(!isSemi()) {
		//Amplicon *ap = static_cast<Amplicon*>(tmpl);
		Amplicon& semiAmp = malbac.getSemiAmplicon(tmplIndx);
		//Fragment *fp = (Fragment*) (ap->getTmpl());
		Fragment& frag = malbac.getFrag(semiAmp.tmplIndx);
		//frag.describe();
		char* fragSeq = frag.getSequence();
		char* semiSeq_o = new char[strlen(fragSeq)+1];
		strcpy(semiSeq_o, fragSeq);
		semiSeq_o = getComplementSeq(semiSeq_o);
		
		AmpError* errs = semiAmp.getErrs();
		semiLength = semiAmp.getLength();
		char* semiSeq = new char[semiLength+1];
		semiSeq[semiLength] = '\0';
		sindx = 0;
		spos = semiAmp.getStartPos();
		m = spos;
		if(errs) {
			k = errs->getAlt();
			while(k != bases.size()) {
				j = errs->getPos();
				while(m < j+spos) {
					semiSeq[sindx++] = semiSeq_o[m++];
				}
				if(errs->getErrType() == 0) {
					m++;
				}
				else {
					semiSeq_o[j+spos] = bases[k];
				}
				errs++;
				k = errs->getAlt();
			}			
		}
		while(sindx < semiLength) {
			semiSeq[sindx++] = semiSeq_o[m++];
		}
		delete[] semiSeq_o;
		
		// convert to 5'->3'
		reverse(semiSeq, semiSeq+strlen(semiSeq));
		char* fullSeq = getComplementSeq(semiSeq);
		
		length = getLength();
		char* ret = new char[length+1];
		ret[length] = '\0';
		sindx = 0;
		spos = getStartPos();
		m = spos;
		errs = ampErrs;
		if(errs) {
			k = errs->getAlt();
			while(k != bases.size()) {
				j = errs->getPos();
				while(m < j+spos) {
					ret[sindx++] = fullSeq[m++];
				}
				if(errs->getErrType() == 0) {
					m++;
				}
				else {
					fullSeq[j+spos] = bases[k];
				}
				errs++;
				k = errs->getAlt();
			}			
		}
		while(sindx < length) {
			ret[sindx++] = fullSeq[m++];
		}
		delete[] fullSeq;
		
		reverse(ret, ret+strlen(ret));
		return ret;
	}
	else {
		//Fragment *fp = (Fragment*) tmpl;
		Fragment& frag = malbac.getFrag(tmplIndx);
		char* fragSeq = frag.getSequence();
		//char* semiSeq_o = getComplementSeq(fragSeq);
		char* semiSeq_o = new char[strlen(fragSeq)+1];
		strcpy(semiSeq_o, fragSeq);
		semiSeq_o = getComplementSeq(semiSeq_o);
		
		length = getLength();
		char* ret = new char[length+1];
		ret[length] = '\0';
		sindx = 0;
		spos = getStartPos();
		m = spos;
		AmpError* errs = ampErrs;
		
		if(errs) {
			k = errs->getAlt();
			while(k != bases.size()) {
				j = errs->getPos();
				while(m < j+spos) {
					ret[sindx++] = semiSeq_o[m++];
				}
				if(errs->getErrType() == 0) {
					m++;
				}
				else {
					semiSeq_o[j+spos] = bases[k];
				}
				errs++;
				k = errs->getAlt();
			}
		}
		while(sindx < length) {
			ret[sindx++] = semiSeq_o[m++];
		}
		/*
		if(errs) {
			k = errs->getAlt();
			while(k != bases.size()) {
				j = errs->getPos();
				if(errs->getErrType() == 0) {
					semiSeq_o[j+spos] = 'X';
				}
				else {
					semiSeq_o[j+spos] = bases[k];
				}
				errs++;
				k = errs->getAlt();
			}
		}
		for(i = 0, j = 0; i < length; j++) {
			if(semiSeq_o[j+spos] != 'X') {
				ret[i++] = semiSeq_o[j+spos];
			}
		}
		*/
		delete[] semiSeq_o;
		
		//cerr << "semi:(" << spos << ") " << ret << endl;
		
		reverse(ret, ret+strlen(ret));
		return ret;
	}	
}

double Amplicon::getWeightedLength() {
	int i, k, spos, gc;
	char c;
	char *s;
	double weight, weightLen = 0;
	unsigned int fragSize = Profile::fragSize;
	char *ampliconSeq = getSequence();
	unsigned int ampliconSize = strlen(ampliconSeq);
	k = ampliconSize/fragSize;
	for(i = 0; i < k; i++) {
		spos = i*fragSize;
		c = ampliconSeq[spos+fragSize];
		ampliconSeq[spos+fragSize] = '\0';
		s = &ampliconSeq[spos];
		gc = calculateGCPercent(s);
		ampliconSeq[spos+fragSize] = c;
		
		weight = profile.getGCFactor(gc)/fragSize;
		weightLen += weight;
	}
	if(k*fragSize < ampliconSize) {
		spos = k*fragSize;
		s = &ampliconSeq[spos];
		gc = calculateGCPercent(s);
		
		weight = profile.getGCFactor(gc)*(ampliconSize-k*fragSize)/(fragSize*fragSize);
		weightLen += weight;
	}
	
	delete[] ampliconSeq;
	
	return weightLen;
}

void* Amplicon::yieldReads(const void* args) {
	unsigned long* indexs = (unsigned long*) args;
	unsigned long sindx = indexs[0];
	unsigned long eindx = indexs[1];
	unsigned long i;
	int j, k, n, fragCount, failCount, seqLen;
	vector<Amplicon>& amplicons = malbac.getFullAmplicons();
	vector<short int>& readNumbers = malbac.getReadNumbers();
	
	char *seq, *ampliconSeq, *fragSeq;
	bool paired = config.isPairedEnd();
	int ampliconLen, readLength = config.getIntPara("readLength");
	
	unsigned long bufferSize = 50000000;
	char buf[10*readLength], buf1[10*readLength], buf2[10*readLength];
	char *results;
	char *outBuffer, *outBuffer1, *outBuffer2;
	unsigned long outIndx = 0, outIndx1 = 0, outIndx2 = 0;
	if(paired) {
		outBuffer1 = new char[bufferSize];
		outBuffer2 = new char[bufferSize];
	}
	else {
		outBuffer = new char[bufferSize];
	}
	
	long pos;
	for(i = sindx; i <= eindx; i++) {
		ampliconSeq = amplicons[i].getSequence();
		ampliconLen = strlen(ampliconSeq);
		if(ampliconLen < readLength) {
			continue;
		}
		n = readNumbers[i];
		fragCount = 0;
		failCount = 0;
		while(n > 0) {
			fragCount++;
			if(!paired) {
				pos = threadPool->randomInteger(0, ampliconLen-readLength+1);
				char c = ampliconSeq[pos+readLength];
				ampliconSeq[pos+readLength] = '\0';
				results = profile.predict(&ampliconSeq[pos], 1);
				ampliconSeq[pos+readLength] = c;
				
				k = 0;
				sprintf(&buf[k], "@%d_%d\n", i, fragCount);
				k = strlen(buf);
				seqLen = strlen(results)/2;
				strncpy(&buf[k], &results[0], seqLen);
				strcpy(&buf[k+seqLen], "\n+\n");
				strncpy(&buf[k+seqLen+3], &results[seqLen], seqLen);
				buf[k+2*seqLen+3] = '\n';
				buf[k+2*seqLen+4] = '\0';
				delete[] results;
				
				if(strlen(buf)+outIndx < bufferSize) {
					strcpy(&outBuffer[outIndx], buf);
					outIndx += strlen(buf);
				}
				else {
					swp->write(outBuffer);
					strcpy(outBuffer, buf);
					outIndx = strlen(buf);
				}
				
				n--;
			}
			else {
				int insertSize = profile.yieldInsertSize();
				if(insertSize < readLength || insertSize > ampliconLen) {
					failCount++;
					if(failCount > 1000) {
						break;
					}
					continue;
				}
				pos = threadPool->randomInteger(0, ampliconLen-insertSize+1);
				char c = ampliconSeq[pos+readLength];
				ampliconSeq[pos+readLength] = '\0';
				results = profile.predict(&ampliconSeq[pos], 1);
				ampliconSeq[pos+readLength] = c;
				
				k = 0;
				sprintf(&buf1[k], "@%d_%d/1\n", i, fragCount);
				k = strlen(buf1);
				seqLen = strlen(results)/2;
				strncpy(&buf1[k], &results[0], seqLen);
				strcpy(&buf1[k+seqLen], "\n+\n");
				strncpy(&buf1[k+seqLen+3], &results[seqLen], seqLen);
				buf1[k+2*seqLen+3] = '\n';
				buf1[k+2*seqLen+4] = '\0';
				delete[] results;
				
				char* seq = &ampliconSeq[pos+insertSize-readLength];
				c = ampliconSeq[pos+insertSize];
				ampliconSeq[pos+insertSize] = '\0';
				seq = getComplementSeq(seq);
				reverse(seq, seq+readLength);
				results = profile.predict(seq, 0);
				ampliconSeq[pos+insertSize] = c;
				
				k = 0;
				sprintf(&buf2[k], "@%d_%d/2\n", i, fragCount);
				k = strlen(buf2);
				seqLen = strlen(results)/2;
				strncpy(&buf2[k], &results[0], seqLen);
				strcpy(&buf2[k+seqLen], "\n+\n");
				strncpy(&buf2[k+seqLen+3], &results[seqLen], seqLen);
				buf2[k+2*seqLen+3] = '\n';
				buf2[k+2*seqLen+4] = '\0';
				delete[] results;
				
				if(strlen(buf1)+outIndx1 < bufferSize && strlen(buf2)+outIndx2 < bufferSize) {
					strcpy(&outBuffer1[outIndx1], buf1);
					strcpy(&outBuffer2[outIndx2], buf2);
					outIndx1 += strlen(buf1);
					outIndx2 += strlen(buf2);
				}
				else {
					swp->write(outBuffer1, outBuffer2);
					strcpy(outBuffer1, buf1);
					strcpy(outBuffer2, buf2);
					outIndx1 = strlen(buf1);
					outIndx2 = strlen(buf2);
				}
				
				n -= 2;
			}
		}
		delete[] ampliconSeq;
	}
	
	if(paired) {
		if(strlen(outBuffer1) > 0) {
			swp->write(outBuffer1, outBuffer2);
		}
		delete[] outBuffer1;
		delete[] outBuffer2;
	}
	else {
		if(strlen(outBuffer) > 0) {
			swp->write(outBuffer);
		}
		delete[] outBuffer;
	}
	
	return NULL;
}

