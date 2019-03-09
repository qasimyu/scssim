// ***************************************************************************
// Amplicon.h (c) 2018 Zhenhua Yu <qasim0208@163.com>
// Health Informatics Lab, Ningxia University
// All rights reserved.

#ifndef _AMPLICON_H
#define _AMPLICON_H

#include <cassert>
#include <vector>
#include <string>

using namespace std;

//enum errorType {ADO, FP};

class AmpError {
	private:
		unsigned char data[4]; // bits: etype(2), pos(27), alt(3)
		
	public:
		AmpError() {}
		AmpError(unsigned char etype, unsigned int pos, unsigned char alt);
		unsigned char getErrType();
		unsigned int getPos();
		unsigned char getAlt();
		void setAlt(unsigned char alt);
		unsigned char* getData() {return data;}
};


class Amplicon {
	private:
		unsigned char data[8]; // bits: isSemi(1), startPos(17), length(17), gcContent(17), primerNum(12)
		//void *tmpl;
		unsigned long tmplIndx;
		AmpError *ampErrs;
		
	public:	
		Amplicon() {ampErrs = NULL;}
		//Amplicon(bool isSemi, void *tmpl, AmpError *ampErrs, unsigned int startPos, unsigned int length, unsigned int gcContent);
		Amplicon(bool isSemi, unsigned long tmplIndx, AmpError *ampErrs, unsigned int startPos, unsigned int length, unsigned int gcContent);
		
		void clear();
		
		bool isSemi();
		unsigned int getStartPos();
		unsigned int getLength();
		void setLength(unsigned int length);
		unsigned int getGCcontent();
		void setGCcontent(unsigned int gcContent);
		void setPrimers(unsigned short int primers);
		unsigned short int getPrimers();
		
		//void* getTmpl() {return tmpl;}
		AmpError* getErrs() {return ampErrs;}
		void setErrs(AmpError *ampErrs) {this->ampErrs = ampErrs;}
		
		//void amplify(vector<Amplicon>& fullAmplicons);
		void amplify(vector<Amplicon>& fullAmplicons, unsigned long tmplIndx);
		static void* batchAmplify(const void* args);
		char* getSequence();
		
		double getWeightedLength();
		static void* yieldReads(const void* args);
};


#endif

