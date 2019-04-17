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

struct AmpliconNode;
typedef AmpliconNode * AmpliconLink;

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
		void *tmpl;
		AmpError *ampErrs;
		char *sequence;
		
	public:	
		Amplicon() {ampErrs = NULL; sequence = NULL;}
		Amplicon(bool isSemi, void *tmpl, AmpError *ampErrs, unsigned int startPos, unsigned int length, unsigned int gcContent);
		
		void clear();
		
		bool isSemi();
		unsigned int getStartPos();
		unsigned int getLength();
		void setLength(unsigned int length);
		unsigned int getGCcontent();
		void setGCcontent(unsigned int gcContent);
		void setPrimers(unsigned short int primers);
		unsigned short int getPrimers();
		
		unsigned int getData();
		void setData(unsigned int n);
		void* getTmpl() {return tmpl;}
		AmpError* getErrs() {return ampErrs;}
		void setErrs(AmpError *ampErrs) {this->ampErrs = ampErrs;}
		void setSequence(char *seq) {sequence = seq;}
		char* getSequence();
		static void* batchGetSequences(const void* args);
		
		void amplify(AmpliconLink& results, AmpliconLink source);
		static void* batchAmplify(const void* args);
		
		double getWeightedLength();
		static void* yieldReads(const void* args);
};

struct AmpliconNode {
	Amplicon amplicon;
	struct AmpliconNode *link;
};

AmpliconLink createNullLink();

void insertLinkList(AmpliconLink& linkList, Amplicon& amplicon);

void printLinkList(AmpliconLink linkList);

#endif

