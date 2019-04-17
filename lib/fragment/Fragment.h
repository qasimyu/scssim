// ***************************************************************************
// Fragment.h (c) 2018 Zhenhua Yu <qasim0208@163.com>
// Health Informatics Lab, Ningxia University
// All rights reserved.

#ifndef _FRAGMENT_H
#define _FRAGMENT_H

#include <cassert>
#include <vector>
#include <string>
#include <pthread.h>

#include "Amplicon.h"

using namespace std;

class Fragment {
	private:
		string chr;
		long startPos;
		int length;
		int strand;
		
		int gcContent;
		int primerNum;
		
		char *sequence;
		
	public:	
		static int minSize, maxSize;
		Fragment(string chr, long startPos, int length, int strand);
		~Fragment();

		static void setMaxSize(unsigned int size) {maxSize = size;}
		static int getMaxSize() {return maxSize;}
		static void setMinSize(unsigned int size) {minSize = size;}
		static int getMinSize() {return minSize;}
		
		string getChr() {return chr;}
		long getStartPos() {return startPos;}
		int getLength() {return length;}
		void setLength(int length) {this->length = length;}
		int getStrand() {return strand;}
		int getGCcontent() {return gcContent;}
		
		void describe();
		
		void setPrimers(int primers) {primerNum = primers;}
		
		void createSequence();
		char* getSequence() {return sequence;}
		
		void amplify(AmpliconLink& results);
		static void* batchAmplify(const void* args);
};


#endif

