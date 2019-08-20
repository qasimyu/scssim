// ***************************************************************************
// Malbac.h (c) 2018 Zhenhua Yu <qasim0208@163.com>
// Health Informatics Lab, Ningxia University
// All rights reserved.

#ifndef _MALBAC_H
#define _MALBAC_H

#include <vector>
#include <string>
#include <pthread.h>

#include "Fragment.h"
#include "Amplicon.h"

using namespace std;

class PrimerIndex {
	public:
		int index;
		long count;
		map<char, PrimerIndex> nextIndexs;
		PrimerIndex() {index = -1;}
};

class Malbac {
	private:
		double refGCcontent;
		unsigned long totalPrimers;
		char **primers;
		int ptypeCount;
		PrimerIndex rIndex;
		
		vector<Fragment> fragments;
		//vector<Amplicon> semiAmplicons;
		//vector<Amplicon> fullAmplicons;
		AmpliconLink semiAmplicons;
		AmpliconLink fullAmplicons;
		unsigned int* readNumbers;
		
		mutable pthread_mutex_t pm_amp, pm_primer;
		
		void createPrimers();
		
		void calGCOfRef();
		void calGCOfProducts();
		
		void amplifyFrags();
		void amplifySemiAmplicons();
		
		void setPrimers(bool onlyFrags);
		void saveFullAmplicons(ofstream& ofs);
		
		void setReadCounts(long reads);
		
	public:
		Malbac();
		~Malbac();
		
		long getPrimerCount(const char *s);
		int updatePrimerCount(const char *s, int n);
		
		vector<Fragment>& getFrags() {return fragments;}
		Fragment& getFrag(unsigned long i) {return fragments[i];}
		//vector<Amplicon>& getSemiAmplicons() {return semiAmplicons;}
		AmpliconLink getSemiAmplicons() {return semiAmplicons;}
		//Amplicon& getSemiAmplicon(unsigned long i) {return semiAmplicons[i];}
		//vector<Amplicon>& getFullAmplicons() {return fullAmplicons;}
		AmpliconLink getFullAmplicons() {return fullAmplicons;}
		unsigned long getAmpliconCount(AmpliconLink linkList);
		
		unsigned int* getReadNumbers() {return readNumbers;}
		
		void extendSemiAmplicons(AmpliconLink amplicons);
		void extendFullAmplicons(AmpliconLink amplicons);
		
		void createFrags();
		void amplify();
		void amplifyAndSaveProducts();
		
		void yieldReads();
};


#endif

