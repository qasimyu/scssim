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

class Malbac {
	private:
		vector<Fragment> fragments;
		vector<Amplicon> semiAmplicons;
		vector<Amplicon> fullAmplicons;
		vector<short int> readNumbers;
		
		mutable pthread_mutex_t pm;
		
		void amplifyFrags();
		void amplifySemiAmplicons();
		
		void setPrimers(bool onlyFrags);
		void saveFullAmplicons(ofstream& ofs);
		
		void setReadCounts(unsigned long reads);
		
	public:
		Malbac();
		~Malbac() {}
		
		vector<Fragment>& getFrags() {return fragments;}
		Fragment& getFrag(unsigned long i) {return fragments[i];}
		vector<Amplicon>& getSemiAmplicons() {return semiAmplicons;}
		Amplicon& getSemiAmplicon(unsigned long i) {return semiAmplicons[i];}
		vector<Amplicon>& getFullAmplicons() {return fullAmplicons;}
		vector<short int>& getReadNumbers() {return readNumbers;}
		
		void extendSemiAmplicons(vector<Amplicon>& amplicons);
		void extendFullAmplicons(vector<Amplicon>& amplicons);
		
		void createFrags();
		void amplify();
		void amplifyAndSaveProducts();
		
		void yieldReads();
};


#endif

