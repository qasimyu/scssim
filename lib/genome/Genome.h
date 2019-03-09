// ***************************************************************************
// Genome.h (c) 2018 Zhenhua Yu <qasim0208@163.com>
// Health Informatics Lab, Ningxia University
// All rights reserved.

#ifndef _GENOME_H
#define _GENOME_H

#include <vector>
#include <map>
#include <string>

#include "Fragment.h"
#include "snp.h"
#include "Fasta.h"

using namespace std;

class CNV {
	public:
		CNV() {}
		CNV(long spos, long epos, float CN, float mCN)
			: spos(spos), epos(epos), CN(CN), mCN(mCN) {}
		~CNV() {}
		long spos;
		long epos;
		float CN;
		float mCN;
};

class SNV {
	public:
		SNV() {}
		SNV(long pos, char ref, char alt, string type)
			: pos(pos), ref(ref), alt(alt), type(type) {}
		~SNV() {}
		long pos;
		char ref;
		char alt;
		string type;
};

class Insert {
	public:
		Insert() {}
		Insert(long pos, string seq) : pos(pos), seq(seq), type(type) {}
		~Insert() {}
		long pos;
		string seq;
		string type;
};

class Del {
	public:
		Del() {}
		Del(long pos, unsigned int len, string type) : pos(pos), len(len), type(type) {}
		~Del() {}
		long pos;
		unsigned int len;
		string type;
};

class Target {
	public:
		Target() {}
		Target(long spos, long epos)
			: spos(spos), epos(epos) {}
		~Target() {}
		long spos;
		long epos;
};

class Genome {
	private:
		vector<string> chromosomes;
		map<string, vector<CNV> > cnvs;
		map<string, vector<SNV> > snvs;
		map<string, vector<Insert> > inserts;
		map<string, vector<Del> > dels;
		
		map<string, vector<Target> > targets;
		
		SNPOnChr sc;	
		FastaReference fr;
		
		//map<string, string> altSequence;
		string altSequence;
		string curChr;
		
		void divideTargets();
		void generateAltSequence(string chr);
		
		void generateSegment(vector<string>& sequences, string chr, long segStartPos, long segEndPos, int CN, int mCN);
		
	public:
		Genome() {curChr = "";}
		~Genome() {}

		void loadData();
		void loadAbers();
		void loadSNPs(bool isvcf);
		void loadRefSeq();
		void loadTargets();
		void saveSequence();

		vector<string>& getChroms() {return chromosomes;}
		
		vector<CNV>& getCNVs(string chr) {return cnvs[chr];}
		vector<SNV>& getSNVs(string chr) {return snvs[chr];}
		vector<SNP>& getSNPs(string chr) {return sc[chr];}
		vector<Insert>& getInserts(string chr) {return inserts[chr];}
		vector<Del>& getDels(string chr) {return dels[chr];}
		vector<Target>& getTargets(string chr) {return targets[chr];}
		
		map<string, vector<Target> >& getTargets() {return targets;}
		
		long getChromLen(string chr);
		long getGenomeLength();
		
		char* getSubSequence(string chr, int startPos, int length);
		char* getSubAltSequence(string chr, int startPos, int length);
		
		void splitToFrags(vector<Fragment>& fragments);
};


#endif

