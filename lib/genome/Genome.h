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
#include "vcfparser.h"
#include "Fasta.h"

using namespace std;

class CNV {
	public:
		CNV() {spos = -1;}
		CNV(long spos, long epos, float CN, float mCN)
			: spos(spos), epos(epos), CN(CN), mCN(mCN) {}
		long getStartPos() {return spos;}
		long getEndPos() {return epos;}
		float getCopyNumber() {return CN;}
		float getMajorCopyNumber() {return mCN;}
		
		long spos;
		long epos;
		float CN;
		float mCN;
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
		map<string, vector<Deletion> > dels;
		
		map<string, vector<Target> > targets;
		
		SNPOnChr sc;
		VcfParser vcfParser;		
		FastaReference fr;
		
		//map<string, string> altSequence;
		string refSequence;
		string altSequence;
		string curChr;
		
		void loadAbers();
		void loadSNPs();
		void loadRefSeq();
		void loadTargets();
		
		void divideTargets();
		void generateChrSequence(string chr);
		
		void generateSegment(vector<string>& sequences, string chr, long segStartPos, long segEndPos, int CN, int mCN);
		
	public:
		Genome() {curChr = "";}
		~Genome() {}

		void loadData();
		void loadTrainData();
		void saveSequence();

		vector<string>& getChroms() {return chromosomes;}
		
		vector<CNV>& getSimuCNVs(string chr) {return cnvs[chr];}
		vector<SNV>& getSimuSNVs(string chr) {return snvs[chr];}
		vector<SNV>& getRealSNVs(string chr) {return vcfParser.getSNVs(chr);}
		vector<SNP>& getSimuSNPs(string chr) {return sc[chr];}
		vector<Insert>& getSimuInserts(string chr) {return inserts[chr];}
		vector<Insert>& getRealInserts(string chr) {return vcfParser.getInserts(chr);}
		vector<Deletion>& getSimuDels(string chr) {return dels[chr];}
		vector<Deletion>& getRealDels(string chr) {return vcfParser.getDels(chr);}

		vector<Target>& getTargets(string chr) {return targets[chr];}
		
		map<string, vector<Target> >& getTargets() {return targets;}
		
		long getChromLen(string chr);
		long getGenomeLength();
		
		char* getSubSequence(string chr, int startPos, int length);
		char* getSubRefSequence(string chr, int startPos, int length);
		char* getSubAltSequence(string chr, int startPos, int length);
		
		void splitToFrags(vector<Fragment>& fragments);
};


#endif

