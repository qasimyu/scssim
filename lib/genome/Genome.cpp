// ***************************************************************************
// Genome.cpp (c) 2018 Zhenhua Yu <qasim0208@163.com>
// Health Informatics Lab, Ningxia University
// All rights reserved.

#include <iostream>
#include <fstream>
#include <sstream>
#include <cassert>
#include <cmath>
#include <unistd.h>
#include <algorithm>

#include "split.h"
#include "MyDefine.h"
#include "Genome.h"

void Genome::loadData() {
	loadAbers();
	loadSNPs();
	loadRefSeq();
	loadTargets(); // WES is currently not supported
	
	curChr = "NA";
}

void Genome::loadTrainData() {
	vcfParser.setVCF(config.getStringPara("vcf"));
	vcfParser.parse();
	loadRefSeq();
	loadTargets();
	curChr = "NA";
}

void Genome::loadAbers() {
	string aberFile = config.getStringPara("var");
	if(aberFile.empty()) {
		return;
	}
	
	ifstream ifs;
	ifs.open(aberFile.c_str());
	if(!ifs.is_open()) {
		cerr << "can not open file " << aberFile << endl;
		exit(-1);
	}
	
	string line;
	int lineNum = 0;
	int cnvCount = 0, snvCount = 0, insertCount = 0, delCount = 0;
	while(getline(ifs,line)) {
		lineNum++;
		if(line.empty() || line[0] == '#') {
			continue;
		}
		//cerr << line << "\n\t";
		vector<string> fields = split(line,'\t');
		
		string aberType = fields[0];
		//CNV
		if(aberType.compare("c") == 0) {
			if(fields.size() != 6) {
				cerr << "ERROR: line " << lineNum << " has wrong number of fields in file " << aberFile << endl;
				cerr << line << endl;
				exit(1);
			}
			string chr = abbrOfChr(fields[1]);
			long spos = atol(fields[2].c_str());// 1-based
			long epos = atol(fields[3].c_str());// 1-based
			float cn = atof(fields[4].c_str());
			float mcn = atof(fields[5].c_str());
			if(cn < mcn) {
				cerr << "ERROR: total copy number should be not lower than major copy number at line " << lineNum << " in file " << aberFile << endl;
				cerr << line << endl;
				exit(1);
			}
			if(cn-mcn > mcn) {
				mcn = cn-mcn;
			}
			CNV cnv(spos, epos, cn, mcn);
			cnvs[chr].push_back(cnv);
			cnvCount++;
		}
		//SNV
		else if(aberType.compare("s") == 0) {
			if(fields.size() != 6) {
				cerr << "ERROR: line " << lineNum << " has wrong number of fields in file " << aberFile << endl;
				cerr << line << endl;
				exit(1);
			}
			string chr = abbrOfChr(fields[1]);
			long pos = atol(fields[2].c_str());// 1-based
			char ref = fields[3].at(0);
			char alt = fields[4].at(0);
			if(ref == alt) {
				cerr << "ERROR: the mutated allele should be not same as the reference allele at line " << lineNum << " in file " << aberFile << endl;
				cerr << line << endl;
				exit(1);
			}
			string typeCode = fields[5];
			if(typeCode.compare("homo") != 0 && typeCode.compare("het") != 0) {
				cerr << "ERROR: unrecognized SNV type at line " << lineNum << " in file " << aberFile << endl;
				cerr << line << endl;
				exit(1);
			}
			varType type = (typeCode.compare("het") == 0)? het : homo;
			SNV snv(pos, ref, alt, type);
			snvs[chr].push_back(snv);
			snvCount++;
		}
		//Insert
		else if(aberType.compare("i") == 0) {
			if(fields.size() != 5) {
				cerr << "ERROR: line " << lineNum << " has wrong number of fields in file " << aberFile << endl;
				cerr << line << endl;
				exit(1);
			}
			string chr = abbrOfChr(fields[1]);
			long pos = atol(fields[2].c_str());// 1-based
			string seq = fields[3];
			string typeCode = fields[4];
			if(typeCode.compare("homo") != 0 && typeCode.compare("het") != 0) {
				cerr << "ERROR: unrecognized insert type at line " << lineNum << " in file " << aberFile << endl;
				cerr << line << endl;
				exit(1);
			}
			varType type = (typeCode.compare("het") == 0)? het : homo;
			Insert insert(pos, seq, type);
			inserts[chr].push_back(insert);
			insertCount++;
		}
		//Del
		else if(aberType.compare("d") == 0) {
			if(fields.size() != 5) {
				cerr << "ERROR: line " << lineNum << " has wrong number of fields in file " << aberFile << endl;
				cerr << line << endl;
				exit(1);
			}
			string chr = abbrOfChr(fields[1]);
			long pos = atol(fields[2].c_str());// 1-based
			int length = atoi(fields[3].c_str());
			string typeCode = fields[4];
			if(typeCode.compare("homo") != 0 && typeCode.compare("het") != 0) {
				cerr << "ERROR: unrecognized deletion type at line " << lineNum << " in file " << aberFile << endl;
				cerr << line << endl;
				exit(1);
			}
			varType type = (typeCode.compare("het") == 0)? het : homo;
			Deletion del(pos, length, type);
			dels[chr].push_back(del);
			delCount++;
		}
		else {
			cerr << "ERROR: unrecognized aberraton type at line " << lineNum << " in file " << aberFile << endl;
			cerr << line << endl;
			exit(1);
		}
	}
	ifs.close();
	cerr << "\nDetails of the aberrations loaded from file " << aberFile << " are as follows:" << endl;
	cerr << "CNV: " << cnvCount << endl;
	cerr << "SNV: " << snvCount << endl;
	cerr << "Insert: " << insertCount << endl;
	cerr << "Deletion: " << delCount << endl;
}

void Genome::loadSNPs() {
	string snpFile = config.getStringPara("snp");
	if(snpFile.empty()) {
		return;
	}
	sc.readSNPs(snpFile);
	cerr << "\n" << sc.SNPNumber() << " SNPs to simulate were loaded from file " << snpFile << endl;
}

void Genome::loadRefSeq() {
	string refFile = config.getStringPara("ref");
	string tmp = refFile;
	if(tmp.empty()) {
		cerr << "reference sequence file not specified!" << endl;
		exit(1);
	}
	if(strcmp(tmp.substr(tmp.length()-3).c_str(), ".gz") == 0) {
		string cmd = "gzip -cd " + tmp + " > "+tmp.substr(0, tmp.length()-3);
		system(cmd.c_str());
		tmp = tmp.substr(0, refFile.length()-3);
	}
	fr.open(tmp, false);
	chromosomes = (*(fr.index)).sequenceNames;
	if(chromosomes.empty()) {
		cerr << "ERROR: reference sequence cannot be empty!" << endl;
		exit(1);
	}
	cerr << "\nReference sequence was loaded from file " << refFile << endl;
}

void Genome::loadTargets() {
	string targetFile = config.getStringPara("target");
	if(targetFile.empty()) {
		return;
	}
	map<string, vector<Target> >::iterator m_it;
	
	ifstream ifs;
	ifs.open(targetFile.c_str());
	if(!ifs.is_open()) {
		cerr << "can not open target file " << targetFile << endl;
		exit(-1);
	}
	
	string line;
	int lineNum = 0;
	int targetNum = 0;
	while(getline(ifs,line)) {
		lineNum++;
		//cerr << line << "\n\t";
		vector<string> fields = split(line,'\t');
		if(fields.size() < 3) {
			cerr << "ERROR: line " << lineNum << " should have at least 3 fields in file " << targetFile << endl;
			cerr << line << endl;
			exit(1);
		}
		
        string chr = abbrOfChr(fields[0]);
		long chrLen = getChromLen(chr);
		if(chrLen <= 0) {
			continue;
		}
		Target target;
		target.spos = max((long) 1, atol(fields[1].c_str())-50+1); // 0-based 
		//target.epos = min(chrLen, atol(fields[2].c_str())+50); // 1-based
		long tmp;
		if(atol(fields[2].c_str()) <= 0) {
			tmp = chrLen-(-atol(fields[2].c_str()))%chrLen;
		}
		else {
			tmp = atol(fields[2].c_str());
		}
		target.epos = min(chrLen, tmp+50); // 1-based
		targets[chr].push_back(target);
		targetNum++;
	}
	ifs.close();
	cerr << "\ntotal " << targetNum << " targets were loaded from file " << targetFile << endl;
	
	divideTargets();
}

long Genome::getChromLen(string chr) {
	FastaIndexEntry fie;
	int i;
	for(i = 0; i < chromosomes.size(); i++) {
		if(chromosomes[i].compare(chr) == 0) {
			break;
		}
	}
	if(i == chromosomes.size()) {
		return 0;
	}
	fie = fr.index->entry(chr);
	return fie.length;
}

long Genome::getGenomeLength() {
	long n = 0;
	for(int i = 0; i < chromosomes.size(); i++) {
		n += getChromLen(chromosomes[i]);
	}
	return n;
}

char* Genome::getSubSequence(string chr, int startPos, int length) {
	string seq = fr.getSubSequence(chr, startPos, length);
	transform(seq.begin(), seq.end(), seq.begin(), (int (*)(int))toupper);
	char* s = new char[seq.length()+1];
	strcpy(s, seq.c_str());
	return s;
}

char* Genome::getSubRefSequence(string chr, int startPos, int length) {
	if(curChr.compare(chr) != 0) {
		generateChrSequence(chr);
	}
	string seq = refSequence.substr(startPos, length);
	char* s = new char[seq.length()+1];
	strcpy(s, seq.c_str());
	return s;
}

char* Genome::getSubAltSequence(string chr, int startPos, int length) {
	if(curChr.compare(chr) != 0) {
		generateChrSequence(chr);
	}
	//string seq = altSequence[chr].substr(startPos, length);
	string seq = altSequence.substr(startPos, length);
	char* s = new char[seq.length()+1];
	strcpy(s, seq.c_str());
	return s;
}

void Genome::generateChrSequence(string chr) {	
	int i, j;
	vector<string>::iterator it = find(chromosomes.begin(), chromosomes.end(), chr);
	if(it == chromosomes.end()) {
		cerr << "ERROR: unrecognized chromosome identifier \"" << chr 
			 << "\" when generating alternative sequence!" << endl;
		exit(1);
	}
	
	curChr = chr;
	vector<SNV>& snvsOfChr = getRealSNVs(chr);
	vector<Insert>& insertsOfChr = getRealInserts(chr);
	vector<Deletion>& delsOfChr = getRealDels(chr);
	
	refSequence = fr.getSubSequence(chr, 0, getChromLen(chr));
	altSequence = refSequence;
	for(j = 0; j < snvsOfChr.size(); j++) {
		SNV snv = snvsOfChr[j];
		long pos = snv.getPosition();
		altSequence[pos-1] = snv.getAlt();
		if(snv.getType() == homo) {
			refSequence[pos-1] = snv.getAlt();
		}
	}
	transform(refSequence.begin(), refSequence.end(), refSequence.begin(), (int (*)(int))toupper);
	transform(altSequence.begin(), altSequence.end(), altSequence.begin(), (int (*)(int))toupper);
}

void Genome::saveSequence() {
	int ploidy = config.getIntPara("ploidy");
	int mCN = (int) ceil((float) ploidy/2);
	int i, j, k;
	
	string outFile = config.getStringPara("output");
	ofstream ofs;
	ofs.open(outFile.c_str());
	if(!ofs.is_open()) {
		cerr << "can not open file " << outFile << endl;
		exit(-1);
	}
	
	for(i = 0; i < chromosomes.size(); i++) {
		string chr = chromosomes[i];
		vector<CNV>& cnvsOfChr = getSimuCNVs(chr);
		long segStartPos = 1, segEndPos = -1;
		
		vector<string> chrSequences(ploidy, "");
		for(j = 0; j < cnvsOfChr.size(); j++) {		
			if(segStartPos > getChromLen(chr)) {
				break;
			}
			CNV cnv = cnvsOfChr[j];
			cnv.epos = min(cnv.epos, getChromLen(chr));
			if(segStartPos < cnv.spos) {
				generateSegment(chrSequences, chr, segStartPos, cnv.spos-1, ploidy, mCN);
			}
			generateSegment(chrSequences, chr, cnv.spos, cnv.epos, cnv.CN, cnv.mCN);
			segStartPos = cnv.epos+1;
		}
		if(segStartPos <= getChromLen(chr)) {
			segEndPos = getChromLen(chr);
			generateSegment(chrSequences, chr, segStartPos, getChromLen(chr), ploidy, mCN);
		}
		
		for(j = 0; j < ploidy; j++) {
			stringstream ss;
			ss << j+1;
			ofs << ">" << chr << "_" << ss.str() << "_" << getChromLen(chr) << endl;
			unsigned int sindx = 0;
			unsigned int length = chrSequences[j].length();
			int width = 100;
			while(sindx < length) {
				if(sindx+width > length) {
					ofs << chrSequences[j].substr(sindx) << endl;
				}
				else {
					ofs << chrSequences[j].substr(sindx, width) << endl;
				}
				sindx += width;
			}
		}
	}
	ofs.close();
}

void Genome::generateSegment(vector<string>& sequences, string chr, long segStartPos, long segEndPos, int CN, int mCN) {
	if(CN == 0) {
		return;
	}

	char* refSeq = getSubSequence(chr, segStartPos-1, segEndPos-segStartPos+1);
	if(refSeq == NULL) {
		return;
	}

	transform(refSeq, refSeq+strlen(refSeq), refSeq, (int (*)(int))toupper);
	
	vector<SNP>& snpsOfChr = getSimuSNPs(chr);
	vector<SNV>& snvsOfChr = getSimuSNVs(chr);
	vector<Insert>& insertsOfChr = getSimuInserts(chr);
	vector<Deletion>& delsOfChr = getSimuDels(chr);
	
	int i, j, k, n;
	int ploidy = config.getIntPara("ploidy");
	unsigned int refSize = strlen(refSeq);
	vector<int>::iterator it;
	
	vector<int> mIndx;
	vector<int> seqReps;
	
	if(CN < ploidy) {
		for(i = 0; i < CN; i++) {
			while(1) {
				j = randomInteger(0, ploidy);
				it = find(seqReps.begin(), seqReps.end(), j);
				if(it == seqReps.end()) {
					seqReps.push_back(j);
					break;
				}
			}
		}
		for(i = 0; i < mCN; i++) {
			mIndx.push_back(seqReps[i]);
		}
	}
	else {
		for(i = 0; i < ploidy; i++) {
			seqReps.push_back(1);
		}
		n = CN-ploidy;
		k = randomInteger(0, ploidy);
		for(i = n; i >= 0; i--) {
			if(seqReps[k]+i == mCN) {
				seqReps[k] += i;
				mIndx.push_back(k);
				break;
			}
			else if(seqReps[k]+i == CN-mCN) {
				seqReps[k] += i;
				for(j = 0; j < ploidy; j++) {
					if(j != k) {
						mIndx.push_back(j);
					}
				}
				break;
			}
		}
		if(i >= 0) {
			n -= i;
			while(n > 0) {
				j = randomInteger(0, ploidy);
				if(j != k) {
					seqReps[j]++;
					n--;
				}
			}
		}
		else {
			while(n > 0) {
				j = randomInteger(0, ploidy);
				seqReps[j]++;
				n--;
			}
			for(i = 0; i < ploidy; i++) {
				mIndx.push_back(i);
			}
		}
	}
	
	vector<string> segSeqs;
	if(CN < ploidy) {
		for(i = 0; i < ploidy; i++) {
			it = find(seqReps.begin(), seqReps.end(), i);
			if(it != seqReps.end()) {
				segSeqs.push_back(refSeq);
			}
			else {
				segSeqs.push_back("");
			}
		}
	}
	else {
		for(i = 0; i < ploidy; i++) {
			string tmp = "";
			for(j = 0; j < seqReps[i]; j++) {
				tmp += refSeq;
			}
			segSeqs.push_back(tmp);
		}
	}
	delete[] refSeq;
	
	//simu SNP
	k = 0;
	for(i = 0; i < snpsOfChr.size(); i++) {
		SNP snp = snpsOfChr[i];
		long pos = snp.getPosition();
		if(pos >= segStartPos && pos <= segEndPos) {
			int sindx = pos-segStartPos;
			for(j = 0; j < ploidy; j++) {
				it = find(mIndx.begin(), mIndx.end(), j);
				if((k == 0 && it == mIndx.end()) || (k == 1 && it != mIndx.end())) {
					continue;
				}
				string& segSeq = segSeqs[j];
				unsigned int segLen = segSeq.length();
				for(int t = 0; t < segLen/refSize; t++) {
					segSeq[sindx+t*refSize] = snp.getNucleotide();
				}
			}
			k = (k+1)%2;
		}
	}
	
	//simu SNV
	k = 0;
	for(i = 0; i < snvsOfChr.size(); i++) {
		SNV snv = snvsOfChr[i];
		long pos = snv.getPosition();
		if(pos >= segStartPos && pos <= segEndPos) {
			int sindx = pos-segStartPos;
			varType type = snv.getType();
			char alt = snv.getAlt();
			if(type == homo) {
				for(j = 0; j < ploidy; j++) {
					string& segSeq = segSeqs[j];
					unsigned int segLen = segSeq.length();
					for(int t = 0; t < segLen/refSize; t++) {
						segSeq[sindx+t*refSize] = alt;
					}
				}
			}
			else {
				for(j = 0; j < ploidy; j++) {
					it = find(mIndx.begin(), mIndx.end(), j);
					if((k == 0 && it == mIndx.end()) || (k == 1 && it != mIndx.end())) {
						continue;
					}
					string& segSeq = segSeqs[j];
					unsigned int segLen = segSeq.length();
					for(int t = 0; t < segLen/refSize; t++) {
						segSeq[sindx+t*refSize] = alt;
					}
				}
				k = (k+1)%2;
			}
		}
	}
	
	//simu Insert
	map<int, map<int, int> > insertsPerhaploidy;
	map<int, int>::iterator m_it;
	int* insertLens = new int[ploidy];
	memset(insertLens, 0, ploidy*sizeof(int));
	k = 0;
	for(i = 0; i < insertsOfChr.size(); i++) {
		Insert insert = insertsOfChr[i];
		long pos = insert.getPosition();
		if(pos >= segStartPos && pos <= segEndPos) {
			int sindx = pos-segStartPos;
			varType type = insert.getType();
			string seq = insert.getSequence();
			if(type == homo) {
				for(j = 0; j < ploidy; j++) {
					int offset = 0;
					map<int, int>& insertedSeq = insertsPerhaploidy[j];
					for(m_it = insertedSeq.begin(); m_it!= insertedSeq.end(); m_it++) {
						if((*m_it).first <= sindx) {
							offset += (*m_it).second;
						}
					}
					string& segSeq = segSeqs[j];
					n = segSeq.length()/(refSize+insertLens[j]);
					int len = seq.length();
					for(int t = 0; t < n; t++) {
						segSeq.insert(sindx+offset+t*(refSize+insertLens[j]+len), seq);
					}
					insertLens[j] += seq.length();
					insertedSeq.insert(make_pair(sindx, seq.length()));
				}
			}
			else {
				for(j = 0; j < ploidy; j++) {
					it = find(mIndx.begin(), mIndx.end(), j);
					if((k == 0 && it == mIndx.end()) || (k == 1 && it != mIndx.end())) {
						continue;
					}
					int offset = 0;
					map<int, int>& insertedSeq = insertsPerhaploidy[j];
					for(m_it = insertedSeq.begin(); m_it!= insertedSeq.end(); m_it++) {
						if((*m_it).first <= sindx) {
							offset += (*m_it).second;
						}
					}
					string& segSeq = segSeqs[j];
					n = segSeq.length()/(refSize+insertLens[j]);
					int len = seq.length();
					for(int t = 0; t < n; t++) {
						segSeq.insert(sindx+offset+t*(refSize+insertLens[j]+len), seq);
					}
					insertLens[j] += seq.length();
					insertedSeq.insert(make_pair(sindx, seq.length()));
				}
				k = (k+1)%2;
			}
		}
	}
	
	//simu Del
	map<int, map<int, int> > delsPerhaploidy;
	int* delLens = new int[ploidy];
	memset(delLens, 0, ploidy*sizeof(int));
	for(i = 0; i < delsOfChr.size(); i++) {
		Deletion del = delsOfChr[i];
		long pos = del.getPosition();
		if(pos >= segStartPos && pos <= segEndPos) {
			int sindx = pos-segStartPos;
			varType type = del.getType();
			int delLen = del.getLength();
			if(type == homo) {
				for(j = 0; j < ploidy; j++) {
					int offset = 0;
					map<int, int>& insertedSeq = insertsPerhaploidy[j];
					for(m_it = insertedSeq.begin(); m_it!= insertedSeq.end(); m_it++) {
						if((*m_it).first <= sindx) {
							offset += (*m_it).second;
						}
					}
					map<int, int>& delSeq = delsPerhaploidy[j];
					for(m_it = delSeq.begin(); m_it!= delSeq.end(); m_it++) {
						if((*m_it).first <= sindx) {
							offset -= (*m_it).second;
						}
					}
					if(sindx+offset < 0) {
						continue;
					}
					string& segSeq = segSeqs[j];
					n = segSeq.length()/(refSize+insertLens[j]-delLens[j]);
					for(int t = 0; t < n; t++) {
						segSeq.erase(sindx+offset+t*(refSize+insertLens[j]-delLens[j]-delLen), delLen);
					}
					delLens[j] += delLen;
					delSeq.insert(make_pair(sindx, delLen));
				}
			}
			else {
				for(j = 0; j < ploidy; j++) {
					it = find(mIndx.begin(), mIndx.end(), j);
					if((k == 0 && it == mIndx.end()) || (k == 1 && it != mIndx.end())) {
						continue;
					}
					int offset = 0;
					map<int, int>& insertedSeq = insertsPerhaploidy[j];
					for(m_it = insertedSeq.begin(); m_it!= insertedSeq.end(); m_it++) {
						if((*m_it).first <= sindx) {
							offset += (*m_it).second;
						}
					}
					map<int, int>& delSeq = delsPerhaploidy[j];
					for(m_it = delSeq.begin(); m_it!= delSeq.end(); m_it++) {
						if((*m_it).first <= sindx) {
							offset -= (*m_it).second;
						}
					}
					if(sindx+offset < 0) {
						continue;
					}
					string& segSeq = segSeqs[j];
					n = segSeq.length()/(refSize+insertLens[j]-delLens[j]);
					for(int t = 0; t < n; t++) {
						segSeq.erase(sindx+offset+t*(refSize+insertLens[j]-delLens[j]-delLen), delLen);
					}
					delLens[j] += delLen;
					delSeq.insert(make_pair(sindx, delLen));
				}
				k = (k+1)%2;
			}
		}
	}
	delete[] insertLens;
	delete[] delLens;
	
	for(i = 0; i < ploidy; i++) {
		if(segSeqs[i].length() > 0) {
			transform(segSeqs[i].begin(), segSeqs[i].end(), segSeqs[i].begin(), (int (*)(int))toupper);
		}
		sequences[i] += segSeqs[i];
	}
	
}

void Genome::divideTargets() {
	if(targets.empty()) {
		return;
	}
	int i, j, k;
	unsigned int targetMaxSize = config.getIntPara("fragSize");
	map<string, vector<Target> > newTargets;
	map<string, vector<Target> >::iterator it, m_it;
	for(it = targets.begin(); it != targets.end(); it++) {
		string chr = it->first;
		vector<Target>& targetsOfChr = it->second;
		for(j = 0; j < targetsOfChr.size(); j++) {
			Target target = targetsOfChr[j];
			long spos = target.spos;
			long tsize = target.epos-target.spos+1;
			int k = tsize/targetMaxSize;
			for(i = 0; i < k; i++) {
				Target newTarget;
				newTarget.spos = spos;
				if(i == k-1) {
					newTarget.epos = target.epos;
				}
				else {
					newTarget.epos = spos+targetMaxSize-1;
				}
				spos = newTarget.epos+1;
				for(m_it = newTargets.begin(); m_it != newTargets.end(); m_it++) {
					if((*m_it).first.compare(chr) == 0) {
						(*m_it).second.push_back(newTarget);
						break;
					}
				}
				if(m_it == newTargets.end()) {
					vector<Target> temp;
					temp.push_back(newTarget);
					newTargets.insert(make_pair(chr, temp));
				}
			}
			if(spos <= target.epos) {
				Target newTarget;
				newTarget.spos = spos;
				newTarget.epos = target.epos;
				for(m_it = newTargets.begin(); m_it != newTargets.end(); m_it++) {
					if((*m_it).first.compare(chr) == 0) {
						(*m_it).second.push_back(newTarget);
						break;
					}
				}
				if(m_it == newTargets.end()) {
					vector<Target> temp;
					temp.push_back(newTarget);
					newTargets.insert(make_pair(chr, temp));
				}	
			}
		}
	}
	targets = newTargets;
	
}

void Genome::splitToFrags(vector<Fragment>& fragments) {
	int i, j, k;
	int fragLen;
	for(i = 0; i < chromosomes.size(); i++) {
		string chr = chromosomes[i];
		long fragStartPos = 1;
		long chrLen = getChromLen(chr);
		while(fragStartPos <= chrLen) {
			fragLen = randomInteger(Fragment::minSize, Fragment::maxSize+1);
			if(fragStartPos+fragLen-1 > chrLen) {
				break;
			}
			Fragment frag1(chr, fragStartPos, fragLen, 1);
			fragments.push_back(frag1);
			//complementary strand
			Fragment frag2(chr, fragStartPos, fragLen, -1);
			fragments.push_back(frag2);
			fragStartPos += fragLen;
		}
		if(fragStartPos <= chrLen) {
			Fragment frag1(chr, fragStartPos, chrLen-fragStartPos+1, 1);
			fragments.push_back(frag1);
			Fragment frag2(chr, fragStartPos, chrLen-fragStartPos+1, 1);
			fragments.push_back(frag2);
		}	
	}
	for(i = 0; i < fragments.size(); i++) {
		fragments[i].createSequence();
	}
}


