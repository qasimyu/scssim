// ***************************************************************************
// MyDefine.cpp (c) 2018 Zhenhua Yu <qasim0208@163.com>
// Health Informatics Lab, Ningxia University
// All rights reserved.

#include <vector>

#include "MyDefine.h"

using namespace std;

MyDefine::MyDefine() {
}

MyDefine::~MyDefine() {
}

string current_version = "1.0";

double ZERO_FINAL = 2.2204e-16;

Config config;
Genome genome;
Malbac malbac;
Profile profile;
ThreadPool* threadPool;
SeqWriter* swp;

//****** normalize matrix ******//
Matrix<double> norm_trans(Matrix<double> &T, double thres) {
	Matrix<double> ret(T);
	double *p1 = ret.getEntrance();
	double ZERO_FINAL = 2.2204e-16;
	int rows = ret.getROWS();
	int cols = ret.getCOLS();
	for(size_t i = 0; i < rows; i++) {
		Matrix<double> temp = ret.Row(i);
		double tmp = temp.sum();
		if(p1[i*cols+i] < thres*tmp) {
			temp.set(0, i, 0);
			Matrix<double> a = temp*((1-thres)/temp.sum());
			a.set(0, i, thres);
			ret.setRow(i, a);
		}
		else {
			Matrix<double> a = temp/(tmp+ZERO_FINAL);
			ret.setRow(i, a);
		}
	}
	return ret;
}

//****** pdf of normal distribution ******//
double normpdf(double x, double mu, double sigma) {
	double PI = 3.1415926;
	return exp(-pow(x-mu,2)/(2*pow(sigma,2)))/(sqrt(2*PI)*sigma);
}

//****** pdf of student't distribution ******//
double stupdf(double x, double mu, double nu, double sigma) {
	double PI = 3.1415926;
	double ret = lgamma((nu+1)*0.5)-(nu+1)*0.5*log(1+pow(x-mu,2)/(nu*sigma*sigma))
				-lgamma(nu*0.5)-0.5*log(nu)-log(sigma)-0.5*log(PI);
	ret = exp(ret);
	return ret;
}

//****** random number from poisson distribution ******//
long poissRand(double lambda)
{
	long x = -1;
	double u;
	double log1 = 0, log2 = -lambda;
	do {
		u = randomDouble(0, 1);
		log1 += log(u);
		x++;
	} while(log1 >= log2);
	return x;
}

//****** produce a matrix containing random intergers according to the probability distribution ******//
Matrix<int> randsrc(int m, int n, Matrix<int> alphabet, Matrix<double> aprob, bool iscdf) {
	Matrix<int> ret;
	if(m <= 0 || n <= 0) {
		return ret;
	}
	
	assert(alphabet.getROWS() == 1 && alphabet.getCOLS() > 0);
	assert(aprob.getROWS() == 1 && aprob.getCOLS() == alphabet.getCOLS());
	
	Matrix<double> prob;
	if(iscdf) {
		prob = aprob;
	}
	else {
		prob = aprob.cumsum();
	}

	srand(time(0));
	
	ret.resize(m, n, false);
	int* p = ret.getEntrance();
	int i, j, k;
	for(i = 0; i < m; i++) {
		for(j = 0; j < n; j++) {
			p[i*n+j] = randsrc(alphabet, prob, true);
		}
	}
	
	return ret;
}

//****** produce a random integer according to the probability distribution ******//
int randsrc(Matrix<int>& alphabet, Matrix<double>& aprob, bool iscdf) {
	assert(alphabet.getROWS() == 1 && alphabet.getCOLS() > 0);
	assert(aprob.getROWS() == 1 && aprob.getCOLS() == alphabet.getCOLS());
	
	int ac = aprob.getCOLS();
	double* p1;
	if(iscdf) {
		p1 = aprob.getEntrance();
	}
	else {
		Matrix<double> prob = aprob.cumsum();
		p1 = prob.getEntrance();
	}
	
	int* p2 = alphabet.getEntrance();
	double r = randomDouble(ZERO_FINAL, 1);
	for(size_t k = 0; k < ac; k++) {
		if(r <= p1[k]) {
			return p2[k];
		}
	}
	return p2[ac-1];
}

//****** produce a matrix containing random index according to the probability distribution ******//
Matrix<int> randIndx(int m, int n, Matrix<double> aprob, bool iscdf) {
	Matrix<int> ret;
	if(m <= 0 || n <= 0) {
		return ret;
	}
	
	assert(aprob.getROWS() == 1 && aprob.getCOLS() >= 1);
	
	Matrix<double> prob;
	if(iscdf) {
		prob = aprob;
	}
	else {
		prob = aprob.cumsum();
	}

	srand(time(0));
	
	ret.resize(m, n, false);
	int* p = ret.getEntrance();
	int i, j, k;
	for(i = 0; i < m; i++) {
		for(j = 0; j < n; j++) {
			p[i*n+j] = randIndx(prob, true);
		}
	}
	
	return ret;
}

int randIndx(Matrix<double>& aprob, bool iscdf) {
	int ac = aprob.getCOLS();
	double* p1;
	Matrix<double> prob;
	if(iscdf) {
		p1 = aprob.getEntrance();
	}
	else {
		prob = aprob.cumsum();
		p1 = prob.getEntrance();
	}
	
	double r = randomDouble(ZERO_FINAL, 1);
	for(size_t k = 0; k < ac; k++) {
		if(r <= p1[k]) {
			return k;
		}
	}
	return ac-1;
}

void* batchSampling(const void* args) {
	double* paras = (double*) args;
	unsigned int i, j;
	unsigned int num_samples = paras[0];
	unsigned int ac = paras[1];
	
	for(i = 0; i < num_samples; i++) {
		j = randIndx(&paras[ac+2], ac);
		paras[j+2]++;
	}
}

unsigned int* randIndx_hp(Matrix<double>& prob, unsigned long n, unsigned int* ret, bool addto) {
	unsigned int i, j = 0, ac = prob.getCOLS();
	double* p = prob.getEntrance();
	unsigned int loadPerThread = max((unsigned int) 1, min((unsigned int) 1000, ac/threadPool->getThreadNumber()));
	unsigned int sindx = 0, eindx;
	unsigned long count = 0;
	vector<double*> threadParas;
	vector<double> totalProbs;
	while(sindx < ac) {
		if(sindx+loadPerThread > ac) {
			eindx = ac-1;
		}
		else {
			eindx = sindx+loadPerThread-1;
		}
		double* paras = new double[(eindx-sindx+1)*2+2];
		memset(&paras[2], 0, sizeof(double)*(eindx-sindx+1));
		double totalProb = 0;
		for(i = sindx; i <= eindx; i++) {
			totalProb += p[i];
		}
		for(i = sindx; i <= eindx; i++) {
			paras[i-sindx+eindx-sindx+3] = paras[i-sindx+eindx-sindx+2] + p[i]/totalProb;
		}
		paras[0] = (unsigned int) (totalProb*n);
		paras[1] = eindx-sindx+1;
		count += paras[0];
		totalProbs.push_back(totalProb);
		threadParas.push_back(paras);
		sindx += loadPerThread;
	}
	n -= count;
	if(n > 0) {
		unsigned int m = totalProbs.size();
		double* probs = new double[m];
		probs[0] = totalProbs[0];
		for(i = 1; i < m; i++) {
			probs[i] = probs[i-1]+totalProbs[i];
		}
		while(n-- > 0) {
			j = randIndx(probs, m);
			threadParas[j][0] += 1;
		}
		delete[] probs;	
	}
	totalProbs.clear();
	
	for(i = 0; i < threadParas.size(); i++) {
		threadPool->pool_add_work(&batchSampling, threadParas[i], i);
	}
	threadPool->wait();
	
	sindx = 0;
	for(i = 0; i < threadParas.size(); i++) {
		double* paras = threadParas[i];
		eindx = sindx+paras[1]-1;
		j = 2;
		if(addto) {
			while(sindx <= eindx) {
				ret[sindx++] += paras[j++];
			}
		}
		else {
			while(sindx <= eindx) {
				ret[sindx++] = paras[j++];
			}
		}
		delete[] paras;
	}
}

unsigned int randIndx(double *cdf, unsigned int ac) {
	double r = threadPool->randomDouble(ZERO_FINAL, 1);
	for(size_t k = 0; k < ac; k++) {
		if(r <= cdf[k]) {
			return k;
		}
	}
	return ac-1;
}

//****** produce random double ******//
double randomDouble(double start, double end) {
	return start+(end-start)*(rand()/(RAND_MAX+1.0));
}

//****** produce random int ******//
long randomInteger(long start, long end) {
	return start+(end-start)*(rand()/(RAND_MAX+1.0));
}

//****** trim string ******//
string trim(const string &str, const char *charlist) {
	string ret(str);
	size_t indx = ret.find_first_not_of(charlist);
	if(indx != string::npos) {
		ret.erase(0, indx);
		indx = ret.find_last_not_of(charlist);
		ret.erase(indx+1);
	}
	else {
		ret.erase();
	}
	return ret;
}

//****** abbreviation of chromosome name ******//
string abbrOfChr(string chr) {
	string abbr_chr = chr;
	size_t i = abbr_chr.find("chrom");
	if(i == string::npos) {
		i = abbr_chr.find("chr");
		if(i != string::npos) {
			abbr_chr = abbr_chr.substr(i+3,abbr_chr.size()-3);
		}
	}
	else {
		abbr_chr = abbr_chr.substr(i+5,abbr_chr.size()-5);
	}
	return abbr_chr;
}

//****** index of a base ******//
int getIndexOfBase(char base) {
	static string bases = config.getStringPara("bases");
	for(size_t i = 0; i < bases.length(); i++) {
		if(base == bases[i]) {
			return i;
		}
	}
	return -1;
}

//****** get next valid line from a file ******//
bool getNextLine(ifstream& ifs, string& line, int& lineNum) {
	line = "";
	while(getline(ifs, line)) {
		lineNum++;
		if(!line.empty() && line.at(0) != '#') {
			break;
		}
	}
	if(line.empty()) {
		return false;
	}
	return true;
}

//****** complementary sequence ******//
char getComplementBase(char base) {
	char c;
	switch(base){
		case 'A': c = 'T'; break;
		case 'T': c = 'A'; break;
		case 'C': c = 'G'; break;
		case 'G': c = 'C'; break;
		case 'a': c = 't'; break;
		case 't': c = 'a'; break;
		case 'c': c = 'g'; break;
		case 'g': c = 'c'; break;
		case 'N': c = 'N'; break;
		default : c = 'N';
	}
	return c;
}

char* getComplementSeq(char* sequence) {
	if(sequence == NULL) {
		return NULL;
	}
	int i, n = strlen(sequence);
	for(i = 0; i < n; i++) {
		sequence[i] = getComplementBase(sequence[i]);
	}
	return sequence;
}

//****** calculate GC percent ******//
int calculateGCPercent(char* sequence) {
	int gcCount = 0;
	int nCount = 0;
	if(sequence == NULL || sequence[0] == '\0') {
		return 0;
	}
	int n = strlen(sequence);
	for(size_t i = 0; i < n; i++) {
		if(sequence[i] == 'G' || sequence[i] == 'C') {
			gcCount++;
		}
		else if(sequence[i] == 'N') {
			nCount++;
		}
	}
	/*
	if(nCount > n/2) {
		return -1;
	}
	*/
	if(nCount > 0) {
		return -1;
	}
	return 100*gcCount/(n-nCount);
}

//****** calculate GC content ******//
double calculateGCContent(char* sequence) {
	int gcCount = 0;
	int nCount = 0;
	if(sequence == NULL || sequence[0] == '\0') {
		return 0;
	}
	int n = strlen(sequence);
	for(size_t i = 0; i < n; i++) {
		if(sequence[i] == 'G' || sequence[i] == 'C') {
			gcCount++;
		}
		else if(sequence[i] == 'N') {
			nCount++;
		}
	}
	/*
	if(nCount > n/2) {
		return -1;
	}
	*/
	if(nCount > 0) {
		return -1;
	}
	return 1.0*gcCount/(n-nCount);
}

int countGC(char* sequence) {
	int gcCount = 0, nCount = 0;
	if(sequence == NULL || sequence[0] == '\0') {
		return 0;
	}
	int i, n = strlen(sequence);
	for(i = 0; i < n; i++) {
		if(sequence[i] == 'G' || sequence[i] == 'C') {
			gcCount++;
		}
		else if(sequence[i] == 'N') {
			nCount++;
		}
	}
	if(nCount > 0) {
		return 0;
	}
	return gcCount;
}


