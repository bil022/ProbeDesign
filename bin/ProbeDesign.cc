#include <iostream>
#include <vector>
#include <map>
#include <set>
#include <cassert>
#include <algorithm>
#include <fstream>
#include <sstream>
#ifdef __unix__
#include <boost/regex.hpp>
using namespace boost;
#else
#include <regex>
#endif

using namespace std;

#define N 17
#define SZ (1L<<(N<<1))

class BitVector : public vector<bool> {
	public:
		BitVector() {
			resize(SZ, false);
		}
		void save(ofstream & ofs) {
			long wc=0;
			ofs << size() << "\n";
			unsigned char base=0, cnt=0, mask=0x3F;
			vector<bool>::iterator itr;
			for (itr=begin(); itr!=end(); itr++) {
				base<<=1;
				if (*itr) base|=(unsigned char)1;
				base&=mask; cnt++;
				if (cnt==6) {
					ofs << (unsigned char)('0'+base);
					cnt=0; wc++;
					if (wc%60==0)
						ofs << "\n";
				}
			}
			if (cnt) {
				while (cnt<6) {
					base<<=1;
					base&=mask; cnt++;
				}
				ofs << (unsigned char)('0'+base);
			}
			ofs << "\n";
		}
		void load(ifstream & ifs, bool verify=false) {
			long sz=0, off=0, wc=0, lc=0;
			ifs >> sz; wc=(sz+5)/6;
			cerr << "sz=" << sz << "\n";
			unsigned char mask=0x20;
			string ln;
			vector<bool>::iterator itr=begin();
			while (ifs >> ln) {
				// cerr << "#Bit_load: ["<<ln<<"]\n";
				string::iterator ptr=ln.begin();
				while (ptr!=ln.end()) {
					if (*ptr<'0')
						cerr << "[" << *ptr << ":" << (unsigned int)*ptr << ":" << sz << ":" << wc << "]\n";
					assert(*ptr>='0');
					unsigned char ch=*ptr-'0';
					for (int i=0; i<6; i++, ch<<=1) {
						if (ch&mask) {
							if (verify) {
								assert(*itr==true);
							}
							*itr=true;
						}
						if (sz) {
							itr++;
							sz--;
						}
						off++;
					}
					ptr++; wc--;
				}
				lc++;
				if (wc<=0)
					break;
			}
			assert(sz==0);
		}
};

#define MAX_CNT 255
/*
   class ByteVector : public vector<unsigned char> {
   public:
   ByteVector() {
   resize(SZ, 0);
   }
   void save(ofstream & ofs) {
   vector<unsigned char>::iterator itr=begin();
   long code=0, last_code=0;
   while (itr!=end()) {
   if (*itr) {
   ofs << (code-last_code) << "\t" << (unsigned int)*itr << "\n";
   last_code=code;
   }
   itr++; code++;
   }
   ofs << "0\t0\n";
   }
   void load(ifstream & ifs, bool verify=false) {
   long code, base=0; unsigned int cnt;
   while (ifs >> code >> cnt) {
   if (!cnt)
   break;
   if (verify) {
   assert((*this)[code+base]==(unsigned char)cnt);
   } else {
   unsigned int origin=(*this)[code+base], update=(unsigned int)cnt;
   if ((origin+update)>MAX_CNT) {
   cnt=(unsigned char)MAX_CNT;
   }
   }

   (*this)[code+base]=(unsigned char)cnt;
   base+=code;
   }
   }
   };*/

class ByteMap : public map<long, unsigned char> {
	public:
		void save(ofstream & ofs) {
			map<long, unsigned char>::iterator itr=begin();
			while (itr!=end()) {
				ofs << itr->first << "\t" << (unsigned int)itr->second << "\n";
				itr++;
			}
			ofs << "0\t0\n";
		}
		void load(ifstream & ifs, bool verify=false) {
			long code; unsigned int cnt;
			while (ifs >> code >> cnt) {
				if (!cnt) {
					break;
				}
				if (verify) {
					assert((*this)[code]==(unsigned char)cnt);
				} else {
					unsigned int origin=(*this)[code], update=(unsigned int)cnt;
					if ((origin+update)>MAX_CNT) {
						cnt=(unsigned char)MAX_CNT;
					}
				}
				(*this)[code]=(unsigned char)cnt;
			}
		}
};

static int pref_gc=18;

#ifndef MAX_REP_SUM
#define MAX_REP_SUM 75
#endif

#define WIN 30
#define MIN_GC 6
#define MAX_GC (WIN-MIN_GC)

#define DIRTY_INIT 1
#define DIRTY_MAX_CNTS 2
#define DIRTY_MAX_REP 4
#define DIRTY_MIN_GC 8
#define DIRTY_MAX_GC 16
#define DIRTY_OVERLAP 32

class Probe {
	char ch; short dirty, qual, winQual; // qual: 17-mer, winQual: probe qual
	int pos, off_cnts, max_cnts, rep, GC, rep_sum, GC_sum;
public:
	Probe(int pos, char ch, int off_cnts, short qual, int max_cnts, int rep_sum, int GC_sum)
		:pos(pos),ch(ch),off_cnts(off_cnts), qual(qual), winQual(qual),
		max_cnts(max_cnts), rep_sum(rep_sum),
		GC_sum(GC_sum), rep(0), GC(0)
	{
		dirty=DIRTY_INIT;
	}
	// update GC & rep
	bool isGood(Probe& probe) {
		GC=GC_sum-probe.GC_sum;
		rep=rep_sum-probe.rep_sum;
		dirty=0;
		if (max_cnts>probe.max_cnts)
			dirty|=DIRTY_MAX_CNTS;
		if (rep>=MAX_REP_SUM)
			dirty|=DIRTY_MAX_REP;
		if (GC<MIN_GC)
			dirty|=DIRTY_MIN_GC;
		if (GC>MAX_GC)
			dirty|=DIRTY_MAX_GC;
		return !dirty;
	}
	inline int getPos() { return pos; }
	inline short getQual() { return qual; }
	inline void updateWinQual(short q) { if (winQual>q) winQual=q; }
	inline char getWinQual() { return winQual; }
	inline int getMaxCnts() { return max_cnts; }
	inline int GCs() { return GC; }
	inline int off() { return rep; }
	inline short isDirty() { return dirty; }
	inline void markDirty() {
        dirty|=DIRTY_OVERLAP;
    }
	void dump(int nth, string& ref) {
		int s=pos+1-WIN, len=WIN; if (s<0) { s=0; len=pos+1; }
		cerr << (nth+1) << "\t" << pos << "\t" << ref.substr(s, len)
			<< "\t" << ch << "\t" << GC << "\t" << rep << "\t"
			<< winQual << "\t" << dirty << "[" << (char)('0'+dirty) << "]\n";
	}
	void print(int nth, string& ref) {
		cout << (nth+1) << "\t" << pos << "\t" << ref.substr(pos+1-WIN, WIN)
			<< "\t" << ch << "\t" << GC << "\t" << rep << "\t" << winQual << "\n";
	}
};

class Regex {
	public:
		// gffread -W -w transcriptW.fa -g ../hg38.fa genes.gtf
		// >NR_046018 gene=DDX11L1 loc:chr1|11874-14409|+
		static bool getGene(string& input, string & id) {
			const char* re=".*gene=([\\S]+).*";
			id="N/A";
			return match(input, re, id);
		}

		static bool getStrand(string& input, char & strand) {
			const char* re=".*loc:\\S+\\|\\d+\\-\\d+\\|([+-]).*";
			strand='+';
			string id;
			if (match(input, re, id)) {
				strand=id[0];
				// cerr << "Get strand " << strand << "@" << input << "\n";
				return true;
			}
			cerr << "Unknown strand " << input << "\n";
			return false;
		}
#ifdef __unix__
		static bool match(string& input, const char* re, string& id) {
			cmatch matched;
			regex expr(re);
			// cerr << "search '" << input << "' for '" << re << "'\n";
			if (regex_match(input.c_str(), matched, expr)) {
				id=string(matched[1].first, matched[1].second);
				return true;
			}
			return false;
		}
#else
		static bool match(string& input, const char* re, string& id) {
			regex e(re);
			smatch sm;
			if (regex_search(input, sm, e)) {
				id=sm[1];
				return true;
			}
			return false;
		}
#endif
};

#define MAX_OVERLAP 10
class ProbeWindow {
	string id, seq;
	vector<Probe> probes;
	int max_cnts, rep_sum, GC;
	public:
	void init(string& header) {
		id=header; seq.clear();
		probes.clear();
		max_cnts=rep_sum=GC=0;
	}
	void append(char ch, int rep, short qual=1) {
		ch=toupper(ch);
		seq+=ch;
		if (rep>=MAX_REP_SUM) {
			max_cnts++;
			rep_sum+=MAX_REP_SUM;
		} else {
			rep_sum+=rep;
		}
		// assert(rep_sum>0);
		if (ch=='G'||ch=='C') {
			GC++;
		}
		Probe probe((int)probes.size(), ch, rep, qual, max_cnts, rep_sum, GC);
		probes.push_back(probe);
	}

	void dump() {
		for (int i=0; i<probes.size(); i++)
			probes[i].dump(i, seq);
	}

	static char strand;
	static bool byProbe (Probe* a,Probe* b) {
		int aQ=a->getWinQual(), bQ=b->getWinQual();
		if (aQ!=bQ)
			return aQ>bQ;	
		int aGC=abs(a->GCs()-pref_gc), bGC=abs(b->GCs()-pref_gc);
		// cerr << aGC << " vs " << bGC << "\n";
		if (aGC!=bGC)
			return aGC<bGC;
		int cntA=a->getMaxCnts();
		int cntB=b->getMaxCnts();
		if (cntA!=cntB) {
			return cntA < cntB;
		}
		int offA=a->off(), offB=b->off();
		if (offA!=offB) {
			return offA<offB;
		}
		int posA=a->getPos(), posB=b->getPos();
		if (posA!=posB) {
			return posA<posB;
		}

		return a<b;
	}

	void printProbes() {
		for (int i=0; i<probes.size(); i++) {
			cerr << (char)('0'+probes[i].isDirty());
		}	cerr << "\n";
	}

	void design() {
		if (!probes.size())
			return;
		strand='+';
		Regex::getStrand(id, strand);
		cout << id << "\n";
		// cerr << probes.size() << " probes\n";
		vector<Probe*> pool;
		//1-17-40-79 ... 0-16-39-78
		for (int i=WIN-1; i<probes.size(); i++) {
			if (probes[i].isGood(probes[i-(WIN-1)])) {
				// update winQual
				for (int j=0; j<(WIN-N); j++) {
					//cerr << "Bug:i/j=" << i<<"/"<<j<<":"<<probes[i-j].getQual()<<"\n";
					assert(probes[i-j].getQual()>0);
					probes[i].updateWinQual(probes[i-j].getQual());
				}
				pool.push_back(&probes[i]);
			}
		}
		cerr << id << "\tMAX_REP_SUM:[" << MAX_REP_SUM << "]:pool["<<pool.size()<<"]\n";
        printProbes();
		//dump();
		//cerr << "Before sort:\n"; printProbes();
		sort(pool.begin(), pool.end(), byProbe);
		//cerr << "After sort:\n"; printProbes();
        //cerr << "\n";

		vector<Probe*>::iterator itr=pool.begin();
		int probes_found=0;
		for (itr=pool.begin(); itr!=pool.end(); itr++) {
			Probe& curr=**itr;

			int i, contig=0, max_contig=0, pos=curr.getPos();
			cerr << "Checking\t" << pos << "\n";
            // printProbes();

			for (i=pos-(WIN-1); i<=pos; i++) {
				if (i<0) i=0; if (i>=probes.size()) break;
				if (probes[i].isDirty())
					contig=0;
				else
					contig++;
				if (contig>max_contig)
					max_contig=contig;
			}
			bool found=false, dirty=false;
			if ((max_contig+MAX_OVERLAP)>=WIN) {
				//cout << pos << "\tOL= " << overlap << "\n";
				curr.print(probes_found, seq);
				found=true;
				probes_found++;
				for (i=pos-(WIN-1); i<=pos; i++) {
					if (i<0) i=0; if (i>=probes.size()) break;
					//cerr << "DIRTY?" << i << "\n";
                    //printProbes();
                    if (!(probes[i].isDirty()&DIRTY_OVERLAP)) {
                        // cerr << probes[i].isDirty() << "?\n";
                        dirty=true;
                        probes[i].markDirty();
                    }
                    //printProbes();
					//cerr << "DIRTY^" << i << "\n";
					//dirty=true;
					//cout << ","<<i;
				}	//cout << "\n";
                if (dirty) {
                    cerr << "Mark\n";
                    printProbes();
                } else {
                    cerr << "No change\n";
                }
            } else {
                cerr << "Short contig\t" << max_contig << "\n";
            }
			//cerr << "contig: " << max_contig << " @ " << curr.getPos() << " "; curr.dump(-1,seq);
			//cerr << (found?'Y':'N') << (dirty?'Y':'N');
			//cerr << "After:\t" << pos << "\t"; printProbes();
		}

		cerr << id << "\n";
		cerr << "#\tpos\tseq\tch\tGC\trep\twinQ\tdirty\n";
		for (int i=0; i<probes.size(); i++)
            		probes[i].dump(i, seq);
	}
};

char ProbeWindow::strand='+';

class GeneCounter {
	// DNA
	BitVector any;
	ByteMap repeats;

	// RNA
	map<string, size_t> genes;
	map<long, set<size_t> > offTargets;

	char decoded[N+1];
	public:
	GeneCounter() {
		// https://www.geeksforgeeks.org/fast-io-for-competitive-programming/
		ios_base::sync_with_stdio(false);
		cin.tie(NULL);
		decoded[N]='\0';
		assert(sizeof(long)==8);
	}

	void decode(long code) {
		long tmp=code;
		for (int i=0; i<N; i++) {
			int c=tmp&3, off=N-1-i; tmp>>=2;
			switch (c) {
				case 0: decoded[off]='A'; break;
				case 1: decoded[off]='C'; break;
				case 2: decoded[off]='T'; break;
				default: decoded[off]='G'; break;
			}
		}
	}

	static void toUpper(string& seq) {
		string::iterator itr;
		for (itr=seq.begin(); itr!=seq.end(); itr++) {
			*itr=toupper(*itr);
		}
	}

	static void revcmp(string& seq) {
		reverse(seq.begin(), seq.end());
		string::iterator itr;
		for (itr=seq.begin(); itr!=seq.end(); itr++) {
			switch (*itr) {
				case 'A': *itr='T'; break;
				case 'C': *itr='G'; break;
				case 'T': *itr='A'; break;
				case 'G': *itr='C'; break;
				default: *itr='N'; break;
			}
		}
	}

	// RNA
	void process(string& seq, string& gid) {
		if (genes.find(gid)==genes.end()) {
			genes[gid]=genes.size()+1;
		}
		toUpper(seq);
		count(seq, genes[gid]);
	}

	// DNA
	void process(string& seq) {
		toUpper(seq);
		count(seq);
	}

	void processRev(string& seq) {
		revcmp(seq);
		count(seq);
	}

	void count(string& seq, size_t gid=0) {
		size_t n=seq.size();
		if (!n)
			return;
		cerr << n << " bp\n";
		long code=0, len=0, cut=SZ-1;
		for (size_t i=0; i<n; i++) {
			char ch=seq[i];
			switch (ch) {
				case 'A': len++; code<<=2; break;
				case 'C': len++; code<<=2; code|=1; break;
				case 'T': len++; code<<=2; code|=2; break;
				case 'G': len++; code<<=2; code|=3; break;
				default: len=0; code=0; break;
			}
			code&=cut;
			if (len>=N) {
				if (gid) {  // RNA
					offTargets[code].insert(gid);
				} else {    // DNA
					if (!any[code]) {
						any[code]=true;
						continue;
					}
					if (repeats[code]<MAX_CNT) {
						repeats[code]++;
					}
				}
			}
		}
	}

	void update() {
		// save gene/gid & off-targets
		map<string, size_t>::iterator gene_iter=genes.begin();
		while (gene_iter!=genes.end()) {
			cerr << gene_iter->first << "\t" << gene_iter->second << "\n";
			gene_iter++;
		}

		// update any/repeats
		map<long, set<size_t> >::iterator iter=offTargets.begin();
		while (iter!=offTargets.end()) {
			long code=iter->first;
			set<size_t>& sets=iter->second;
			set<size_t>::iterator set_iter=sets.begin();
			cerr << code;
			while (set_iter!=sets.end()) {
				cerr << "\t" << *set_iter;
				set_iter++;
			}   cerr << "\n";

			any[code]=true;
			size_t sz=sets.size();
			if (sz>1) {
				repeats[code]=sz-1;
			}
			iter++;
		}
	}

	void save(const char* file) {
		ofstream ofs(file);
		cerr << "Saving any ...\n";
		any.save(ofs);
		cerr << "Saving repeats ...\n";
		repeats.save(ofs);
		ofs.close();

		ifstream ifs(file);
		cerr << "Checking any ...\n";
		any.load(ifs, true);
		cerr << "Checking repeats ...\n";
		repeats.load(ifs, true);
		ifs.close();
	}

	void load(const char* file) {
		ifstream ifs(file);
		any.load(ifs);
		repeats.load(ifs);
		ifs.close();
	}

	void scan(string& seq, ofstream& ofs) {
		size_t n=seq.size();
		if (!n)
			return;
		cerr << n << " bp\n";
		toUpper(seq);
		vector<unsigned char> rep_cnts(n, 0);

		cerr << "Begin rep_cnts scan\n";
		long code=0, len=0, cut=SZ-1;
		size_t i;
		for (i=0; i<n; i++) {
			char ch=seq[i];
			switch (ch) {
				case 'A': len++; code<<=2; break;
				case 'C': len++; code<<=2; code|=1; break;
				case 'T': len++; code<<=2; code|=2; break;
				case 'G': len++; code<<=2; code|=3; break;
				default: len=0; code=0; break;
			}
			code&=cut;
			if (len>=N) {
				// assert(any[code]);
				rep_cnts[i]=repeats[code];
			}
			ofs << ch << "\t" << (unsigned int)rep_cnts[i] << "\n";
		}
		cerr << "End rep_cnts scan\n";
	}

	void scan(string& seq, string& qual, ofstream& ofs) {
		size_t n=seq.size();
		if (!n)
			return;
		cerr << n << " bp\n";
		toUpper(seq);
		vector<unsigned char> rep_cnts(n, 0);

		cerr << "Begin rep_cnts scan\n";
		long code=0, len=0, cut=SZ-1;
		size_t i;
		for (i=0; i<n; i++) {
			char ch=seq[i];
			switch (ch) {
				case 'A': len++; code<<=2; break;
				case 'C': len++; code<<=2; code|=1; break;
				case 'T': len++; code<<=2; code|=2; break;
				case 'G': len++; code<<=2; code|=3; break;
				default: len=0; code=0; break;
			}
			code&=cut;
			if (len>=N) {
				// assert(any[code]);
				rep_cnts[i]=repeats[code];
			}
			char q='0';
			if (i>=(N-1)) {
				q=qual[i];
				for (size_t j=0; j<N; j++) {
					if (q<qual[i-j])
						q=qual[i-j];
				}
			}
			ofs << ch << "\t" << (unsigned int)rep_cnts[i] << "\t" << short(q-'0') << "\n";
		}
		cerr << "End rep_cnts scan\n";
	}

	void markRep() {
		vector<bool>::iterator itr=any.begin();
		long code=0;
		while (itr!=any.end()) {
			if (*itr) {
				repeats[code]=(unsigned char)MAX_CNT;
			}
			itr++; code++;
		}
	}

	static void design(ifstream& ifs) {
		// cerr << "Begin design\n";
		ProbeWindow probeWindow;
		string ln;
		char ch; short qual; int rep;
		bool less=true;
		while (getline(ifs, ln)) {
			ch=ln[0];
			if (ch=='>' || ch=='@') {
				probeWindow.design();
				probeWindow.init(ln);
				if (ch=='@')
					less=false;
				continue;
			}
			stringstream ss(ln);
			if (less) {
				ss >> ch >> rep;
				probeWindow.append(ch, rep);
			} else {
				ss >> ch >> rep >> qual;
				probeWindow.append(ch, rep, qual);
			}
		}
		probeWindow.design();
	}

	static size_t search(string& chr, string& seq, string& query) {
		if (!seq.size())
			return 0;
		toUpper(seq); toUpper(query);
		size_t found=0, off=0;
		cerr << chr << "\n";
		// cout << "Search " << query << " in " << chr << "\n";
		while ((off=seq.find(query, off)) != string::npos) {
			cout << chr << ":" << off << "\t+\n";
			off++; found++;
		}
		off=0;
		revcmp(seq);
		while ((off=seq.find(query, off)) != string::npos) {
			cout << chr << ":" << off << "\t-\n";
			off++; found++;
		}
		// cout << found << " found\n";
		return found;
	}
};

int usage(string& prog) {
	cerr << "Usage:\n\t" << prog << " build <in:ref.fa> <out:ref.idx>\nor\n\t" << prog << " design <in:ref.idx> <in:rep.idx> <in:input.fa> <out:design.txt>\n";
	return -1;
}

int buildUsage(string& prog) {
	cerr << "Usage:\n\t" << prog << " build <in:ref.fa> <out:ref.idx>\nor\n";
	return -1;
}

int scanUsage(string& prog) {
	cerr << "Usage:\n\t" << prog << " scan <in:ref.idx> <in:rep.idx> <in:input.fa> <out:input.cnt>\n";
	return -1;
}

int designUsage(string& prog) {
	cerr << "Usage:\n\t" << prog << " design <in:input.cnt>\n";
	return -1;
}

int main(int argc, char* argv[]) {
	string prog=argv[0];
	if (argc<2) {
		return usage(prog);
	}
	string task=argv[1];
	bool DNA_PROBE=true;
	if (task.compare("build_rna")==0) {
		DNA_PROBE=false;
		task="build";
	}
	if (task.compare("build")==0) {
		if (argc!=4) { return buildUsage(prog); }
		GeneCounter counter;
		string ref=argv[2], idx=argv[3];
		cerr << "Build " << ref << "\n";
		ifstream ifs(ref.c_str());
		string ln, seq, curr_gene;
		while (getline(ifs, ln)) {
			if (ln[0]=='>') {
				// process
				if (seq.size()) {
					if (DNA_PROBE) {
						counter.process(seq);
						counter.processRev(seq);
					} else {
						assert(curr_gene.size());
						counter.process(seq, curr_gene);
					}
					seq.clear();
				}

				if (!DNA_PROBE) {
					Regex::getGene(ln, curr_gene);
				}
				cerr << ln << "\n";
				continue;
			}
			seq+=ln;
		}

		// process
		if (seq.size()) {
			if (DNA_PROBE) {
				counter.process(seq);
				counter.processRev(seq);
			} else {
				assert(curr_gene.size());
				counter.process(seq, curr_gene);
			}
			seq.clear();
		}

		cerr << "Building " << ref << "\n";
		if (!DNA_PROBE) {
			counter.update();
		}
		counter.save(idx.c_str());
	} else if (task.compare("scan")==0) {
		if (argc!=6) { return scanUsage(prog); }
		string ref=argv[2], rep=argv[3], input=argv[4], cnts=argv[5];

		GeneCounter counter;
		cerr << "Loading repeats " << rep << "\n";
		counter.load(rep.c_str());
		counter.markRep();

		cerr << "Loading " << ref << "\n";
		counter.load(ref.c_str());

		cerr << "Scanning " << input << "\n";
		ifstream ifs(input.c_str());
		assert(ifs.good());
		ofstream ofs(cnts.c_str());
		string ln, seq;
		getline(ifs, ln);
		if (!ln.size())
			return scanUsage(prog);
		if (ln[0]=='>') {
			ofs << ln << "\n";
			while (getline(ifs, ln)) {
				if (ln[0]=='>') {
					counter.scan(seq, ofs);
					seq.clear();
					ofs << ln << "\n";
					continue;
				}
				seq+=ln;
			}
			if (seq.size())
				counter.scan(seq, ofs);
		} else if (ln[0]=='@') { // fastq format
			string seq, plus, qual;
			do {
				assert(ln[0]=='@');
				ofs << ln << "\n";
				getline(ifs, seq);
				getline(ifs, plus); assert(plus[0]=='+');
				getline(ifs, qual);
				counter.scan(seq, qual, ofs);
			} while (getline(ifs, ln));
		} else {
			return scanUsage(prog);
		}
		ofs.close();
	} else if (task.compare("design")==0) {
		pref_gc=atoi(argv[2]);
		// cerr << "GC: " << pref_gc << "\n";
		assert(pref_gc>0 && pref_gc<WIN);
		string cnts=argv[3];
		ifstream ifs(cnts.c_str());
		GeneCounter::design(ifs);
		ifs.close();
	} else if (task.compare("search")==0){
		assert(argc==4);
		ifstream ref(argv[2]);
		string query=argv[3], ln, chr, seq;
		while (getline(ref, ln)) {
			if (ln[0]=='>') {
				GeneCounter::search(chr, seq, query);
				seq.clear();
				chr=ln.substr(1);
				continue;
			}
			seq+=ln;
		}
		GeneCounter::search(chr, seq, query);
		ref.close();
	} else {
		return usage(prog);
	}

	return 0;
}
