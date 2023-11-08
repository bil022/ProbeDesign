#include <iostream>
#include <fstream>
#include <sstream>
#include <cstdlib>
#include <cstring>
#include <set>
#include <map>
#include <vector>
#include <cassert>
#include <algorithm>

using namespace std;

class Bed {
public:
	Bed(int h, int t, char s):head(h), tail(t), strand(s) {
	}
	Bed(int h, int t):head(h), tail(t), strand('+') {
		if (head>=tail) {
			int tmp=head;
			head=tail;
			tail=tmp;
			strand='-';
		}
        }
	int head, tail; char strand;
	bool operator< (const Bed& b) const {
		if (strand=='+') {
			if (head!=b.head)
				return head<b.head;
			return tail<b.tail;
		}

		if (head!=b.head)
			return head>b.head;
		return tail>b.tail;
	}
};

class BedGraph : public Bed {
public:
	BedGraph(Bed& b, int v):Bed(b),v(v) {
	}
	int v;
	bool operator< (const BedGraph& b) const {
		if (v!=b.v)
			return v>b.v;
		return (Bed)(*this)<(Bed)b;
	}
};

class Pileup {
public:
  Pileup() {
    reset();
  }
  void reset() {
    segs.clear();
    pts.clear();
    strand='?';
  }
  bool add(Bed b) {
    segs.push_back(b);
    if (segs.size()==1) {
      strand=b.strand;
    } else {
      if (b.head==b.tail) {
        b.strand=strand;
      }
      if (strand!=b.strand) {
	cerr << "?[" << strand << "] vs [" << b.head << "-" << b.tail << b.strand << "]\n";
        return false;
      }
    }
    pts.insert(b.head); pts.insert(b.tail);
    return true;
  }
  void report(string& chr) {
    // cout << "Begin " << chr << "\n";
    map<Bed, int> sum;
    vector<Bed>::iterator itr_seg=segs.begin();
    char strand=itr_seg->strand;
    for ( ;itr_seg!=segs.end(); itr_seg++) {
      Bed& b=*itr_seg;
      set<int>::iterator low=pts.lower_bound(b.head);
      set<int>::iterator high=pts.upper_bound(b.tail);
      set<int>::iterator itr=low;
      int left=*itr; itr++;
      while (itr!=high) {
        int right=*itr;
       	sum[Bed(left, right, strand)]++;
        left=right; itr++;
      }
    }
 
    set<BedGraph> bg; 
    map<Bed, int>::iterator itr_bed=sum.begin();
    while (itr_bed!=sum.end()) {
      Bed b=itr_bed->first;
      int v=itr_bed->second;
      bg.insert(BedGraph(b, v));
      itr_bed++;
    }

    set<BedGraph>::iterator itr_bg=bg.begin();
    while (itr_bg!=bg.end()) {
      BedGraph bg=*itr_bg;
      cout << chr << "\t" << bg.head << "\t" << bg.tail << "\t" << bg.v << endl;
      itr_bg++;
    }

    reset();
    // cout << "End " << chr << "\n";
  }
private:
  vector<Bed> segs;
  set<int> pts;
  char strand;
} pileup;

void pileup_bed() {
  string ln, chr, prev; int h, t;
  while (getline(cin, ln)) {
    stringstream ss(ln);
    ss >> chr >> h >> t;
    assert(!ss.fail());
    if (chr.compare(prev)) {
      if (prev.size()) {
        pileup.report(prev);
      }
      prev=chr;
    }
    if (!pileup.add(Bed(h,t))) {
      cerr << ln << endl;
      assert(false);
    }
  }
  if (chr.size()) {
    pileup.report(chr);
  }
}

map<string, string> fa;
void load_fa(const char* file) {
  ifstream ifs(file); assert(ifs.good());
  string ln, chr, seq;
  while (getline(ifs, ln)) {
    assert(ln.size());
    if (ln[0]=='>') {
      if (seq.size()) {
         assert(chr.size());
         fa[chr]=seq;
         seq.clear();
         cerr << chr << "\t" << fa[chr].size() << "\n";
      }
      chr=ln.substr(1);
    } else {
      seq+=ln;
    }
  }
  if (seq.size()) {
    assert(chr.size());
    fa[chr]=seq;
    seq.clear();
    cerr << chr << "\t" << fa[chr].size() << "\n";
  }
}

// https://codereview.stackexchange.com/questions/201990/dna-reverse-complement-as-fast-as-possible
class DNA {
public:
  static void reverse_complement(string& seq) {
    string::iterator front = seq.begin();
    string::iterator back = front + seq.size() - 1;
    while (back > front) {
        if (*front == '\n') {
            ++front;
        } else if (*back == '\n') {
            --back;
        } else {
            *back = complement(*back);
            *front = complement(*front);
            std::swap(*back, *front);
            ++front;
            --back;
        }
    }
  } 
private:
  static char complement(char ch) {
    switch (toupper(ch)) {
      case 'A': ch = 'T'; break; 
      case 'C': ch ='G'; break;  
      case 'T': ch ='A'; break;  
      case 'G': ch ='C'; break;  
      case 'N': ch ='N'; break;  
      case 'U': ch ='A'; break;  
      case 'M': ch ='K'; break;  
      case 'R': ch ='Y'; break;  
      case 'W': ch ='W'; break;  
      case 'S': ch ='S'; break;  
      case 'Y': ch ='R'; break;  
      case 'K': ch ='M'; break;  
      case 'V': ch ='B'; break;  
      case 'H': ch ='D'; break;  
      case 'D': ch ='H'; break;  
      case 'B': ch ='V'; break;
      default: assert(false);
    }
    return ch;
  }
};

char op_replace(char a) { if (a==':') return ' '; return a; }
char UC(char a) { return toupper(a); }

void pileup_fa() {
  string ln, chr, gid; char strand; int h, t, v;
  while (getline(cin, ln)) {
    transform(ln.begin(), ln.end(), ln.begin(), op_replace);
    stringstream ss(ln);
    ss >> chr >> strand >> gid >> h >> t >> v;
    if (ss.fail()) { cerr << ln << "\n"; }
    assert(!ss.fail());
    if (fa.find(chr)==fa.end()) cerr << chr << endl;
    assert(fa.find(chr)!=fa.end());
    assert(h && t>h);
    string& seq=fa[chr]; assert(seq.size()>t);
    string exon=seq.substr(h-1, t+1-h);
    cout << ">" << chr << ":" << h << "-" << t << ":" << exon.size() ;
    cout  << ":" << gid << ":" << v << ":" << strand << "\n";
    transform(exon.begin(), exon.end(), exon.begin(), UC);
    if (strand=='-') {
      DNA::reverse_complement(exon);
    }
    cout << exon << "\n";
  }
}

int main(int argc, char*argv[]) {
  string task;
  if (argc>1)
    task=argv[1];
  if (task.compare("pileup")==0) {
    pileup_bed();
  } else if (task.compare("pileup_fa")==0 && argc>2) {
    load_fa(argv[2]);
    pileup_fa();
  } else {
    cerr << "Usage:\n\t" << argv[0] << " pileup\n";
    cerr << "\t" << argv[0] << " pileup_fa\n";
  }

  return 0;
}
