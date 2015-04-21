/**** system headers ****/
#include <iostream>
#include <iomanip>
#include <fstream>
#include <math.h>
#include <complex>
#include <algorithm>
#include <map>
#include <string>
#include <vector>
#include <unistd.h>
#include <zlib.h>
//#include <sys/types.h>
//#include <cstdlib>
using namespace std;

/**** user headers ****/
#include "functions.h"
#include "kseq.h"
KSEQ_INIT(gzFile, gzread)

/*ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*
 * DNA:
 * Nucleotide Code:  Base:
 * ----------------  -----
 * A.................Adenine
 * C.................Cytosine
 * G.................Guanine
 * T (or U)..........Thymine (or Uracil)
 * R.................A or G
 * Y.................C or T
 * S.................G or C
 * W.................A or T
 * K.................G or T
 * M.................A or C
 * B.................C or G or T
 * D.................A or G or T
 * H.................A or C or T
 * V.................A or C or G
 * N.................any base
 * . or -............gap
 * 
 * Note: b1 is supposed to be from read, b2 from primer
 */
bool base_equal(const char b1, const char b2) {
  
  static std::map <char, string> iupac;
  iupac['A'] = "A";
  iupac['C'] = "C";
  iupac['G'] = "G";
  iupac['T'] = "T";
  iupac['R'] = "AG";
  iupac['Y'] = "CT";
  iupac['S'] = "CG";
  iupac['W'] = "AT";
  iupac['K'] = "GT";
  iupac['M'] = "AC";
  iupac['B'] = "CGT";
  iupac['D'] = "AGT";
  iupac['H'] = "ACT";
  iupac['V'] = "ACG";
  iupac['N'] = "ACGT";
  
  if ( b1=='N' ) return(false); // no call base
  
  if ( b1==b2 ) return(true); // base equal
  
  if ( b2=='?' || b2=='.' || b2=='N') return(true); // wildcard chracters
  
  if ( iupac.find(b2) == iupac.end() ) {
    cerr << "Non IUPAC code found: " << b2 << endl;
    exit(1);
  }
  
  return( iupac[b2].find(b1) != std::string::npos );
}

/*ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/* map primer or linker at the front of a read
 * @parameters:
 * s1: the front part of a read, must be longer than primer or linker
 * s2: primer or linker sequence
 * @returns:
 * qend: the right most position where primer or linker map, 0 based
 *       if qend<0, qend is not calculated
 * return: edit_distance
 */
void map_with_edit_distance_3( const std::string& s1, const std::string& s2, 
			       const int threshold,
			       int& q0, int& qend, int& edit_distance )
{
  q0=-1;
  qend=-1;
  edit_distance=s1.length();
  bool goodMap=true;
  
  //int maxindel=max(1, (int)s2.length()/10);
  int maxindel=1;
  
  int q_cur=0, ed=0;
  
  // direct character comparison 
  vector<bool> basematch(s2.length(), false);
  vector<int> baseidx(s2.length(), 0);
  for(size_t i=0; i<s2.length(); ++i) {
    baseidx[i]=i;
    basematch[i]=base_equal(s1[i], s2[i]); 
    ed += (!base_equal(s1[i], s2[i])) ; 
  }
  
  // find longest mismatch and try deletion or insertion
  bool isthereindel=true;
  int m0=0, mend=0;
  for(size_t i=0; i<s2.length(); ++i) {
    if ( basematch[i] ) continue;
    size_t j=i; 
    for(j=i; j<s2.length(); ++j) if ( basematch[j] ) { j-=1; break;}
    if ( j-i > mend-m0 ) {m0=i; mend=j;}
  }
  if ( mend-m0+1 < 4 ) isthereindel=false;
  
  if ( ! isthereindel ) {
    edit_distance=ed;
    q0=0; while( ! basematch[q0] && q0<s2.length()-1 ) q0+=1;
    qend=s2.length()-1; while( ! basematch[qend] && qend>1 ) qend-=1;
    return;
  }
  
  
  vector<bool> matched(s2.length()+1, true);
  vector<int> pos(s2.length()+1, 0);
  for(size_t i=0; i<s2.length(); ++i) {
    matched[i+1]=false;
    if ( s2[i]=='N' ) {
      matched[i+1]=true;
      q_cur = i+1;
    }
    int dx=0;
    for( dx = matched[i] ? 0 : -maxindel; dx<maxindel; ++dx) {
      if ( q_cur+dx <0 || q_cur+dx>=s1.length() ) continue;
      if ( base_equal(s1[q_cur+dx], s2[i]) ) { 
	matched[i+1]=true;
	break; 
      }
    }
    if ( matched[i+1] ) {
      ed += abs(dx)*matched[i]*(dx>0);
      pos[i] = q_cur+dx;
      q_cur += dx+1; 
    }
    else { ed+=1; pos[i] = q_cur; q_cur+=1; }
    
  }
  
  q0=pos[0];
  qend=q_cur-1;
  edit_distance=ed;
  
  // cerr << s1 << "\t" << s1.length() << "\n" 
  //    << s2 << "\t" << s2.length() << "\t" << qend << "\t" << ed << endl;
  return;
  
  //cerr << s1 << "\t" << s1.length() << "\n" 
  //    << s2 << "\t" << s2.length() << "\t" << qend << "\t" << p[n2] << endl;
  //for( i=0; i<=n1; ++i ) cerr << c3[i] << " "; cerr << "\n-----" << endl;
  
}

/*ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/* map primer or linker at the front of a read
 * using Levenshtein's Edit Distance
 * @parameters:
 * s1: the front part of a read, must be 2+ bases longer than primer or linker
 * s2: primer or linker sequence
 * @returns:
 * q0:   where primer start to match, 0 based
 * qend: the right most position where primer or linker map, 0 based
 *       if q0<0 || qend<0, mapping failed
 * edit_distance
 */
void map_with_edit_distance(const std::string& s1, const std::string& s2, 
			    const int threshold,
			    int& q0, int& qend, int& edit_distance )
{
  q0=-1;
  qend=-1;
  edit_distance=s2.length();
  bool goodMap=true;
  
  const int cost_del = 1;
  const int cost_ins = 1;
  const int cost_sub = 1;
  
  int c1[4096];
  int c2[4096];
  int c3[4096];
  int* p = c1;
  int* q = c2;
  int* r;
  
  int i,j;
  int n1 = s1.length();
  int n2 = s2.length();
  if ( n1<n2 ) {
    cerr << "string1 should be equal or longer than string2" << endl;
    exit(1);
  }
  if ( n1 > 4096 || n2 > 4096 ) {
    cerr << "map_with_edit_distance(), change size of array and make\n";
    exit(1);
  }
  
  // just check a little bit longer on S1
  n1=n2+2;
  
  for(i=0; i<=n1; ++i) c3[i]=n1;
  c3[0]=0;
  
  p[0] = 0;
  for( j = 1; j <= n2; ++j ) p[j] = p[j-1] + cost_ins;
  
  for( i = 1; i <= n1; ++i ){
    int edMin=n1;
    q[0] = p[0] + cost_del;
    for( j = 1; j <= n2; ++j ) {
      int d_del = p[j] + cost_del;
      int d_ins = q[j-1] + cost_ins;
      // int d_sub = p[j-1] + ( s1[i-1] == s2[j-1] ? 0 :cost_sub );
      int d_sub = p[j-1] + ( base_equal(s1[i-1], s2[j-1]) ? 0 : cost_sub );
      q[j] = std::min( std::min( d_del, d_ins ), d_sub );
      if ( q[j] < edMin ) edMin=q[j];
    }
    c3[i]=edMin;
    r = p;
    p = q;
    q = r;
    
    if ( i<n2 && edMin>threshold+1 ) { goodMap=false; break; }
  }
  
  if ( goodMap ) {
    // return values
    edit_distance=p[n2];
    qend=n1;
    for( qend=n1; qend>1; --qend) if ( c3[qend]==c3[qend-1] ) break;
    edit_distance=c3[qend];
    if ( qend>0 ) qend-=1; // convert to 0 based index
    for( q0=0; q0<n2; ++q0) if ( c3[q0]==c3[q0+1] ) break;
  }
  
  return;
  
  //cerr << s1 << "\t" << s1.length() << "\n" 
  //    << s2 << "\t" << s2.length() << "\t" << qend << "\t" << p[n2] << endl;
  //for( i=0; i<=n1; ++i ) cerr << c3[i] << " "; cerr << "\n-----" << endl;
  
}

void map_with_edit_distance_2( const std::string& s1, const std::string& s2, 
			       const int threshold,
			       int& q0, int& qend, int& edit_distance )
{
  q0=-1;
  qend=-1;
  edit_distance=s1.length();
  bool goodMap=true;
  
  //int maxindel=max(1, (int)s2.length()/10);
  int maxindel=1;
  
  int q_cur=0, ed=0;
  
  // direct character comparison 
  vector<bool> basematch(s2.length(), false);
  vector<int> baseidx(s2.length(), 0);
  vector<int> edidx(s2.length(), 0);
  for(size_t i=0; i<s2.length(); ++i) {
    baseidx[i]=i;
    basematch[i]=base_equal(s1[i], s2[i]); 
    ed += ! base_equal(s1[i], s2[i]) ; 
    edidx[i]=ed;
  }
  
  // return if match well
  if ( ed<=2 ) {
    edit_distance=ed;
    q0=0; while( ! basematch[q0] && q0<s2.length()-1 ) q0+=1;
    qend=s2.length()-1; while( ! basematch[qend] && qend>1 ) qend-=1;
    return;
  }
  
  for(size_t i=0; i<s2.length(); ++i) {
    if ( basematch[i] ) continue;
    
    // try del on S1, S1 missing s2[i]
    size_t ed_del=edidx[i];
    for(size_t j=i+1; j<s2.length(); ++j) {
      ed_del += ! base_equal(s1[j-1], s2[j]); 
      cerr << j << "\t" << s1[j-1] << "\t" << s2[j] << "\t" << ed_del << endl;
      if (ed_del>threshold) break;
    }
    cerr << "---------" << endl;
    // try insertion on S1, S1 add one more at s2[i]
    size_t ed_ins=edidx[i];
    for(size_t j=i; j<s2.length() && j+1<s1.length(); ++j) {
      ed_ins += ! base_equal(s1[j+1], s2[j]); 
      cerr << j << "\t" << s1[j+1] << "\t" << s2[j] << "\t" << ed_ins << endl;
      if (ed_ins>threshold) break;
    }
    if ( min(ed_del, ed_ins)<threshold ) {
      cerr << i << "\t" << edidx[i] << "\t" << ed << endl;
      cerr << ed_del << "\t" << ed_ins << endl;
      edit_distance=min(ed_del, ed_ins);
      q0=0;
      qend = ed_del<ed_ins ? s2.length()-1-1 : s2.length()-1+1; 
      return;
    }
  }
  
  // single base indel failed, now use edit_distance
  map_with_edit_distance(s1, s2, threshold, q0, qend, edit_distance );
  
  // cerr << s1 << "\t" << s1.length() << "\n" 
  //    << s2 << "\t" << s2.length() << "\t" << qend << "\t" << ed << endl;
  return;
  
  //cerr << s1 << "\t" << s1.length() << "\n" 
  //    << s2 << "\t" << s2.length() << "\t" << qend << "\t" << p[n2] << endl;
  //for( i=0; i<=n1; ++i ) cerr << c3[i] << " "; cerr << "\n-----" << endl;
  
}


/*ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
int main_usage() 
{
  cerr << "Usage:\n"
       << "   findprimer -f fastq -p primer -m n -o outFile\n"
       << "\nOptions:\n"
       << "   fastq    fastq file\n"
       << "   primer   primer or linker sequence\n"
       << "   n        allowed base mismatches, optional, default=1+primer/20\n"
       << "   outFile  output folder, optional, default=.\n"
       << "\nOutput: rows of the following columns\n"
       << "   sequence_name\n"
       << "   sequence_name_comment\n"
       << "   primer/linker\n"
       << "   primer/linker_length\n"
       << "   primer/linker_on_sequence\n"
       << "   edit distance\n"
       << "   primer_start\t#0 based\n"
       << "   primer_end\t#0 based\n"
       << "\nNote:\n"
       << "   1. primer is converted to upper cases.\n"
       << "   2. sequence is not converted.\n"
       << "   3. N, ?, . in primer matches all.\n"
       << "   4. N in sequence does not matche any.\n"
       << endl;
  return 0;
}

/*ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
int main(int argc, char* argv[])
{
  if ( argc < 2 ) return main_usage();
  
  int res;
  string mycommand="";
  string primer="", R1File="", outFile="";
  int nMismatch=1; 
  long int BUFFERSIZE=(int)600E6; // buffer size to hold reads, two buffers used
  
  /*  read in parameters
   *  R1File: read1 fastq file
   *  primer: primer or linker sequence
   *  nMismatch: allowed base mismatches 
   *  outFile: output file, default to current .
   */
  vector<string> ARGV(0);
  for(int i=0;i<argc;++i) ARGV.push_back(string(argv[i]));
  for(int i=1;i<(int)ARGV.size();++i) {
#define _next2 ARGV[i]=""; ARGV[i+1]=""; continue;
#define _next1 ARGV[i]=""; continue;
    if ( ARGV[i]=="-p" ) { 
      primer=ARGV[i+1]; 
      primer=to_upper(primer);
      nMismatch=1+primer.length()/10;
      _next2; 
    }
    if ( ARGV[i]=="-f" ) { 
      R1File=ARGV[i+1]; 
      outFile=R1File+".primer";
      _next2; 
    }
    if ( ARGV[i]=="-m" ) { nMismatch=atoi(ARGV[i+1].c_str()); _next2; }
    if ( ARGV[i]=="-o" ) { outFile=ARGV[i+1]; _next2; }
  }
  
  cerr << "Reads:\t" << R1File << "\n"
       << "Primer:\t" << primer << "\n"
       << "misMatch:\t" << nMismatch << "\n"
       << "Output:\t" << outFile << "\n"
       << endl;
  
  //string cmd="mkdir -p " + folder;
  //res = system( cmd.c_str() );
  
  /*  1. read in barcode and sample id from barcodeFile 
   *     barcode: vector
   *     sampleid: vector
   *     barcode and sampleid are row matched
   *  2. check edit distance between barcodes
   *  3. check conflicks
   */
  /*
  string s1= "GAAAATCTCTAGCAGT";
  string s2= "GAAATCTCTAGCAGT";
  string s3= "GAAAAGTCTCTAGCAGT";
  string s4= "AGAAAATCTCTAGCAGT";
  string s5= "TAAAATCTCTAGCAGT";
  string p1= "GAAAATCTCTAGCA";
  
  int q0, qend, ed;
  map_with_edit_distance_2(s1, p1, 2, q0, qend, ed);
  map_with_edit_distance_2(s2, p1, 2, q0, qend, ed);
  map_with_edit_distance_2(s3, p1, 2, q0, qend, ed);
  map_with_edit_distance_2(s4, p1, 2, q0, qend, ed);
  map_with_edit_distance_2(s5, p1, 2, q0, qend, ed);
  //return 0;
  */

  ofstream FOUT(outFile.c_str());

  gzFile fp;
  kseq_t *seq;
  int l;
  fp = gzopen(R1File.c_str(), "r");
  if(fp==Z_NULL) { cerr << "Can't open " << R1File << endl; exit(1); }
  seq = kseq_init(fp);
  size_t icount=0;
  while ((l = kseq_read(seq)) >= 0) {
    
    string read = string(seq->seq.s);
    
    int p0, pend,editDistance;
    map_with_edit_distance(read, primer, nMismatch, 
			   p0, pend, editDistance);
    
    //cerr << read.substr(0, primer.length()+4 )  << "\n"
    //	 << primer << "\t" << p0 << "\t" << pend << "\t" << editDistance
    //	 << "\n-------" << icount << endl;
    
    string tmps = "@"+string(seq->name.s);
    if (seq->comment.l) tmps += " " + string(seq->comment.s);
    tmps += "\n";
    tmps += string(seq->seq.s)+"\n+\n";
    if (seq->qual.l) tmps += string(seq->qual.s)+"\n";
    
    FOUT << "@"+string(seq->name.s) << "\t"
	 << string(seq->comment.s) << "\t"
	 << primer << "\t"
	 << p0 << "\t" 
	 << pend << "\t" 
	 << editDistance
	 << "\n";
    
    
    ++icount;
    if ( icount%1000000==0 ) cerr << icount << endl;
    // if ( icount > 100 ) break;
  }
  
  gzclose(fp); //close input
  
  
  kseq_destroy(seq);
  return 0;
} 
/*ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/

