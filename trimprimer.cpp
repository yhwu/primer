/**** system headers ****/
#include <iostream>
#include <iomanip>
#include <fstream>
#include <math.h>
#include <complex>
#include <algorithm>
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
bool base_equal(const char b1, const char b2) {
  if ( b1=='N' ) return(false);
  if ( b2=='N' ) return(true);
  else return( b1==b2 );
}


/*ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
int main_usage() 
{
  cerr << "Usage:\n"
       << "   trimprimer -f fastq -t trimfile -o outFile\n"
       << "\nOptions:\n"
       << "   fastq    fastq file\n"
       << "   trimfile primer or linker sequence\n"
       << "   outFile  output folder, optional, default=.\n"
       << "\nNote: the trimfile should contain rows with the following field\n"
       << "   sequence_name  start with @\n"
       << "   start          ignored, always cut from 0\n"
       << "   end            end is the last chracter removed, 0 based\n"
       << "\nNote:\n"
       << "   sequence_name does not contain sequence comment, so be careful\n"
       << "   with paired end reads.\n" 
       << "   sequence not in the trim file will be discarded\n"
       << "   sequence names must be in the same order as in the fastq file\n"
       << endl;
  return 0;
}


/*ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
int main(int argc, char* argv[])
{
  if ( argc < 2 ) return main_usage();
  
  int res;
  string mycommand="";
  string R1File="", trimFile="", outFile="";
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
    if ( ARGV[i]=="-f" ) { R1File=ARGV[i+1]; 
      outFile=R1File+".trimmed";
      _next2; 
    }
    if ( ARGV[i]=="-t" ) { trimFile=ARGV[i+1]; _next2; } 
    if ( ARGV[i]=="-o" ) { outFile=ARGV[i+1]; _next2; }
  }
  
  cerr << "Reads:\t" << R1File << "\n"
       << "trimFile:\t" << trimFile << "\n"
       << "Output:\t" << outFile << "\n"
       << endl;
  
  // initialize to read fastq file
  gzFile fp;
  kseq_t *seq;
  int l;
  fp = gzopen(R1File.c_str(), "r");
  if(fp==Z_NULL) { cerr << "Can't open " << R1File << endl; exit(1); }
  seq = kseq_init(fp);
  
  // initialize zout 
  gzFile zout;
  zout = gzopen(outFile.c_str(), "wb");
  if (!zout) { cerr << "Can't open " << outFile << endl; exit(0); }
  
  ifstream FIN(trimFile.c_str());
  if (!FIN) { cerr << "Can't open " << trimFile << endl; exit(0); }
  string trim_name;
  int q0;
  int qend;
  int ed;
  size_t icount=0;
  string outBuffer="";
  while( FIN >> trim_name >> q0 >> qend >> ed ) {
    if ( trim_name[0]=='#' ) continue;
    bool ifdrop=false;
    if ( q0<0 || qend<0 ) ifdrop=true;
    if ( ed> qend/2 ) ifdrop=true;
    if ( ifdrop ) {
      cerr << "read is dropped:\t" << trim_name << "\t" << q0 << "\t" << qend << "\t" << ed << endl;
      continue;
    }
    
    // read fastq file to find read
    bool found=false;
    while ((l = kseq_read(seq)) >= 0) {
      string seq_name = "@"+string(seq->name.s);
      if ( seq_name != trim_name ) continue;
      
      found=true;
      string tmps = seq_name;
      if (seq->comment.l) tmps += " " + string(seq->comment.s);
      tmps += "\n";
      tmps += string(seq->seq.s).substr(qend+1) + "\n+\n";
      if (seq->qual.l) tmps += string(seq->qual.s).substr(qend+1) +"\n";
      outBuffer+=tmps;
      if ( outBuffer.length()>60000000 ) {
	gzwrite(zout, &outBuffer[0], outBuffer.size() );
	outBuffer="";
      }
      break;
    }
    
    if ( !found ) cerr << trim_name << "  NOT FOUND" << endl;
    if ( ++icount%1000000 == 0 ) cerr << icount << " processed" << endl;
  }
  if ( outBuffer.length()>1 ) {
    gzwrite(zout, &outBuffer[0], outBuffer.size() );
    outBuffer="";
  }
  gzclose(zout); //close input
  gzclose(fp); //close input
  
  
  kseq_destroy(seq);
  return 0;
} 
/*ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/

