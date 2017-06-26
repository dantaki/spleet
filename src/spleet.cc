#include "aux.h"
#include <iostream>
#include <fstream>
using namespace BamTools;
using namespace std;
int main(int argc, char *argv[])
{
	string splash=		"\neeeee   eeeee   e         eeee   eeee   eeeee\n" 
				"8   \"   8   8   8         8      8        8   \n"
				"8eeee   8eee8   8e        8eee   8eee     8e  \n"
				"   88   88      88        88     88       88  \n"
				"8ee88   88      88eee     88ee   88ee     88  \n\n"
				"spleet         Validate SV with Split-Reads\n"
				"Version: 1.0	Author: Danny Antaki <dantaki@ucsd.edu>\n\n"
				"Usage: spleet -i <in.bam> -r <sv.bed> -x <FLOAT> -q <INT> -o <output.txt>\n\nOptions:\n"
				"    -i        Input: BAM filename\n"
				"    -r        Input: SV bed file\n"
				"    -x        Minimum reciprocal overlap [0.5]\n" 
				"    -q        Mapping quality threshold [10]\n"
				"    -o        Output: filename\n";
	string ifh; string bed; string ofh; int Q=10; double OVR=0.5;
	if ( (argc==1) ||
	     (argc==2 && string(argv[1]) == "-h") ||
	     (argc==2 && string(argv[1]) == "-help") ||
	     (argc==2 && string(argv[1]) == "--help"))
		{ cerr << splash << endl; return 1;}
	for(int i=1; i<argc; i++){
		if(string(argv[i]) == "-i" || string(argv[i]) == "--in" || string(argv[i]) == "-in"){ ifh = string(argv[i+1]); i++; continue; }
		if(string(argv[i]) == "-l" || string(argv[i]) == "-r" || string(argv[i]) == "--l" || string(argv[i])=="--r" || string(argv[i])=="-bed" || string(argv[i])=="--bed"){ bed=string(argv[i+1]); i++; continue; }
		if(string(argv[i]) == "-o" || string(argv[i]) == "--out" || string(argv[i]) == "-out"){ ofh = string(argv[i+1]); i++; continue; }
		if(string(argv[i]) == "-q" || string(argv[i]) == "--q"){ Q = atoi(argv[i+1]); i++; continue; }
		if(string(argv[i]) == "-x" || string(argv[i]) == "--x"){ OVR = ::atof(argv[i+1]); i++; continue; }
		cerr << "ERROR: Unknown option "<<string(argv[i])<< '\n' << splash <<endl;
		return 1;
	}
	if(ifh == "") { cerr << "ERROR: No BAM file given"<< endl; return 1; }
	if(bed == "") { cerr << "ERROR: No BED file given"<< endl; return 1; }
	if(ofh == ""){ ofh = ifh.substr(0,ifh.find_last_of('.'))+"_spleet.txt";  }		
	ifstream fin(bed.c_str());
	BamReader bam;
	if(!fin) { cerr << "ERROR: Cannot open " << ifh << endl; return 1; } 
	if (!bam.Open(ifh)){cerr << "ERROR: " << ifh << " could not be opened!"<< endl; return 1; }
	BamAlignment al;
	const RefVector refs = bam.GetReferenceData();
	map<string, int> chrom = hashRef(refs);
	bool chrFlag=0; if (!refs[0].RefName.compare(0, 3, "chr")){ chrFlag=1; } //chrFlag=true when reference is prefixed with "chr"  
	ofstream out(ofh.c_str());
	out << "CHR\tSTART\tEND\tLENGTH\tOVERLAP\tREADNAME\tALIGNMENT\tSTRANDS\tSV\tTYPE" << endl;
	string line; while(getline(fin,line,'\n'))
	{
		vector<string> pos = split(line,'\t');
		if(pos.size()<4){ cerr << "ERROR: Malformed BED record " << line << endl; } 
		int32_t start = atoi(pos[1].c_str()); int32_t end = atoi(pos[2].c_str()); 	
		string CHR=pos[0]; string TYPE=pos[3];
		if(chrFlag==false && !CHR.compare(0,3,"chr")) { CHR.erase(0,3); }
		if(chrom.count(CHR)==0) { cerr << "ERROR: " << CHR << " not found in BAM file reference" << endl; return 1; }
		if(!bam.LocateIndex()) { cerr << "ERROR: Cannot find index for BAM" << endl; return 1; }
		if(!bam.SetRegion(chrom[CHR],start-10001,chrom[CHR],end+10001)) { cerr << "ERROR: Is your BAM file indexed?" << endl; return 1; } 
		char SV = '-'; // DEL == -; DUP == + 
		if(TYPE.find("DEL")){ SV='+';}
		while (bam.GetNextAlignmentCore(al))
		{
			if(al.MapQuality < Q || 
				al.IsDuplicate()==true || 
				al.IsFailedQC()==true || 
				al.IsMapped()==false )
				{ continue; } 	
			al.BuildCharData();
			if (!al.HasTag("SA")){ continue; } 
			string aln = "PRIMARY";
			if(al.AlignmentFlag & 2048) { aln = "SUPPLEMENTARY"; }
			if(al.AlignmentFlag & 256) { aln = "SECONDARY"; }
			string primChrom = refs[al.RefID].RefName;
			char primOri='+'; if(al.IsReverseStrand()) { primOri='-'; }
			int32_t primLeft = al.Position;
			int32_t primRight = rightPos(primLeft,al.CigarData);
			string tags; al.GetTag("SA",tags);
			vector<string> splitTag = split(tags,';');
			for(unsigned int i=0; i!=splitTag.size(); ++i){	
				vector<string> splitAln= split(splitTag[i],',');
				string splitChrom = splitAln[0];
				if(splitChrom != primChrom) { continue; } 
				int32_t splitLeft = atoi(splitAln[1].c_str());
				char splitOri = '+';  if(splitAln[2] == "-") { splitOri='-'; }
				if((SV=='-' || SV == '+') && primOri!=splitOri) { continue; }
				string splitCig = splitAln[3]; vector<CigarOp> splitCigar = rollCig(splitCig);
				int32_t splitRight = rightPos(splitLeft,splitCigar);
				primLeft++; primRight++; //convert to 1-base
				int32_t brks[4] = {primLeft, primRight, splitLeft, splitRight};
				sort(brks,brks+4);
				int32_t brkStart=0; int32_t brkEnd=0;
				if(SV == '-' && brks[2]-brks[1]+1 > 2) { brkStart=brks[1]; brkEnd=brks[2]; brkStart++; brkEnd--; } 
				else if(SV == '-' && brks[2]-brks[1]+1 <=2) { continue; }
				else if(SV == '+' && primRight <= splitLeft) { brkStart=brks[0]; brkEnd=brks[3]; }
				else if(SV == '+' && primRight > splitLeft) { brkStart=splitLeft; brkEnd=primRight; }
				if(brkStart==0 || brkEnd==0){continue; }
				float ovr = overlap(start,end,brkStart,brkEnd);
				if(ovr < OVR) { continue; }
				out << primChrom << '\t' << brkStart << '\t' << brkEnd << '\t' << brkEnd-brkStart+1 << '\t' << ovr
					<< '\t' << al.Name << '\t' << aln << '\t' << primOri << '|' << splitOri << '\t' 
					<< CHR << ':' << start << '-' << end << '\t' << TYPE << endl; 
			}
		}	
		
	}
	bam.Close(); fin.close(); out.close();
	return 0;
}
