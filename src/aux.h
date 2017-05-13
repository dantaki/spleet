#ifndef AUX_H
#define AUX_H
#include "api/BamReader.h"
#include <string>
#include <map>
#include <vector>
#include <algorithm>
#include <sstream>
inline std::vector<std::string> split(std::string &s, char delim)
{
  std::vector<std::string> v;
  std::stringstream ss(s); std::string i;
  while(std::getline(ss, i, delim)) { v.push_back(i); }
  return v;
}
inline int32_t rightPos(int32_t &pos, std::vector<BamTools::CigarOp> &cigar)
{
        int len=0;
        for(std::vector<BamTools::CigarOp>::iterator it=cigar.begin(); it != cigar.end(); ++it){
                if(it->Type == 'D' || it->Type == 'M' || it->Type == '=' || it->Type == 'X') { len+=it->Length; }
        }
        len+=pos;
        return len;
}
inline std::vector<BamTools::CigarOp> rollCig(std::string &cigar)
{
        std::vector<BamTools::CigarOp> splitCigar;
        std::istringstream parser(cigar);
        char sC; uint32_t sL;
        while(parser >> sL >> sC) { splitCigar.push_back(BamTools::CigarOp(sC, sL)); }
        return splitCigar;
}
inline float overlap(int32_t &s1, int32_t &e1, int32_t &s2, int32_t &e2)
{
        int s[2]={s1,s2}; int e[2] = {e1,e2};
        std::sort(s,s+2); std::sort(e,e+2);
        int ovr = e[0]-s[1]+1;
        float o[2] = { (float)ovr/(e2-s2+1), (float)ovr/(e1-s1+1)};
        std::sort(o,o+2);
        return o[0];
}
std::map <std::string, int> hashRef (const BamTools::RefVector &ref);
#endif
