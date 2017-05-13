#include "aux.h"
using namespace BamTools;
using namespace std;
map <string, int> hashRef (const RefVector &ref)
{
        map<string, int> hash;
        for(uint16_t i=0; i!= ref.size(); ++i){ hash[ref[i].RefName]=i;}
        return hash;
}
