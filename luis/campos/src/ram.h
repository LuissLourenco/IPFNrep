#ifndef _RAM_H_
#define _RAM_H_

#include <iostream>
#include <fstream>
#include <sstream>
#include <cstring>
#include <vector>
#include "sys/types.h"
#include "sys/sysinfo.h"

using namespace std;


double situacao(){
	struct sysinfo memInfo;
	sysinfo (&memInfo);
	long long totalVirtualMem = (memInfo.totalram + memInfo.totalswap) * memInfo.mem_unit;
	long long virtualMemUsed = (memInfo.totalram - memInfo.freeram) * memInfo.mem_unit;
	virtualMemUsed += (memInfo.totalswap - memInfo.freeswap) * memInfo.mem_unit;
	long long totalPhysMem = memInfo.totalram * memInfo.mem_unit;
	long long physMemUsed = (memInfo.totalram - memInfo.freeram) * memInfo.mem_unit;
	return (double)physMemUsed/(double)totalPhysMem *100.;
}


#endif