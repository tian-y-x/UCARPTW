//
// Created by 段心童 on 2021/12/21.
//

#ifndef RESEARCH_READ_UCARP_TW_H
#define RESEARCH_READ_UCARP_TW_H

#include <string>
#include "instance.h"

using namespace std;

string& trim_tw(string& s);

Instance* data_process_macOS_tw(string dir_path);

Instance* data_process_win_tw(string dir_path);


#endif //RESEARCH_READ_UCARP_TW_H
