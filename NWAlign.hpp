//
//  NWAlign.hpp
//  mapris
//
//  Created by mark enstrom on 4/16/18.
//  Copyright © 2018 mark enstrom. All rights reserved.
//

#ifndef nwAlign_hpp
#define nwAlign_hpp


#include <stdio.h>
#include <iostream>
#include <fstream>
#include <string.h>
#include <vector>
#include <map>
#include <unordered_map>
#include <mutex>
#include <condition_variable>
#include <queue>
#include <list>
#include <stdlib.h>
#include <assert.h>
#include <ctime>
#include <unistd.h>
#include <algorithm>
#include <iomanip>
using namespace std;


class NWAlign {
public:
    NWAlign() {};
    int _alignStart{0};
    int _alignLength{0};
    int alignWithLeadingGap(string &s1, string &s2);
    int align(string s1, string s2);
};

#endif /* NWAlign_hpp */
