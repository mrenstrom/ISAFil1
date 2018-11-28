//
//  readFasta.hpp
//  mapris
//
//  Created by mark enstrom on 3/26/18.
//  Copyright © 2018 mark enstrom. All rights reserved.
//

#ifndef readFasta_hpp
#define readFasta_hpp

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
using namespace std;





#define KMER 8

class Capture {
public:
    string  _seq;
    string  _id;
    size_t  _count;
    unordered_map<std::string,int> _famSeq;
    unordered_map<std::string,size_t> _seqMap;
    bool    _valid;
    int     _span;
    Capture(string seq,string id,int i) {
        _seq = seq;
        _id = id;
        _count = i;
        _valid=true;
        _famSeq = unordered_map<std::string,int>();
        _seqMap = unordered_map<std::string,size_t>();
        //
        // !!! duplicate just to keep counts accurate
        // !!! can delete later
        //
        _famSeq[seq] = i;
        //
        // build seq map (8) out to max of 40
        //
        size_t l = seq.length() - KMER;
        if (l > (40-KMER)) l = (40-KMER);
        for (int i = 0; i < (l-KMER);i++) {
            string s = seq.substr(i,KMER);
            ++_seqMap[s];
        }
    };
    void addSeq(string seq) {
        ++_famSeq[seq];
        _count++;
    }
    void merge(Capture *p) {
        for (auto pr:p->_famSeq) {
            string s = pr.first;
            int count = pr.second;
            _famSeq[s] += count;
        }
    }
    
    bool similar(string s2) {
        for (int i = 0; i < (s2.length()-KMER); i += KMER) {
            if  (_seqMap.find(s2.substr(i,KMER)) != _seqMap.end()) {
                return true;
            }
        }
        return false;
    }
    
    void dump(string s) {
        for (int i = 0; i < (s.length()-8); i += 8) {
            cout << "comp " << s.substr(i,8) << endl;
        }

        for (auto pr:_seqMap) {
            cout << pr.first << "\t" << pr.second << endl;
        }
    }
};


typedef Capture *PCapture;
typedef vector<PCapture> CapVec;
typedef CapVec *PCapVec;


typedef std::unordered_map<std::string,size_t> RisMap;


class ReadFasta {
public:
    vector<string> _sourceFiles;
    string _dstDir;
    string _test;
    string _fileID;
    //
    // read stats
    //
    int _totalSequence   = 0;
    int _failQuality     = 0;
    int _lengthError1    = 0;
    int _prLtrFound      = 0;
    int _goodSeqFullLtr  = 0;
    int _goodSeqShortLtr = 0;
    int _prLtrError      = 0;
    int _lengthError2    = 0;
    int _vectorError     = 0;
    int _plasmidError    = 0;
    
    vector<Capture*> _capture = vector<Capture*>();
    vector<Capture*> _short   = vector<Capture*>();
    vector<Capture*> _primer  = vector<Capture*>();
    vector<Capture*> _noPrimer = vector<Capture*>();
    vector<Capture*> _Ltr = vector<Capture*>();
    vector<Capture*> _noVec = vector<Capture*>();
    vector<Capture*> _hasVec = vector<Capture*>();
    vector<Capture*> _fullSeq = vector<Capture*>();
    RisMap _risMap;

    ReadFasta(string dstDir,string id,string test,string vectorName) {
        _dstDir = dstDir;
        _test = test;
        _risMap = RisMap();
        _fileID = id;
        _sourceFiles = vector<string>();
        
        if (vectorName == "HIV") {
            _primerSeq  = "TGTGACTCTGGTAACTAGAGATCCCTC";
            _pr_ltr     = "TGTGACTCTGGTAACTAGAGATCCCTCAGACCCTTTTAGTCAGTGTGGAAAATCTCTAGC";
            _pr_ltr_s   = "TGTGACTCTGGTAACTAGAGATCCCTCAGACCCTTTTAGTCAGTGTGGAAAA";
            _pr_vec     = "GAAAGTAAAGCCAGAGGAGATCTC";
            _pr_vec2    = "GAAAGTAAAGCCAGAGGAGATCTC";
            cout << vectorName << endl;
        } else if (vectorName == "LentiLong") {
                         //0         0         0         0         0         0         0         0         0         0         0
            _pr_ltr     = "AGTAGTGTGTGCCCGTCTGTTGTGTGACTCTGGTAACTAGAGATCCCTCAGACCCTTTTAGTCAGTGTGGAAAATCTCTAGCA";
            _pr_ltr_s   = "AGTAGTGTGTGCCCGTCTGTTGTGTGACTCTGGTAACTAGAGATCCCTCAGACCCTTTTAGTCAGTGTGGAAAA";
            _pr_vec     = "TGGCGCCCGAACAGGGACTTGAAAGCGA";
            _pr_vec2    = "AACAGGGACTTGAAAGCGAAAGGGAAAC";
            cout << vectorName << endl;
        } else if (vectorName == "MiseqLentiLong") {
                         //0         0         0         0         0         0         0         0         0         0         0
            _pr_ltr     = "AGCTTGCCTTGAGTGCTTCAAGTAGTGTGTGCCCGTCTGTTGTGTGACTCTGGTAACTAGAGATCCCTCAGACCCTTTTAGTCAGTGTGGAAAATCTCTAGCA";
            _pr_ltr_s   = "AGCTTGCCTTGAGTGCTTCAAGTAGTGTGTGCCCGTCTGTTGTGTGACTCTGGTAACTAGAGATCCCTCAGACCCTTTTAGTCAGTGTGGAAAA";
            _pr_vec     = "TGGCGCCCGAACAGGGACTTGAAAGCGA";
            _pr_vec2    = "AACAGGGACTTGAAAGCGAAAGGGAAAC";
            cout << vectorName << endl;
        } else if (vectorName == "MiseqLentiLong115") {
                         //0         0         0         0         0         0         0         0         0         0         0
            _pr_ltr     = "AGCTTGCCTTGAGTGCTTCAAGTAGTGTGTGCCCGTCTGTTGTGTGACTCTGGTAACTAGAGATCCCTCAGACCCTTTTAGTCAGTGTGGAAAATCTCTAGCAGT";
            _pr_ltr_s   = "AGCTTGCCTTGAGTGCTTCAAGTAGTGTGTGCCCGTCTGTTGTGTGACTCTGGTAACTAGAGATCCCTCAGACCCTTTTAGTCAGTGTGGAAAA";
            _pr_vec     = "TCATGTCATCTTATTATTCAGTATTTATAACTTGCAAA";
            _pr_vec2    = "AACAGGGACTTGAAAGCGAAAGGGAAAC";
            _pr_plasmid = "TGGCGCCCGAACAGGGACCTGAAAGCGAAAGGGAAACC"; 
            cout << vectorName << endl;
        } else {
            cout << "Bad vector name" << endl;
        }
    }
    void readFastA(string filename);
    void readFastQ(string filename);
    bool findSequence(string seq, string tag, string target);
    bool read(string filename);
    bool process();
    bool readQ(string filename,string outName);
    void findAndTrimPrimer();
    bool findAndTrimPrimerSeq(string& seq);
    void trimLTR();
    void lookForVector();
    bool addSequence(string s);
    void errorCorrection();
    //
    // vector sequences must all be defined
    //
    string _primerSeq  = "";
    string _pr_ltr     = "";
    string _pr_ltr_s   = "";
    string _pr_vec     = "";
    string _pr_vec2    = "";
    string _pr_plasmid = ""; 
    //
    // 454-Modified Linker Cassette Oligos:
    // 5' - 3'
    //
    string _LCTLS1 = "CCTAACTGCTGTGCCACTGAATTCAGATC";
    string _LCTLS2 =     "ACTGCTGTGCCACTGAATTCAGATC";
    string _LCTLS3 =         "CTGTGCCACTGAATTCAGATC";
    string _LCTLS4 =             "GCCACTGAATTCAGATC";
    string _LCTU =  "GACCCGGGAGATCTGAATTCAGTGGCACAGCAGTTAGG";
    
};

#endif /* readFasta_hpp */
