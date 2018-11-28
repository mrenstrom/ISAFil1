//
//  readFasta.cpp
//  mapris
//
//  Created by mark enstrom on 3/26/18.
//  Copyright © 2018 mark enstrom. All rights reserved.
//

#include "readFasta.hpp"
#include <cstdio>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <string>
#include <array>
#include <iomanip>
#include <numeric>
#include <ctime>
#include "NWAlign.hpp"
using namespace std;
//
// vector sort classes
//
struct myclassS {
    bool operator() (Capture *pa1,Capture *pa2) { return (pa1->_count > pa2->_count);}
} smallCompObj;


struct myclassL {
    bool operator() (Capture *pa1,Capture *pa2) { return (pa1->_count < pa2->_count);}
} largeCompObj;


struct myclassLen {
    bool operator() (Capture *pa1,Capture *pa2) { return (pa1->_seq.length() > pa2->_seq.length());}
} smallLenCompObj;


struct myclassString {
    bool operator() (Capture *pa1,Capture *pa2) { return (pa1->_seq > pa2->_seq);}
} smallStrCompObj;

//  readFile.cpp
//  Created by mark enstrom on 12/3/17.
//  Copyright © 2017 mark enstrom. All rights reserved.
//
/*--------------------------------------------------------------------------------------------
 *
 * like strncmp - exact string comp
 * compare a dna(char) sequence and return # differences
 *
 *
 *--------------------------------------------------------------------------------------------*/
int dnancmp(char *p1, char *p2,size_t n) {
    int error = 0;
    for (int i = 0; i < n; i++) {
        if (p1[i] != p2[i]) {
            error += 1;
            if (error > 4) return error;
        }
    }
    return error;
}
/*--------------------------------------------------------------------------------------------
 *
 * addSequence 
 *
 *	check for ggod primer and ltr ...and no vector, then add to map
 *
 *
 *--------------------------------------------------------------------------------------------*/
bool ReadFasta::addSequence(string s) {
    bool bRet = findAndTrimPrimerSeq(s);
    if (bRet) {
        ++_risMap[s];
    }
    return bRet;
}
/*--------------------------------------------------------------------------------------------
 *
 * Read for fasta
 *
 *
 *>BMTW6:02212:01662
 *AGTAGTATGTGGCCCGCTCTGTTGTGTGACTCTGGTAACTAGAGATCCCTCAGACCCTTT
 *TAGTCAGTGCAAGGCAGAGATATACTGGCCTTCCTGCTAGGGGAGAGTGTGAGCCACACT
 *CACGCTCTCTCGTTCTCCACTGGAGGGGACACCCGCACACAAGGAGCTTGCAGCCCAGGA
 *GGACCAGGAAACCGG
 *>BMTW6:03334:00712
 *CATAGTGGTGA
 *
 *
 *--------------------------------------------------------------------------------------------*/
void
ReadFasta::readFastA(string filename)
{
    std::ifstream inFile;
    inFile.open(filename);
    char hdr[512];
    char seq[512];
    bool bSequenceComplete = false;
    //
    // read first header line
    //
    inFile.getline(hdr, 512);
    //
    // read rest of file
    //
    string s = "";
    int iCount  = 0;
    int badPrimer  = 0;
    do {
        inFile.getline(seq, 512);
        if (!inFile) {
            break;
        }
        if (seq[0] == '>') {
            bSequenceComplete = true;
        } else {
            s += string(seq);
        }
        
        if (bSequenceComplete) {
            _totalSequence++;
            //
            // try to add
            //
            if (!addSequence(s)) badPrimer++;
            ++iCount;
            s = "";
            bSequenceComplete = false;
            
            if ((iCount % 10000) == 0) {
                cout << "read " << iCount << " seq: good = " << _goodSeqShortLtr + _goodSeqFullLtr << endl;
                
            }
            //if (iCount > 1000) break;
        }
    } while (true);
    //
    // finish
    //
    if (s != "") {
        if (!addSequence(s)) badPrimer++;
    }
    cout << "read " << iCount << " Sequence strings \n";
    cout << "risMap size = " <<  _risMap.size() << endl;
}
    
    

/*--------------------------------------------------------------------------------------------
 *
 * Read for fastq
 *
 *
 *
 * @M03100:56:000000000-AH66G:1:1101:15597:1533 1:N:0:1
 * AGTAGTGTGTGCCCGTCTGTCCCGTGGGGCCTAACTGCTGTGCCACTGAATTCAGATC
 * +
 * 8ICIIFIFIEIFFGIAFFIGE@@BCB,777IIIIIIIIIIIFIIIIIIIIGIIIIAGG
 *
 *
 *--------------------------------------------------------------------------------------------*/
void
ReadFasta::readFastQ(string filename)
{
    
    std::ifstream inFile;
    inFile.open(filename);
    if (!inFile) {
        cout << "can't open " << filename << "\n";
        return;
    }
    char hdr[4096];
    char seq[4096];
    char spacer[10];
    char score[4096];
    int iCount = 0;
    do {
        inFile.getline(hdr, 4096);
        if (!inFile) {
            break;
        }
        
        inFile.getline(seq, 4096);
        if (!inFile) {
            break;
        }
        
        inFile.getline(spacer, 10);
        if (!inFile) {
            break;
        }
        
        inFile.getline(score, 4096);
        if (!inFile) {
            break;
        }
        string s_seq(seq);
        //
        // quality check
        //
        _totalSequence++;
        int badBases = 0;
        for (int i = 0; i < (strlen(score)-1); i++)
        {
            int qScore = (int)score[i] - 33;
            if (qScore < 20)
            {
                //seq[i] = 'N';
		// don't change to N, clustal wont work
                badBases++;
            }
        }
        if (badBases < 4) {
            addSequence(s_seq);
        } else {
            _failQuality++;
        }
        ++iCount;
        if ((iCount % 10000) == 0) {
            cout << "read " << iCount << " seq: good = " << _goodSeqShortLtr + _goodSeqFullLtr << endl;
        }
        
    } while (true);
    
    cout << "read " << iCount << " Sequence strings \n";
    cout << "risMap size = " <<  _risMap.size() << endl;
}



/*--------------------------------------------------------------------------------------------
 * ReadFasta::findAndTrimPrimerSeq
 *
 *
 *  look for defined primer and ltr sequence while allowing errors. Trim if found.
 *
 *
 *--------------------------------------------------------------------------------------------*/
bool ReadFasta::findAndTrimPrimerSeq(string& seq)
{
    NWAlign nwAlign =  NWAlign();
    bool _matchFullLTR = false;
    bool _matchShortLTR = false;
    //
    // calc edit (lev) distance threholds for primer/ltr and vector
    // more lenient on vector
    //
    int ltrError   = (int)_pr_ltr.length()/8;
    int ltrError_s = (int)_pr_ltr_s.length()/8;
    int vecError   = (int)_pr_vec.length()/5;
    if (seq.length() >= (_pr_ltr_s.length() + 30)) {
	//
        // to keep string lengths as short as possible: 
        //
	// primer start must be within ?? of seq start ...
        // we know some LentiLong and MiseqLentiLong overlap so must cover dist (30) between these
	//
	string sTmp = seq.substr(0,min(seq.length(),_pr_ltr.length()+30));
        //
        // look for close match of primer + start of LTR
        //
        int editDist = nwAlign.alignWithLeadingGap(_pr_ltr,sTmp);
        //
        // if full ltr not found try smaller version
        //
        if ((editDist > ltrError)) {
            editDist = nwAlign.alignWithLeadingGap(_pr_ltr_s,sTmp);
            if ((editDist <= ltrError_s)) {
                _matchShortLTR = true;
            }
        } else {
            _matchFullLTR = true;
        }
        //
        // if primer and ltr found then process
        //
        if (_matchFullLTR || _matchShortLTR) {
            _prLtrFound++;
            //
            // trim primer and start of LTR off of sequence
            //
            //cout << "before : " << seq << " " << nwAlign._alignStart << " " << nwAlign._alignLength << endl;
            seq = seq.substr(nwAlign._alignStart+nwAlign._alignLength);
            //cout << "after :  " << seq << endl;
            //
            // look for vector
            //
            NWAlign nwVec = NWAlign();
            int vecDist = nwVec.alignWithLeadingGap(_pr_vec, seq);
            if (vecDist <= vecError) {
                //
                // vector found
                //
                _vectorError++;
                return false;
            }
            //
            // 2nd vec sequence
            //
            vecDist = nwVec.alignWithLeadingGap(_pr_vec2, seq);
            if (vecDist <= vecError) {
                //
                // vector found
                //
                _vectorError++;
                return false;
            }
            //
            // look for plasmid if defined (3' of final LTR)
            //
            if (_pr_plasmid != "") {
                int plasmidDist = nwVec.alignWithLeadingGap(_pr_plasmid,seq);
                if (plasmidDist <= 6) {
                    // plasmid found
                    cout << "plasmid found " << plasmidDist << endl;
                    _plasmidError++;
                    return false;
                }
            }
            //
            // look for linker
            //
            NWAlign nwLink = NWAlign();
            int linkDist = nwLink.alignWithLeadingGap(_LCTLS1,seq);
            if (linkDist <= 7) {
                //
                // trim linker off
                //
                seq = seq.substr(0,nwLink._alignStart);
            }
            //
            // add to vector of sequences
            //
            if (seq.length() >= 30) {
                if (_matchFullLTR) {
                    _goodSeqFullLtr++;
                } else {
                    _goodSeqShortLtr++;
                }
                return true;
            } else {
                _lengthError2++;
                return false;
            }
        } else {
            //cout << "0         1         2         3         4         5         6         7         8         9         0         1         2         3         4         5         6         7         8         9\n";

            //cout << seq << endl;
            //cout << _pr_ltr << endl;
            //cout << "pt_ltr: " << editDist << "\t" << nwAlign._alignStart << "\t" << nwAlign._alignLength << endl;
            _prLtrError++;
            return false;
        }
    } else {
//        cout << "bad init length\n";
//        cout << seq << endl;
//        cout << pr_ltr << endl;
        _lengthError1++;
        return false;
    }
}

/*--------------------------------------------------------------------------------------------
 *
 * Read all source files and merge
 *
 *
 *
 *
 *--------------------------------------------------------------------------------------------*/

bool ReadFasta::read(string source)
{
    _sourceFiles.push_back(source);
    cout << "Read "  << source << "\n";
    //
    // read fastA or fastQ
    //
    if (source.find("fastq") != string::npos) {
        readFastQ(source);
    } else {
        readFastA(source);
    }
    if (_risMap.size() > 0) return true;
    return false;
}

/*--------------------------------------------------------------------------------------------
 *
 * errorCorrection
 *
 *
 *
 *
 *--------------------------------------------------------------------------------------------*/
void ReadFasta::errorCorrection() {
    //
    // copy map to vector
    //
    int falsePositive = 0;
    int id = 0;
    for (auto pr:_risMap) {
        PCapture pcap =  new Capture(pr.first,to_string(id++),(int)pr.second);
        _capture.push_back(pcap);
    }
    //
    // sort on size
    //
    std::sort(_capture.begin(),_capture.end(),smallCompObj);
    int comps=0;
    NWAlign nw = NWAlign();
    //
    // search largest first, merge close
    //
    vector<int> prelimMatch = vector<int>();
    vector<int> finalMatch = vector<int>();
    
    clock_t begin = clock();
    int matchPrelim = 0;
    int matchComplete = 0;
    for (int i = 0; i < _capture.size(); i++) {
        PCapture pi = _capture[i];
        if ((pi->_valid) && (pi->_count > 1)) {
            comps = 0;
            matchPrelim = 0;
            matchComplete = 0;
            //
            // merge from smallest to largest
            //
            for (int j = (int)_capture.size()-1;j > i ;j--) {
                PCapture pj = _capture[j];
                if (pj->_valid) {
                    comps++;
                    //
                    //  how about match first 8 within 1???
		            //
                    bool bSim = pi->similar(pj->_seq);
                    if (!bSim) {
                        continue;
                    }
                    matchPrelim++;
                    
                    //int matchCount = pi->similarCap(pj);
                    //if (matchCount < 16) continue;
                    //prelimMatch.push_back(matchCount);
                    
                    int d = (int)min(pi->_seq.length(),pj->_seq.length());
                    if (d > 60) d = 60;
                    
                    string s1 = pi->_seq.substr(0,d);
                    string s2 = pj->_seq.substr(0,d);
                    
                    int editDist = nw.align(s1, s2);
//                    int lDist = levenshtein_distance(s1, s2);
//                    cout << "\n  Similar =  " << bSim << endl;
//                    cout << " nwalign = " << editDist << endl;
//                    cout << " lev     = " << lDist << endl;
//                    cout << pi->_seq << endl;
//                    cout << pj->_seq << endl;
//                    cout << editDist << endl;
                    
                    
                    int max_error = d/15;
                    if (editDist <= max_error) {
                        //
                        // could do a larger check here
                        //
                        pi->merge(pj);
                        pj->_valid = false;
                        matchComplete++;
                    } else {
                        falsePositive++;
                    }
                }
            }
        }
        if ((i % 100) == 0) {
            
            clock_t end = clock();
            double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
            begin = end;
            
            cout << "combine seq count = " << i << " j still valid = " << comps << " fp " << falsePositive << " time " << elapsed_secs << endl;
            cout << comps << "  " << matchPrelim << " " << matchComplete << endl;
        }
    }
    for (auto p:_capture) {
        if (p->_valid) {
            _short.push_back(p);
        }
    }
    cout << "starting seq " << _capture.size() << endl;
    cout << "first merge  " << _short.size() << endl;
    
    std::sort(_short.begin(),_short.end(),smallCompObj);
    //
    // log output
    //
    if (true) {
        int n1 = 0;
        int n2 = 0;
        int n10 = 0;
        int n100 = 0;
        int n1000 = 0;
        int nBig = 0;
        for (auto p: _short) {
            if (p->_count <= 1) {
                n1++;
            } else if (p->_count <= 2) {
                n2++;
            } else if (p->_count <= 10) {
                n10++;
            } else if (p->_count <= 100) {
                n100++;
            } else if (p->_count <= 1000) {
                n1000++;
            } else nBig++;
        }
        cout << " big = " << nBig << endl;
        cout << "1000 = " << n1000 << endl;
        cout << " 100 = " << n100 << endl;
        cout << "  10 = " << n10 << endl;
        cout << "   2 = " << n2 << endl;
        cout << "   1 = " << n1 << endl;
    }
    
    for (auto p : _short) {
        if (p->_valid)  {
            _noVec.push_back(p);
        }
    }
    
    cout << "final unique seq size = " << _noVec.size() << endl;
    //
    // sort by count
    //
    std::sort(_noVec.begin(),_noVec.end(),smallCompObj);
    //
    // calc span
    //
    for (auto p:_noVec) {
        std::map<int,int> spanMap = std::map<int,int>();
        for (auto pr: p->_famSeq) {
            string s = pr.first;
            ++spanMap[(int)s.length()];
        }
        p->_span = (int)spanMap.size();
    }
}


/*--------------------------------------------------------------------------------------------
 *
 * find good fasta sequences:
 *      good primer and ltr
 *      no vector
 *
 *
 *
 *--------------------------------------------------------------------------------------------*/
bool ReadFasta::process()
{
    //
    // copy map to vector
    //
    int id = 0;
    for (auto pr:_risMap) {
        PCapture pcap =  new Capture(pr.first,to_string(id++),(int)pr.second);
        _capture.push_back(pcap);
    }
    //
    // sort on size
    //
    std::sort(_capture.begin(),_capture.end(),smallCompObj);//
    //
    // >mapris:count
    // refseq
    //
    std::ofstream myFasta;
    std::string sOut = _dstDir + _test + ".fasta";
    cout << "Writing output to " << sOut << "\n";
    //
    // 
    //
    myFasta.open (sOut);
    for (auto *pCap:_capture) {
        myFasta << ">MAPRIS:" << _fileID << ':' << setfill('0') << setw(5) << pCap->_id << ":" << pCap->_count << endl;
        myFasta << pCap->_seq << endl;
    }
    myFasta.close();
    


    ofstream myStats;
    string sStats =  _dstDir + _test +  "_stat.txt";
    myStats.open(sStats);
    for (auto s:_sourceFiles) {
        myStats << "Source:" << s << "\n";
    }
    myStats << "Total Sequences           = " << setw(8) << _totalSequence << endl;
    myStats << "Fail Quality              = " << setw(8) << _failQuality << endl;
    myStats << "Fail Lenth 1              = " << setw(8) << _lengthError1 << endl;
    myStats << "Primer & LTR not found    = " << setw(8) << _prLtrError << endl;
    myStats << "Primer and LTR found      = " << setw(8) << _prLtrFound << endl;
    myStats << "Contains vector           = " << setw(8) << _vectorError << endl;
    myStats << "Contains plasmid          = " << setw(8) << _plasmidError << endl;
    myStats << "Fail Lenth 2              = " << setw(8) << _lengthError2 << endl;
    myStats << "good full primer and ltr  = " << setw(8) << _goodSeqFullLtr << endl;
    myStats << "good short primer and ltr = " << setw(8) << _goodSeqShortLtr << endl;
    myStats << "Total good seq            = " << setw(8) << _goodSeqShortLtr + _goodSeqFullLtr << endl;
    
    myStats.close();
    
    for (auto pc: _capture) {
        delete pc;
    }
    _capture.clear();
    return true;
}
