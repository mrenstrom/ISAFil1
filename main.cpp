//
//  main.cpp
//  mapris
//
//  Created by mark enstrom on 3/26/18.
//  Copyright © 2018 mark enstrom. All rights reserved.
//

#include "readFasta.hpp"
#include "NWAlign.hpp"
#include <cstdio>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <string>
#include <array>
#include <iomanip>
#include <unistd.h>
using namespace std;

/*--------------------------------------------------------------------------------------------
 * main routine for filtering fasta for ISA processing
 *
 * baseDir:     base of ISA tree:
 * subject:     ie Z09132
 * test:        ie 32DPT_CD3
 * vectorName:  identifier of vector for Primer/LTR/Vector/Linker sequenecs
 * sourceFiles: IT14_KiemLab47KGH.fasta IT15_KiemLab47KGH.fasta IT16_KiemLab47KGH.fasta
 * --- based on parameters program assumes:
 *      source:       base/Z09132/fasta_files/<filename>
 *      destination:  base/Z09132/filtered_files/subject_test_rep.fasta
 *--------------------------------------------------------------------------------------------*/
int main(int argc, const char * argv[]) {
    //string master =  "/fh/fast/kiem_h/grp/kiemlab/PROJECTS/RIS/Master Sequence Files/";
    string master = "/Users/mark_enstrom/Barcode/RIS/";
    //
    // initial read of source
    //
    
    std::cout << "Read FASTA\n";
    
    if (argc < 6) {
        cout << "Usage: ISAfilter  destDir fileID test vectorName sourceFile1 sourceFile2...\n";
        exit(0);
    }
    string dstDir = argv[1];
    string fileID = argv[2];
    string test = argv[3];
    string vectorName = argv[4];
    
    ReadFasta rf = ReadFasta(dstDir,fileID,test,vectorName);
    //
    // read in sequences form each file and do quality filter
    //
    int current = 5;
    while (current < argc) {
        string source = argv[current++];
        
        // server
        bool bRet = rf.read(master + source);
        // local
        //bool bRet = rf.read(source);
        
        
        
        if (!bRet) {
            cout << "error reading file " << source << endl;
            exit(0);
        }
    }
    //
    // combine do error correction on sequences
    //
    rf.process();
    
    return 0;
}
