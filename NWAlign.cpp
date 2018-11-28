//
//  NWAlign.cpp
//  mapris
//
//  Created by mark enstrom on 4/16/18.
//  Copyright © 2018 mark enstrom. All rights reserved.
//

#include "NWAlign.hpp"


/*--------------------------------------------------------------------------------------------
 * "Dynamic Programming" sequence alignement. The first sequence is the query sequence and
 * is being matched to the longer seq. The first seq must not be longer than the second.
 * Initial unaligned spaces have no penalty but there is
 * a penalty for non-complete match
 *
 *
 * s1 is assmued to be a constant and s2 is vaiable. _alignStart is the location on s2 where
 * s1 matches. _alignLength is the length s1 aligns to s2 (insertions in s2 increase this length)
 *--------------------------------------------------------------------------------------------*/

int NWAlign::alignWithLeadingGap(string &s1, string &s2) {
    //
    // pad shorter string with "U"
    //
    if (s1.length() > s2.length()) {
        //
        // looking for 21 in s2....broken if s2 shorter
        //
        return 1000;
    }
 
    //cout << "s1 = " << s1 << endl;
    //cout << "L(s1) = " << s1.length() << endl;
    //cout << "s2 = " << s2 << endl;
    //cout << "L(s2) = " << s2.length() << endl;
    char *p1 = (char *)s1.c_str();
    int m1   = (int)s1.length();
    char *p2 = (char*)s2.c_str();
    int m2   = (int)s2.length();
    //
    // starting
    //
    //        G    C    A    T    G
    //
    //    0   0    0    0    0    0
    //
    // G -1   1    0   -1   -2   -3
    //
    // A -2   0    0    1    0   -1
    //
    // T -3  -1   -1    0    2    1
    //
    // T -4  -2   -2   -1    1    1
    //
    // A -5  -3   -3   -1    0    0
    //
    //
    // row by row build matrix
    //
    int ned[m2+1][m1+1];
    //
    // init score array to zero
    //
    for (int j = 0; j <= m2; j++) {
        for (int i = 0; i <= m1; i++) {
            ned[j][i] = 0;
        }
    }
    //
    // set first row to - count
    //
    int val = 0;
    //
    //
    //
    //    for (int j = 0; j <= m2; j++) {
    //        ned[j][0] = val--;
    //    }
    //
    // set first col to minus count
    //
    val = 0;
    for (int i = 0; i <= m1; i++) {
        ned[0][i] = val--;
    }
    //
    // build score matrix - names of row, col backwards
    //
    for (int row = 0; row < m2; row++) {
        for (int col = 0; col < m1; col++) {
            char dna1 = p1[col];
            char dna2 = p2[row];
            int nedRow = row+1;
            int nedCol = col+1;
            int scoreLeft = ned[nedRow-1][nedCol] -1;
            int scoreUp   = ned[nedRow][nedCol-1] -1;
            int scoreDiag = ned[nedRow-1][nedCol-1];
            if (dna1 == dna2) {
                //scoreDiag += 1;
            } else {
                scoreDiag -= 1;
            }
            int score = std::max({scoreLeft,scoreUp,scoreDiag});
            ned[nedRow][nedCol] = score;
        }
    }
    //
    // add unaligned end error
    //
    for (int i = 0; i <= m1; i++) {
        int lengthError = m1 - i;
        if (lengthError > 0) {
            ned[m2][i] -= lengthError;
        }
    }
    
    //
    // dbg print
    //
#if 0
    for (int row = 0; row < m2+1; row++) {
        for (int col = 0; col < m1+1; col++) {
            cout << setw(3) << ned[row][col] << " ";
        }
        cout << endl;
    }
#endif
    //    for (int row = 0; row <= m2; row++) {
    //        for (int col = 0; col <= m1; col++) {
    //            cout << setw(4) << ned[row][col];
    //        }
    //        cout << "\n";
    //    }
    //
    // search last row and column for minimum edit distance
    //
    int iDist = -1000000;
    int iMin = 0;
    
    int jDist = -1000000;
    int jMin = 0;
    
    for (int i = 0; i <= m1; i++) {
        int score = ned[m2][i];
        if (score > iDist) {
            iDist = score;
            iMin = i;
        }
    }
    
    for (int j = 0; j <= m2; j++) {
        int score = ned[j][m1];
        if (score > jDist) {
            jDist = score;
            jMin = j;
        }
    }
    //
    // walk back to find alignment
    //
    int iEnd = 0;
    int jEnd = 0;
    int retDist = 0;
    
    if (iDist >= jDist) {
        iEnd = iMin;
        jEnd = m2;
        retDist = iDist;
    } else {
        iEnd = m1;
        jEnd = jMin;
        retDist = jDist;
    }
    //cout << "Start walk back: iEnd = " << iEnd << " jEnd = " << jEnd << "retDist = " << retDist << endl;
    //cout << "iDist = " << iDist << endl;
    //cout << "jDist = " << jDist << endl;
    //
    // find best path back
    //
    int aLength = 0;
    while ((iEnd > 0) && (jEnd > 0)) {
        //
        // walk back...diag = match
        //          ...i = insert
        //          ...j - del
        int aDiag = ned[jEnd-1][iEnd-1];
        int aIns  = ned[jEnd-1][iEnd];
        int aDel  = ned[jEnd][iEnd-1];
        //cout << "comp[" << jEnd<<"]["<<iEnd<<"] = " << aDiag << "  " << aIns << "  " << aDel << endl;
        if (aDiag >= aIns) {
            if (aDiag >= aDel) {
                iEnd--;
                jEnd--;
                aLength++;
                //cout << "match " << endl;
            } else {
                iEnd--;
                //cout << "del " << endl;
                
            }
        } else {
            if (aDel >= aIns) {
                iEnd--;
                //cout << "del " << endl;
            } else {
                jEnd--;
                //cout << "ins " << endl;
                aLength++;
            }
        }
        //cout << "aLength = " << aLength << endl;
        //cout << "jEnd    = " << jEnd << endl;
    }
    
    _alignLength = aLength;
    _alignStart  = jEnd;
    
    
    //cout << "return = " << retDist << "\n";
    return -retDist;
    
}

/*--------------------------------------------------------------------------------------------
 *
 * Helper max function
 *
 *
 *--------------------------------------------------------------------------------------------*/
int max3 (int a, int b, int c) {
    if (a > b) {
        if ( a > c) {
            return a;
        } else {
            return c;
        }
    } else {
        if (b > c) {
            return b;
        } else {
            return c;
        }
    }
}

/*--------------------------------------------------------------------------------------------
 *
 * old alignment
 *
 *--------------------------------------------------------------------------------------------*/
int NWAlign::align(string s1, string s2) {
    //
    // pad shorter string with "U"
    //
    if (s1.length() > s2.length()) {
        while (s2.length() < s1.length()) {
            s2 = s2 + "U";
        }
    } else if (s2.length() > s1.length()) {
        while(s1.length() < s2.length()) {
            s1 = s1 + "U";
        }
    }
    
    char *p1 = (char *)s1.c_str();
    int m1   = (int)s1.length();
    char *p2 = (char*)s2.c_str();
    int m2   = (int)s2.length();
    
    
    
    //
    // starting
    //
    //        G    C    A    T    G
    //
    //    0  -1   -2   -3   -4   -5
    //
    // G -1   1    0   -1   -2   -3
    //
    // A -2   0    0    1    0   -1
    //
    // T -3  -1   -1    0    2    1
    //
    // T -4  -2   -2   -1    1    1
    //
    // A -5  -3   -3   -1    0    0
    //
    //
    // row by row build matrix
    //
    int ned[m2+1][m1+1];
    //
    // init score array to zero
    //
    for (int j = 0; j <= m2; j++) {
        for (int i = 0; i <= m1; i++) {
            ned[j][i] = 0;
        }
    }
    //
    // set first column to - count
    //
    int val = 0;
    for (int j = 0; j <= m2; j++) {
        ned[j][0] = val--;
    }
    //
    // set first row to minus count
    //
    val = 0;
    for (int i = 0; i <= m1; i++) {
        ned[0][i] = val--;
    }
    //
    // build score matrix
    //
    for (int row = 0; row < m2; row++) {
        for (int col = 0; col < m1; col++) {
            char dna1 = p1[col];
            char dna2 = p2[row];
            int nedRow = row+1;
            int nedCol = col+1;
            int scoreLeft = ned[nedRow-1][nedCol] -1;
            int scoreUp   = ned[nedRow][nedCol-1] -1;
            int scoreDiag = ned[nedRow-1][nedCol-1];
            if (dna1 == dna2) {
                //scoreDiag += 1;
            } else {
                scoreDiag -= 1;
            }
            int score = std::max({scoreLeft,scoreUp,scoreDiag});
            ned[nedRow][nedCol] = score;
        }
    }
    
    //    for (int row = 0; row <= m2; row++) {
    //        for (int col = 0; col <= m1; col++) {
    //            cout << setw(4) << ned[row][col];
    //        }
    //        cout << "\n";
    //    }
    //
    // search last row and column for minimum edit distance
    //
    int editDist = -1000000;
    for (int i = 0; i <= m1; i++) {
        int score = ned[m2][i];
        if (score > editDist) editDist = score;
    }
    for (int j = 0; j <= m2; j++) {
        int score = ned[j][m1];
        if (score > editDist) editDist = score;
    }
    //cout << "return = " << editDist << "\n";
    //
    // not set
    //
    _alignStart = -1;
    _alignLength = -1;
    return -editDist;
    
}
