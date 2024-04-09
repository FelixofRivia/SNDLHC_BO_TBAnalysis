#include <iostream>
#include <TFile.h>
#include <TKey.h>
#include <TList.h>

void testM(int runN, bool isTB) {
    //std::vector<const char*> fileNames;
    char* outputFileName;

    if (isTB) {
        outputFileName = Form("/eos/user/f/fmei/snd_analysis/TB/run_%06d/skimmed_%06d.root", runN, runN);
    }
    else {
        outputFileName = Form("/eos/user/f/fmei/snd_analysis/TI18/run_%06d/skimmed_%06d.root", runN, runN);
    }

    // Check if the output file exists
    if (gSystem->AccessPathName(outputFileName) == 0) {
        std::cout<<"WARNING: Output file already exists, overwriting...\n";
        // Output file exists; remove it before running hadd
        gSystem->Unlink(outputFileName);
    }

    // Construct the hadd command
    TString haddCommand = "hadd ";
    haddCommand += outputFileName;
    if (isTB) {
        haddCommand += Form(" /eos/user/f/fmei/snd_analysis/TB/run_%06d/sndsw_raw*", runN);
    }
    else {
        haddCommand += Form(" /eos/user/f/fmei/snd_analysis/TI18/run_%06d/sndsw_raw*", runN);
    }

    // Execute the hadd command
    gSystem->Exec(haddCommand);
    std::cout<<"DONE\n";
}