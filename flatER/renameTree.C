#include "TFile.h"
#include "TTree.h"
#include "TObject.h"
#include "TSystem.h"
#include <cstdio>


// root -l -b -q 'renameTree.C("PATH to file","OLD tree","NEW tree")'

void renameTree(const char* fname,
                const char* oldName,
                const char* newName)
{
    TFile f(fname, "UPDATE");
    if (!f.IsOpen()) {
        printf("Cannot open file: %s\n", fname);
        return;
    }

    TTree* t = (TTree*) f.Get(oldName);
    if (!t) {
        printf("Tree '%s' not found in %s\n", oldName, fname);
        f.Close();
        return;
    }

    t->SetName(newName);
    t->Write("", TObject::kOverwrite);

    printf("Renamed tree '%s' -> '%s' in %s\n", oldName, newName, fname);

    f.Close();
}
