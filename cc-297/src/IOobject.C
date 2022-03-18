#include "IOobject.H"

// constructors and destructors


// member functions
void IOobject::readSetup( )
{
    ifstream readFile;
    readFile.open( INP_FILE_PATH + "setup.inp" );
    assert( readFile.is_open() );

    string meof;
    while ( readFile >> meof )
    {
        cout << meof << "\n" ;
    }

    readFile.close();
}
