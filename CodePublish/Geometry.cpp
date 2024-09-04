
#include <fstream>
#include <cmath>
#include <cassert>
using namespace std;

#include "Geometry.h"
#include "Species.h"

// Constructor
// -----------
// in the constructor, we save the species reference, then read from
// the input file.
Geometry::Geometry(Species & S, fstream &fin, fstream &fout) :
    Facilitator(fin, fout), theSpecies(S)
{
    Initialize((CharString) "Geometry", fin);
    NumberOfGeometries = Read(fin);
//    theCoordinates.ResetToFront();
//    while (theCoordinates)
//    {
//      Coordinates thisCoord = theCoordinates.Get();
//      cout << thisCoord.SpecName << " " << thisCoord.MetalCoord << " " << thisCoord.NumAtoms << endl;
//      ++theCoordinates;
//    }
//    exit(1);
    (void)NumberOfGeometries;
}
// Read()
// ------
// We read from the input file.
// The previous call to Initialize("Geometry") put us just after the
// word GEOMETRY in the input file.
//
int Geometry::Read(fstream &fin)
{
    CharString End = "END";
    CharString Next, SubNext;
    double SurfCharge;
    double Pos[3];
    Eshort MMFF;
    float Charge;
    
    EPosition<double> NewCoordinates;
    
    Coordinates NewList;
    
    int Num;
    EPosition<double> Basis[3];
    fin >> Next;
    int i = 0;
    while(Next.lower() != End.lower())
    {
        // Find the species named
        Num = 0;
        fin >> SurfCharge;
        while (fin >> Next && (Next.lower() != "surfacecharge" && Next.lower() != End.lower()))
        {
          NewList.SurfaceCharge = SurfCharge;
          NewList.SpeciesName = theSpecies.Getid(Next);
          NewList.SpecName = theSpecies.getName(theSpecies.Getid(Next));
          // get the coordination number of the species
          fin >> NewList.MetalCoord;
          NewList.NumAtoms = 0;
          
          // The first coordinates listed is of grid site. The other two are any points in the plane of surface
          for(int i = 0; i < 3; ++i)
          {
              fin >> SubNext >> Pos[0] >> Pos[1] >> Pos[2];
              Basis[i].SetDistance(Pos[0], Pos[1], Pos[2]);
          }
          
          // Now, for each atom listed...
          while(fin >> SubNext && (SubNext.lower() != "finish"))
          {
              // Get the x, y, z, MMFF and Charge for that atom
              fin >> Pos[0] >> Pos[1] >> Pos[2] >> MMFF >> Charge;
              assert(NewList.NumAtoms < MAX_ATOMS);
              
              // Save the x, y, and z
              NewCoordinates.SetDistance(Pos[0], Pos[1], Pos[2]);
              
              
              // use 'SubNext' to determine the id.  Then set the rest
              // of the data.
              NewList.AtomNames[NewList.NumAtoms] = theSpecies.Getid(SubNext);
              NewList.AtomName[NewList.NumAtoms] = theSpecies.getName(theSpecies.Getid(SubNext));
              NewList.Position[NewList.NumAtoms] = NewCoordinates;
              NewList.MMFF94_Type[NewList.NumAtoms] = MMFF;
              NewList.Charge[NewList.NumAtoms] = Charge;
              ++NewList.NumAtoms;
          }
          theSpecies.get(NewList.SpeciesName)->NumAtoms = NewList.NumAtoms;
          for ( int i = 0; i < theSpecies.get(NewList.SpeciesName)->NumAtoms; ++i)
          {  
            theSpecies.get(NewList.SpeciesName)->AtomName[i] = NewList.AtomName[i];
          }
          // Now make them into a standard orientation
          ReOrient(Basis, NewList.Position, NewList.NumAtoms);
          theCoordinates.Add(NewList);
          ++Num;
        }  
    }
    return Num;
}



// Reorient()
// ----------
// This function reorients the coordinates
// Basically, we set the first metal atom to be the origin.  
void Geometry::ReOrient(EPosition<double> Basis[], EPosition<double> Coords[],
    Eshort Num)
{
    double basvec[3][3] = {{1.0, 0.0, 0.0}, {0.0, 1.0, 0.0}, {0.0, 0.0, 1.0}};
    double t[3][3] = {{1.0, 0.0, 0.0}, {0.0, 1.0, 0.0}, {0.0, 0.0, 1.0}};
    
    double MetalVecA[3];
    double MetalVecB[3];
    
    // Finding theta, phi, and r
    //double xo, yo, zo;
    //Basis[0].GetDistance(xo, yo, zo);
    //double x, y, z;
    
    // DW 8/26/04: This is my code to replace Eric's code (below).
    MetalVecA[0] = Basis[1].GetXDistance() - Basis[0].GetXDistance();
    MetalVecA[1] = Basis[1].GetYDistance() - Basis[0].GetYDistance();
    MetalVecA[2] = Basis[1].GetZDistance() - Basis[0].GetZDistance();
    
    MetalVecB[0] = Basis[2].GetXDistance() - Basis[0].GetXDistance();
    MetalVecB[1] = Basis[2].GetYDistance() - Basis[0].GetYDistance();
    MetalVecB[2] = Basis[2].GetZDistance() - Basis[0].GetZDistance();
    // translate the adsorbate coordinates to metal atom one
    // calculation of the rotation matrix
    rotation_matrix_calc(MetalVecA, MetalVecB, basvec, t);
    double xo, yo, zo;
    Basis[0].GetDistance(xo, yo, zo);   
    for (int i=0; i<3; i++) {
        double xm[3];
        Basis[i].GetDistance(xm[0], xm[1], xm[2]);
        xm[0] = xm[0] - xo;
        xm[1] = xm[1] - yo;
        xm[2] = xm[2] - zo;
        for(int j = 0; j < 3; ++j)
        {
            xm[j] = t[0][j] * xm[0] +
                t[1][j] * xm[1] +
                t[2][j] * xm[2];
            
        }
        Basis[i].SetDistance(xm[0], xm[1], xm[2]);
    } 
    for(int i = 0; i < Num; ++i)
    {
        // DW 8/26/04: I could have left this outside of the loop, but
        // since it is only run a few times, it isn't crucial that it
        // be super fast.  It is more clear having it in here.
        // double xo, yo, zo;
        // Basis[0].GetDistance(xo, yo, zo);
        
        
        double x, y, z;
        Coords[i].GetDistance(x, y, z);
        
        // DW 8/26/04: This is not needed because we set the distance
        //             later on.  
        //Coords[i].SetDistance(x - xo, y - yo, z - zo); // is this needed?
        
        double r[3];
        double rp[3];

        
        r[0] = x - xo;
        r[1] = y - yo;
        r[2] = z - zo;
        
        // Now do the transformation to the new basis set.
        for(int j = 0; j < 3; ++j)
        {
            rp[j] =
                t[0][j] * r[0] +
                t[1][j] * r[1] +
                t[2][j] * r[2];
            
        }
        Coords[i].SetDistance(rp[0], rp[1], rp[2]);
    }
}


// rotation_matrix_calc()
// ----------------------
// Calculate the Rotation Matrix to transform from our coordinates to
// a normal set of coords.
//
// 'va' and 'vb' are the basis vectors that describe the plane that
// the metal atoms are in.
//
// 'basvec' is what we multiply the resulting matrix by.  Since we
// only pass in the unit matrix, this has no effect.
//
// 't' is the output matrix
//
void Geometry::rotation_matrix_calc(double va[], double vb[],
    double basvec[3][3], double t[][3])
{
    double vap[3], vbp[3], norm[3];
    double u[3][3];
    
    // First, make va and vb unit vectors
    double lena = sqrt(va[0] * va[0] + va[1] * va[1] + va[2] * va[2]);
    double lenb = sqrt(vb[0] * vb[0] + vb[1] * vb[1] + vb[2] * vb[2]);
    for(int i = 0; i < 3; ++i)
    {
        va[i] = va[i] / lena;
        vb[i] = vb[i] / lenb;
    }
    
    // <vap> = <va> [[basvec]]
    // where <> denotes a vector and [[ ]] denotes a matrix.
    // basvec[j] is a column vector.
    for(int j = 0; j < 3; ++j)
    {
        vap[j] =
            va[0] * basvec[j][0] +
            va[1] * basvec[j][1] +
            va[2] * basvec[j][2];
        
        vbp[j] =
            vb[0] * basvec[j][0] +
            vb[1] * basvec[j][1] +
            vb[2] * basvec[j][2];
    }
    
    // Find cross product of input metal vectors
    // This makes norm[] perpendicular to the metal surface.
    norm[0] = vap[1] * vbp[2] - vap[2] * vbp[1];
    norm[1] = vap[2] * vbp[0] - vap[0] * vbp[2];
    norm[2] = vap[0] * vbp[1] - vap[1] * vbp[0];
    
    double len =
        sqrt(norm[0] * norm[0] + norm[1] * norm[1] + norm[2] * norm[2]);
    
    norm[0] = norm[0] / len;
    norm[1] = norm[1] / len;
    norm[2] = norm[2] / len;
    
    // this is the first basis ( z direction )
    u[2][0] = norm[0];
    u[2][1] = norm[1];
    u[2][2] = norm[2];
    
    // the second basis, we will take as the vector
    // to the second metal atom
    // unit length was already calculated
    // Also, we know that va[] is perpendicular to norm[]
    u[1][0] = va[0];
    u[1][1] = va[1];
    u[1][2] = va[2];
    
    // the third basis is the cross product of the first two
    // This ensures that all three vectors are orthogonal.
    // Also, since norm[] and va[] are normal and perpendicular, this
    // guarantees that this third vector is normal as well.  Thus,
    // these are an orthonormal basis.
    u[0][0] = norm[1] * va[2] - va[1] * norm[2];
    u[0][1] = norm[2] * va[0] - va[2] * norm[0];
    u[0][2] = norm[0] * va[1] - va[0] * norm[1];
    
    
    // ensure unit length
    len = sqrt(u[0][0] * u[0][0] + u[0][1] * u[0][1] + u[0][2] * u[0][2]);
    for(int k = 0; k < 3; ++k)
    {
        u[0][k] = u[0][k] / len;
    }
    
    // since the basis is orthogonal, the inverse is the
    // transpose of the matrix
    for(int m = 0; m < 3; ++m)
    {
        for(int n = 0; n < 3; ++n)
        {
            t[m][n] = u[n][m];
        }
    }
    
    
    // Now we multiply [[basvec]] and [[t]] and store the result in [[temp]].
    double temp[3][3];
    for(int p = 0; p < 3; ++p)
    {
        for(int q = 0; q < 3; ++q)
        {
            temp[p][q] = 0.;
            for(int r = 0; r < 3; ++r)
            {
                temp[p][q] += basvec[r][q] * t[p][r];
            }
        }
    }
    
    // Finally, copy the result to [[ t ]]
    for(int ii = 0; ii < 3; ++ii)
    {
        for(int jj = 0; jj < 3; ++jj)
        {
            t[ii][jj] = temp[ii][jj];
        }
    }


    return;
}
