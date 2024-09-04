#include <stdio.h>
#include <fstream>
#include <iostream>
#include <cassert>
#include <string>
#include <string.h>
#include <sstream>
using namespace std;

int main(int argc, char *argv[])
{ 	
    fstream fin;
    fstream fout;
    char Output_File[100];
    strcpy(Output_File, "Coordinates_short0.xyz");
    fin.open(argv[1], ios::in /*|ios::nocreate*/);
    fout.open(Output_File, /*fstream::in |*/ fstream::out);
    string start_frame = "100";
    string end_frame = "250";
    string line1;
    string line2;
    int start = 0;
    //Going to the start frame
    while (getline(fin, line1))
    {
      getline(fin, line2);
      if (line2 == start_frame)
      {
        start = 1;
        
      }
      if (line2 == end_frame)
      {
        start = 0;
	break;
      }
      if (start == 1)
      {
        fout << line1 << endl;
        fout << line2 << endl;
      }
      stringstream ss;
      ss << line1;
      int num;
      ss >> num;
      for (int i = 0; i < num; i++)
      {
        string line3;
        getline(fin, line3);
        if (start == 1)
          fout << line3 << endl;
      }
    }
    fin.close();
    
//    fstream fout2;
//    char Output_File2[100];
//    strcpy(Output_File2, "Coordinates_short_Han");
//    fin.open(argv[1], ios::in /*|ios::nocreate*/);
//    fout2.open(Output_File2, /*fstream::in |*/ fstream::out);
//    start = 0;
//    //Going to the start frame
//    while (getline(fin, line1))
//    {
//      getline(fin, line2);
//      if (line2 == start_frame)
//      {
//        start = 1;
//      }
//      if (line2 == end_frame)
//      {
//        start = 0;
//        break;
//      }
//      stringstream ss;
//      ss << line1;
//      int num;
//      ss >> num;
//      for (int i = 0; i < num; i++)
//      {
//        string line3;
//        fin >> line3;
//        if (start == 1)
//        {
//            fout2 << line3 << endl;
//            fout2 << num << endl;
//        }
//        for (int j = 0; j< 3; j++)
//        {
//          fin >> line3;
//          if (start == 1)
//            fout2 << line3 << endl;
//        }
//        fout2 << endl;
//        getline(fin, line3);
//        
//      }
//      fout2 << "END" << endl;
//    }
    return 0;
}
