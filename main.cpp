#include <iostream>
#include <vector>
#include <fstream>  // in- and output of files
#include <string>
#include <sstream>  // stringstream
#include <memory>   // smart pointers

#include <Readfile.h>
#include <Einsteinvisco.h>
#include <GKvisco.h>

int main(int argc, char *argv[]) // int main NOT void main
{
 
  std::string filename;
  std::string grofilename;
  if (argc > 2) 
  { 
    filename = argv[1]; 
    grofilename = argv[2];
    std::cout << filename << "  " << grofilename << std::endl;
  }
  else 
  {  
    std::cout << "Something went wrong with the arguments" << std::endl; 
    std::cout << "Usage: ./main.o $filename $filename_of_gro" << std::endl;
    exit(0);
  }

  std::unique_ptr<Readfile> ptrReadfile;
  ptrReadfile = std::make_unique<Readfile>(filename, grofilename);

  ptrReadfile->data_to_vector();
  ptrReadfile->get_volume_from_gro();
  std::cout << "time difference between data is:  " << ptrReadfile->calc_timestep() << " ps" << std::endl;

  //std::unique_ptr<EinsteinVisco> ptrEinsteinvisco;
  //ptrEinsteinvisco = std::make_unique<EinsteinVisco>(ptrReadfile);

  std::unique_ptr<GKvisco> ptrGKvisco;
  ptrGKvisco = std::make_unique<GKvisco>(ptrReadfile);

  ptrGKvisco->calcvisco();

  //ptrEinsteinvisco->calcevisco();

  return 0; // everything went right.
} 
