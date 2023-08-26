#include <iostream>
#include <vector>
#include <fstream>  // in- and output of files
#include <string>
#include <sstream>  // stringstream
#include <stdlib.h> // atof (string to double conversion)
#include <algorithm> // transform
#include <cmath>     // exponent for volune

#include <Readfile.h>

//Constructor
Readfile::Readfile (std::string filename, std::string grofilename)
{
  std::vector<double> m_time;
  double              m_volume;
  std::vector<double> m_temperature;
  std::vector<std::vector<double>> m_alloffdiagonals;
  m_filename = filename;
  m_grofilename = grofilename;

  std::cout << "Readfile Constructed" << std::endl;
}

void Readfile::data_to_vector () // int main NOT void main
{
  std::ifstream f(m_filename);
  std::string line;
  std::cout << m_filename << std::endl;

  std::size_t offdiagsize = 9;

  while (std::getline(f, line)) 
  {
    std::istringstream ss(line);
    std::string test;
    ss >> test;

    std::string testfirstchar;
    testfirstchar = test[0];


    if (testfirstchar != "@" && testfirstchar != "#") 
    {
      std::cout << "Testfirstchar" << std::endl;
      double temperature;
      //double volume;
      std::vector<std::string> offdiag(offdiagsize);
      std::string empty;
      std::vector<double> offdiagdouble(offdiagsize);
 
      m_time.push_back(std::atof(test.c_str()));

      // store as Pxy Pxz Pyx Pyz Pzx Pzy
      ss >> temperature >> empty >> offdiag[6] >> offdiag[0] >> offdiag[1] >> offdiag[2] >> offdiag[7] >> offdiag[3] >> offdiag[4] >> offdiag[5] >> offdiag[8];

      std::cout << offdiag[8] << std::endl;

      std::transform(offdiag.begin(), offdiag.end(), offdiagdouble.begin(), [](const std::string& val)
      {
        std::cout << val << std::endl;
        return std::stod(val);
      });

      //m_volume.push_back(volume);
      m_temperature.push_back(temperature);
      m_alloffdiagonals.push_back(offdiagdouble);

      std::cout << "end of first block" << std::endl;

      break;
    }
  }

  while (std::getline(f, line))
  {
    std::istringstream ss(line);

    std::vector<std::string> offdiag(offdiagsize);

    std::string empty;
    std::vector<double> offdiagdouble(offdiag.size());
    double time;
    //double volume;
    double temperature;
    
    ss >> time >> temperature >> empty >> offdiag[6] >> offdiag[0] >> offdiag[1] >> offdiag[2] >> offdiag[7] >> offdiag[3] >> offdiag[4] >> offdiag[5] >> offdiag[8];

    //transform string vector to double vector
    std::transform(offdiag.begin(), offdiag.end(), offdiagdouble.begin(), [](const std::string& val)
    {
      return std::stod(val);
    });

    m_time.push_back(time);
    //m_volume.push_back(volume);
    m_temperature.push_back(temperature);

    m_alloffdiagonals.push_back(offdiagdouble);
  }

  std::cout << "Reading of Data complete. Number of Lines read and stored:  " << m_alloffdiagonals.size() << std::endl;
}

void Readfile::get_volume_from_gro()
{
    std::ifstream fin;
    fin.open(m_grofilename);
    if(fin.is_open()) {
        fin.seekg(-2,std::ios_base::end);                // go to one spot before the EOF

        bool keepLooping = true;
        while(keepLooping) {
            char ch;
            fin.get(ch);                            // Get current byte's data

            if((int)fin.tellg() <= 1) {             // If the data was at or before the 0th byte
                fin.seekg(0);                       // The first line is the last line
                keepLooping = false;                // So stop there
            }
            else if(ch == '\n') {                   // If the data was a newline
                keepLooping = false;                // Stop at the current position.
            }
            else {                                  // If the data was neither a newline nor at the 0 byte
                fin.seekg(-2,std::ios_base::cur);        // Move to the front of that data, then to the front of the data before it
            }
        }

        std::string lastLine;
        getline(fin,lastLine);                      // Read the current line
        std::cout << "Result: " << lastLine << std::endl;     // Display it
        std::istringstream ss(lastLine);
        std::string edgelength;
        ss >> edgelength;
        std::cout << "Edgelength:   " << edgelength << std::endl;
        m_volume = pow(std::stod(edgelength), 3);
        std::cout << "Volume =     " << m_volume << "nm^3" << std::endl;
}
}


void Readfile::m_transpose_offdiag()
{
  for (int i=0; i<9; ++i)
  {
    std::vector<double> temptranspose;
    for (std::size_t j=0; j < m_alloffdiagonals.size(); ++j)
    {
      temptranspose.push_back(m_alloffdiagonals[j][i]);
    }
    m_transposedoffdiag.push_back(temptranspose);
  }
  //for (int i=0; i<6; ++i)
  //{
  //  for (std::size_t j=0; j < m_transposedoffdiag[0].size(); ++j)
  //  {
  //    std::cout << m_transposedoffdiag[i][j] << std::endl;
  //  }
  //}
}


std::vector<std::vector<double>> Readfile::get_offdiags() 
{
  return m_alloffdiagonals;
}

std::vector<std::vector<double>> Readfile::get_transposedoffdiag()
{
  m_transpose_offdiag();
  std::cout << m_transposedoffdiag.size() << "     " << m_transposedoffdiag[0].size() << std::endl;
  return m_transposedoffdiag;
}

std::vector<double> Readfile::get_time()
{
  return m_time;
}

double Readfile::calc_timestep()
{
  double timestep = m_time[1] - m_time[0];
  return timestep;
}

std::vector<double> Readfile::get_temperature()
{
  return m_temperature;
}

double Readfile::get_volume()
{
  return m_volume;
}
