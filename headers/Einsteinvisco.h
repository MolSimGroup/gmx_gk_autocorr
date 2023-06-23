#ifndef EINSTEINVISCO_H
#define EINSTEINVISCO_H

#include <Readfile.h>

class EinsteinVisco
{
  private:
    //member variables
    std::vector<std::vector<double>> m_offdiagoriginal, m_alloffdiagint;
    std::vector<std::vector<double>> m_transposesym;
    std::vector<double>              m_temperature, m_volume, m_time;

    double m_volavg, m_tempavg, m_timestep;

    //member functions
    void m_averagevolume();
    void m_averagetemperature();

    void m_transpose_symmetrized(std::vector<std::vector<double>>);
    void m_transpose_integral();
    void m_symmetrize();
    void m_integrate();
    void m_viscosity();

  public:
    //constructor
    EinsteinVisco(std::unique_ptr<Readfile> const &ptrReadfile);


    //functions
    void calcevisco();
};

#endif
