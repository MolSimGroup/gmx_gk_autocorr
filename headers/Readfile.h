#ifndef READFILE_H
#define READFILE_H

class Readfile
{
  private:

    //member variables
    std::string m_filename, m_grofilename;
    std::vector<std::vector<double>> m_alloffdiagonals, m_transposedoffdiag;
    std::vector<double> m_time, m_temperature;
    double m_volume;
    //member functions
    void m_transpose_offdiag();

  public:
    Readfile (std::string filename, std::string grofilename);
    
    void data_to_vector();
    void get_volume_from_gro();
    std::vector<std::vector<double>> get_offdiags();
    std::vector<std::vector<double>> get_transposedoffdiag();
    std::vector<double> get_time();
    std::vector<double> get_temperature();
    double get_volume();

    double calc_timestep();
};

#endif
