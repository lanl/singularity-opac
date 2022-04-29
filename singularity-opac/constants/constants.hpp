// ======================================================================
// © 2021. Triad National Security, LLC. All rights reserved.  This
// program was produced under U.S. Government contract
// 89233218CNA000001 for Los Alamos National Laboratory (LANL), which
// is operated by Triad National Security, LLC for the U.S.
// Department of Energy/National Nuclear Security Administration. All
// rights in the program are reserved by Triad National Security, LLC,
// and the U.S. Department of Energy/National Nuclear Security
// Administration. The Government is granted for itself and others
// acting on its behalf a nonexclusive, paid-up, irrevocable worldwide
// license in this material to reproduce, prepare derivative works,
// distribute copies to the public, perform publicly and display
// publicly, and to permit others to do so.
// ======================================================================

#ifndef SINGULARITY_OPAC_CONSTANTS_CONSTANTS_
#define SINGULARITY_OPAC_CONSTANTS_CONSTANTS_

#ifdef SINGULARITY_ENABLE_EXCEPTIONS
#include <stdexcept>
#define CONSTANTS_ERROR(x) (throw std::runtime_error(x))
#else
#define CONSTANTS_ERROR(x) printf("%s", x)
#endif
#define UNDEFINED_ERROR CONSTANTS_ERROR("DEFINE ME\n")

namespace singularity {

struct BaseUnity {
  static constexpr Real avogadro = 1.;
  static constexpr Real fine_structure = 1.;
  static constexpr Real planck = 1.;
  static constexpr Real boltzmann = 1.;
  static constexpr Real elementary_charge = 1.;
  static constexpr Real speed_of_light = 1.;
  static constexpr Real gravitational_constant = 1.;
  static constexpr Real acceleration_from_gravity = 1.;
  static constexpr Real electron_mass = 1.;
  static constexpr Real proton_mass = 1.;
  static constexpr Real neutron_mass = 1.;
  static constexpr Real atomic_mass_unit = 1.;
  static constexpr Real vacuum_permeability = 1.;
  static constexpr Real fermi_coupling_constant = 1.;
  static constexpr Real axial_vector_coupling_constant = -1.;
};

struct BaseSICODATA2010 {
  // Avogadro constant
  static constexpr Real avogadro = 6.02214129e23;
  // Fine structure constant
  static constexpr Real fine_structure = 7.2973525698e-3;
  // Planck constant
  static constexpr Real planck = 6.62606957e-34;
  // Boltzmann constant
  static constexpr Real boltzmann = 1.380648800e-23;
  // Elementary charge
  static constexpr Real elementary_charge = 1.602176565e-19;
  // Speed of light (exact)
  static constexpr Real speed_of_light = 2.99792458e8;
  // Gravitational constant
  static constexpr Real gravitational_constant = 6.67384e-11;
  // Standard acceleration of gravity
  static constexpr Real acceleration_from_gravity = 9.80665;
  // Electron rest mass
  static constexpr Real electron_mass = 9.10938291e-31;
  // Proton rest mass
  static constexpr Real proton_mass = 1.672621777e-27;
  // Neutron rest mass
  static constexpr Real neutron_mass = 1.674927351e-27;
  // Atomic mass unit
  static constexpr Real atomic_mass_unit = 1.660538921e-27;
  // Permeability of free space
  static constexpr Real vacuum_permeability =
      4.0 * M_PI * 1.0e-7;
  // Fermi coupling constant (Units of erg^{-2})
  static constexpr Real fermi_coupling_constant = 4.543791885043567014;
  // Axial-vector coupling
  static constexpr Real axial_vector_coupling_constant = -1.23;
};

struct UnitConversionDefault {
  static constexpr Real length = 1.;
  static constexpr Real mass = 1.;
  static constexpr Real time = 1.;
  static constexpr Real temperature = 1.;
  static constexpr Real current = 1.;
  static constexpr Real charge = 1.;
};

struct UnitConversionSIToCGS {
  static constexpr Real length = 1.e2;                    // centimeter
  static constexpr Real mass = 1.e3;                      // gram
  static constexpr Real time = 1.;                        // second
  static constexpr Real temperature = 1.;                 // Kelvin
  static constexpr Real current = 1.e-1;                  // Biot
  static constexpr Real charge = 2.997924580e9;           // Statcoulomb
};

template <typename BASE=BaseSICODATA2010, typename CONVERSION=UnitConversionSIToCGS>
struct PhysicalConstants {
 protected:
  static constexpr Real length = CONVERSION::length;
  static constexpr Real mass = CONVERSION::mass;
  static constexpr Real time = CONVERSION::time;
  static constexpr Real temperature = CONVERSION::temperature;
  static constexpr Real current = CONVERSION::current;
  static constexpr Real charge = CONVERSION::charge;

  // Derived unit conversions
  static constexpr Real force = mass * length / (time * time);
  static constexpr Real energy = force * length;
  static constexpr Real power = energy / time;

 public:
  static constexpr Real avogadro = BASE::avogadro;
  static constexpr Real na = avogadro;

  static constexpr Real fine_structure = BASE::fine_structore;
  static constexpr Real alpha = fine_structure;

  static constexpr Real planck = BASE::planck * energy * time;
  static constexpr Real h = planck;

  // Reduced Planck constant (derived)
  static constexpr Real reduced_planck = planck / (2.0 * M_PI);
  static constexpr Real hbar = reduced_planck;

  static constexpr Real boltzmann = BASE::boltzmann * energy / temperature;
  static constexpr Real kb = boltzmann;

  // Molar gas constant (derived)
  static constexpr Real gas_constant = avogadro * boltzmann;
  static constexpr Real r_gas = gas_constant;

  static constexpr Real elementary_charge = BASE::elementary_charge * charge;
  static constexpr Real qe = elementary_charge;

  static constexpr Real speed_of_light = BASE::speed_of_light * length / time;
  static constexpr Real c = speed_of_light;

  static constexpr Real gravitational_constant =
      BASE::gravitational_constant * length * length * length / (mass * time * time);
  static constexpr Real g_newt = gravitational_constant;

  static constexpr Real acceleration_from_gravity =
      BASE::acceleration_from_gravity * length / (time * time);
  static constexpr Real g_accel = acceleration_from_gravity;

  static constexpr Real electron_mass = BASE::electron_mass * mass;
  static constexpr Real me = electron_mass;

  static constexpr Real proton_mass = BASE::proton_mass * mass;
  static constexpr Real mp = proton_mass;

  static constexpr Real neutron_mass = BASE::neutron_mass * mass;
  static constexpr Real mn = neutron_mass;

  static constexpr Real atomic_mass_unit = BASE::atomic_mass_unit * mass;
  static constexpr Real amu = atomic_mass_unit;

  // Stefan-Boltzmann constant (derived)
  static constexpr Real stefan_boltzmann = 2.0 * M_PI * M_PI * M_PI * M_PI *
                                             M_PI * kb * kb * kb * kb /
                                             (15.0 * h * h * h * c * c);
  static constexpr Real sb = stefan_boltzmann;

  // Radiation constant (derived)
  static constexpr Real radiation_constant = 4.0 * sb / c;
  static constexpr Real ar = radiation_constant;

  // Faraday constant (derived)
  static constexpr Real faraday_constant = na * qe;
  static constexpr Real faraday = faraday_constant;

  // Permeability of free space
  static constexpr Real vacuum_permeability =
      BASE::vacuum_permeability * force / (current * current);
  static constexpr Real mu0 = vacuum_permeability;

  // Permittivity of free space (derived)
  static constexpr Real vacuum_permittivity = 1.0 / (mu0 * c * c);
  static constexpr Real eps0 = vacuum_permittivity;

  // Electron volt (derived)
  static constexpr Real electron_volt = qe * energy;
  static constexpr Real eV = electron_volt;

  static constexpr Real fermi_coupling_constant =
      BASE::fermi_coupling_constant / (energy * energy);
  static constexpr Real Fc = fermi_coupling_constant;

  // Fiducial weak cross section (derived; units of cm^2)
  static constexpr Real fiducial_weak_cross_section =
      4. * length * length * Fc * Fc * c * c * hbar * hbar * me * me * c * c *
      c * c / M_PI;
  static constexpr Real nu_sigma0 = fiducial_weak_cross_section;

  static constexpr Real axial_vector_coupling_constant = -1.23;
  static constexpr Real gA = axial_vector_coupling_constant;
};

using PhysicalConstantsUnity = PhysicalConstants<BaseUnity, UnitConversionDefault>;
using PhysicalConstantsSI = PhysicalConstants<BaseSICODATA2010, UnitConversionDefault>;
using PhysicalConstantsCGS = PhysicalConstants<BaseSICODATA2010, UnitConversionSIToCGS>;

/*
// Unit system conversion factors

struct SI {
  static constexpr double length = 1.;      // meter
  static constexpr double mass = 1.;        // kilogram
  static constexpr double time = 1.;        // second
  static constexpr double temperature = 1.; // Kelvin
  static constexpr double current = 1.;     // Amp
  static constexpr double charge = 1.;      // Coulomb
  static constexpr double capacitance = 1.; // Farad
  static constexpr double angle = 1.;       // Radian
};

struct CGS {
  static constexpr double length = 1.e2;                    // centimeter
  static constexpr double mass = 1.e3;                      // gram
  static constexpr double time = 1.;                        // second
  static constexpr double temperature = 1.;                 // Kelvin
  static constexpr double current = 1.e-1;                  // Biot
  static constexpr double charge = 2.997924580e9;           // Statcoulomb
  static constexpr double capacitance = 8.9831483395497e11; // Statfarad
  static constexpr double angle = 1.;                       // Radian
};

struct HVal1 {
  static constexpr double h = 1.;
};

template <typename H=HVal1>
class TrivialConstants {
  public:
  constexpr TrivialConstants() {}

  static constexpr double avogadro = 1.;
  static constexpr double na = avogadro;
  static constexpr double fine_structure = 1.;
  static constexpr double alpha = fine_structure;
  static constexpr double planck = 1.;//H::h;
  static constexpr double h = planck;
  static constexpr double reduced_planck = planck / (2.0 * M_PI);
  static constexpr double hbar = reduced_planck;
  static constexpr double gas_constant = 1.;
  static constexpr double r_gas = gas_constant;
  static constexpr double boltzmann = 1.;
  static constexpr double kb = boltzmann;
  static constexpr double electron_charge = 1.;
  static constexpr double qe = electron_charge;
  static constexpr double speed_of_light = 1.;
  static constexpr double c = speed_of_light;
  static constexpr double gravitational_constant = 1.;
  static constexpr double g_newt = gravitational_constant;
  static constexpr double acceleration_from_gravity = 1.;
  static constexpr double g_accel = acceleration_from_gravity;
  static constexpr double electron_mass = 1.;
  static constexpr double me = electron_mass;
  static constexpr double proton_mass = 1.;
  static constexpr double mp = proton_mass;
  static constexpr double stefan_boltzmann = 1.;
  static constexpr double sb = stefan_boltzmann;
  static constexpr double faraday_constant = 1.;
  static constexpr double faraday = faraday_constant;
  static constexpr double vacuum_permeability = 1.;
  static constexpr double mu0 = vacuum_permeability;
  static constexpr double vacuum_permittivity = 1.;
  static constexpr double eps0 = vacuum_permittivity;
  static constexpr double classical_electron_radius = 1.;
  static constexpr double re = classical_electron_radius;
  static constexpr double electron_volt = 1.;
  static constexpr double eV = electron_volt;
  static constexpr double atomic_mass_unit = 1.;
  static constexpr double amu = atomic_mass_unit;
  static constexpr double fermi_coupling_constant = 1.;
  static constexpr double Fc = fermi_coupling_constant;
  static constexpr double fiducial_weak_cross_section =
      4. * Fc * Fc * c * c * hbar * hbar * me * me * c * c *
      c * c / M_PI;
  static constexpr double nu_sigma0 = fiducial_weak_cross_section;
  static constexpr double axial_vector_coupling_constant = 1.;
  static constexpr double gA = axial_vector_coupling_constant;

};

// Defines and encapsulates physical and mathematical constants in a purely
// constexpr way, with support for multiple unit systems.
template <typename UNITSYSTEM>
class PhysicalConstants {
 protected:
  static constexpr double length = UNITSYSTEM::length;
  static constexpr double mass = UNITSYSTEM::mass;
  static constexpr double time = UNITSYSTEM::time;
  static constexpr double temperature = UNITSYSTEM::temperature;
  static constexpr double current = UNITSYSTEM::current;
  static constexpr double charge = UNITSYSTEM::charge;
  static constexpr double capacitance = UNITSYSTEM::capacitance;
  static constexpr double angle = UNITSYSTEM::angle;

  // Derived unit conversions
  static constexpr double force = mass * length / (time * time);
  static constexpr double energy = force * length;
  static constexpr double power = energy / time;

 public:
  constexpr PhysicalConstants() {}

  // Avogadro constant (CODATA 2010 value)
  static constexpr double avogadro = 6.02214129e23;
  static constexpr double na = avogadro;

  // Fine structure constant (CODATA 2010 value)
  static constexpr double fine_structure = 7.2973525698e-3;
  static constexpr double alpha = fine_structure;

  // Planck constant (CODATA 2010 value)
  static constexpr double planck = 6.62606957e-34 * energy * time;
  static constexpr double h = planck;

  // Reduced Planck constant
  static constexpr double reduced_planck = planck / (2.0 * M_PI);
  static constexpr double hbar = reduced_planck;

  // Molar gas constant (CODATA 2010 value)
  static constexpr double gas_constant = 8.3144621 * energy / temperature;
  static constexpr double r_gas = gas_constant;

  // Boltzmann constant (CODATA 2010 value)
  static constexpr double boltzmann = 1.380648800e-23 * energy / temperature;
  static constexpr double kb = boltzmann;

  // Electron charge (CODATA 2018 exact value)
  static constexpr double electron_charge = 1.602176565e-19 * charge;
  static constexpr double qe = electron_charge;

  // Speed of light (CODATA 2018 exact value)
  static constexpr double speed_of_light = 2.99792458e8 * length / time;
  static constexpr double c = speed_of_light;

  // Gravitational constant (CODATA 2010 value)
  static constexpr double gravitational_constant =
      6.67384e-11 * length * length * length / (mass * time * time);
  static constexpr double g_newt = gravitational_constant;

  // Standard acceleration of gravity (CODATA 2010 value)
  static constexpr double acceleration_from_gravity =
      9.80665 * length / (time * time);
  static constexpr double g_accel = acceleration_from_gravity;

  // Electron rest mass (CODATA 2010 value)
  static constexpr double electron_mass = 9.10938291e-31 * mass;
  static constexpr double me = electron_mass;

  // Proton rest mass (CODATA 2010 value)
  static constexpr double proton_mass = 1.672621777e-27 * mass;
  static constexpr double mp = proton_mass;

  // Stefan-Boltzmann constant
  static constexpr double stefan_boltzmann = 2.0 * M_PI * M_PI * M_PI * M_PI *
                                             M_PI * kb * kb * kb * kb /
                                             (15.0 * h * h * h * c * c);
  static constexpr double sb = stefan_boltzmann;

  // Faraday constant
  static constexpr double faraday_constant = 96485.33645957 * capacitance;
  static constexpr double faraday = faraday_constant;

  // Permeability of free space
  static constexpr double vacuum_permeability =
      4.0 * M_PI * 1.0e-7 * force / (current * current);
  static constexpr double mu0 = vacuum_permeability;

  // Permittivity of free space
  static constexpr double vacuum_permittivity =
      8.85418782e-12 * capacitance / length;
  static constexpr double eps0 = vacuum_permittivity;

  // Classical electron radius
  static constexpr double classical_electron_radius = 2.81794033e-15 * length;
  static constexpr double re = classical_electron_radius;

  // Electron volt
  static constexpr double electron_volt = 1.602176565e-19 * energy;
  static constexpr double eV = electron_volt;

  // Atomic mass unit (CODATA 2010 value)
  static constexpr double atomic_mass_unit = 1.660538921e-27 * mass;
  static constexpr double amu = atomic_mass_unit;

  // Fermi coupling constant (Units of erg^{-2})
  static constexpr double fermi_coupling_constant =
      4.543791885043567014 / (energy * energy);
  static constexpr double Fc = fermi_coupling_constant;

  // Fiducial weak cross section (Units of cm^2)
  static constexpr double fiducial_weak_cross_section =
      4. * length * length * Fc * Fc * c * c * hbar * hbar * me * me * c * c *
      c * c / M_PI;
  static constexpr double nu_sigma0 = fiducial_weak_cross_section;

  // Axial-vector coupling
  static constexpr double axial_vector_coupling_constant = -1.23;
  static constexpr double gA = axial_vector_coupling_constant;
};

// These can be removed for C++17 and up
template <typename T>
constexpr double PhysicalConstants<T>::avogadro;
template <typename T>
constexpr double PhysicalConstants<T>::na;
template <typename T>
constexpr double PhysicalConstants<T>::fine_structure;
template <typename T>
constexpr double PhysicalConstants<T>::alpha;
template <typename T>
constexpr double PhysicalConstants<T>::planck;
template <typename T>
constexpr double PhysicalConstants<T>::h;
template <typename T>
constexpr double PhysicalConstants<T>::reduced_planck;
template <typename T>
constexpr double PhysicalConstants<T>::hbar;
template <typename T>
constexpr double PhysicalConstants<T>::gas_constant;
template <typename T>
constexpr double PhysicalConstants<T>::r_gas;
template <typename T>
constexpr double PhysicalConstants<T>::boltzmann;
template <typename T>
constexpr double PhysicalConstants<T>::kb;
template <typename T>
constexpr double PhysicalConstants<T>::electron_charge;
template <typename T>
constexpr double PhysicalConstants<T>::qe;
template <typename T>
constexpr double PhysicalConstants<T>::speed_of_light;
template <typename T>
constexpr double PhysicalConstants<T>::c;
template <typename T>
constexpr double PhysicalConstants<T>::gravitational_constant;
template <typename T>
constexpr double PhysicalConstants<T>::g_newt;
template <typename T>
constexpr double PhysicalConstants<T>::acceleration_from_gravity;
template <typename T>
constexpr double PhysicalConstants<T>::g_accel;
template <typename T>
constexpr double PhysicalConstants<T>::electron_mass;
template <typename T>
constexpr double PhysicalConstants<T>::me;
template <typename T>
constexpr double PhysicalConstants<T>::proton_mass;
template <typename T>
constexpr double PhysicalConstants<T>::mp;
template <typename T>
constexpr double PhysicalConstants<T>::stefan_boltzmann;
template <typename T>
constexpr double PhysicalConstants<T>::sb;
template <typename T>
constexpr double PhysicalConstants<T>::faraday_constant;
template <typename T>
constexpr double PhysicalConstants<T>::faraday;
template <typename T>
constexpr double PhysicalConstants<T>::vacuum_permeability;
template <typename T>
constexpr double PhysicalConstants<T>::mu0;
template <typename T>
constexpr double PhysicalConstants<T>::vacuum_permittivity;
template <typename T>
constexpr double PhysicalConstants<T>::eps0;
template <typename T>
constexpr double PhysicalConstants<T>::classical_electron_radius;
template <typename T>
constexpr double PhysicalConstants<T>::re;
template <typename T>
constexpr double PhysicalConstants<T>::electron_volt;
template <typename T>
constexpr double PhysicalConstants<T>::eV;
template <typename T>
constexpr double PhysicalConstants<T>::atomic_mass_unit;
template <typename T>
constexpr double PhysicalConstants<T>::amu;
template <typename T>
constexpr double PhysicalConstants<T>::fermi_coupling_constant;
template <typename T>
constexpr double PhysicalConstants<T>::Fc;
template <typename T>
constexpr double PhysicalConstants<T>::fiducial_weak_cross_section;
template <typename T>
constexpr double PhysicalConstants<T>::nu_sigma0;
template <typename T>
constexpr double PhysicalConstants<T>::axial_vector_coupling_constant;
template <typename T>
constexpr double PhysicalConstants<T>::gA;
*/
} // namespace singularity

#undef CONSTANTS_ERROR
#undef UNDEFINED_ERROR

#endif // SINGULARITY_OPAC_CONSTANTS_CONSTANTS_
