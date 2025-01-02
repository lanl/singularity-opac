// ======================================================================
// © 2021-2024. Triad National Security, LLC. All rights reserved.  This
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

#include <cstring>
#include <type_traits>

#include <ports-of-call/portability.hpp>

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
  static constexpr Real vacuum_permeability = 4.0 * M_PI * 1.0e-7;
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
  static constexpr Real length = 1.e2;          // centimeter
  static constexpr Real mass = 1.e3;            // gram
  static constexpr Real time = 1.;              // second
  static constexpr Real temperature = 1.;       // Kelvin
  static constexpr Real current = 1.e-1;        // Biot
  static constexpr Real charge = 2.997924580e9; // Statcoulomb
};

template <typename BASE = BaseSICODATA2010,
          typename CONVERSION = UnitConversionSIToCGS>
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

  static constexpr Real fine_structure = BASE::fine_structure;
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

  static constexpr Real gravitational_constant = BASE::gravitational_constant *
                                                 length * length * length /
                                                 (mass * time * time);
  static constexpr Real g_newt = gravitational_constant;

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

  // Electron volt (derived)
  static constexpr Real electron_volt = BASE::elementary_charge * energy;
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

struct RuntimePhysicalConstants {
  template <typename T>
  PORTABLE_INLINE_FUNCTION RuntimePhysicalConstants(T pc)
      : length(1.), time(1.), mass(1.), temp(1.), na(pc.na), alpha(pc.alpha),
        h(pc.h), hbar(pc.hbar), kb(pc.kb), r_gas(pc.r_gas), qe(pc.qe), c(pc.c),
        g_newt(pc.g_newt), me(pc.me), mp(pc.mp), mn(pc.mn), amu(pc.amu),
        sb(pc.sb), ar(pc.ar), eV(pc.eV), Fc(pc.Fc), nu_sigma0(pc.nu_sigma0),
        gA(pc.gA) {}

  template <typename T>
  PORTABLE_INLINE_FUNCTION
  RuntimePhysicalConstants(T pc, const Real time_, const Real mass_,
                           const Real length_, const Real temp_)
      : time(time_), mass(mass_), length(length_), temp(temp_), na(pc.na),
        alpha(pc.alpha), h(pc.h * time / mass * std::pow(length, -2)),
        hbar(pc.hbar * time / mass * std::pow(length, -2)),
        kb(pc.kb * std::pow(time, 2) / mass * std::pow(length, -2) * temp),
        r_gas(pc.r_gas * std::pow(time, 2) / mass * std::pow(length, -2)),
        qe(pc.qe), c(pc.c * time / length),
        g_newt(pc.g_newt * std::pow(length, -3) / mass * std::pow(time, 2)),
        me(pc.me / mass), mp(pc.mp / mass), mn(pc.mn / mass),
        amu(pc.amu / mass),
        sb(pc.sb * std::pow(time, 3) / mass * std::pow(temp, 4)),
        ar(pc.ar * std::pow(time, 2) * length / mass * std::pow(temp, 4)),
        eV(pc.eV * std::pow(time, 2) / mass * std::pow(length, -2)),
        Fc(pc.Fc * std::pow(time, 2) * length / mass),
        nu_sigma0(pc.nu_sigma0 / std::pow(length, 2)), gA(pc.gA) {}

  // Conversion factors to the T pc unit system from code units
  const Real time;
  const Real mass;
  const Real length;
  const Real temp;

  const Real na;
  const Real alpha;
  const Real h;
  const Real hbar;
  const Real kb;
  const Real r_gas;
  const Real qe;
  const Real c;
  const Real g_newt;
  const Real me;
  const Real mp;
  const Real mn;
  const Real amu;
  const Real sb;
  const Real ar;
  const Real eV;
  const Real Fc;
  const Real nu_sigma0;
  const Real gA;

  // Runtime physical constants are trivially copyable
  bool operator==(const RuntimePhysicalConstants &other) const {
    return std::memcmp(this, &other, sizeof(*this));
  }
};

using PhysicalConstantsUnity =
    PhysicalConstants<BaseUnity, UnitConversionDefault>;
using PhysicalConstantsSI =
    PhysicalConstants<BaseSICODATA2010, UnitConversionDefault>;
using PhysicalConstantsCGS =
    PhysicalConstants<BaseSICODATA2010, UnitConversionSIToCGS>;

} // namespace singularity

#endif // SINGULARITY_OPAC_CONSTANTS_CONSTANTS_
