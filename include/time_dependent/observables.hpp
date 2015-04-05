#pragma once
#include <petsc_cpp/Petsc.hpp>

namespace Erwin
{

/*
  So, We are interetested in developing a type that contains a list of observers.  Each observer is then called at each timestep (we should give time information to the observers?), and it does its thing.

  Notes:

  //this is run each at each timestep:
  //Examples:  dipole observable,
  //           energy observable,
  //           population observable,
  //           zeros of field observable (using TimeStepper)
  //           eigenstates observable
  //non-constant observables:
  //           absorber mask
  //           eigenstate mask
  void
  Observable::operator()(Matrix& A, Vector& U, TimeStepper& ts);

  //For monitoring (run during the monitor stage:)
  string
  Observable::last();

  //static method for getting the name of the observable:
  string
  Observable::name();

  Inside our propagation code:

  for(auto& O: Observers)
      O(A, U, ts);

  Inside our Monitor code:

  for(auto& O: Observers)
      cout << O.name << ": " << O.last() << ", ";
 */

struct Observable {
    virtual void operator()(const petsc::Matrix& A, const petsc::Vector& U, petsc::TimeStepper& ts) = 0;
    virtual void modify(petsc::Matrix& A, petsc::Vector& U, petsc::Vector& F, petsc::TimeStepper& ts) = 0;
    virtual std::string last() const = 0;
    virtual std::string name() const = 0;
    virtual ~Observable();
};

//to kill the weak vtable error... need to figure out why that happens
inline Observable::~Observable() {}

}
