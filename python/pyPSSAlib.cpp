
// Build pSSAlib using:
// CC=$HOME/bin/gcc CXX=$HOME/bin/g++ CFLAGS="-fPIC" CXXFLAGS="-fPIC" ./configure --prefix=$HOME/.local --disable-sbml --disable-optimisation
// make clean; make -j && make install


// To compile:
// TODO

#include <pssalib/PSSA.h>
#include <pssalib/util/Timing.h>

#include <time.h>       /* time_t, struct tm, time, localtime */

#include <vector>
#include <string>
#include <iostream>

#include <boost/config.hpp>
#include <boost/current_function.hpp>

#include <boost/python.hpp>
#include <boost/python/numpy.hpp>
#include <boost/python/stl_iterator.hpp>
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>

#include <boost/algorithm/string/join.hpp>
#include <boost/algorithm/string/case_conv.hpp>

#ifdef BOOST_NO_CXX11_HDR_THREAD
#  include <boost/thread/thread_only.hpp>
#  include <boost/chrono/time_point.hpp>
#else
#  include <thread>
#endif

// inspired by https://stackoverflow.com/questions/30670909/boostpython-extract-c-value-from-numpy-ndarray-return-value
#define BOOST_PY_ASSERT(expr) { if(!(expr)) { \
  PyErr_SetString(PyExc_TypeError, (boost::format("BOOST_PY_ASSERT(%1%:%2%): !(%3%) '%4%'") % (__FILE__) % (__LINE__) % (expr) % (#expr)).str().c_str()); \
  boost::python::throw_error_already_set(); \
}; };

#define BOOST_PY_ERRMSG(msg) \
  { PyErr_SetString(PyExc_TypeError, (boost::format("BOOST_PY_ERRMSG(%1%:%2%): %3%") % (__FILE__) % (__LINE__) % (msg)).str().c_str()); };

#include "../include/tools.h"
#include "../include/models.h"

class pSSAlibWrapper
{
protected:
  enum tagTestCase
  {
    tcCLC,
    tcCA,
    tcHomo,
    tcSBD,
    tcTCS,
    tcED,
    tcNone
  };

  enum tagTestCase testCase;

#ifdef BOOST_NO_CXX11_HDR_ATOMIC
  //! Atomic flag to interrupt execution
  boost::atomic<bool>           bProcessing;
#else
  //! Atomic flag to interrupt execution
  std::atomic<bool>             bProcessing;
#endif

  //! Storage for simulation results
  std::unique_ptr<unsigned int> arRawTrajectories;
  size_t                        szTimePointsCount,
                                szSpeciesCount;

  /**
   * Initialize the data model for a test case
   * 
   * @param model reference to a @link pssalib::datamodel::detail::Model object
   * @return @true if initialization succeeds, @false otherwise
   */
  bool initializeModel(pssalib::datamodel::detail::Model &model) const
  {
    boost::python::stl_input_iterator<double> begin(lstParams), end;
    std::vector<double> vParams(begin, end);

    double * arParams = vParams.data();
    size_t   szParams = vParams.size();

    switch(testCase)
    {
    case tcCLC:
      generateCyclicLinearChain(model, arParams, szParams);
      break;
    case tcCA:
      generateColloidalAggregation(model, arParams, szParams);
      break;
    case tcHomo:
      generateHomoreaction(model, arParams, szParams);
      break;
    case tcSBD:
      generateSingleBirthDeath(model, arParams, szParams);
      break;
    case tcTCS:
      generateTwoComponentSystem(model, arParams, szParams);
      break;
    case tcED:
      generateEnzymaticDegradation(model, arParams, szParams);
      break;
    default:
      BOOST_PY_ERRMSG("Invalid test case identifier")
      return false;
    }

    return true;
  }

  bool initializeSimulation(pssalib::datamodel::SimulationInfo &simInfo)
  {
    if(!initializeModel(simInfo.getModel()))
      return false;

    // setup the simulation params
    simInfo.unSamplesTotal = 1;
    simInfo.dTimeStart = dTimeStart;
    simInfo.dTimeEnd = dTimeEnd;
    simInfo.dTimeStep = dTimeStep;

    // determine the dimensions for simulation results (# species x # time points)
    szSpeciesCount = simInfo.getModel().getSpeciesCount();
    szTimePointsCount = pssalib::timing::getNumTimePoints(simInfo.dTimeStart, simInfo.dTimeEnd, simInfo.dTimeStep);

    // storage for trajectories
    arRawTrajectories.reset(new unsigned int [szTimePointsCount*szSpeciesCount]);
//     std::fill(arRawTrajectories.begin(), arRawTrajectories.end(), 0);
    simInfo.ptrarRawTrajectory = arRawTrajectories.get();

    // setup output
    simInfo.unOutputFlags = pssalib::datamodel::SimulationInfo::ofNone
      | pssalib::datamodel::SimulationInfo::ofRawTrajectory
      | pssalib::datamodel::SimulationInfo::ofTrajectory
//       | pssalib::datamodel::SimulationInfo::ofStatus
      | pssalib::datamodel::SimulationInfo::ofLog
//       | pssalib::datamodel::SimulationInfo::ofTrace
//       | pssalib::datamodel::SimulationInfo::ofInfo
//       | pssalib::datamodel::SimulationInfo::ofWarning
      | pssalib::datamodel::SimulationInfo::ofError
//       | pssalib::datamodel::SimulationInfo::eofModuleGrouping
//       | pssalib::datamodel::SimulationInfo::eofModuleSampling
//       | pssalib::datamodel::SimulationInfo::eofModuleUpdate
      ;
    simInfo.resetOutput();

    simInfo.setOutputStreamBuf(pssalib::datamodel::SimulationInfo::ofLog, std::cerr.rdbuf());

    return true;
  }

  void runSimulation(pssalib::datamodel::SimulationInfo *ptrSimInfo)
  {
    pssalib::PSSA simEngine;

    simEngine.setMethod(pssalib::PSSA::M_PDM);

    simEngine.run(ptrSimInfo);

    bProcessing = false;
  }

//
public:

  //! list of simulation parameters
  boost::python::list lstParams;
  //! start time of the simulation output (any results before this time point are excluded from the output)
  double dTimeStart;
  //! time step of the simulation output
  double dTimeStep;
  //! end time of the simulation output (simulation will end at this time point)
  double dTimeEnd;

  //! Default constructor
  pSSAlibWrapper():
    bProcessing(false),
    szTimePointsCount(0),
    szSpeciesCount(0),
    testCase(tcNone),
    dTimeStart(0.0),
    dTimeStep(1e-1),
    dTimeEnd(1.0)
  {
    // do nothing
  }

  //! Copy constructor
  pSSAlibWrapper(const pSSAlibWrapper &other):
    bProcessing(other.bProcessing.load()),
    szTimePointsCount(other.szTimePointsCount),
    szSpeciesCount(other.szSpeciesCount),
    testCase(other.testCase),
    lstParams(other.lstParams),
    dTimeStart(other.dTimeStart),
    dTimeStep(other.dTimeStep),
    dTimeEnd(other.dTimeEnd)
  {
    // do nothing
    arRawTrajectories.reset(new unsigned int [szTimePointsCount*szSpeciesCount]);
  }

  std::string str() const
  {
    namespace bpy = boost::python;

    std::string params_str = bpy::extract<std::string>(bpy::str(lstParams));

    return (boost::format("tc='%1%'; params=%2%; start=%3%; step=%4%; end=%5%") % getTestCase() % params_str.c_str() % dTimeStart % dTimeStep % dTimeEnd).str();
  }

  std::string repr() const
  {
    namespace bpy = boost::python;

    std::string params_str = bpy::extract<std::string>(bpy::str(lstParams));

    return (boost::format("pSSAlib(tc='%1%'; params=%2%; start=%3%; step=%4%; end=%5%)") % getTestCase() % params_str.c_str() % dTimeStart % dTimeStep % dTimeEnd).str();
  }

  /**
   * Run the simulation
   * @return @true if simulation run succeeds, @false otherwise
   */
  bool run()
  {
    if(bProcessing) return false;
    bProcessing = true;

    pssalib::datamodel::SimulationInfo simInfo;

    if(!initializeSimulation(simInfo))
    {
      bProcessing = false;
      return false;
    }

    // run the simulation
#ifdef BOOST_NO_CXX11_HDR_THREAD
    boost::thread
#else
    std::thread
#endif
      simThread(&pSSAlibWrapper::runSimulation, this, &simInfo);

    while(bProcessing)
    {
#ifdef BOOST_NO_CXX11_HDR_THREAD
      boost::this_thread::sleep_for(boost::chrono::seconds(1));
#else      
      std::this_thread::sleep_for(std::chrono::seconds(1));
#endif
      if(PyErr_CheckSignals() == -1)
      {
        BOOST_PY_ERRMSG("simulation cancelled by keyboard interrupt")

        simInfo.bInterruptRequested = true;
      }
    }

    simThread.join();

    // done
    bProcessing = false;
    return true ^ simInfo.bInterruptRequested;
  }

  /**
   * Set test case name
   * 
   * @param strTCId id of the test case, must be one of 'clc', 'ca', 'homo', 'sbd', 'tcs' or 'ed'.
   */
  void setTestCase(std::string strTCId)
  {
    boost::to_lower(strTCId);
    if(strTCId == "clc")
      testCase = tcCLC;
    else if(strTCId == "ca")
      testCase = tcCA;
    else if(strTCId == "homo")
      testCase = tcHomo;
    else if(strTCId == "sbd")
      testCase = tcSBD;
    else if(strTCId == "tcs")
      testCase = tcTCS;
    else if(strTCId == "ed")
      testCase = tcED;
    else
      testCase = tcNone;
  }

  /**
   * Get test case name
   * 
   * @return id of the test case, one of 'clc', 'ca', 'homo', 'sbd', 'tcs' or 'ed'.
   */
  std::string getTestCase(void) const
  {
    switch(testCase)
    {
      case tcCLC:
        return "clc";
        break;
      case tcCA:
        return "ca";
        break;
      case tcHomo:
        return "homo";
        break;
      case tcSBD:
        return "sbd";
        break;
      case tcTCS:
        return "tcs";
        break;
      case tcED:
        return "ed";
        break;
      case tcNone:
      default:
        return "";
        break;
    }
  }

  /**
   * Get simulation results
   * 
   * @return numpy.ndarray with number of rows equal to number of time points and number of column equal to number of species in the simulation. An empty array object if no simulation was performed.
   */
  boost::python::numpy::ndarray getResult() const
  {
    namespace bpy = boost::python;
    return bpy::numpy::from_data(
            arRawTrajectories.get(),
            bpy::numpy::dtype::get_builtin<unsigned int>(),
            bpy::make_tuple(szTimePointsCount,szSpeciesCount),
            bpy::make_tuple(sizeof(unsigned int)*szSpeciesCount,sizeof(unsigned int)),
            bpy::object());
  }

  /**
   * Get number of simulation time points
   * 
   * @return number of time points between start time and end time with time step intervals.
   */
  std::size_t getTimePointsCount() const
  {
    if(szTimePointsCount == 0)
      return pssalib::timing::getNumTimePoints(dTimeStart, dTimeEnd, dTimeStep);
    return szTimePointsCount;
  };

  /**
   * Get number species
   * 
   * @return number of species in the model.
   */
  std::size_t getSpeciesCount() const
  {
//     namespace bpy = boost::python;
// 
//     pssalib::datamodel::detail::Model model;
//     if(!initializeModel(model))
//     {
//       bpy::throw_error_already_set();
//     }
// 
//     return model.getSpeciesCount();
    return szSpeciesCount;
  }

  boost::python::numpy::ndarray getODEs(boost::python::list &lstPopulation)
  {
    namespace bpy = boost::python;

    pssalib::datamodel::detail::Model model;
    if(!initializeModel(model))
    {
      bpy::throw_error_already_set();
    }
    // normalize the model
    model.normalize();

    bpy::stl_input_iterator<double> begin(lstPopulation), end;
    std::vector<double> vPopulation(begin, end);

    if(vPopulation.size() < model.getSpeciesCount())
    {
      BOOST_PY_ERRMSG("invalid arguments: population array size is less than number of species in the model")
      bpy::throw_error_already_set();
    }

    // compute the ODEs
    boost::scoped_array<double> ardF(new double[model.getReactionsCount()]);
    bpy::numpy::ndarray result = bpy::numpy::zeros(bpy::make_tuple(model.getSpeciesCount()), bpy::numpy::dtype::get_builtin<double>());

    computeF(model, vPopulation.data(), ardF.get());
    computeODEs(model, ardF.get(), reinterpret_cast<double*>(result.get_data()));

    return result;
  }

  boost::python::numpy::ndarray getJacobian(boost::python::list &lstPopulation)
  {
    namespace bpy = boost::python;

    pssalib::datamodel::detail::Model model;
    if(!initializeModel(model))
    {
      bpy::throw_error_already_set();
    }
    // normalize the model
    model.normalize();

    bpy::stl_input_iterator<double> begin(lstPopulation), end;
    std::vector<double> vPopulation(begin, end);

    if(vPopulation.size() < model.getSpeciesCount())
    {
      BOOST_PY_ERRMSG("invalid arguments: population array size is less than number of species in the model")
      bpy::throw_error_already_set();
    }

    // compute the Jacobian matrix
    boost::scoped_array<double> ardF(new double[model.getReactionsCount()]);
    bpy::numpy::ndarray result = bpy::numpy::zeros(bpy::make_tuple(model.getSpeciesCount(),model.getSpeciesCount()), bpy::numpy::dtype::get_builtin<double>());

    computeJ(model, vPopulation.data(), reinterpret_cast<double*>(result.get_data()));

    return result;
  }

  boost::python::numpy::ndarray getLyapunovQ(boost::python::list &lstPopulation)
  {
    namespace bpy = boost::python;

    pssalib::datamodel::detail::Model model;
    if(!initializeModel(model))
    {
      bpy::throw_error_already_set();
    }
    // normalize the model
    model.normalize();

    bpy::stl_input_iterator<double> begin(lstPopulation), end;
    std::vector<double> vPopulation(begin, end);

    if(vPopulation.size() < model.getSpeciesCount())
    {
      BOOST_PY_ERRMSG("invalid arguments: population array size is less than number of species in the model")
      bpy::throw_error_already_set();
    }

    // compute the Jacobian matrix
    boost::scoped_array<double> ardF(new double[model.getReactionsCount()]);
    bpy::numpy::ndarray result = bpy::numpy::zeros(bpy::make_tuple(model.getSpeciesCount(),model.getSpeciesCount()), bpy::numpy::dtype::get_builtin<double>());

    computeQ(model, vPopulation.data(), reinterpret_cast<double*>(result.get_data()));

    return result;
  }

  boost::python::numpy::ndarray getStoichiometry(boost::python::list &lstPopulation)
  {
    namespace bpy = boost::python;

    pssalib::datamodel::detail::Model model;
    if(!initializeModel(model))
    {
      bpy::throw_error_already_set();
    }
    // normalize the model
    model.normalize();

    bpy::stl_input_iterator<double> begin(lstPopulation), end;
    std::vector<double> vPopulation(begin, end);

    if(vPopulation.size() < model.getSpeciesCount())
    {
      BOOST_PY_ERRMSG("invalid arguments: population array size is less than number of species in the model")
      bpy::throw_error_already_set();
    }

    // compute the stoichiometry matrix
    bpy::numpy::ndarray result = bpy::numpy::zeros(bpy::make_tuple(model.getReactionsCount(),model.getSpeciesCount()), bpy::numpy::dtype::get_builtin<double>());

    computeS(model, vPopulation.data(), reinterpret_cast<double*>(result.get_data()));

    return result;
  }
};


boost::python::numpy::ndarray moments(boost::python::numpy::ndarray &input, boost::python::list &order)
{
  namespace bpy = boost::python;

  BOOST_PY_ASSERT((bpy::len(order) != 0))
  BOOST_PY_ASSERT((bpy::len(input) != 0))
  BOOST_PY_ASSERT((input.get_nd() == 2))

  bpy::stl_input_iterator<unsigned int> begin(order), end;
  std::vector<unsigned int> vOrder(begin, end);

  size_t szRowCount = input.shape(0);
  size_t szColCount = input.shape(1);

  bpy::numpy::ndarray result = bpy::numpy::zeros(bpy::make_tuple(vOrder.size(),szColCount), bpy::numpy::dtype::get_builtin<double>());

  if(bpy::numpy::equivalent(input.get_dtype(), bpy::numpy::dtype::get_builtin<double>()))
    computeMoments<double>(reinterpret_cast<double *>(input.get_data()), szRowCount, szColCount, vOrder.data(), vOrder.size(), reinterpret_cast<double*>(result.get_data()));
  else if(bpy::numpy::equivalent(input.get_dtype(), bpy::numpy::dtype::get_builtin<unsigned int>()))
    computeMoments<unsigned int>(reinterpret_cast<unsigned int *>(input.get_data()), szRowCount, szColCount, vOrder.data(), vOrder.size(), reinterpret_cast<double*>(result.get_data()));
  else if(bpy::numpy::equivalent(input.get_dtype(), bpy::numpy::dtype::get_builtin<int>()))
    computeMoments<int>(reinterpret_cast<int *>(input.get_data()), szRowCount, szColCount, vOrder.data(), vOrder.size(), reinterpret_cast<double*>(result.get_data()));
  else
  {
    BOOST_PY_ERRMSG("unsupported data type")
    bpy::throw_error_already_set();
  }

  return result;
};

struct pSSAlibWrapper_pickle_suite : boost::python::pickle_suite
{
  static
  boost::python::tuple
  getinitargs(const pSSAlibWrapper& o)
  {
    return boost::python::tuple();
  }

  static
  boost::python::tuple
  getstate(const pSSAlibWrapper& o)
  {
    return boost::python::make_tuple(o.getTestCase().c_str(), o.lstParams, o.dTimeStart, o.dTimeStep, o.dTimeEnd);
  }

  static
  void
  setstate(pSSAlibWrapper& o, boost::python::tuple state)
  {
    namespace bpy = boost::python;
    o.setTestCase(std::string(bpy::extract<std::string::const_pointer>(state[0])));
    o.lstParams = bpy::extract<bpy::list>(state[1]);
    o.dTimeStart = bpy::extract<double>(state[2]);
    o.dTimeStep = bpy::extract<double>(state[3]);
    o.dTimeEnd = bpy::extract<double>(state[4]);
  }
};


BOOST_PYTHON_MODULE(pyPSSAlib)
{
  Py_Initialize();
  boost::python::numpy::initialize();

  def("moments", &moments);

  boost::python::class_<pSSAlibWrapper>("pSSAlib")
    .def_pickle(pSSAlibWrapper_pickle_suite())
    .def("__str__", &pSSAlibWrapper::str)
    .def("__repr__", &pSSAlibWrapper::repr)
    .def("run", &pSSAlibWrapper::run)
    .def("odes", &pSSAlibWrapper::getODEs)
    .def("jacobian", &pSSAlibWrapper::getJacobian)
    .def("lyapunovQ", &pSSAlibWrapper::getLyapunovQ)
    .def("stoichiometry", &pSSAlibWrapper::getStoichiometry)
    .def_readonly("num_timepoints", &pSSAlibWrapper::getTimePointsCount)
    .def_readonly("num_species", &pSSAlibWrapper::getSpeciesCount)
    .def_readwrite("params", &pSSAlibWrapper::lstParams)
    .def_readwrite("time_start", &pSSAlibWrapper::dTimeStart)
    .def_readwrite("time_step", &pSSAlibWrapper::dTimeStep)
    .def_readwrite("time_end", &pSSAlibWrapper::dTimeEnd)
    .add_property("test_case", &pSSAlibWrapper::getTestCase, &pSSAlibWrapper::setTestCase)
    .add_property("result", &pSSAlibWrapper::getResult)
  ;
};
