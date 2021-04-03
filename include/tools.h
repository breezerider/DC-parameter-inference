#include <pssalib/datamodel/detail/Model.h>
#include <pssalib/datamodel/detail/Species.h>
#include <pssalib/datamodel/detail/Reaction.h>
#include <pssalib/datamodel/detail/SpeciesReference.h>

#ifndef PSSA_TOOLS_H_
#define PSSA_TOOLS_H_

#ifdef BOOST_PY_ERRMSG

/*
 * Test computeS, computeF, computeJ and computeQ:
 * eng = pyPSSAlib.pSSAlib()
 * eng.test_case = 'ca'
 * eng.params = [2, 2.1, 0.1, 1.0, 0.01, 0.1, 15, 0, 0]
 * eng.time_start = 3000.0
 * eng.time_step = 0.1
 * eng.time_end = 4000.0
 * 
 * # steady state solution
 * k = eng.params
 * S_ss = [0, 0];
 * S_ss[0] = (np.sqrt(k[4]**2 + 8 * k[1] * k[2] * k[5] / (k[3] + k[5])) - k[4])*(k[3] + k[5])/(4.0 * k[2] * k[5])
 * S_ss[1] = k[2] * S_ss[0]**2 / (k[3] + k[5])
 * print(f's_ss = {S_ss}')
 * 
 * # Jacobian
 * J = [[-4.0*k[2]*S_ss[0]-k[4], 2.0*k[3]], [2.0*k[2]*S_ss[0], -k[5]-k[3]]]
 * print(f'J_anl = {J}')
 * 
 * # Rate functions
 * F = [k[1], k[2]*S_ss[0]**2, k[3]*S_ss[1], k[4]*S_ss[0], k[5]*S_ss[1]]
 * 
 * # Steady state molecule #
 * S_ss = [s * k[6] for s in S_ss]
 * print(f'S_ss = {S_ss}')
 * 
 * J = eng.jacobian(S_ss)
 * print(f'J_comp = {J}')
 * # Stoichiometry matrix
 * S = eng.stoichiometry(S_ss)
 * print(f'S = {S}')
 * 
 * # Lyapunov Q
 * print(f'Q_anl = {k[6] * S.T.dot(np.diag(F)).dot(S)}')
 * print(f'Q_comp = {eng.lyapunovQ(S_ss)}')
 */

/*
 * Compute the stoichiometry matrix
 */
void
computeS(
  pssalib::datamodel::detail::Model & model,
  const double * pardPopulation,
  double * pardOutput
)
{
  for(UINTEGER ri = 0; ri < model.getReactionsCount(); ++ri)
  {
    pssalib::datamodel::detail::Reaction * r = model.getReaction(ri);

    for(UINTEGER sri = 0; sri < r->getSpeciesReferencesCount(); ++sri)
    {
      const pssalib::datamodel::detail::SpeciesReference * sr = r->getSpeciesReferenceAt(sri);

      // fill the stoichiometry matrix S
      if(!sr->isReservoir())
        pardOutput[ri*model.getSpeciesCount()+sr->getIndex()] = (sri < r->getReactantsCount()) ? double(-1.0 * sr->getStoichiometryAbs()) : double(sr->getStoichiometryAbs());
    }
  }
}

/*
 * Compute the reaction rate for mass-action kinetics
 */
void
computeF(
  const pssalib::datamodel::detail::Model & model,
  const double * pardPopulation,
  double * pardF
)
{
  double Vinv = model.getCompartmentVolume();
  if(Vinv <= 0.0)
    Vinv = 1.0;
  else
    Vinv = 1.0 / Vinv;

  for(UINTEGER ri = 0; ri < model.getReactionsCount(); ++ri)
  {
    const pssalib::datamodel::detail::Reaction * r = model.getReaction(ri);

    pardF[ri] = r->getForwardRate(); // fill the reaction rates vector F

    for(UINTEGER sri = 0; sri < r->getReactantsCount(); ++sri)
    {
      const pssalib::datamodel::detail::SpeciesReference * sr = r->getSpeciesReferenceAt(sri);

      // fill the reaction rates vector F
      if(!sr->isReservoir() && (sri < r->getReactantsCount()))
        pardF[ri] *= std::pow(pardPopulation[sr->getIndex()] * Vinv, sr->getStoichiometryAbs());
    }
  }
}

/*
 * Compute the value of mass-action system of ODEs at a given point
 */
void
computeODEs(
  const pssalib::datamodel::detail::Model & model,
  const double * pardF,
  double * pardODEs
)
{
  for(UINTEGER ri = 0; ri < model.getReactionsCount(); ++ri)
  {
    const pssalib::datamodel::detail::Reaction * r = model.getReaction(ri);

    for(UINTEGER sri = 0; sri < r->getSpeciesReferencesCount(); ++sri)
    {
      const pssalib::datamodel::detail::SpeciesReference * sr = r->getSpeciesReferenceAt(sri);

      // fill the ODEs vector
      if(!sr->isReservoir())
      {
        if(sri < r->getReactantsCount())
          pardODEs[sr->getIndex()] -= sr->getStoichiometryAbs() * pardF[ri];
        else
          pardODEs[sr->getIndex()] += sr->getStoichiometryAbs() * pardF[ri];
      }
    }
  }
}

/*
 * Compute the Jacobian of mass-action system of ODEs at a given point
 */
void
computeJ(
  pssalib::datamodel::detail::Model & model,
  const double * pardPopulation,
  double * pardJ
)
{
  double Vinv = model.getCompartmentVolume();
  if(Vinv <= 0.0)
    Vinv = 1.0;
  else
    Vinv = 1.0 / Vinv;

  for(UINTEGER ri = 0; ri < model.getReactionsCount(); ++ri)
  {
    const pssalib::datamodel::detail::Reaction * r = model.getReaction(ri);

    for(UINTEGER sri = 0; sri < r->getReactantsCount(); ++sri)
    {
      const pssalib::datamodel::detail::SpeciesReference * srI = r->getSpeciesReferenceAt(sri);

      double jri = r->getForwardRate();

      if(srI->isReservoir()) break;

      for(UINTEGER srj = 0; srj < r->getReactantsCount(); ++srj)
      {
        const pssalib::datamodel::detail::SpeciesReference * srJ = r->getSpeciesReferenceAt(srj);

        if(srI->getIndex() != srJ->getIndex())
          jri *= std::pow(pardPopulation[srJ->getIndex()] * Vinv, srJ->getStoichiometryAbs());
        else if(srJ->getStoichiometryAbs() > 1)
          jri *= double(srJ->getStoichiometryAbs()) * std::pow(pardPopulation[srJ->getIndex()] * Vinv, srJ->getStoichiometryAbs() - 1);
      }

      for(UINTEGER srj = 0; srj < r->getSpeciesReferencesCount(); ++srj)
      {
        const pssalib::datamodel::detail::SpeciesReference * srJ = r->getSpeciesReferenceAt(srj);

        if(srJ->isReservoir()) continue;

        if(srj < r->getReactantsCount())
          pardJ[srI->getIndex() + srJ->getIndex() * model.getSpeciesCount()] -= jri * double(srJ->getStoichiometryAbs());
        else
          pardJ[srI->getIndex() + srJ->getIndex() * model.getSpeciesCount()] += jri * double(srJ->getStoichiometryAbs());
      }
    }
  }
}

/*
 * Compute the Q-term for Lyapunov Equation
 */
void
computeQ(
  pssalib::datamodel::detail::Model & model,
  const double * pardPopulation,
  double * pardQ
)
{
  boost::scoped_array<double> ardF(new double[model.getReactionsCount()]);
  computeF(model, pardPopulation, ardF.get());

  for(UINTEGER ri = 0; ri < model.getReactionsCount(); ++ri)
  {
    pssalib::datamodel::detail::Reaction * r = model.getReaction(ri);

    // Scale with compartment volume
    ardF[ri] *= model.getCompartmentVolume();

    for(UINTEGER sri = 0; sri < r->getSpeciesReferencesCount(); ++sri)
    {
      const pssalib::datamodel::detail::SpeciesReference * srI = r->getSpeciesReferenceAt(sri);

      for(UINTEGER srj = 0; srj < r->getSpeciesReferencesCount(); ++srj)
      {
        const pssalib::datamodel::detail::SpeciesReference * srJ = r->getSpeciesReferenceAt(srj);

        // fill the stoichiometry matrix S
        if(!srI->isReservoir()&&!srJ->isReservoir())
        {
          double dq = ardF[ri] * srI->getStoichiometryAbs() * srJ->getStoichiometryAbs();
          if((sri < r->getReactantsCount())^(srj < r->getReactantsCount()))
            pardQ[srI->getIndex() + srJ->getIndex() * model.getSpeciesCount()] -= dq;
          else
            pardQ[srI->getIndex() + srJ->getIndex() * model.getSpeciesCount()] += dq;
        }
      }
    }
  }
}

/*
 * Compute given statistical moments of a data array
 */
template<typename T>
void computeMoments(T *arData, size_t szDataRows, size_t szDataCols, unsigned int* arOrder, size_t szOrder, double * arOutput)
{
  for(size_t orderIdx = 0; orderIdx < szOrder; ++orderIdx)
  {
    size_t order = arOrder[orderIdx];

    for(size_t colIdx = 0; colIdx < szDataCols; ++colIdx)
    {
      double mu = 0.0, sigma = 0.0, mu_i = 0.0;

      if(0 == order)
      {
        arOutput[orderIdx*szDataCols+colIdx] = 0.0;
        continue;
      }

      for(size_t rowIdx = 0; rowIdx < szDataRows; ++rowIdx)
        mu += (double)arData[rowIdx*szDataCols+colIdx];
      mu /= (double)szDataRows;

      if(1 < order)
      {
        for(size_t rowIdx = 0; rowIdx < szDataRows; ++rowIdx)
          sigma += std::pow((double)arData[rowIdx*szDataCols+colIdx] - mu, 2.0);
        sigma = std::sqrt(sigma / (double)szDataRows);

        if(2 < order)
        {
          if(sigma > 0.0)
          {
            for(size_t rowIdx = 0; rowIdx < szDataRows; ++rowIdx)
              mu_i += std::pow((double)arData[rowIdx*szDataCols+colIdx] - mu, (double)order);
            mu_i = std::pow(std::abs(mu_i / (double)szDataRows), 1.0 / ((double)order)) / sigma;
          }
          else
            mu_i = 0.0;
        }
        else
          mu_i = sigma;
      }
      else
        mu_i = mu;

      arOutput[orderIdx*szDataCols+colIdx] = mu_i;
    }
  }
}

#else /* MATLAB version */

void
generateS(
  pssalib::datamodel::detail::Model & model,
  const double * pardPopulation,
  mxArray ** parS
)
{
  if(!mxIsSparse(*parS)) return;

  double  * prs = mxGetPr(*parS);
  mwIndex * irs = mxGetIr(*parS);
  mwIndex * jcs = mxGetJc(*parS);

  mwIndex nnzi = 0;
  for(UINTEGER ri = 0; ri < model.getReactionsCount(); ++ri)
  {
    pssalib::datamodel::detail::Reaction * r = model.getReaction(ri);

    jcs[ri] = nnzi;
    for(UINTEGER sri = 0; sri < r->getSpeciesReferencesCount(); ++sri)
    {
      const pssalib::datamodel::detail::SpeciesReference * sr = r->getSpeciesReferenceAt(sri);

      // fill the stoichiometry matrix S
      if(!sr->isReservoir())
      {
        irs[nnzi] = sr->getIndex();
        prs[nnzi] = sr->getStoichiometryAbs();
        if(sri < r->getReactantsCount())
          prs[nnzi] *= -1;
        ++nnzi;
      }
    }
  }
  jcs[model.getReactionsCount()] = nnzi;
}

void
computeF(
  const pssalib::datamodel::detail::Model & model,
  const double * pardPopulation,
  double * pardF
)
{
  double Vinv = model.getCompartmentVolume();
  if(Vinv <= 0.0)
    Vinv = 1.0;
  else
    Vinv = 1.0 / Vinv;

  for(UINTEGER ri = 0; ri < model.getReactionsCount(); ++ri)
  {
    const pssalib::datamodel::detail::Reaction * r = model.getReaction(ri);

    pardF[ri] = r->getForwardRate(); // fill the reaction rates vector F

    for(UINTEGER sri = 0; sri < r->getReactantsCount(); ++sri)
    {
      const pssalib::datamodel::detail::SpeciesReference * sr = r->getSpeciesReferenceAt(sri);

      // fill the reaction rates vector F
      if(!sr->isReservoir() && (sri < r->getReactantsCount()))
        pardF[ri] *= std::pow(pardPopulation[sr->getIndex()] * Vinv, sr->getStoichiometryAbs());
    }
  }
}

void
computeODEs(
  const pssalib::datamodel::detail::Model & model,
  const double * pardF,
  double * pardODEs
)
{
  for(UINTEGER ri = 0; ri < model.getReactionsCount(); ++ri)
  {
    const pssalib::datamodel::detail::Reaction * r = model.getReaction(ri);

    for(UINTEGER sri = 0; sri < r->getSpeciesReferencesCount(); ++sri)
    {
      const pssalib::datamodel::detail::SpeciesReference * sr = r->getSpeciesReferenceAt(sri);

      // fill the ODEs vector
      if(!sr->isReservoir())
      {
        if(sri < r->getReactantsCount())
          pardODEs[sr->getIndex()] -= sr->getStoichiometryAbs() * pardF[ri];
        else
          pardODEs[sr->getIndex()] += sr->getStoichiometryAbs() * pardF[ri];
      }
    }
  }
}

void
computeJ(
  pssalib::datamodel::detail::Model & model,
  const double * pardPopulation,
  double * pardJ
)
{
  double Vinv = model.getCompartmentVolume();
  if(Vinv <= 0.0)
    Vinv = 1.0;
  else
    Vinv = 1.0 / Vinv;

  for(UINTEGER ri = 0; ri < model.getReactionsCount(); ++ri)
  {
    const pssalib::datamodel::detail::Reaction * r = model.getReaction(ri);

    for(UINTEGER sri = 0; sri < r->getReactantsCount(); ++sri)
    {
      const pssalib::datamodel::detail::SpeciesReference * srI = r->getSpeciesReferenceAt(sri);

      double jri = r->getForwardRate();

      if(srI->isReservoir()) break;

      for(UINTEGER srj = 0; srj < r->getReactantsCount(); ++srj)
      {
        const pssalib::datamodel::detail::SpeciesReference * srJ = r->getSpeciesReferenceAt(srj);

        if(srI->getIndex() != srJ->getIndex())
          jri *= std::pow(pardPopulation[srJ->getIndex()] * Vinv, srJ->getStoichiometryAbs());
        else if(srJ->getStoichiometryAbs() > 1)
          jri *= double(srJ->getStoichiometryAbs()) * std::pow(pardPopulation[srJ->getIndex()] * Vinv, srJ->getStoichiometryAbs() - 1);
      }

      for(UINTEGER srj = 0; srj < r->getSpeciesReferencesCount(); ++srj)
      {
        const pssalib::datamodel::detail::SpeciesReference * srJ = r->getSpeciesReferenceAt(srj);

        if(srJ->isReservoir()) continue;

        if(srj < r->getReactantsCount())
          pardJ[srJ->getIndex() + srI->getIndex() * model.getSpeciesCount()] -= jri * double(srJ->getStoichiometryAbs());
        else
          pardJ[srJ->getIndex() + srI->getIndex() * model.getSpeciesCount()] += jri * double(srJ->getStoichiometryAbs());
      }
    }
  }
}

#endif /* BOOST_PY_ERRMSG */

#endif /* PSSA_TOOLS_H_ */
