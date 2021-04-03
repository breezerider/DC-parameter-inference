#include <pssalib/datamodel/detail/Model.h>
#include <pssalib/datamodel/detail/Species.h>
#include <pssalib/datamodel/detail/Reaction.h>
#include <pssalib/datamodel/detail/SpeciesReference.h>

#ifndef PSSA_MODELS_H_
#define PSSA_MODELS_H_

#ifdef BOOST_PY_ERRMSG
#define MODEL_ERROR(txt) BOOST_PY_ERRMSG(txt)
#else
#define MODEL_ERROR(txt) mexErrMsgIdAndTxt("MATLAB:mexpssa:internal", txt)
#endif

void
generateCyclicLinearChain(
    pssalib::datamodel::detail::Model & model,
    double * pdParams,
    size_t   szParams
)
{
  size_t szPopulation = pdParams[0];
  if((1 + szPopulation + 1 + szPopulation) != szParams) // system size, reactions rates, reactor volume and initial species population
    MODEL_ERROR("Invalid parameters vector for the CyclicLinearChain test case, must contain system size, reactions rates, reactor volume and initial species population.");
  double * pdPopulation = pdParams + szPopulation + 1;

  // clear any previous model definitions
  model.free();

  // initialize the model
  model.setId("CyclicLinearChain");
  model.setCompartmentVolume(pdParams[szPopulation]);
  model.allocSpecies(szPopulation);
  model.allocReactions(szPopulation);

  // Species
  BOOSTFORMAT fmtSpeciesName("S%d");
  for(size_t s = 0; s < szPopulation; s++)
  {
    pssalib::datamodel::detail::Species * species =
      model.getSpecies(s);
    species->setId((fmtSpeciesName % s).str());
    species->setInitialAmount(std::floor(pdPopulation[s]));
  }

  // Reactions
  BOOSTFORMAT fmtReactionName("R%d");
  pssalib::datamodel::detail::Reaction * reaction = NULL;
  pssalib::datamodel::detail::SpeciesReference * spRef = NULL;
  for(size_t r = 0; r < szPopulation; r++)
  {
    reaction = model.getReaction(r);
    reaction->setId((fmtReactionName % r).str());
    reaction->setReversible(false);

    // rate constant
    reaction->setForwardRate(pdParams[r]);

    // species refs
    reaction->allocSpeciesRefs(1, 1);

    // reactant
    spRef = reaction->getReactantsListAt(0);
    spRef->setStoichiometry(1);
    spRef->setIndex(r);

    // product
    spRef = reaction->getProductsListAt(0);
    spRef->setStoichiometry(1);
    if(r != (szPopulation - 1))
      spRef->setIndex(r + 1); // propagate down the chain
    else
      spRef->setIndex(0); // close the loop
  }
}

void
generateColloidalAggregation(
    pssalib::datamodel::detail::Model & model,
    double * pdParams,
    size_t   szParams
)
{
  size_t szPopulation = pdParams[0];
  if((1 + szPopulation*szPopulation/2 + szPopulation + 2 + szPopulation) != szParams) // system size, reactions rates, reactor volume and initial species population
    MODEL_ERROR("Invalid parameters vector for the ColloidalAggregation test case, must contain exactly (1 + (1 + floor(n^2/2) + n) + 1 + n) elements, where n is the number of species in the system. It must contain the system size, reactions rates, reactor volume and initial species population");
  double * pdPopulation = pdParams + (szPopulation*szPopulation/2 + szPopulation + 2) + 1;

  // clear any previous model definitions
  model.free();

  // initialize the model
  model.setId("ColloidalAggregation");
  model.setCompartmentVolume(pdParams[szPopulation*szPopulation/2 + szPopulation + 1 + 1]);
  model.allocSpecies(szPopulation);
  model.allocReactions(szPopulation*szPopulation/2 + szPopulation + 1);

  // Species
  BOOSTFORMAT fmtSpeciesName("S%d");
  for(size_t s = 0; s < szPopulation; ++s)
  {
    pssalib::datamodel::detail::Species * species =
      model.getSpecies(s);
    species->setId((fmtSpeciesName % s).str());
    species->setInitialAmount(std::floor(pdPopulation[s]));
  }

  // Reactions
  size_t r = 0;
  BOOSTFORMAT fmtReactionName("R%d");
  pssalib::datamodel::detail::Reaction * reaction = NULL;
  pssalib::datamodel::detail::SpeciesReference * spRef = NULL;

  ////////////////////////
  // 0 -> S_1
  reaction = model.getReaction(r);
  reaction->setId((fmtReactionName % r).str());
  reaction->setReversible(false);

  // rate constant
  reaction->setForwardRate(pdParams[r+1]);

  // species refs
  reaction->allocSpeciesRefs(1, 1);

  // reactant
  spRef = reaction->getReactantsListAt(0);
  spRef->makeReservoir();

  // product
  spRef = reaction->getProductsListAt(0);
  spRef->setStoichiometry(1);
  spRef->setIndex(0);

  r++;

  ////////////////////////
  // S_n + S_m --> S_(n+m)
  for(size_t s1 = 0; s1 < (szPopulation / 2); s1++)
  {
    for(size_t s2 = s1; s2 < (szPopulation - s1 - 1); s2++, r++)
    {
      reaction = model.getReaction(r);
      reaction->setId((fmtReactionName % r).str());
      reaction->setReversible(false);

      // rate constant
      reaction->setForwardRate(pdParams[r+1]);

      // species refs
      reaction->allocSpeciesRefs(2, 1);

      // reactant 1
      spRef = reaction->getReactantsListAt(0);
      spRef->setStoichiometry(1);
      spRef->setIndex(s1);

      // reactant 2
      spRef = reaction->getReactantsListAt(1);
      spRef->setStoichiometry(1);
      spRef->setIndex(s2);

      // product
      spRef = reaction->getProductsListAt(0);
      spRef->setStoichiometry(1);
      spRef->setIndex(s1 + s2 + 1);
    }
  }
  ////////////////////////
  // S_p --> S_q + S_(p-q)
  for(size_t s1 = 1; s1 < szPopulation; s1++)
  {
    for(size_t s2 = 0; s2 < ((s1 - 1) / 2 + 1); s2++, r++)
    {
      reaction = model.getReaction(r);
      reaction->setId((fmtReactionName % r).str());
      reaction->setReversible(false);

      // rate constant
      reaction->setForwardRate(pdParams[r+1]);

      // species refs
      reaction->allocSpeciesRefs(1, 2);

      // reactant
      spRef = reaction->getReactantsListAt(0);
      spRef->setStoichiometry(1);
      spRef->setIndex(s1);

      // product 1
      spRef = reaction->getProductsListAt(0);
      spRef->setStoichiometry(1);
      spRef->setIndex(s2);

      // product 2
      spRef = reaction->getProductsListAt(1);
      spRef->setStoichiometry(1);
      spRef->setIndex(s1 - s2 - 1);
    }
  }

  ////////////////////////
  // S_i -> 0
  for(size_t s = 0; s < szPopulation; ++s, ++r)
  {
    reaction = model.getReaction(r);
    reaction->setId((fmtReactionName % r).str());
    reaction->setReversible(false);

    // rate constant
    reaction->setForwardRate(pdParams[r+1]);

    // species refs
    reaction->allocSpeciesRefs(1, 1);

    // reactant
    spRef = reaction->getReactantsListAt(0);
    spRef->setStoichiometry(1);
    spRef->setIndex(s);

    // product
    spRef = reaction->getProductsListAt(0);
    spRef->makeReservoir();
  }
}

void
generateSingleBirthDeath(
    pssalib::datamodel::detail::Model & model,
    double * pdParams,
    size_t   szParams
)
{
  if(4 != szParams)
    MODEL_ERROR("Invalid parameters vector or species number for the SingleBirthDeath test case. One chemical specie and three parameters (two reaction rates and the reactor volume) are expected.");
  double * pdPopulation = pdParams + 3;

  // clear any previous model definitions
  model.free();

  // initialize the model
  model.setId("SingleBirthDeath");
  model.setCompartmentVolume(pdParams[2]);
  model.allocSpecies(1);
  model.allocReactions(2);

  ////////////////////////
  // Species
  pssalib::datamodel::detail::Species * species = NULL;

  // Species S
  species = model.getSpecies(0);
  species->setId(STRING("S"));
  species->setInitialAmount(std::floor(pdPopulation[0]));

  ////////////////////////
  // Reactions
  pssalib::datamodel::detail::Reaction * reaction = NULL;
  pssalib::datamodel::detail::SpeciesReference * spRef = NULL;

  ////////////////////////
  // 0 -> S
  reaction = model.getReaction(0);
  reaction->setId("S_Generation");
  reaction->setReversible(false);

  // rate constant
  reaction->setForwardRate(pdParams[0]);

  // species refs
  reaction->allocSpeciesRefs(1, 1);

  // reactant
  spRef = reaction->getReactantsListAt(0);
  spRef->makeReservoir();

  // product
  spRef = reaction->getProductsListAt(0);
  spRef->setStoichiometry(1);
  spRef->setIndex(0);
  ////////////////////////
  // S -> 0
  reaction = model.getReaction(1);
  reaction->setId("S_Degradation");
  reaction->setReversible(false);

  // rate constant
  reaction->setForwardRate(pdParams[1]);

  // species refs
  reaction->allocSpeciesRefs(1, 1);

  // reactant
  spRef = reaction->getReactantsListAt(0);
  spRef->setStoichiometry(1);
  spRef->setIndex(0);

  // product
  spRef = reaction->getProductsListAt(0);
  spRef->makeReservoir();
}

void
generateHomoreaction(
    pssalib::datamodel::detail::Model & model,
    double * pdParams,
    size_t   szParams
)
{
  if(4 != szParams)
    MODEL_ERROR("Invalid parameters vector or species number for the Homoreaction test case. One specie number and three parameters (two reaction rates and the reactor volume) are expected.");
  double * pdPopulation = pdParams + 3;

  // clear any previous model definitions
  model.free();

  // initialize the model
  model.setId("Homoreaction");
  model.setCompartmentVolume(pdParams[2]);
  model.allocSpecies(1);
  model.allocReactions(2);

  ////////////////////////
  // Species
  pssalib::datamodel::detail::Species * species = NULL;

  // Species A
  species = model.getSpecies(0);
  species->setId(STRING("A"));
  species->setInitialAmount(std::floor(pdPopulation[0]));

  ////////////////////////
  // Reactions
  pssalib::datamodel::detail::Reaction * reaction = NULL;
  pssalib::datamodel::detail::SpeciesReference * spRef = NULL;

  ////////////////////////
  // 0 -> A
  reaction = model.getReaction(0);
  reaction->setId("A_Generation");
  reaction->setReversible(false);

  // rate constant
  reaction->setForwardRate(pdParams[0]);

  // species refs
  reaction->allocSpeciesRefs(1, 1);

  // reactant
  spRef = reaction->getReactantsListAt(0);
  spRef->makeReservoir();

  // product
  spRef = reaction->getProductsListAt(0);
  spRef->setStoichiometry(1);
  spRef->setIndex(0);
  ////////////////////////
  // A + A -> 0
  reaction = model.getReaction(1);
  reaction->setId("A_Dimerization");
  reaction->setReversible(false);

  // rate constant
  reaction->setForwardRate(pdParams[1]);

  // species refs
  reaction->allocSpeciesRefs(1, 1);

  // reactant
  spRef = reaction->getReactantsListAt(0);
  spRef->setStoichiometry(2);
  spRef->setIndex(0);

  // product
  spRef = reaction->getProductsListAt(0);
  spRef->makeReservoir();
}

void
generateTwoComponentSystem(
    pssalib::datamodel::detail::Model & model,
    double * pdParams,
    size_t   szParams
)
{
  if(11 != szParams)
    MODEL_ERROR("Invalid parameters vector or species number for the Two Component System test case. Five specie numbers and six parameters (five reaction rates and the reactor volume) are expected.");
  double * pdPopulation = pdParams + 6;

  // clear any previous model definitions
  model.free();

  // initialize the model
  model.setId("TwoComponentSystem");
  model.setCompartmentVolume(pdParams[5]);
  model.allocSpecies(5);
  model.allocReactions(5);

  ////////////////////////
  // Species
  pssalib::datamodel::detail::Species * species = NULL;

  // Species I - initiator
  species = model.getSpecies(0);
  species->setId(STRING("I"));
  species->setInitialAmount(std::floor(pdPopulation[0]));

  // Species HK - histidine kinase
  species = model.getSpecies(1);
  species->setId(STRING("HK"));
  species->setInitialAmount(std::floor(pdPopulation[1]));

  // Species HK - phosphorylated histidine kinase
  species = model.getSpecies(2);
  species->setId(STRING("HKp"));
  species->setInitialAmount(std::floor(pdPopulation[2]));

  // Species RR - response regulator
  species = model.getSpecies(3);
  species->setId(STRING("RR"));
  species->setInitialAmount(std::floor(pdPopulation[3]));

  // Species RR - phosphorylated response regulator
  species = model.getSpecies(4);
  species->setId(STRING("RRp"));
  species->setInitialAmount(std::floor(pdPopulation[4]));

  ////////////////////////
  // Reactions
  pssalib::datamodel::detail::Reaction * reaction = NULL;
  pssalib::datamodel::detail::SpeciesReference * spRef = NULL;

  ////////////////////////
  // I + HK -k1-> I + HKp
  reaction = model.getReaction(0);
  reaction->setId("HK_Phosporylation");
  reaction->setReversible(false);

  // rate constant
  reaction->setForwardRate(pdParams[0]);

  // species refs
  reaction->allocSpeciesRefs(2, 2);

  // reactant 1
  spRef = reaction->getReactantsListAt(0);
  spRef->setStoichiometry(1);
  spRef->setIndex(0);

  // reactant 2
  spRef = reaction->getReactantsListAt(1);
  spRef->setStoichiometry(1);
  spRef->setIndex(1);

  // product 1
  spRef = reaction->getProductsListAt(0);
  spRef->setStoichiometry(1);
  spRef->setIndex(0);

  // product 2
  spRef = reaction->getProductsListAt(1);
  spRef->setStoichiometry(1);
  spRef->setIndex(2);
  ////////////////////////
  // HKp + RR -k2-> HK + RRp
  reaction = model.getReaction(1);
  reaction->setId("RR_Phosporylation");
  reaction->setReversible(false);

  // rate constant
  reaction->setForwardRate(pdParams[1]);

  // species refs
  reaction->allocSpeciesRefs(2, 2);

  // reactant 1
  spRef = reaction->getReactantsListAt(0);
  spRef->setStoichiometry(1);
  spRef->setIndex(2);

  // reactant 2
  spRef = reaction->getReactantsListAt(1);
  spRef->setStoichiometry(1);
  spRef->setIndex(3);

  // product 1
  spRef = reaction->getProductsListAt(0);
  spRef->setStoichiometry(1);
  spRef->setIndex(1);

  // product 2
  spRef = reaction->getProductsListAt(1);
  spRef->setStoichiometry(1);
  spRef->setIndex(4);
  ////////////////////////
  // HKp -k3-> HK
  reaction = model.getReaction(2);
  reaction->setId("HK_Autodephosporylation");
  reaction->setReversible(false);

  // rate constant
  reaction->setForwardRate(pdParams[2]);

  // species refs
  reaction->allocSpeciesRefs(1, 1);

  // reactant
  spRef = reaction->getReactantsListAt(0);
  spRef->setStoichiometry(1);
  spRef->setIndex(2);

  // product
  spRef = reaction->getProductsListAt(0);
  spRef->setStoichiometry(1);
  spRef->setIndex(1);
  ////////////////////////
  // HK + RRp -k4-> HKp + RR
  reaction = model.getReaction(3);
  reaction->setId("RR_Dephosporylation");
  reaction->setReversible(false);

  // rate constant
  reaction->setForwardRate(pdParams[3]);

  // species refs
  reaction->allocSpeciesRefs(2, 2);

  // reactant 1
  spRef = reaction->getReactantsListAt(0);
  spRef->setStoichiometry(1);
  spRef->setIndex(1);

  // reactant 2
  spRef = reaction->getReactantsListAt(1);
  spRef->setStoichiometry(1);
  spRef->setIndex(4);

  // product 1
  spRef = reaction->getProductsListAt(0);
  spRef->setStoichiometry(1);
  spRef->setIndex(2);

  // product 2
  spRef = reaction->getProductsListAt(1);
  spRef->setStoichiometry(1);
  spRef->setIndex(3);
  ////////////////////////
  // RRp -k5-> RR
  reaction = model.getReaction(4);
  reaction->setId("RR_Autodephosporylation");
  reaction->setReversible(false);

  // rate constant
  reaction->setForwardRate(pdParams[4]);

  // species refs
  reaction->allocSpeciesRefs(1, 1);

  // reactant
  spRef = reaction->getReactantsListAt(0);
  spRef->setStoichiometry(1);
  spRef->setIndex(4);

  // product
  spRef = reaction->getProductsListAt(0);
  spRef->setStoichiometry(1);
  spRef->setIndex(3);
}

void
generateEnzymaticDegradation(
    pssalib::datamodel::detail::Model & model,
    double * pdParams,
    size_t   szParams
)
{
  if(11 != szParams)
    MODEL_ERROR("Invalid parameters vector or species number for the Enzymatic Degration test case. Seven parameters (six reaction rates and the reactor volume) and initial populations of four species are expected.");
  double * pdPopulation = pdParams + 7;

  // clear any previous model definitions
  model.free();

  // initialize the model
  model.setId("EnzymaticDegration");
  model.setCompartmentVolume(pdParams[6]);
  model.allocSpecies(4);
  model.allocReactions(6);

  ////////////////////////
  // Species
  pssalib::datamodel::detail::Species * species = NULL;

  // Species mRNA
  species = model.getSpecies(0);
  species->setId(STRING("mRNA"));
  species->setInitialAmount(std::floor(pdPopulation[0]));

  // Species Complex = Enzyme + Protein
  species = model.getSpecies(1);
  species->setId(STRING("Complex"));
  species->setInitialAmount(std::floor(pdPopulation[1]));

  // Species Enzyme
  species = model.getSpecies(2);
  species->setId(STRING("Enzyme"));
  species->setInitialAmount(std::floor(pdPopulation[2]));

  // Species Protein
  species = model.getSpecies(3);
  species->setId(STRING("Protein"));
  species->setInitialAmount(std::floor(pdPopulation[3]));

  ////////////////////////
  // Reactions
  pssalib::datamodel::detail::Reaction * reaction = NULL;
  pssalib::datamodel::detail::SpeciesReference * spRef = NULL;

  ////////////////////////
  // 0 -k0-> mRNA
  reaction = model.getReaction(0);
  reaction->setId("mRNA_Production");
  reaction->setReversible(false);

  // rate constant
  reaction->setForwardRate(pdParams[0]);

  // species refs
  reaction->allocSpeciesRefs(1, 1);

  // reactant
  spRef = reaction->getReactantsListAt(0);
  spRef->makeReservoir();

  // product
  spRef = reaction->getProductsListAt(0);
  spRef->setStoichiometry(1);
  spRef->setIndex(0);

  ////////////////////////
  // mRNA -k1-> mRNA + Protein
  reaction = model.getReaction(1);
  reaction->setId("Protein_Translation");
  reaction->setReversible(false);

  // rate constant
  reaction->setForwardRate(pdParams[1]);

  // species refs
  reaction->allocSpeciesRefs(1, 2);

  // reactant
  spRef = reaction->getReactantsListAt(0);
  spRef->setStoichiometry(1);
  spRef->setIndex(0);

  // product 1
  spRef = reaction->getProductsListAt(0);
  spRef->setStoichiometry(1);
  spRef->setIndex(0);

  // product 2
  spRef = reaction->getProductsListAt(1);
  spRef->setStoichiometry(1);
  spRef->setIndex(3);

  ////////////////////////
  // Protein + Enzyme -k2-> Complex
  reaction = model.getReaction(2);
  reaction->setId("Complex_Formation");
  reaction->setReversible(false);

  // rate constant
  reaction->setForwardRate(pdParams[2]);

  // species refs
  reaction->allocSpeciesRefs(2, 1);

  // reactant 1
  spRef = reaction->getReactantsListAt(0);
  spRef->setStoichiometry(1);
  spRef->setIndex(2);

  // reactant 2
  spRef = reaction->getReactantsListAt(1);
  spRef->setStoichiometry(1);
  spRef->setIndex(3);

  // product
  spRef = reaction->getProductsListAt(0);
  spRef->setStoichiometry(1);
  spRef->setIndex(1);

  ////////////////////////
  // Complex -k3-> Protein + Enzyme
  reaction = model.getReaction(3);
  reaction->setId("Complex_Detachment");
  reaction->setReversible(false);

  // rate constant
  reaction->setForwardRate(pdParams[3]);

  // species refs
  reaction->allocSpeciesRefs(1, 2);

  // reactant
  spRef = reaction->getReactantsListAt(0);
  spRef->setStoichiometry(1);
  spRef->setIndex(1);

  // product 1
  spRef = reaction->getProductsListAt(0);
  spRef->setStoichiometry(1);
  spRef->setIndex(2);

  // product 2
  spRef = reaction->getProductsListAt(1);
  spRef->setStoichiometry(1);
  spRef->setIndex(3);

  ////////////////////////
  // Complex -k4-> Enzyme
  reaction = model.getReaction(4);
  reaction->setId("Complex_Degraded");
  reaction->setReversible(false);

  // rate constant
  reaction->setForwardRate(pdParams[4]);

  // species refs
  reaction->allocSpeciesRefs(1, 1);

  // reactant
  spRef = reaction->getReactantsListAt(0);
  spRef->setStoichiometry(1);
  spRef->setIndex(3);

  // product
  spRef = reaction->getProductsListAt(0);
  spRef->setStoichiometry(1);
  spRef->setIndex(2);

  ////////////////////////
  // mRNA -k4-> 0
  reaction = model.getReaction(5);
  reaction->setId("mRNA_Degradation");
  reaction->setReversible(false);

  // rate constant
  reaction->setForwardRate(pdParams[5]);

  // species refs
  reaction->allocSpeciesRefs(1, 1);

  // reactant
  spRef = reaction->getReactantsListAt(0);
  spRef->setStoichiometry(1);
  spRef->setIndex(0);

  // product
  spRef = reaction->getProductsListAt(0);
  spRef->makeReservoir();
}

#endif /* PSSA_MODELS_H_ */
