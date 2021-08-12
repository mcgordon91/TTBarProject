#ifndef MUON_MATCHER
#define MUON_MATCHER

#include <iostream>
#include <string>
#include <vector>

typedef ROOT::VecOps::RVec<Float_t>                        RVec_f;
typedef ROOT::VecOps::RVec<Float_t>::const_iterator        RVec_f_iter;
typedef ROOT::VecOps::RVec<Int_t>                          RVec_i;
typedef ROOT::VecOps::RVec<Int_t>::const_iterator          RVec_i_iter;
typedef ROOT::VecOps::RVec<std::string>                    RVec_str;
typedef ROOT::VecOps::RVec<std::string>::const_iterator    RVec_str_iter;

using namespace std;

namespace MuonMatcher
{
    class MatchOppositelyChargedMuons
    {
        private:
            RVec_f Muon_pt;
            RVec_f Muon_eta;
            RVec_f Muon_phi;
            RVec_f Muon_mass;
            RVec_i Muon_charge;


        public:
            MatchOppositelyChargedMuons(RVec_f Muon_pt, RVec_f Muon_eta, RVec_f Muon_phi, RVec_f Muon_mass, RVec_i Muon_charge);
            int Test();
            RVec_f MuonSorter();
            RVec_f InvariantMassCalculator();
    };

    MatchOppositelyChargedMuons::MatchOppositelyChargedMuons(RVec_f Muon_pt, RVec_f Muon_eta, RVec_f Muon_phi, RVec_f Muon_mass, RVec_i Muon_charge)
    {
        this->Muon_pt = Muon_pt;
        this->Muon_eta = Muon_eta;
        this->Muon_phi = Muon_phi;
        this->Muon_mass = Muon_mass;
        this->Muon_charge = Muon_charge;
    }

    int MatchOppositelyChargedMuons::Test()
    {
        return 0;
    }


    /* This function matches each muon with oppositely charged muons. */
    RVec_f MatchOppositelyChargedMuons::InvariantMassCalculator()
    { 
        RVec_f pt {};
        RVec_f eta {};
        RVec_f phi {};
        RVec_f mass {};
        RVec_f InvariantMasses {};
        int i = 0;
        float im = 0;

        /* Loop over the set of muons to determine which muons have +1 charge, then match them with all the ones with -1 charge. */
        for(int charge : this->Muon_charge)
        {
            if(charge == 1)
            {
                /* If charges are opposite, calculate the invariant mass of them */
                for(int j = i+1; j < this->Muon_charge.size(); j++)
                {
                    if(charge == -1)
                    {
                        pt.push_back(Muon_pt[i]);
                        eta.push_back(Muon_eta[i]);
                        phi.push_back(Muon_phi[i]);
                        mass.push_back(Muon_mass[i]);

                        im = ROOT::VecOps::InvariantMass(pt, eta, phi, mass);

                        InvariantMasses.push_back(im);

                        pt.clear();
                        eta.clear();
                        phi.clear();
                        mass.clear();
                    }
                }
            }

            i++;
        }

        return InvariantMasses;
    }
}

#endif