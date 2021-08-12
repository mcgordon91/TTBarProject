#ifndef MUON_MATCHER
#define MUON_MATCHER

#include <iostream>

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
        RVec_f rvf;
        RVec_i rvi;
        float sum;
        
    public:
        MatchOppositelyChargedMuons();
        MatchOppositelyChargedMuons(RVec_f rvf, RVec_i rvi);
        float Test();
    };
    
    MatchOppositelyChargedMuons::MatchOppositelyChargedMuons()
    {
        this->rvf = {};
        this->rvi = {};
        cout << "Hi" << endl;
    }
    
    MatchOppositelyChargedMuons::MatchOppositelyChargedMuons(RVec_f rvf, RVec_i rvi)
    {
        this->rvf = rvf;
        this->rvi = rvi;
        cout << "Hi" << endl;
    }
    
    float MatchOppositelyChargedMuons::Test()
    {
        sum = ROOT::VecOps::Sum(this->rvf) + ROOT::VecOps::Sum(this->rvi);
        return sum;
    }
}

ROOT::RDF::RNode something(ROOT::RDF::RNode rdf)
{
    auto ret = rdf;
    ret = ret.Define("Muon_Matcher_Class", MuonMatcher::MatchOppositelyChargedMuons, {"Muon_pt", "Muon_charge"})\
        .Define("Muon_Matcher_Result", "Muon_Matcher_Class.Test()");
}

#endif