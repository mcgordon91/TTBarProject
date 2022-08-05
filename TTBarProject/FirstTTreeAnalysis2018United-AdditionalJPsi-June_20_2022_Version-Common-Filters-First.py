#!/usr/bin/env python
# coding: utf-8

# In[1]:


import ROOT
import glob
import numpy as np
import os
ROOT.gROOT.ProcessLine(".L FTFunctions.cpp")
ROOT.gInterpreter.Declare("""
    const UInt_t barWidth = 60;
    ULong64_t processed = 0, totalEvents = 0;
    std::string progressBar;
    std::mutex barMutex; 
    auto registerEvents = [](ULong64_t nIncrement) {totalEvents += nIncrement;};
    ROOT::RDF::RResultPtr<ULong64_t> AddProgressBar(ROOT::RDF::RNode df, int everyN=10000, int totalN=100000) {
        registerEvents(totalN);
        auto c = df.Count();
        c.OnPartialResultSlot(everyN, [everyN] (unsigned int slot, ULong64_t &cnt){
            std::lock_guard<std::mutex> l(barMutex);
            processed += everyN; //everyN captured by value for this lambda
            progressBar = "[";
            for(UInt_t i = 0; i < static_cast<UInt_t>(static_cast<Float_t>(processed)/totalEvents*barWidth); ++i){
                progressBar.push_back('|');
            }
            // escape the '\' when defined in python string
            std::cout << "\\r" << std::left << std::setw(barWidth) << progressBar << "] " << processed << "/" << totalEvents << std::flush;
        });
        return c;
    }
""")
# Enables multithreading
useRange = False
if not useRange:
    nThreads = 8
    ROOT.ROOT.EnableImplicitMT(nThreads)


# In[2]:


listOfFilesData = glob.glob("/eos/user/m/migordon/SWAN_projects/JPsiSkims/SingleMuon/*/*.root")
listOfFilesMonteCarloTTToSemiLeptonic = glob.glob("/eos/user/m/migordon/SWAN_projects/JPsiSkims/TTToSemiLeptonic_TuneCP5_13TeV-powheg-pythia8/*/*.root")
listOfFilesMonteCarloTTTo2L2Nu = glob.glob("/eos/user/m/migordon/SWAN_projects/JPsiSkims/TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8/*/*.root")
listOfFilesMonteCarloWJetsToLNu = glob.glob("/eos/user/m/migordon/SWAN_projects/JPsiSkims/WJetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8/*/*.root")
listOfFilesMonteCarloST_tW_top = glob.glob("/eos/user/m/migordon/SWAN_projects/JPsiSkims/ST_tW_top_5f_inclusiveDecays_TuneCP5_13TeV-powheg-pythia8/*/*.root")
listOfFilesMonteCarloST_tchannel_top = glob.glob("/eos/user/m/migordon/SWAN_projects/JPsiSkims/ST_t-channel_top_4f_InclusiveDecays_TuneCP5_13TeV-powheg-madspin-pythia8/*/*.root")

listOfFilesMonteCarloTTToSemiLeptonic = listOfFilesMonteCarloTTToSemiLeptonic[0:1]
listOfFilesMonteCarloTTTo2L2Nu = listOfFilesMonteCarloTTTo2L2Nu[0:1]
listOfFilesMonteCarloWJetsToLNu = listOfFilesMonteCarloWJetsToLNu[0:1]
listOfFilesMonteCarloST_tW_top = listOfFilesMonteCarloST_tW_top[0:1]
listOfFilesMonteCarloST_tchannel_top = listOfFilesMonteCarloST_tchannel_top[0:1]

dictOfListOfFiles = {"Data" : listOfFilesData, "MonteCarloTTToSemiLeptonic" : listOfFilesMonteCarloTTToSemiLeptonic, "MonteCarloTTTo2L2Nu" : listOfFilesMonteCarloTTTo2L2Nu, "MonteCarloWJetsToLNu" : listOfFilesMonteCarloWJetsToLNu, "MonteCarloST_tW_top" : listOfFilesMonteCarloST_tW_top, "MonteCarloST_tchannel_top" : listOfFilesMonteCarloST_tchannel_top} 


# In[3]:


chain = {}
meta = {}
rdf = {}
mrdf = {}
nevents = {}
sumweight = {}
neventsVal = {}

mureport = {}

rdfPassedIsolatedLeptonTrigger = {}
rdfIsolatedLeptonNoHighWeights = {}
rdfIsolatedMuonNoHighWeights = {}
rdfIsolatedElectronNoHighWeights = {}
rdfIsolatedMuonAfterMETCut = {}
rdfIsolatedElectronAfterMETCut = {}

rdfAnyNumberOfJets = {}
rdfTemp = {}
rdfJetAndIsolatedLeptonFiltered = {}

rdfJPsiMuons = {}


hist = {}
report = {}


nparray = {}
nparraynode = {}

# Muon_pfIsoId is PFIso ID from miniAOD selector (1=PFIsoVeryLoose, 2=PFIsoLoose, 3=PFIsoMedium, 4=PFIsoTight, 5=PFIsoVeryTight, 6=PFIsoVeryVeryTight)
LeadingIsolatedMuonMask = "Muon_pt > 30 && abs(Muon_eta) < 2.4 && Muon_mediumId == true && Muon_pfIsoId >= 4"
# SubleadingIsolatedMuonMask = ""
LeadingIsolatedElectronMask = "Electron_pt > 30 && abs(Electron_eta < 2.4) && Electron_cutBased == 4"
# SubleadingIsolatedElectronMask = ""
#DimuonMask = ""
#DielectronMask = ""
#MuonElectronMask = ""
JetMask = "ROOT::VecOps::RVec<Int_t> jmask = (Jet_pt >= 30 && abs(Jet_eta) <= 2.5 && Jet_jetId >= 2); "                          "for(int i=0; i < LeadingIsolatedMuon_pt.size(); ++i){"                              "ROOT::VecOps::RVec<Float_t> dr;"                              "for(int j=0; j < jmask.size(); ++j){"                                  "dr.push_back(ROOT::VecOps::DeltaR(Jet_eta.at(j), LeadingIsolatedMuon_eta.at(i), Jet_phi.at(j), LeadingIsolatedMuon_phi.at(i)));}"                                  "jmask = jmask && dr >= 0.4;"                                  "dr.clear();}"                          "return jmask;"
JPsiCandidateMask = "Muon_pt > 3 && abs(Muon_eta) <= 2.4 && Muon_mediumId == true && leading_isolated_muon_mask == false"


# In[4]:


for sample, fileList in dictOfListOfFiles.items():
    
    if sample == 'Data':
        vecList = ROOT.std.vector(str)()

        for element in dictOfListOfFiles['Data']:
            vecList.push_back(element)

        rdf['Data'] = ROOT.ROOT.RDataFrame("Events", vecList)    
        mureport['Data'] = rdf['Data'].Report()
        
        print(rdf['Data'].Report())
        mrdf['Data'] = ROOT.ROOT.RDataFrame("Runs", vecList)
        
        if useRange:
            rdf['Data'] = rdf['Data'].Range(5000)
            nrange = 5000
            printcode = ' if(rdfentry_ % 5000 == 0) { std::cout << "Processed entry " << rdfentry_ << "/' + str(nrange) + '" << std::endl; } return rdfentry_;'
            print("Data run")
            rdf['Data'] = rdf['Data'].Define("my_rdfentry", printcode)
            

        
    else:
        chain[sample] = ROOT.TChain("Events")
        meta[sample] = ROOT.TChain("Runs")

        for file in fileList:
            
            chain[sample].Add(file)
            meta[sample].Add(file)

        rdf[sample] = ROOT.ROOT.RDataFrame(chain[sample])
        mureport[sample] = rdf[sample].Report()
        mrdf[sample] = ROOT.ROOT.RDataFrame(meta[sample])
        
        if useRange:
            rdf[sample] = rdf[sample].Range(5000)
            
            printcode = ' if(rdfentry_ % 5000 == 0) { std::cout << "Processed entry " << rdfentry_ << " " << rdfslot_ << std::endl; } return rdfentry_;'
            print("MC run")
            rdf[sample] = rdf[sample].Define("my_rdfentry", printcode)
        
        nevents[sample] = mrdf[sample].Sum("genEventCount")
        sumweight[sample] = mrdf[sample].Sum("genEventSumw")


# In[5]:


#Semileptonic ttbar xsection: 364.3109
#Single mu trigger for 2017 (B,C,D,E,F): "HLT_IsoMu27"
#"HLT_Ele35_WPTight_Gsf"
#lumiDict = {"2017": 41.53, "2018": 59.97}
wgtFormula = {}

# wgtFormula used to weight each event
# XS = Literature Cross section of the process of interest (in picobarnes; the 1000 converts to femotobarnes), lumi = presumed luminosity of the data one is normalizing against; XS * lumi = # of expected events;
# genWeight = quantity stored in every event which comes from the Monte Carlo generator telling you what the value of the generated event is (usually close to 1); it can be + or -; it also contains matching
#     effeciency; tells you the Monte Carlo defined value of the event
# sW = sum of weights; normalizes the genWeight
lumiDict = {"2018": 59.97} #  brilcalc lumi --normtag /cvmfs/cms-bril.cern.ch/cms-lumi-pog/Normtags/normtag_PHYSICS.json -u /fb --begin 302031 --end 302663 --hltpath "HLT_IsoMu27*"
wgtFormula['Data'] = "1"
wgtFormula['MonteCarloTTToSemiLeptonic'] = "{XS:f} * {lumi:f} * 1000 * genWeight / {sW:f}".format(XS=364.31, lumi=lumiDict["2018"], sW=float(sumweight['MonteCarloTTToSemiLeptonic'].GetValue()))
wgtFormula['MonteCarloTTTo2L2Nu'] = "{XS:f} * {lumi:f} * 1000 * genWeight / {sW:f}".format(XS=87.33, lumi=lumiDict["2018"], sW=float(sumweight['MonteCarloTTTo2L2Nu'].GetValue()))
wgtFormula['MonteCarloWJetsToLNu'] = "{XS:f} * {lumi:f} * 1000 * genWeight / {sW:f}".format(XS=61526.7, lumi=lumiDict["2018"], sW=float(sumweight['MonteCarloWJetsToLNu'].GetValue()))
wgtFormula['MonteCarloST_tW_top'] = "{XS:f} * {lumi:f} * 1000 * genWeight / {sW:f}".format(XS=71.7, lumi=lumiDict["2018"], sW=float(sumweight['MonteCarloST_tW_top'].GetValue()))
wgtFormula['MonteCarloST_tchannel_top'] = "{XS:f} * {lumi:f} * 1000 * genWeight / {sW:f}".format(XS=130, lumi=lumiDict["2018"], sW=float(sumweight['MonteCarloST_tchannel_top'].GetValue()))


# In[6]:


cpp_code = """
typedef ROOT::VecOps::RVec<Float_t>                        RVec_f;
typedef ROOT::VecOps::RVec<Int_t>                          RVec_i;
typedef ROOT::VecOps::RVec<Long_t>                         RVec_l;

class MatchOppositelyChargedMuons
{
    private:
        RVec_f Paired_Muon_pt;
        RVec_f Paired_Muon_eta;
        RVec_f Paired_Muon_phi;
        RVec_f Paired_Muon_mass;
        RVec_i Paired_Muon_charge;
        long Paired_EventNumber;
        int FlagPair;
        RVec_f Isolated_Muon_pt;
        RVec_f Isolated_Muon_eta;
        RVec_f Isolated_Muon_phi;
        RVec_f Isolated_Muon_mass;
        RVec_i Isolated_Muon_charge;


    public:
    MatchOppositelyChargedMuons(RVec_f Paired_Muon_pt, RVec_f Paired_Muon_eta, RVec_f Paired_Muon_phi, RVec_f Paired_Muon_mass, RVec_i Paired_Muon_charge, long Paired_EventNumber, int FlagPair);
        MatchOppositelyChargedMuons(RVec_f Paired_Muon_pt, RVec_f Paired_Muon_eta, RVec_f Paired_Muon_phi, RVec_f Paired_Muon_mass, RVec_i Paired_Muon_charge, long Paired_EventNumber, int FlagPair, RVec_f Isolated_Muon_pt, RVec_f Isolated_Muon_eta, RVec_f Isolated_Muon_phi, RVec_f Isolated_Muon_mass, RVec_i Isolated_Muon_charge);
        RVec_f InvariantMassCalculator();
        RVec_f OppositelyChargedMuonInvariantMassCalculator();
        RVec_f ThreeMuonInvariantMassCalculator();
        RVec_f DeltaEtaCalculator();
        RVec_f DeltaPhiCalculator();
        RVec_f DeltaRCalculator();
        int ReturnFlagPair();
};

MatchOppositelyChargedMuons::MatchOppositelyChargedMuons(RVec_f Paired_Muon_pt, RVec_f Paired_Muon_eta, RVec_f Paired_Muon_phi, RVec_f Paired_Muon_mass, RVec_i Paired_Muon_charge, long Paired_EventNumber, int FlagPair)
{
    this->Paired_Muon_pt = Paired_Muon_pt;
    this->Paired_Muon_eta = Paired_Muon_eta;
    this->Paired_Muon_phi = Paired_Muon_phi;
    this->Paired_Muon_mass = Paired_Muon_mass;
    this->Paired_Muon_charge = Paired_Muon_charge;
    this->Paired_EventNumber = Paired_EventNumber;
    this->FlagPair = FlagPair;
    this->Isolated_Muon_pt = {};
    this->Isolated_Muon_eta = {};
    this->Isolated_Muon_phi = {};
    this->Isolated_Muon_mass = {};
    this->Isolated_Muon_charge = {};
}

MatchOppositelyChargedMuons::MatchOppositelyChargedMuons(RVec_f Paired_Muon_pt, RVec_f Paired_Muon_eta, RVec_f Paired_Muon_phi, RVec_f Paired_Muon_mass, RVec_i Paired_Muon_charge, long Paired_EventNumber, int FlagPair, RVec_f Isolated_Muon_pt, RVec_f Isolated_Muon_eta, RVec_f Isolated_Muon_phi, RVec_f Isolated_Muon_mass, RVec_i Isolated_Muon_charge)
{
    this->Paired_Muon_pt = Paired_Muon_pt;
    this->Paired_Muon_eta = Paired_Muon_eta;
    this->Paired_Muon_phi = Paired_Muon_phi;
    this->Paired_Muon_mass = Paired_Muon_mass;
    this->Paired_Muon_charge = Paired_Muon_charge;
    this->Paired_EventNumber = Paired_EventNumber;
    this->FlagPair = FlagPair;
    this->Isolated_Muon_pt = Isolated_Muon_pt;
    this->Isolated_Muon_eta = Isolated_Muon_eta;
    this->Isolated_Muon_phi = Isolated_Muon_phi;
    this->Isolated_Muon_mass = Isolated_Muon_mass;
    this->Isolated_Muon_charge = Isolated_Muon_charge;
}


/* This function matches each muon with oppositely charged muons. */
RVec_f MatchOppositelyChargedMuons::OppositelyChargedMuonInvariantMassCalculator()
{ 
    RVec_f pt {};
    RVec_f eta {};
    RVec_f phi {};
    RVec_f mass {};
    
    int FirstMuonCharge = 0;
    int SecondMuonCharge = 0;
    
    RVec_f InvariantMasses {};
    
    float im = 0;
 
    /* Loop over the set of muons to determine which muons have +1 charge, then match them with all the ones with -1 charge. */
    for(int i = 0; i < this->Paired_Muon_charge.size(); i++)
    {    

        FirstMuonCharge = this->Paired_Muon_charge[i];
        
        /* If charges are opposite, calculate the invariant mass of them */
        for(int j = i+1; j < this->Paired_Muon_charge.size(); j++)
        {
            SecondMuonCharge = this->Paired_Muon_charge[j];
            
            if(FirstMuonCharge * SecondMuonCharge == -1)
            {
                pt.push_back(Paired_Muon_pt[i]);
                eta.push_back(Paired_Muon_eta[i]);
                phi.push_back(Paired_Muon_phi[i]);
                mass.push_back(Paired_Muon_mass[i]);
                    
                pt.push_back(Paired_Muon_pt[j]);
                eta.push_back(Paired_Muon_eta[j]);
                phi.push_back(Paired_Muon_phi[j]);
                mass.push_back(Paired_Muon_mass[j]);
                    
                im = ROOT::VecOps::InvariantMass(pt, eta, phi, mass);
                
                InvariantMasses.push_back(im);
                
                this->FlagPair = 1;
                    
                pt.clear();
                eta.clear();
                phi.clear();
                mass.clear();
            }
        }
    }
        
    return InvariantMasses;
}


RVec_f MatchOppositelyChargedMuons::ThreeMuonInvariantMassCalculator()
{ 
    RVec_f pt {};
    RVec_f eta {};
    RVec_f phi {};
    RVec_f mass {};
    
    int FirstMuonCharge = 0;
    int SecondMuonCharge = 0;
    
    RVec_f InvariantMasses {};
    
    float im = 0;
 
    for(int k = 0; k < this->Isolated_Muon_pt.size(); k++)
    {
    
        /* Loop over the set of muons to determine which muons have +1 charge, then match them with all the ones with -1 charge. */
        for(int i = 0; i < this->Paired_Muon_charge.size(); i++)
        {    

            FirstMuonCharge = this->Paired_Muon_charge[i];

            /* If charges are opposite, calculate the invariant mass of them */
            for(int j = i+1; j < this->Paired_Muon_charge.size(); j++)
            {
                SecondMuonCharge = this->Paired_Muon_charge[j];

                if(FirstMuonCharge * SecondMuonCharge == -1)
                {
                    pt.push_back(Paired_Muon_pt[i]);
                    eta.push_back(Paired_Muon_eta[i]);
                    phi.push_back(Paired_Muon_phi[i]);
                    mass.push_back(Paired_Muon_mass[i]);

                    pt.push_back(Paired_Muon_pt[j]);
                    eta.push_back(Paired_Muon_eta[j]);
                    phi.push_back(Paired_Muon_phi[j]);
                    mass.push_back(Paired_Muon_mass[j]);
                    
                    pt.push_back(Isolated_Muon_pt[k]);
                    eta.push_back(Isolated_Muon_eta[k]);
                    phi.push_back(Isolated_Muon_phi[k]);
                    mass.push_back(Isolated_Muon_mass[k]);

                    im = ROOT::VecOps::InvariantMass(pt, eta, phi, mass);

                    InvariantMasses.push_back(im);

                    this->FlagPair = 1;

                    pt.clear();
                    eta.clear();
                    phi.clear();
                    mass.clear();
                }
            }
        }
    }
        
    return InvariantMasses;
}


RVec_f MatchOppositelyChargedMuons::DeltaEtaCalculator()
{ 
    float DeltaEtaIndividual = 0;
    
    RVec_f DeltaEtaRVec {};
    
    for(int i = 0; i < this->Isolated_Muon_eta.size(); i++)
    {
    
        for(int j = 0; j < this->Paired_Muon_eta.size(); j++)
        {    
                    DeltaEtaIndividual = this->Isolated_Muon_eta[i] - this->Paired_Muon_eta[j];

                    DeltaEtaRVec.push_back(DeltaEtaIndividual);
        }
    }
        
    return DeltaEtaRVec;
}


RVec_f MatchOppositelyChargedMuons::DeltaPhiCalculator()
{ 
    RVec_f PhiAll {};
    RVec_f PhiJPsi {};
    
    RVec_f DeltaPhiIndividual {};
    
    RVec_f DeltaPhiRVec {};
    
    for(int i = 0; i < this->Isolated_Muon_phi.size(); i++)
    {
    
        for(int j = 0; j < this->Paired_Muon_phi.size(); j++)
        {    

                    PhiAll.push_back(this->Isolated_Muon_phi[i]);

                    PhiJPsi.push_back(this->Paired_Muon_phi[j]);

                    DeltaPhiIndividual = ROOT::VecOps::DeltaPhi(PhiAll, PhiJPsi);

                    for(int k = 0; k < DeltaPhiIndividual.size(); k++)
                    {
                        DeltaPhiRVec.push_back(DeltaPhiIndividual[k]);
                    }

                    PhiAll.clear();
                    PhiJPsi.clear();
        }
    }
        
    return DeltaPhiRVec;
}


RVec_f MatchOppositelyChargedMuons::DeltaRCalculator()
{ 
    RVec_f EtaAll {};
    RVec_f PhiAll {};
    RVec_f EtaJPsi {};
    RVec_f PhiJPsi {};
    
    RVec_f DeltaRIndividual {};
    
    RVec_f DeltaRRVec {};
    
    for(int i = 0; i < this->Isolated_Muon_phi.size(); i++)
    {
    
        for(int j = 0; j < this->Paired_Muon_phi.size(); j++)
        {    
                    EtaAll.push_back(this->Isolated_Muon_eta[i]);
                    PhiAll.push_back(this->Isolated_Muon_phi[i]);
                    
                    EtaJPsi.push_back(this->Paired_Muon_eta[j]);
                    PhiJPsi.push_back(this->Paired_Muon_phi[j]);

                    DeltaRIndividual = ROOT::VecOps::DeltaR(EtaAll, EtaJPsi, PhiAll, PhiJPsi);

                    for(int k = 0; k < DeltaRIndividual.size(); k++)
                    {
                        DeltaRRVec.push_back(DeltaRIndividual[k]);
                    }

                    EtaAll.clear();
                    PhiAll.clear();
                    EtaJPsi.clear();
                    PhiJPsi.clear();
        }
    }
        
    return DeltaRRVec;
}


int MatchOppositelyChargedMuons::ReturnFlagPair()
{
    return this->FlagPair;
}
"""

ROOT.gInterpreter.Declare(cpp_code)


# In[7]:


def IsolatedLeptonSelection():

    for sample in dictOfListOfFiles:       

        rdfPassedIsolatedLeptonTrigger[sample] = rdf[sample].Filter("HLT_IsoMu24 == true", "HLTLeptonTrigger")            .Define("LumiXS",wgtFormula[sample])            .Define("leading_isolated_muon_mask", LeadingIsolatedMuonMask)            .Define("leading_isolated_electron_mask", LeadingIsolatedElectronMask)            .Define("jpsi_mu_candidate_mask", JPsiCandidateMask)
            #^ HLT_Ele32_WPTight_Gsf == true", "HLTLeptonTrigger")\
        if sample == 'Data':
            rdfIsolatedMuonNoHighWeights[sample] = rdfPassedIsolatedLeptonTrigger[sample].Filter("Sum(leading_isolated_muon_mask) == 1 && Sum(leading_isolated_electron_mask) == 0", "Exactly one isolated muon and exactly zero isolated electrons")                .Define("LeadingIsolatedMuon_pt", "Muon_pt[leading_isolated_muon_mask].at(0, -10)")                .Define("LeadingIsolatedMuon_eta", "Muon_eta[leading_isolated_muon_mask].at(0, -10)")                .Define("LeadingIsolatedMuon_phi", "Muon_phi[leading_isolated_muon_mask].at(0, -5)")                .Define("LeadingIsolatedMuon_mass", "Muon_mass[leading_isolated_muon_mask].at(0, -2.71)")                .Define("LeadingIsolatedMuon_charge", "Muon_charge[leading_isolated_muon_mask].at(0, -5)")                .Define("IsolatedMuon_pt", "Muon_pt[leading_isolated_muon_mask]")                .Define("IsolatedMuon_eta", "Muon_eta[leading_isolated_muon_mask]")                .Define("IsolatedMuon_phi", "Muon_phi[leading_isolated_muon_mask]")                .Define("IsolatedMuon_mass", "Muon_mass[leading_isolated_muon_mask]")                .Define("IsolatedMuon_charge", "Muon_charge[leading_isolated_muon_mask]")                .Define("IsolatedMuon_pdgId", "Muon_pdgId[leading_isolated_muon_mask]")                .Define("IsolatedElectron_pt", "Electron_pt[leading_isolated_electron_mask]")                .Define("IsolatedElectron_eta", "Electron_eta[leading_isolated_electron_mask]")                .Define("IsolatedElectron_phi", "Electron_phi[leading_isolated_electron_mask]")                .Define("IsolatedElectron_mass", "Electron_mass[leading_isolated_electron_mask]")                .Define("IsolatedElectron_charge", "Electron_charge[leading_isolated_electron_mask]")                .Define("IsolatedElectron_pdgId", "Electron_pdgId[leading_isolated_electron_mask]")                .Define("IsolatedLepton_pt", "Concatenate(IsolatedMuon_pt, IsolatedElectron_pt)")                .Define("IsolatedLepton_eta", "Concatenate(IsolatedMuon_eta, IsolatedElectron_eta)")                .Define("IsolatedLepton_phi", "Concatenate(IsolatedMuon_phi, IsolatedElectron_phi)")                .Define("IsolatedLepton_mass", "Concatenate(IsolatedMuon_mass, IsolatedElectron_mass)")                .Define("IsolatedLepton_charge", "Concatenate(IsolatedMuon_charge, IsolatedElectron_charge)")                .Define("IsolatedLepton_pdgid", "Concatenate(IsolatedMuon_pdgId, IsolatedElectron_pdgId)")                .Define("METBeforeMETCut", "MET_pt")
            
            rdfIsolatedElectronNoHighWeights[sample] = rdfPassedIsolatedLeptonTrigger[sample].Filter("Sum(leading_isolated_electron_mask) == 1 && Sum(leading_isolated_muon_mask) == 0", "Exactly one isolated electron and exactly zero isolated muon")                .Define("LeadingIsolatedElectron_pt", "Electron_pt[leading_isolated_electron_mask].at(0, -10)")                .Define("LeadingIsolatedElectron_eta", "Electron_eta[leading_isolated_electron_mask].at(0, -10)")                .Define("LeadingIsolatedElectron_phi", "Electron_phi[leading_isolated_electron_mask].at(0, -5)")                .Define("LeadingIsolatedElectron_mass", "Electron_mass[leading_isolated_electron_mask].at(0, -2.71)")                .Define("LeadingIsolatedElectron_charge", "Electron_charge[leading_isolated_electron_mask].at(0, -5)")                .Define("IsolatedMuon_pt", "Muon_pt[leading_isolated_muon_mask]")                .Define("IsolatedMuon_eta", "Muon_eta[leading_isolated_muon_mask]")                .Define("IsolatedMuon_phi", "Muon_phi[leading_isolated_muon_mask]")                .Define("IsolatedMuon_mass", "Muon_mass[leading_isolated_muon_mask]")                .Define("IsolatedMuon_charge", "Muon_charge[leading_isolated_muon_mask]")                .Define("IsolatedMuon_pdgId", "Muon_pdgId[leading_isolated_muon_mask]")                .Define("IsolatedElectron_pt", "Electron_pt[leading_isolated_electron_mask]")                .Define("IsolatedElectron_eta", "Electron_eta[leading_isolated_electron_mask]")                .Define("IsolatedElectron_phi", "Electron_phi[leading_isolated_electron_mask]")                .Define("IsolatedElectron_mass", "Electron_mass[leading_isolated_electron_mask]")                .Define("IsolatedElectron_charge", "Electron_charge[leading_isolated_electron_mask]")                .Define("IsolatedElectron_pdgId", "Electron_pdgId[leading_isolated_electron_mask]")                .Define("IsolatedLepton_pt", "Concatenate(IsolatedMuon_pt, IsolatedElectron_pt)")                .Define("IsolatedLepton_eta", "Concatenate(IsolatedMuon_eta, IsolatedElectron_eta)")                .Define("IsolatedLepton_phi", "Concatenate(IsolatedMuon_phi, IsolatedElectron_phi)")                .Define("IsolatedLepton_mass", "Concatenate(IsolatedMuon_mass, IsolatedElectron_mass)")                .Define("IsolatedLepton_charge", "Concatenate(IsolatedMuon_charge, IsolatedElectron_charge)")                .Define("IsolatedLepton_pdgid", "Concatenate(IsolatedMuon_pdgId, IsolatedElectron_pdgId)")                .Define("METBeforeMETCut", "MET_pt")


        # Change criteria to nothing below or above 4 standard deviations from the mean
        else:
            rdfIsolatedMuonNoHighWeights[sample] = rdfPassedIsolatedLeptonTrigger[sample].Filter("genWeight < 1000 && Sum(leading_isolated_muon_mask) == 1 && Sum(leading_isolated_electron_mask) == 0", "Exactly one isolated muon and exactly zero isolated electrons")                .Define("LeadingIsolatedMuon_pt", "Muon_pt[leading_isolated_muon_mask].at(0, -10)")                .Define("LeadingIsolatedMuon_eta", "Muon_eta[leading_isolated_muon_mask].at(0, -10)")                .Define("LeadingIsolatedMuon_phi", "Muon_phi[leading_isolated_muon_mask].at(0, -5)")                .Define("LeadingIsolatedMuon_mass", "Muon_mass[leading_isolated_muon_mask].at(0, -2.71)")                .Define("LeadingIsolatedMuon_charge", "Muon_charge[leading_isolated_muon_mask].at(0, -5)")                .Define("IsolatedMuon_pt", "Muon_pt[leading_isolated_muon_mask]")                .Define("IsolatedMuon_eta", "Muon_eta[leading_isolated_muon_mask]")                .Define("IsolatedMuon_phi", "Muon_phi[leading_isolated_muon_mask]")                .Define("IsolatedMuon_mass", "Muon_mass[leading_isolated_muon_mask]")                .Define("IsolatedMuon_charge", "Muon_charge[leading_isolated_muon_mask]")                .Define("IsolatedMuon_pdgId", "Muon_pdgId[leading_isolated_muon_mask]")                .Define("IsolatedElectron_pt", "Electron_pt[leading_isolated_electron_mask]")                .Define("IsolatedElectron_eta", "Electron_eta[leading_isolated_electron_mask]")                .Define("IsolatedElectron_phi", "Electron_phi[leading_isolated_electron_mask]")                .Define("IsolatedElectron_mass", "Electron_mass[leading_isolated_electron_mask]")                .Define("IsolatedElectron_charge", "Electron_charge[leading_isolated_electron_mask]")                .Define("IsolatedElectron_pdgId", "Electron_pdgId[leading_isolated_electron_mask]")                .Define("IsolatedLepton_pt", "Concatenate(IsolatedMuon_pt, IsolatedElectron_pt)")                .Define("IsolatedLepton_eta", "Concatenate(IsolatedMuon_eta, IsolatedElectron_eta)")                .Define("IsolatedLepton_phi", "Concatenate(IsolatedMuon_phi, IsolatedElectron_phi)")                .Define("IsolatedLepton_mass", "Concatenate(IsolatedMuon_mass, IsolatedElectron_mass)")                .Define("IsolatedLepton_charge", "Concatenate(IsolatedMuon_charge, IsolatedElectron_charge)")                .Define("IsolatedLepton_pdgid", "Concatenate(IsolatedMuon_pdgId, IsolatedElectron_pdgId)")                .Define("METBeforeMETCut", "MET_pt")

            rdfIsolatedElectronNoHighWeights[sample] = rdfPassedIsolatedLeptonTrigger[sample].Filter("genWeight < 1000 && Sum(leading_isolated_electron_mask) == 1 && Sum(leading_isolated_muon_mask) == 0", "Exactly one isolated electron and exactly zero isolated muon")                .Define("LeadingIsolatedElectron_pt", "Electron_pt[leading_isolated_electron_mask].at(0, -10)")                .Define("LeadingIsolatedElectron_eta", "Electron_eta[leading_isolated_electron_mask].at(0, -10)")                .Define("LeadingIsolatedElectron_phi", "Electron_phi[leading_isolated_electron_mask].at(0, -5)")                .Define("LeadingIsolatedElectron_mass", "Electron_mass[leading_isolated_electron_mask].at(0, -2.71)")                .Define("LeadingIsolatedElectron_charge", "Electron_charge[leading_isolated_electron_mask].at(0, -5)")                .Define("IsolatedMuon_pt", "Muon_pt[leading_isolated_muon_mask]")                .Define("IsolatedMuon_eta", "Muon_eta[leading_isolated_muon_mask]")                .Define("IsolatedMuon_phi", "Muon_phi[leading_isolated_muon_mask]")                .Define("IsolatedMuon_mass", "Muon_mass[leading_isolated_muon_mask]")                .Define("IsolatedMuon_charge", "Muon_charge[leading_isolated_muon_mask]")                .Define("IsolatedMuon_pdgId", "Muon_pdgId[leading_isolated_muon_mask]")                .Define("IsolatedElectron_pt", "Electron_pt[leading_isolated_electron_mask]")                .Define("IsolatedElectron_eta", "Electron_eta[leading_isolated_electron_mask]")                .Define("IsolatedElectron_phi", "Electron_phi[leading_isolated_electron_mask]")                .Define("IsolatedElectron_mass", "Electron_mass[leading_isolated_electron_mask]")                .Define("IsolatedElectron_charge", "Electron_charge[leading_isolated_electron_mask]")                .Define("IsolatedElectron_pdgId", "Electron_pdgId[leading_isolated_electron_mask]")                .Define("IsolatedLepton_pt", "Concatenate(IsolatedMuon_pt, IsolatedElectron_pt)")                .Define("IsolatedLepton_eta", "Concatenate(IsolatedMuon_eta, IsolatedElectron_eta)")                .Define("IsolatedLepton_phi", "Concatenate(IsolatedMuon_phi, IsolatedElectron_phi)")                .Define("IsolatedLepton_mass", "Concatenate(IsolatedMuon_mass, IsolatedElectron_mass)")                .Define("IsolatedLepton_charge", "Concatenate(IsolatedMuon_charge, IsolatedElectron_charge)")                .Define("IsolatedLepton_pdgid", "Concatenate(IsolatedMuon_pdgId, IsolatedElectron_pdgId)")                .Define("METBeforeMETCut", "MET_pt")

        rdfIsolatedMuonAfterMETCut[sample] = rdfIsolatedMuonNoHighWeights[sample].Filter("MET_pt > 30", "Muon MET Greater than 30 GeV")            .Define("METAfterMETCut", "MET_pt")
        
        rdfIsolatedElectronAfterMETCut[sample] = rdfIsolatedElectronNoHighWeights[sample].Filter("MET_pt > 30", "Electron MET Greater than 30 GeV")            .Define("METAfterMETCut", "MET_pt")
        
    return rdfIsolatedMuonAfterMETCut, rdfIsolatedElectronAfterMETCut


# In[8]:


rdfIsolatedMuonAfterMETCut, rdfIsolatedElectronAfterMETCut = IsolatedLeptonSelection()


# In[9]:


def JetSelection():

    for sample in dictOfListOfFiles:

        rdfAnyNumberOfJets[sample] = rdfIsolatedMuonAfterMETCut[sample]            .Define("jet_mask", "ROOT::VecOps::RVec<Int_t> jmask = (Jet_pt >= 30 && abs(Jet_eta) <= 2.5 && Jet_jetId >= 2); "                          "for(int i=0; i < IsolatedLepton_pt.size(); ++i){"                              "ROOT::VecOps::RVec<Float_t> dr;"                              "for(int j=0; j < jmask.size(); ++j){"                                  "dr.push_back(ROOT::VecOps::DeltaR(Jet_eta.at(j), IsolatedLepton_eta.at(i), Jet_phi.at(j), IsolatedLepton_phi.at(i)));}"                                  "jmask = jmask && dr >= 0.4;"                                  "dr.clear();}"                          "return jmask;")            .Define("Num_Jets", "Jet_pt[jet_mask].size()")
        
        rdfJetAndIsolatedLeptonFiltered[sample] = rdfAnyNumberOfJets[sample].Filter("Num_Jets >= 2", "At Least Two Jets")            .Define("SJet1_pt", "Jet_pt[jet_mask].size() > 0 ? Jet_pt[jet_mask].at(0) : -500")            .Define("SJet2_pt", "Jet_pt[jet_mask].size() > 1 ? Jet_pt[jet_mask].at(1) : -500")            .Define("SJet1_eta", "Jet_eta[jet_mask].size() > 0 ? Jet_eta[jet_mask].at(0) : 500")            .Define("SJet2_eta", "Jet_eta[jet_mask].size() > 1 ? Jet_eta[jet_mask].at(1) : 500")            .Define("SJet1_phi", "Jet_phi[jet_mask].size() > 0 ? Jet_phi[jet_mask].at(0) : 500")            .Define("SJet2_phi", "Jet_phi[jet_mask].size() > 1 ? Jet_phi[jet_mask].at(1) : 500")            .Define("MediumBJetMask", "Jet_btagDeepFlavB > 0.3033 && jet_mask" )            .Define("MTofMETandMu", "FTA::transverseMassMET(IsolatedMuon_pt, IsolatedMuon_phi, IsolatedMuon_mass, MET_pt, MET_phi)")            .Define("SJet_btagDeepFlavB", "Jet_btagDeepFlavB[jet_mask]")            .Define("Num_BTaggedJets", "Sum(MediumBJetMask)")            .Define("Ht", "Sum(Jet_pt[jet_mask])")
        
    return rdfJetAndIsolatedLeptonFiltered


# In[10]:


rdfJetAndIsolatedLeptonFiltered = JetSelection()


# In[11]:


def JPsiSelection():

    for sample in dictOfListOfFiles:

        rdfJPsiMuons[sample] = rdfJetAndIsolatedLeptonFiltered[sample].Filter("Sum(jpsi_mu_candidate_mask) >= 2", "JPsi Candidate")            .Define("Muon_FourMomentum", "auto fourVecs = ROOT::VecOps::Construct<ROOT::Math::PtEtaPhiMVector>(Muon_pt, Muon_eta, Muon_phi, Muon_mass); return fourVecs;")            .Define("InvariantMassJPsiMuons", "std::cout << rdfentry_ << std::endl; int FlagPair = -1; auto c = MatchOppositelyChargedMuons(Muon_pt, Muon_eta, Muon_phi, Muon_mass, Muon_charge, event, FlagPair); return c.OppositelyChargedMuonInvariantMassCalculator()")            .Define("JPsiCandidate_packedIndices", "ROOT::VecOps::Combinations(Muon_pt[jpsi_mu_candidate_mask], 2)")            .Define("JPsiCandidate_idx0", "JPsiCandidate_packedIndices[0]")            .Define("JPsiCandidate_idx1",  "JPsiCandidate_packedIndices[1]")            .Define("JPsiCandidate_four_momentum", "ROOT::VecOps::Take(Muon_FourMomentum[jpsi_mu_candidate_mask], JPsiCandidate_idx0) + ROOT::VecOps::Take(Muon_FourMomentum[jpsi_mu_candidate_mask], JPsiCandidate_idx1)")            .Define("JPsiCandidate_pt", "ROOT::VecOps::RVec<double> temp; for (auto p: JPsiCandidate_four_momentum) {temp.push_back(p.Pt());} return temp;")            .Define("JPsiCandidate_eta", "ROOT::VecOps::RVec<double> temp; for (auto p: JPsiCandidate_four_momentum) {temp.push_back(p.Eta());} return temp;")            .Define("JPsiCandidate_phi", "ROOT::VecOps::RVec<double> temp; for (auto p: JPsiCandidate_four_momentum) {temp.push_back(p.Phi());} return temp;")            .Define("JPsiCandidate_mass", "ROOT::VecOps::RVec<double> temp; for (auto p: JPsiCandidate_four_momentum) {temp.push_back(p.M());} return temp;")            .Define("JPsiCandidate_charge", "ROOT::VecOps::Take(Muon_charge[jpsi_mu_candidate_mask], JPsiCandidate_idx0) + ROOT::VecOps::Take(Muon_charge[jpsi_mu_candidate_mask], JPsiCandidate_idx1)")            .Define("DeltaEta", "int FlagPair = -1; auto c = MatchOppositelyChargedMuons(IsolatedMuon_pt, IsolatedMuon_eta, IsolatedMuon_phi, IsolatedMuon_mass, IsolatedMuon_charge, event, FlagPair, JPsiCandidate_pt, JPsiCandidate_eta, JPsiCandidate_phi, JPsiCandidate_mass, JPsiCandidate_charge); return c.DeltaEtaCalculator()")            .Define("DeltaPhi", "int FlagPair = -1; auto c = MatchOppositelyChargedMuons(IsolatedMuon_pt, IsolatedMuon_eta, IsolatedMuon_phi, IsolatedMuon_mass, IsolatedMuon_charge, event, FlagPair, JPsiCandidate_pt, JPsiCandidate_eta, JPsiCandidate_phi, JPsiCandidate_mass, JPsiCandidate_charge); return c.DeltaPhiCalculator()")            .Define("DeltaR", "int FlagPair = -1; auto c = MatchOppositelyChargedMuons(IsolatedMuon_pt, IsolatedMuon_eta, IsolatedMuon_phi, IsolatedMuon_mass, IsolatedMuon_charge, event, FlagPair, JPsiCandidate_pt, JPsiCandidate_eta, JPsiCandidate_phi, JPsiCandidate_mass, JPsiCandidate_charge); return c.DeltaRCalculator()")

    return rdfJPsiMuons


# In[12]:


rdfJPsiMuons = JPsiSelection()


# In[13]:


# Make this block a function with the inputs as the nodes we want to attach histograms to (2-17-22)
# rdfWithFourMomentum needs to be changed to whatever our final RDataFrame is (5-23-22)

for sample in dictOfListOfFiles:
    
    if sample not in hist.keys():
        hist[sample] = {}
        report[sample] = rdf[sample].Report()
        #nparraynode[sample] = rdfLeadingMuon[sample]
        
    if sample == 'Data':
        
        hist['Data']["leading_isolated_muon_pt"] = rdfIsolatedMuonNoHighWeights['Data'].Histo1D(("leading_isolated_muon_pt","Leading Isolated Muon Transverse Momentum; Pt (GeV);Events",100,20,220),"LeadingIsolatedMuon_pt")
        hist['Data']["leading_isolated_muon_eta"] = rdfIsolatedMuonNoHighWeights['Data'].Histo1D(("leading_isolated_muon_eta", "Leading Isolated Muon Pseudorapidity; Eta; Events",100,-3,3),"LeadingIsolatedMuon_eta")
        hist['Data']["leading_isolated_muon_phi"] = rdfIsolatedMuonNoHighWeights['Data'].Histo1D(("leading_isolated_muon_phi", "Leading Isolated Muon Angle; Phi (Radians); Events",100,-3.5,3.5),"LeadingIsolatedMuon_phi")
        hist['Data']["leading_isolated_muon_mass"] = rdfIsolatedMuonNoHighWeights['Data'].Histo1D(("leading_isolated_muon_mass", "Leading Isolated Muon Mass; Mass(Gev); Events",10,.1,.2),"LeadingIsolatedMuon_mass")
        hist['Data']["leading_isolated_muon_charge"] = rdfIsolatedMuonNoHighWeights['Data'].Histo1D(("leading_isolated_muon_charge", "Leading Isolated Muon Charge; Charge; Events",5,-2,2),"LeadingIsolatedMuon_charge")
        hist['Data']["leading_isolated_electron_pt"] = rdfIsolatedElectronNoHighWeights['Data'].Histo1D(("leading_isolated_electron_pt","Leading Isolated Electron Transverse Momentum; Pt (GeV);Events",100,20,220),"LeadingIsolatedElectron_pt")
        hist['Data']["leading_isolated_electron_eta"] = rdfIsolatedElectronNoHighWeights['Data'].Histo1D(("leading_isolated_electron_eta", "Leading Isolated Electron Pseudorapidity; Eta; Events",100,-3,3),"LeadingIsolatedElectron_eta")
        hist['Data']["leading_isolated_electron_phi"] = rdfIsolatedElectronNoHighWeights['Data'].Histo1D(("leading_isolated_electron_phi", "Leading Isolated Electron Angle; Phi (Radians); Events",100,-3.5,3.5),"LeadingIsolatedElectron_phi")
        hist['Data']["leading_isolated_electron_mass"] = rdfIsolatedElectronNoHighWeights['Data'].Histo1D(("leading_isolated_electron_mass", "Leading Isolated Electron Mass; Mass(Gev); Events",10,.1,.2),"LeadingIsolatedElectron_mass")
        hist['Data']["leading_isolated_electron_charge"] = rdfIsolatedElectronNoHighWeights['Data'].Histo1D(("leading_isolated_electron_charge", "Leading Isolated Electron Charge; Charge; Events",5,-2,2),"LeadingIsolatedElectron_charge")
        hist['Data']["met_before_met_cut"] = rdfIsolatedMuonNoHighWeights['Data'].Histo1D(("met_before_met_cut", "MET Before MET Cut; Pt (GeV); Events",100,0,250),"METBeforeMETCut")
        hist['Data']["met_after_met_cut"] = rdfIsolatedMuonAfterMETCut['Data'].Histo1D(("met_after_met_cut", "MET After MET Cut; Pt (GeV); Events",100,0,250),"METAfterMETCut")
        hist['Data']["leading_isolated_muon_pt_after_met_cut"] = rdfJPsiMuons['Data'].Histo1D(("leading_isolated_muon_pt_after_met_cut","Leading Isolated Muon Transverse Momentum After MET Cut; Pt (GeV);Events",100,20,220),"LeadingIsolatedMuon_pt")
        hist['Data']["leading_isolated_muon_eta_after_met_cut"] = rdfJPsiMuons['Data'].Histo1D(("leading_isolated_muon_eta_after_met_cut", "Leading Isolated Muon Pseudorapidity After MET Cut; Eta; Events",100,-3,3),"LeadingIsolatedMuon_eta")
        hist['Data']["leading_isolated_muon_phi_after_met_cut"] = rdfJPsiMuons['Data'].Histo1D(("leading_isolated_muon_phi_after_met_cut", "Leading Isolated Muon Angle After MET Cut; Phi (Radians); Events",100,-3.5,3.5),"LeadingIsolatedMuon_phi")
        hist['Data']["leading_isolated_muon_mass_after_met_cut"] = rdfJPsiMuons['Data'].Histo1D(("leading_isolated_muon_mass_after_met_cut", "Leading Isolated Muon Mass After MET Cut; Mass(Gev); Events",10,.1,.2),"LeadingIsolatedMuon_mass")
        hist['Data']["leading_isolated_muon_charge_after_met_cut"] = rdfJPsiMuons['Data'].Histo1D(("leading_isolated_muon_charge_after_met_cut", "Leading Isolated Muon Charge After MET Cut; Charge; Events",5,-2,2),"LeadingIsolatedMuon_charge")
        hist['Data']["leading_isolated_electron_pt_after_met_cut"] = rdfIsolatedElectronAfterMETCut['Data'].Histo1D(("leading_isolated_electron_pt_after_met_cut","Leading Isolated Electron Transverse Momentum After MET Cut; Pt (GeV);Events",100,20,220),"LeadingIsolatedElectron_pt")
        hist['Data']["leading_isolated_electron_eta_after_met_cut"] = rdfIsolatedElectronAfterMETCut['Data'].Histo1D(("leading_isolated_electron_eta_after_met_cut", "Leading Isolated Electron Pseudorapidity After MET Cut; Eta; Events",100,-3,3),"LeadingIsolatedElectron_eta")
        hist['Data']["leading_isolated_electron_phi_after_met_cut"] = rdfIsolatedElectronAfterMETCut['Data'].Histo1D(("leading_isolated_electron_phi_after_met_cut", "Leading Isolated Electron Angle After MET Cut; Phi (Radians); Events",100,-3.5,3.5),"LeadingIsolatedElectron_phi")
        hist['Data']["leading_isolated_electron_mass_after_met_cut"] = rdfIsolatedElectronAfterMETCut['Data'].Histo1D(("leading_isolated_electron_mass_after_met_cut", "Leading Isolated Electron Mass After MET Cut; Mass(Gev); Events",10,.1,.2),"LeadingIsolatedElectron_mass")
        hist['Data']["leading_isolated_electron_charge_after_met_cut"] = rdfIsolatedElectronAfterMETCut['Data'].Histo1D(("leading_isolated_electron_charge_after_met_cut", "Leading Isolated Electron Charge After MET Cut; Charge; Events",5,-2,2),"LeadingIsolatedElectron_charge")

        
        hist['Data']["number_of_jets"] = rdfJPsiMuons['Data'].Histo1D(("number_of_jets", "Number of Jets; Number Of Jets; Events", 20, 0, 20), "Num_Jets") 
        hist['Data']["jet1_pt"] = rdfJPsiMuons['Data'].Histo1D(("jet1_pt", "Jet Transverse Momentum for Leading Jet; Pt (GeV); Events", 100, 20, 200), "SJet1_pt")         
        hist['Data']["jet2_pt"] = rdfJPsiMuons['Data'].Histo1D(("jet2_pt", "Jet Transverse Momentum for Subleading Jet; Pt (GeV); Events", 100, 20, 200), "SJet2_pt")
        hist['Data']["jet1_eta"] = rdfJPsiMuons['Data'].Histo1D(("jet1_eta", "Jet Pseudorapidity for Leading Jet; Eta; Events", 100, -3, 3), "SJet1_eta")
        hist['Data']["jet2_eta"] = rdfJPsiMuons['Data'].Histo1D(("jet2_eta", "Jet Pseudorapidity for Subleading Jet; Eta; Events", 100, -3, 3), "SJet2_eta")
        hist['Data']["jet1_phi"] = rdfJPsiMuons['Data'].Histo1D(("jet1_phi", "Jet Angle for Leading Jet; Phi (Radians); Events", 100, -3.5, 3.5), "SJet1_phi")
        hist['Data']["jet2_phi"] = rdfJPsiMuons['Data'].Histo1D(("jet2_phi", "Jet Angle for Subleading Jet; Phi (Radians); Events", 100, -3.5, 3.5), "SJet2_phi")
        hist['Data']["transverse_mass"] = rdfJPsiMuons['Data'].Histo1D(("transverse_mass", "Transverse Mass; Transverse Mass (GeV); Events", 150, 0, 150), "MTofMETandMu")  
        hist['Data']["ht"] = rdfJPsiMuons['Data'].Histo1D(("ht", "Ht; Ht; Events", 300, 0, 1500), "Ht")
        
        
        hist['Data']["invariant_mass_jpsi_muons"] = rdfJPsiMuons['Data'].Histo1D(("invariant_mass_jpsi_muons", "Invariant Masses for J/Psi Candidate Muons (Oppositely Charged); Invariant Masses; Events", 50, .5, 12), "InvariantMassJPsiMuons")
        hist['Data']["invariant_mass_jpsi_muons"] = rdfJPsiMuons['Data'].Histo1D(("invariant_mass_jpsi_muons", "Invariant Masses for J/Psi Candidate Muons (Oppositely Charged); Invariant Masses; Events", 300, .5, 12), "InvariantMassJPsiMuons")
        hist['Data']["invariant_masses_zoomed"] = rdfJPsiMuons['Data'].Histo1D(("invariant_masses_zoomed", "Invariant Masses for Two Muons (Oppositely Charged); Invariant Masses; Events", 50, 2.8, 3.4), "InvariantMassJPsiMuons")
        hist['Data']["jpsi_muons_pt"] = rdfJPsiMuons['Data'].Histo1D(("Pt for JPsi Muons", "Transverse Momentum for JPsi Muons; Pt; Events", 150, 2.8, 75), "JPsiCandidate_pt")
        hist['Data']["jpsi_muons_eta"] = rdfJPsiMuons['Data'].Histo1D(("Eta for JPsi Muons", "Pseudorapidity for JPsi Muons; Eta; Events", 50, -3, 3), "JPsiCandidate_eta")
        hist['Data']["jpsi_muons_phi"] = rdfJPsiMuons['Data'].Histo1D(("Phi for JPsi Muons", "Angle for JPsi Muons; Phi; Events", 50, -3.5, 3.5), "JPsiCandidate_phi")
#        hist['Data']["jpsi_muons_multiplicity"] = rdfJPsiMuons['Data'].Histo1D(("jpsi_muon_multiplicity", "Number of J/Psi Muons; Number of J/Psi Muons; Events", 8, 0, 8), "JPsiMuons_multiplicity")
        hist['Data']["delta_eta"] = rdfJPsiMuons['Data'].Histo1D(("Delta Eta for Isolated Muon - JPsi Muons", "Delta Eta for Isolated Muon - JPsi Muons; Delta Eta; Events", 50, 0, 6), "DeltaEta")
        hist['Data']["delta_phi"] = rdfJPsiMuons['Data'].Histo1D(("Delta Phi for Isolated Muon - JPsi Muons", "Delta Phi for Isolated Muon - JPsi Muons; Delta Phi; Events", 50, -3.5, 3.5), "DeltaPhi")
        hist['Data']["delta_r"] = rdfJPsiMuons['Data'].Histo1D(("Delta R for Isolated and JPsi Muons", "Delta R for Isolated and JPsi Muons; Delta R; Events", 50, 0, 6), "DeltaR")
        
        

    else:
        
        hist[sample]["leading_isolated_muon_pt"] = rdfIsolatedMuonNoHighWeights[sample].Histo1D((sample + "_" + "leading_isolated_muon_pt", "Monte Carlo " + sample + ";Leading Isolated Muon Transverse Momentum; Pt (GeV);Events",100,20,220),"LeadingIsolatedMuon_pt","LumiXS")
        hist[sample]["leading_isolated_muon_eta"] = rdfIsolatedMuonNoHighWeights[sample].Histo1D((sample + "_" + "leading_isolated_muon_eta", "Monte Carlo " + sample + ";Leading Isolated Muon Pseudorapidity; Eta; Events",100,-3,3),"LeadingIsolatedMuon_eta","LumiXS")
        hist[sample]["leading_isolated_muon_phi"] = rdfIsolatedMuonNoHighWeights[sample].Histo1D((sample + "_" + "leading_isolated_muon_phi", "Monte Carlo " + sample + ";Leading Isolated Muon Angle; Phi (Radians); Events",100,-3.5,3.5),"LeadingIsolatedMuon_phi","LumiXS")
        hist[sample]["leading_isolated_muon_mass"] = rdfIsolatedMuonNoHighWeights[sample].Histo1D((sample + "_" + "leading_isolated_muon_mass", "Monte Carlo " + sample + ";Leading Isolated Muon Mass; Mass(Gev); Events",10,.1,.2),"LeadingIsolatedMuon_mass","LumiXS")
        hist[sample]["leading_isolated_muon_charge"] = rdfIsolatedMuonNoHighWeights[sample].Histo1D((sample + "_" + "leading_isolated_muon_charge", "Monte Carlo " + sample + ";Leading Isolated Muon Charge; Charge; Events",5,-2,2),"LeadingIsolatedMuon_charge","LumiXS")
        hist[sample]["leading_isolated_electron_pt"] = rdfIsolatedElectronNoHighWeights[sample].Histo1D((sample + "_" + "leading_isolated_electron_pt", "Monte Carlo " + sample + ";Leading Isolated Electron Transverse Momentum; Pt (GeV);Events",100,20,220),"LeadingIsolatedElectron_pt","LumiXS")
        hist[sample]["leading_isolated_electron_eta"] = rdfIsolatedElectronNoHighWeights[sample].Histo1D((sample + "_" + "leading_isolated_electron_eta", "Monte Carlo " + sample + ";Leading Isolated Electron Pseudorapidity; Eta; Events",100,-3,3),"LeadingIsolatedElectron_eta","LumiXS")
        hist[sample]["leading_isolated_electron_phi"] = rdfIsolatedElectronNoHighWeights[sample].Histo1D((sample + "_" + "leading_isolated_electron_phi", "Monte Carlo " + sample + ";Leading Isolated Electron Angle; Phi (Radians); Events",100,-3.5,3.5),"LeadingIsolatedElectron_phi","LumiXS")
        hist[sample]["leading_isolated_electron_mass"] = rdfIsolatedElectronNoHighWeights[sample].Histo1D((sample + "_" + "leading_isolated_electron_mass", "Monte Carlo " + sample + ";Leading Isolated Electron Mass; Mass(Gev); Events",10,.1,.2),"LeadingIsolatedElectron_mass","LumiXS")
        hist[sample]["leading_isolated_electron_charge"] = rdfIsolatedElectronNoHighWeights[sample].Histo1D((sample + "_" + "leading_isolated_electron_charge", "Monte Carlo " + sample + ";Leading Isolated Electron Charge; Charge; Events",5,-2,2),"LeadingIsolatedElectron_charge","LumiXS")
        hist[sample]["met_before_met_cut"] = rdfIsolatedMuonNoHighWeights[sample].Histo1D((sample + "_" + "met_before_met_cut", "Monte Carlo " + sample + "; MET Before MET Cut; Pt (GeV); Events",100,0,250), "METBeforeMETCut", "LumiXS")
        hist[sample]["met_after_met_cut"] = rdfIsolatedMuonAfterMETCut[sample].Histo1D((sample + "_" + "met_after_met_cut", "Monte Carlo " + sample + "; MET After MET Cut; Pt (GeV); Events",100,0,250), "METAfterMETCut", "LumiXS")
        hist[sample]["leading_isolated_muon_pt_after_met_cut"] = rdfJPsiMuons[sample].Histo1D((sample + "_" + "leading_isolated_muon_pt_after_met_cut", "Monte Carlo " + sample + ";Leading Isolated Muon Transverse Momentum After MET Cut; Pt (GeV);Events",100,20,220),"LeadingIsolatedMuon_pt","LumiXS")
        hist[sample]["leading_isolated_muon_eta_after_met_cut"] = rdfJPsiMuons[sample].Histo1D((sample + "_" + "leading_isolated_muon_eta_after_met_cut", "Monte Carlo " + sample + ";Leading Isolated Muon Pseudorapidity After MET Cut; Eta; Events",100,-3,3),"LeadingIsolatedMuon_eta","LumiXS")
        hist[sample]["leading_isolated_muon_phi_after_met_cut"] = rdfJPsiMuons[sample].Histo1D((sample + "_" + "leading_isolated_muon_phi_after_met_cut", "Monte Carlo " + sample + ";Leading Isolated Muon Angle After MET Cut; Phi (Radians); Events",100,-3.5,3.5),"LeadingIsolatedMuon_phi","LumiXS")
        hist[sample]["leading_isolated_muon_mass_after_met_cut"] = rdfJPsiMuons[sample].Histo1D((sample + "_" + "leading_isolated_muon_mass_after_met_cut", "Monte Carlo " + sample + ";Leading Isolated Muon Mass After MET Cut; Mass(Gev); Events",10,.1,.2),"LeadingIsolatedMuon_mass","LumiXS")
        hist[sample]["leading_isolated_muon_charge_after_met_cut"] = rdfJPsiMuons[sample].Histo1D((sample + "_" + "leading_isolated_muon_charge_after_met_cut", "Monte Carlo " + sample + ";Leading Isolated Muon Charge After MET Cut; Charge; Events",5,-2,2),"LeadingIsolatedMuon_charge","LumiXS")
        hist[sample]["leading_isolated_electron_pt_after_met_cut"] = rdfIsolatedElectronAfterMETCut[sample].Histo1D((sample + "_" + "leading_isolated_electron_pt_after_met_cut", "Monte Carlo " + sample + ";Leading Isolated Electron Transverse Momentum After MET Cut; Pt (GeV);Events",100,20,220),"LeadingIsolatedElectron_pt","LumiXS")
        hist[sample]["leading_isolated_electron_eta_after_met_cut"] = rdfIsolatedElectronAfterMETCut[sample].Histo1D((sample + "_" + "leading_isolated_electron_eta_after_met_cut", "Monte Carlo " + sample + ";Leading Isolated Electron Pseudorapidity After MET Cut; Eta; Events",100,-3,3),"LeadingIsolatedElectron_eta","LumiXS")
        hist[sample]["leading_isolated_electron_phi_after_met_cut"] = rdfIsolatedElectronAfterMETCut[sample].Histo1D((sample + "_" + "leading_isolated_electron_phi_after_met_cut", "Monte Carlo " + sample + ";Leading Isolated Electron Angle After MET Cut; Phi (Radians); Events",100,-3.5,3.5),"LeadingIsolatedElectron_phi","LumiXS")
        hist[sample]["leading_isolated_electron_mass_after_met_cut"] = rdfIsolatedElectronAfterMETCut[sample].Histo1D((sample + "_" + "leading_isolated_electron_mass_after_met_cut", "Monte Carlo " + sample + ";Leading Isolated Electron Mass After MET Cut; Mass(Gev); Events",10,.1,.2),"LeadingIsolatedElectron_mass","LumiXS")
        hist[sample]["leading_isolated_electron_charge_after_met_cut"] = rdfIsolatedElectronAfterMETCut[sample].Histo1D((sample + "_" + "leading_isolated_electron_charge_after_met_cut", "Monte Carlo " + sample + ";Leading Isolated Electron Charge After MET Cut; Charge; Events",5,-2,2),"LeadingIsolatedElectron_charge","LumiXS")
        
        
        hist[sample]["number_of_jets"] = rdfJPsiMuons[sample].Histo1D((sample + "_" + "number_of_jets", "Monte Carlo " + sample +"; Number Of Jets; Events", 20, 0, 20), "Num_Jets", "LumiXS")
        hist[sample]["jet1_pt"] = rdfJPsiMuons[sample].Histo1D((sample + "_" + "jet1_pt", "Monte Carlo " + sample + "; Pt (GeV); Events", 100, 20, 200), "SJet1_pt", "LumiXS")
        hist[sample]["jet2_pt"] = rdfJPsiMuons[sample].Histo1D((sample + "_" + "jet2_pt", "Monte Carlo " + sample + "; Pt (GeV); Events", 100, 20, 200), "SJet2_pt", "LumiXS")
        hist[sample]["jet1_eta"] = rdfJPsiMuons[sample].Histo1D((sample + "_" + "jet1_eta", "Monte Carlo " + sample + "; Eta; Events", 100, -3, 3), "SJet1_eta", "LumiXS")
        hist[sample]["jet2_eta"] = rdfJPsiMuons[sample].Histo1D((sample + "_" + "jet2_eta", "Monte Carlo " + sample +"; Eta; Events", 100, -3, 3), "SJet2_eta", "LumiXS")
        hist[sample]["jet1_phi"] = rdfJPsiMuons[sample].Histo1D((sample + "_" + "jet1_phi", "Monte Carlo " + sample + "; Phi (Radians); Events", 100, -3.5, 3.5), "SJet1_phi", "LumiXS")
        hist[sample]["jet2_phi"] = rdfJPsiMuons[sample].Histo1D((sample + "_" + "jet2_phi", "Monte Carlo " + sample +"; Phi (Radians); Events", 100, -3.5, 3.5), "SJet2_phi", "LumiXS")
        hist[sample]["transverse_mass"] = rdfJPsiMuons[sample].Histo1D((sample + "_" + "transverse_mass", "Monte Carlo " + sample +"; Transverse Mass (GeV); Events", 150, 0, 150), "MTofMETandMu", "LumiXS")
        hist[sample]["ht"] = rdfJPsiMuons[sample].Histo1D((sample + "_" + "ht", "Monte Carlo " + sample + "; Ht; Events", 300, 0, 1500), "Ht", "LumiXS")
        
        
        hist[sample]["invariant_mass_jpsi_muons"] = rdfJPsiMuons[sample].Histo1D((sample + "_" + "invariant_mass_jpsi_muons", "Monte Carlo " + sample + "; Invariant Masses for J/Psi Candidate Muons (Oppositely Charged); Invariant Masses; Events", 50, .5, 12), "InvariantMassJPsiMuons", "LumiXS")
        hist[sample]["invariant_mass_jpsi_muons"] = rdfJPsiMuons[sample].Histo1D((sample + "_" + "invariant_mass_jpsi_muons", "Monte Carlo " + sample + "; Invariant Masses for J/Psi Candidate Muons (Oppositely Charged); Invariant Masses; Events", 300, .5, 12), "InvariantMassJPsiMuons", "LumiXS")
        hist[sample]["invariant_masses_zoomed"] = rdfJPsiMuons[sample].Histo1D((sample + "_" + "invariant_masses_zoomed", "Monte Carlo " + sample + "; Invariant Masses for Two Muons (Oppositely Charged); Invariant Masses; Events", 50, 2.8, 3.4), "InvariantMassJPsiMuons", "LumiXS")
        hist[sample]["jpsi_muons_pt"] = rdfJPsiMuons[sample].Histo1D((sample + "_" + "jpsi_muons_pt", "Monte Carlo " + sample + "; Transverse Momentum for JPsi Muons; Pt; Events", 150, 2.8, 75), "JPsiCandidate_pt", "LumiXS")
        hist[sample]["jpsi_muons_eta"] = rdfJPsiMuons[sample].Histo1D((sample + "_" + "jpsi_muons_eta", "Monte Carlo " + sample + "; Pseudorapidity for JPsi Muons; Eta; Events", 50, -3, 3), "JPsiCandidate_eta", "LumiXS")
        hist[sample]["jpsi_muons_phi"] = rdfJPsiMuons[sample].Histo1D((sample + "_" + "jpsi_muons_phi", "Monte Carlo " + sample + "; Angle for JPsi Muons; Phi; Events", 50, -3.5, 3.5), "JPsiCandidate_phi", "LumiXS")
#        hist[sample]["jpsi_muons_multiplicity"] = rdfJPsiMuons[sample].Histo1D((sample + "_" + "jpsi_muon_multiplicity", "Monte Carlo" + sample + "; Number of J/Psi Muons; Number of J/Psi Muons; Events", 8, 0, 8), "JPsiMuons_multiplicity", "LumiXS")
        hist[sample]["delta_eta"] = rdfJPsiMuons[sample].Histo1D((sample + "_" + "delta_eta", "Monte Carlo " + sample + "; Delta Eta for Isolated Muon - JPsi Muons; Delta Eta; Events", 50, 0, 6), "DeltaEta", "LumiXS")
        hist[sample]["delta_phi"] = rdfJPsiMuons[sample].Histo1D((sample + "_" + "delta_phi", "Monte Carlo " + sample + "; Delta Phi for Isolated Muon - JPsi Muons; Delta Phi; Events", 50, -3.5, 3.5), "DeltaPhi", "LumiXS")
        hist[sample]["delta_r"] = rdfJPsiMuons[sample].Histo1D((sample + "_" + "delta_r", "Monte Carlo " + sample + "; Delta R for Isolated and JPsi Muons; Delta R; Events", 50, 0, 6), "DeltaR", "LumiXS")


# In[14]:


# def BTagSelection():
#     rdfmuOneBTaggedJet[sample] = rdfTwoPlusJets[sample].Filter("Num_BTaggedJets >= 1", "At Least One B-Tagged Jet")
#     # This fifth filter keeps only events with no vetoed muons
#     rdfmu[sample] = rdfmuOneBTaggedJet[sample].Filter("Sum(isolated_mu_veto) == 0", "No Vetoed Muons")


# In[15]:


# for sample in dictOfListOfFiles:
    
#     if sample not in hist.keys():
#         hist[sample] = {}
#         stats[sample] = {}
#         report[sample] = rdf[sample].Report()
        
#     #stats[sample]["countercode"] = rdf[sample].Min('my_rdfentry')
        
#     if sample == 'Data':
        
#         hist['Data']["isolated_mu_pt"] = rdfmu['Data'].Histo1D(("mu_pt","Muon Transverse Momentum; Pt (GeV);Events",100,20,220),"IsolatedMuon_pt")
#         hist['Data']["isolated_mu_eta"] = rdfmu['Data'].Histo1D(("mu_eta", "Muon Pseudorapidity; Eta; Events",100,-3,3),"IsolatedMuon_eta")
#         hist['Data']["isolated_mu_phi"] = rdfmu['Data'].Histo1D(("mu_phi", "Muon Angle; Phi (Radians); Events",100,-3.5,3.5),"IsolatedMuon_phi")
#         hist['Data']["jet1_pt"] = rdfmu['Data'].Histo1D(("jet1_pt", "Jet Transverse Momentum for Leading Jet; Pt (GeV); Events", 100, 20, 200), "SJet1_pt")
#         hist['Data']["jet2_pt"] = rdfmu['Data'].Histo1D(("jet2_pt", "Jet Transverse Momentum for Subleading Jet; Pt (GeV); Events", 100, 20, 200), "SJet2_pt")
#         hist['Data']["jet1_eta"] = rdfmu['Data'].Histo1D(("jet1_eta", "Jet Pseudorapidity for Leading Jet; Eta; Events", 100, -3, 3), "SJet1_eta")
#         hist['Data']["jet2_eta"] = rdfmu['Data'].Histo1D(("jet2_eta", "Jet Pseudorapidity for Subleading Jet; Eta; Events", 100, -3, 3), "SJet2_eta")
#         hist['Data']["jet1_phi"] = rdfmu['Data'].Histo1D(("jet1_phi", "Jet Angle for Leading Jet; Phi (Radians); Events", 100, -3.5, 3.5), "SJet1_phi")
#         hist['Data']["jet2_phi"] = rdfmu['Data'].Histo1D(("jet2_phi", "Jet Angle for Subleading Jet; Phi (Radians); Events", 100, -3.5, 3.5), "SJet2_phi")
#         hist['Data']["jet_deep"] = rdfmu['Data'].Histo1D(("jet_deep", "Deep Jet B Discriminator; Discriminant Value; Events", 100, 0, 1), "DeepJetB")
#         hist['Data']["number_of_jets"] = rdfmu['Data'].Histo1D(("number_of_jets", "Number of Jets; Number Of Jets; Events", 20, 0, 20), "Num_Jets") 
#         hist['Data']["number_of_muons"] = rdfmu['Data'].Histo1D(("number_of_muons", "Number of Muons; Number of Muons; Events", 5, 0, 5), "Num_Muons")
#         hist['Data']["transverse_mass"] = rdfmu['Data'].Histo1D(("transverse_mass", "Transverse Mass; Transverse Mass (GeV); Events", 150, 0, 150), "MTofMETandMu")
#         hist['Data']["missing_transverse_momentum"] = rdfmu['Data'].Histo1D(("missing_transverse_momentum", "Missing Transverse Momentum; Missing Transverse Momentum(GeV); Events", 150, 0, 300), "MET_pt")
#         hist['Data']["ht"] = rdfmu['Data'].Histo1D(("ht", "Ht; Ht; Events", 300, 0, 1500), "Ht")
#         hist['Data']["isolated_muon_pfRelIso03_all"] = rdfmu['Data'].Histo1D(("muon_pfRelIso03_all", "Muon Pf Rel Iso 03 (All); Muon Pf Rel Iso 03 (All); Events", 60, 0, .3), "IsolatedMuon_pfRelIso03_all")
#         hist['Data']["isolated_muon_pfRelIso03_chg"] = rdfmu['Data'].Histo1D(("muon_pfRelIso03_chg", "Muon Pf Rel Iso 03 (Chg); Muon Pf Rel Iso 03 (Chg); Events", 60, 0, .3), "IsolatedMuon_pfRelIso03_chg")
#         hist['Data']["isolated_muon_pfRelIso04_all"] = rdfmu['Data'].Histo1D(("muon_pfRelIso04_all", "Muon Pf Rel Iso 04 (All); Muon Pf Rel Iso 04 (All); Events", 60, 0, .3), "IsolatedMuon_pfRelIso04_all")
#         hist['Data']["invariant_masses"] = rdfOtherMuons['Data'].Histo1D(("invariant_masses", "Invariant Masses for Two Muons (Oppositely Charged); Invariant Masses; Events", 50, .5, 12), "InvariantMasses")        
#         hist['Data']["invariant_masses_zoomed"] = rdfOtherMuons['Data'].Histo1D(("invariant_masses_zoomed", "Invariant Masses for Two Muons (Oppositely Charged); Invariant Masses; Events", 50, 2.8, 3.4), "InvariantMasses")
#         hist['Data']["jpsi_muons_pt"] = rdfOtherMuons['Data'].Histo1D(("Pt for JPsi Muons", "Transverse Momentum for JPsi Muons; Pt; Events", 150, 2.8, 75), "JPsiMuons_pt")
#         hist['Data']["jpsi_muons_eta"] = rdfOtherMuons['Data'].Histo1D(("Eta for JPsi Muons", "Pseudorapidity for JPsi Muons; Eta; Events", 50, -3, 3), "JPsiMuons_eta")
#         hist['Data']["jpsi_muons_phi"] = rdfOtherMuons['Data'].Histo1D(("Phi for JPsi Muons", "Angle for JPsi Muons; Phi; Events", 50, -3.5, 3.5), "JPsiMuons_phi")
#         hist['Data']["jpsi_muons_multiplicity"] = rdfOtherMuons['Data'].Histo1D(("jpsi_muon_multiplicity", "Number of J/Psi Muons; Number of J/Psi Muons; Events", 8, 0, 8), "JPsiMuons_multiplicity")
#         hist['Data']["delta_eta"] = rdfOtherMuons['Data'].Histo1D(("Delta Eta for Isolated Muon - JPsi Muons", "Delta Eta for Isolated Muon - JPsi Muons; Delta Eta; Events", 50, 0, 6), "DeltaEta")
#         hist['Data']["delta_phi"] = rdfOtherMuons['Data'].Histo1D(("Delta Phi for Isolated Muon - JPsi Muons", "Delta Phi for Isolated Muon - JPsi Muons; Delta Phi; Events", 50, -3.5, 3.5), "DeltaPhi")
#         hist['Data']["delta_r"] = rdfOtherMuons['Data'].Histo1D(("Delta R for Isolated and JPsi Muons", "Delta R for Isolated and JPsi Muons; Delta R; Events", 50, 0, 6), "DeltaR")
#         hist['Data']["invariant_masses_all_muons"] = rdfOtherMuons['Data'].Histo1D(("invariant_masses_all_muons", "Invariant Masses for Three Muons (Isolated and Paired, Oppositely Charged); Invariant Masses; Events", 50, .5, 200), "InvariantMassesAllMuons")

#     else:
        
#         hist[sample]["isolated_mu_pt"] = rdfmu[sample].Histo1D((sample + "_" + "mu_pt","Monte Carlo " + sample + ";Pt (GeV);Events",100,20,220),"IsolatedMuon_pt","LumiXS")
#         hist[sample]["isolated_mu_eta"] = rdfmu[sample].Histo1D((sample + "_" + "mu_eta", "Monte Carlo " + sample + "; Eta; Events",100,-3,3),"IsolatedMuon_eta","LumiXS")
#         hist[sample]["isolated_mu_phi"] = rdfmu[sample].Histo1D((sample + "_" + "mu_phi", "Monte Carlo " + sample + "; Phi (Radians); Events",100,-3.5,3.5),"IsolatedMuon_phi","LumiXS")
#         hist[sample]["jet1_pt"] = rdfmu[sample].Histo1D((sample + "_" + "jet1_pt", "Monte Carlo " + sample + "; Pt (GeV); Events", 100, 20, 200), "SJet1_pt", "LumiXS")
#         hist[sample]["jet2_pt"] = rdfmu[sample].Histo1D((sample + "_" + "jet2_pt", "Monte Carlo " + sample + "; Pt (GeV); Events", 100, 20, 200), "SJet2_pt", "LumiXS")
#         hist[sample]["jet1_eta"] = rdfmu[sample].Histo1D((sample + "_" + "jet1_eta", "Monte Carlo " + sample + "; Eta; Events", 100, -3, 3), "SJet1_eta", "LumiXS")
#         hist[sample]["jet2_eta"] = rdfmu[sample].Histo1D((sample + "_" + "jet2_eta", "Monte Carlo " + sample +"; Eta; Events", 100, -3, 3), "SJet2_eta", "LumiXS")
#         hist[sample]["jet1_phi"] = rdfmu[sample].Histo1D((sample + "_" + "jet1_phi", "Monte Carlo " + sample + "; Phi (Radians); Events", 100, -3.5, 3.5), "SJet1_phi", "LumiXS")
#         hist[sample]["jet2_phi"] = rdfmu[sample].Histo1D((sample + "_" + "jet2_phi", "Monte Carlo " + sample +"; Phi (Radians); Events", 100, -3.5, 3.5), "SJet2_phi", "LumiXS")
#         hist[sample]["jet_deep"] = rdfmu[sample].Histo1D((sample + "_" + "jet_deep", "Monte Carlo " + sample +"; Discriminant Value; Events", 100, 0, 1), "DeepJetB", "LumiXS")
#         hist[sample]["number_of_jets"] = rdfmu[sample].Histo1D((sample + "_" + "number_of_jets", "Monte Carlo " + sample +"; Number Of Jets; Events", 20, 0, 20), "Num_Jets", "LumiXS") 
#         hist[sample]["number_of_muons"] = rdfmu[sample].Histo1D((sample + "_" + "number_of_muons", "Monte Carlo " + sample +"; Number of Muons; Events", 5, 0, 5), "Num_Muons", "LumiXS")
#         hist[sample]["transverse_mass"] = rdfmu[sample].Histo1D((sample + "_" + "transverse_mass", "Monte Carlo " + sample +"; Transverse Mass (GeV); Events", 150, 0, 150), "MTofMETandMu", "LumiXS")
#         hist[sample]["missing_transverse_momentum"] = rdfmu[sample].Histo1D((sample + "_" + "missing_transverse_momentum", "Monte Carlo" + sample + "; Missing Transverse Momentum(GeV); Events", 150, 0, 300), "MET_pt", "LumiXS")
#         hist[sample]["ht"] = rdfmu[sample].Histo1D((sample + "_" + "ht", "Monte Carlo " + sample + "; Ht; Events", 300, 0, 1500), "Ht", "LumiXS")
#         hist[sample]["isolated_muon_pfRelIso03_all"] = rdfmu[sample].Histo1D((sample + "_" + "muon_pfRelIso03_all", "Monte Carlo " + sample + "; Muon Pf Rel Iso 03 (All); Events", 60, 0, .3), "IsolatedMuon_pfRelIso03_all", "LumiXS")
#         hist[sample]["isolated_muon_pfRelIso03_chg"] = rdfmu[sample].Histo1D((sample + "_" + "muon_pfRelIso03_chg", "Monte Carlo " + sample + "; Muon Pf Rel Iso 03 (Chg); Events", 60, 0, .3), "IsolatedMuon_pfRelIso03_chg", "LumiXS")
#         hist[sample]["isolated_muon_pfRelIso04_all"] = rdfmu[sample].Histo1D((sample + "_" + "muon_pfRelIso04_all", "Monte Carlo " + sample + "; Muon Pf Rel Iso 04 (All); Events", 60, 0, .3), "IsolatedMuon_pfRelIso04_all", "LumiXS")
#         hist[sample]["invariant_masses"] = rdfOtherMuons[sample].Histo1D((sample + "_" + "invariant_masses", "Monte Carlo " + sample + "; Invariant Masses for Two Muons (Oppositely Charged); Invariant Masses; Events", 50, .5, 12), "InvariantMasses", "LumiXS")
#         hist[sample]["invariant_masses_zoomed"] = rdfOtherMuons[sample].Histo1D((sample + "_" + "invariant_masses_zoomed", "Monte Carlo " + sample + "; Invariant Masses for Two Muons (Oppositely Charged); Events", 50, 2.8, 3.4), "InvariantMasses", "LumiXS")
#         hist[sample]["jpsi_muons_pt"] = rdfOtherMuons[sample].Histo1D((sample + "_" + "jpsi_muons_pt", "Monte Carlo " + sample + "; Transverse Momentum for JPsi Muons; Pt; Events", 150, 2.8, 75), "JPsiMuons_pt", "LumiXS")
#         hist[sample]["jpsi_muons_eta"] = rdfOtherMuons[sample].Histo1D((sample + "_" + "jpsi_muons_eta", "Monte Carlo " + sample + "; Pseudorapidity for JPsi Muons; Eta; Events", 50, -3, 3), "JPsiMuons_eta", "LumiXS")
#         hist[sample]["jpsi_muons_phi"] = rdfOtherMuons[sample].Histo1D((sample + "_" + "jpsi_muons_phi", "Monte Carlo " + sample + "; Angle for JPsi Muons; Phi; Events", 50, -3.5, 3.5), "JPsiMuons_phi", "LumiXS")
#         hist[sample]["jpsi_muons_multiplicity"] = rdfOtherMuons[sample].Histo1D((sample + "_" + "jpsi_muon_multiplicity", "Monte Carlo" + sample + "; Number of J/Psi Muons; Number of J/Psi Muons; Events", 8, 0, 8), "JPsiMuons_multiplicity", "LumiXS")
#         hist[sample]["delta_eta"] = rdfOtherMuons[sample].Histo1D((sample + "_" + "delta_eta", "Monte Carlo " + sample + "; Delta Eta for Isolated Muon - JPsi Muons; Delta Eta; Events", 50, 0, 6), "DeltaEta", "LumiXS")
#         hist[sample]["delta_phi"] = rdfOtherMuons[sample].Histo1D((sample + "_" + "delta_phi", "Monte Carlo " + sample + "; Delta Phi for Isolated Muon - JPsi Muons; Delta Phi; Events", 50, -3.5, 3.5), "DeltaPhi", "LumiXS")
#         hist[sample]["delta_r"] = rdfOtherMuons[sample].Histo1D((sample + "_" + "delta_r", "Monte Carlo " + sample + "; Delta R for Isolated and JPsi Muons; Delta R; Events", 50, 0, 6), "DeltaR", "LumiXS")
#         hist[sample]["invariant_masses_all_muons"] = rdfOtherMuons[sample].Histo1D((sample + "_" + "invariant_masses_all_muons", "Monte Carlo" + sample + "; Invariant Masses for Three Muons (Isolated and Paired, Oppositely Charged); Invariant Masses; Events", 50, .5, 200), "InvariantMassesAllMuons", "LumiXS")


# In[16]:


for sample in dictOfListOfFiles:

    ROOT.RDF.SaveGraph(rdf[sample], str(sample) + ".dot")
    os.system("dot -Tpdf " + str(sample) + ".dot > " + str(sample) + "_graph.pdf" )


# In[ ]:


for sample in dictOfListOfFiles:
    #nparray[sample] = nparraynode[sample].AsNumpy(["genWeight", "event", "run"])
    
    cutflow = report[sample].GetValue()
    cutflow.Print()
    
    era = "2018"
    process = sample
    channel = "Mu"
    syst = "nominal"

    outFile = ROOT.TFile.Open("{}_{}_{}.root".format(era, channel, process), "RECREATE")
    for name, hist_pointer in hist[sample].items():
        print(hist_pointer)
        hist_value = hist_pointer.GetValue()
        
        hist_value.SetName("{}___{}___{}___{}___{}".format(era, channel, process, name, syst))
        print(hist_value)
        hist_value.Write()

    outFile.Close()


# In[ ]:


for item in mureport:
    print(item)
    mureport[item].Print()
#print(nparray['MonteCarloWJetsToLNu'].keys())

#a = np.sort(nparray['MonteCarloWJetsToLNu']['genWeight'])

#print(a[:10:-1])


# In[ ]:




