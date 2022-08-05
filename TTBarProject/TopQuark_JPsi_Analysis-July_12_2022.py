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

rdfJetAndIsolatedLeptonFiltered = {}

rdfJPsiMuons = {}

rdfJPsisInPeak = {}

rdfRemainingIsolatedMuonMuonIdMediumIsoIdTight = {}
rdfRemainingIsolatedMuonMuonIdLooseIsoIdTight = {}
rdfRemainingIsolatedMuonMuonIdMediumIsoIdMedium = {}
rdfRemainingIsolatedMuonMuonIdLooseIsoIdMedium = {}

hist = {}
report = {}


nparray = {}
nparraynode = {}

# Muon_pfIsoId is PFIso ID from miniAOD selector (1=PFIsoVeryLoose, 2=PFIsoLoose, 3=PFIsoMedium, 4=PFIsoTight, 5=PFIsoVeryTight, 6=PFIsoVeryVeryTight)
LeadingIsolatedMuonMask = "Muon_pt > 30 && abs(Muon_eta) < 2.4 && Muon_mediumId == true && Muon_pfIsoId >= 4"
LeadingIsolatedElectronMask = "Electron_pt > 30 && abs(Electron_eta) < 2.4 && Electron_cutBased == 4"
JetMask = "ROOT::VecOps::RVec<Int_t> jmask = (Jet_pt >= 30 && abs(Jet_eta) <= 2.5 && Jet_jetId >= 2); "                          "for(int i=0; i < LeadingIsolatedMuon_pt.size(); ++i){"                              "ROOT::VecOps::RVec<Float_t> dr;"                              "for(int j=0; j < jmask.size(); ++j){"                                  "dr.push_back(ROOT::VecOps::DeltaR(Jet_eta.at(j), LeadingIsolatedMuon_eta.at(i), Jet_phi.at(j), LeadingIsolatedMuon_phi.at(i)));}"                                  "jmask = jmask && dr >= 0.4;"                                  "dr.clear();}"                          "return jmask;"
JPsiCandidateMask = "Muon_pt > 6 && abs(Muon_eta) <= 2.4 && Muon_mediumId == true && !(leading_isolated_muon_mask)"

# Diagnostic Masks
# LeadingIsolatedMuonMaskPfIsoMedium = "Muon_pt > 30 && abs(Muon_eta) < 2.4 && Muon_mediumId == true && Muon_pfIsoId >= 4 && !(Muon_pt > 30 && abs(Muon_eta) < 2.4 && Muon_tightId == true && Muon_pfIsoId >= 4)"

# LeadingIsolatedMuonMaskMuonIdMediumIsoTight = "Muon_pt > 30 && abs(Muon_eta) < 2.4 && Muon_mediumId == true && Muon_pfIsoId >= 4 && !(Muon_pt > 30 && abs(Muon_eta) < 2.4 && Muon_tightId == true && Muon_pfIsoId >= 4)"
# LeadingIsolatedMuonMaskMuonIdLooseIsoIdTight = "Muon_pt > 30 && abs(Muon_eta) < 2.4 && Muon_looseId == true && Muon_pfIsoId >= 4 && !(Muon_pt > 30 && abs(Muon_eta) < 2.4 && Muon_tightId == true && Muon_pfIsoId >= 4)"
# LeadingIsolatedMuonMaskMuonIdMediumIsoIdMedium = "Muon_pt > 30 && abs(Muon_eta) < 2.4 && Muon_mediumId == true && Muon_pfIsoId >= 3 && !(Muon_pt > 30 && abs(Muon_eta) < 2.4 && Muon_tightId == true && Muon_pfIsoId >= 4)"
# LeadingIsolatedMuonMaskMuonIdLooseIsoIdMedium = "Muon_pt > 30 && abs(Muon_eta) < 2.4 && Muon_looseId == true && Muon_pfIsoId >= 3 && !(Muon_pt > 30 && abs(Muon_eta) < 2.4 && Muon_tightId == true && Muon_pfIsoId >= 4)"


# JPsiMuonPositiveCandidateMask = "Muon_pt > 3 && abs(Muon_eta) <= 2.4 && Muon_mediumId == true && !(Muon_pt > 30 && abs(Muon_eta) < 2.4 && Muon_tightId == true && Muon_pfIsoId >= 4) && Muon_charge == 1"
# JPsiMuonNegativeCandidateMask = "Muon_pt > 3 && abs(Muon_eta) <= 2.4 && Muon_mediumId == true && !(Muon_pt > 30 && abs(Muon_eta) < 2.4 && Muon_tightId == true && Muon_pfIsoId >= 4) && Muon_charge == -1"
# JPsiMuonPositiveCandidatePeakOnlyMask = "Muon_pt > 3.0 && Muon_pt < 3.2 && abs(Muon_eta) <= 2.4 && Muon_mediumId == true && !(Muon_pt > 30 && abs(Muon_eta) < 2.4 && Muon_tightId == true && Muon_pfIsoId >= 4) && Muon_charge == 1"
# JPsiMuonNegativeCandidatePeakOnlyMask = "Muon_pt > 3.0 && Muon_pt < 3.2 && abs(Muon_eta) <= 2.4 && Muon_mediumId == true && !(Muon_pt > 30 && abs(Muon_eta) < 2.4 && Muon_tightId == true && Muon_pfIsoId >= 4) && Muon_charge == -1"
# JPsiMuonPositiveCandidatePfIsoMediumMask = "Muon_pt > 3 && abs(Muon_eta) <= 2.4 && Muon_mediumId == true && !(Muon_pt > 30 && abs(Muon_eta) < 2.4 && Muon_tightId == true && Muon_pfIsoId >= 4) && Muon_charge == 1"
# JPsiMuonNegativeCandidatePfIsoMediumMask = "Muon_pt > 3 && abs(Muon_eta) <= 2.4 && Muon_mediumId == true && !(Muon_pt > 30 && abs(Muon_eta) < 2.4 && Muon_tightId == true && Muon_pfIsoId >= 4) && Muon_charge == -1"


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
typedef ROOT::VecOps::RVec<std::tuple<Float_t, Float_t>>   RVec_ff;

typedef ROOT::VecOps::RVec<ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double>>> RVec_FourVector;

class MuonAndJPsiStatisticsAndKinematics
{
    private:
        RVec_f JPsi_Muon_pt;
        RVec_f JPsi_Muon_eta;
        RVec_f JPsi_Muon_phi;
        RVec_f JPsi_Muon_mass;
        RVec_i JPsi_Muon_charge;
        RVec_f Isolated_Muon_pt;
        RVec_f Isolated_Muon_eta;
        RVec_f Isolated_Muon_phi;
        RVec_f Isolated_Muon_mass;
        RVec_i Isolated_Muon_charge;
        RVec_f Isolated_Muon_PfRelIso03_all;
        RVec_f JPsi_Muon_PfRelIso03_all;


    public:
        MuonAndJPsiStatisticsAndKinematics(RVec_f JPsi_Muon_pt, RVec_f JPsi_Muon_eta, RVec_f JPsi_Muon_phi, RVec_f JPsi_Muon_mass, RVec_i JPsi_Muon_charge);
        MuonAndJPsiStatisticsAndKinematics(RVec_f JPsi_Muon_pt, RVec_f JPsi_Muon_eta, RVec_f JPsi_Muon_phi, RVec_f JPsi_Muon_mass, RVec_i JPsi_Muon_charge, RVec_f Isolated_Muon_pt, RVec_f Isolated_Muon_eta, RVec_f Isolated_Muon_phi, RVec_f Isolated_Muon_mass, RVec_i Isolated_Muon_charge);
        MuonAndJPsiStatisticsAndKinematics(RVec_f JPsi_Muon_pt, RVec_f JPsi_Muon_eta, RVec_f JPsi_Muon_phi, RVec_f JPsi_Muon_mass, RVec_i JPsi_Muon_charge, RVec_f Isolated_Muon_pt, RVec_f Isolated_Muon_eta, RVec_f Isolated_Muon_phi, RVec_f Isolated_Muon_mass, RVec_i Isolated_Muon_charge, RVec_f Isolated_Muon_PfRelIso03_all, RVec_f JPsi_Muon_PfRelIso03_all);
        RVec_FourVector JPsiFourVectorCalculator();
        RVec_FourVector JPsiForMuonsInPeakFourVectorCalculator();
        RVec_FourVector TopQuarkInPeakFourVectorCalculator();
        RVec_f ReturnJPsiPt();
        RVec_f ReturnJPsiEta();
        RVec_f ReturnJPsiPhi();
        RVec_f JPsiMuonInvariantMassCalculator();
        RVec_f IsolatedAndJPsiMuonInvariantMassCalculator();
        RVec_f ReturnMuonPfRelIso03AllForIsolatedMuonInInvariantMass();
        RVec_f ReturnMuonPfRelIso03AllForJPsiMuonsInInvariantMass();
        RVec_f ReturnJPsiMassforJPsiInInvariantMass();
        RVec_f DeltaEtaBetweenIsolatedAndJPsiMuonCalculator();
        RVec_f DeltaPhiBetweenIsolatedAndJPsiMuonCalculator();
        RVec_f DeltaRBetweenIsolatedAndJPsiMuonCalculator();
        RVec_f DeltaEtaBetweenJPsiMuonsCalculator();
        RVec_f DeltaPhiBetweenJPsiMuonsCalculator();
        RVec_f DeltaRBetweenJPsiMuonsCalculator();
        RVec_f ReturnDeltaEtaBetweenJPsiMuonsInPeak();
        RVec_f ReturnDeltaPhiBetweenJPsiMuonsInPeak();
        RVec_f ReturnDeltaRBetweenJPsiMuonsInPeak();
        RVec_f DeltaEtaBetweenIsolatedMuonAndJPsiCalculator();
        RVec_f DeltaPhiBetweenIsolatedMuonAndJPsiCalculator();
        RVec_f DeltaRBetweenIsolatedMuonAndJPsiCalculator();
        RVec_f ReturnLowDeltaRIsolatedMuonAndJPsiMuonPt(bool);
        RVec_f ReturnJPsiInPeakMass();
        RVec_f ReturnJPsiMuonsInPeakCharge();
        RVec_f ReturnJPsiInPeakPt();
        RVec_f ReturnJPsiMuonPtInPeak(bool);
        RVec_f ReturnJPsiMuonInPeakPfRelIso03All();
        RVec_f ReturnTopQuarkInPeakMass();
        RVec_f DeltaEtaBetweenIsolatedMuonAndJPsiInPeakCalculator();
        RVec_f DeltaPhiBetweenIsolatedMuonAndJPsiInPeakCalculator();
        RVec_f DeltaRBetweenIsolatedMuonAndJPsiInPeakCalculator();

};

MuonAndJPsiStatisticsAndKinematics::MuonAndJPsiStatisticsAndKinematics(RVec_f JPsi_Muon_pt, RVec_f JPsi_Muon_eta, RVec_f JPsi_Muon_phi, RVec_f JPsi_Muon_mass, RVec_i JPsi_Muon_charge)
{
    this->JPsi_Muon_pt = JPsi_Muon_pt;
    this->JPsi_Muon_eta = JPsi_Muon_eta;
    this->JPsi_Muon_phi = JPsi_Muon_phi;
    this->JPsi_Muon_mass = JPsi_Muon_mass;
    this->JPsi_Muon_charge = JPsi_Muon_charge;
    this->Isolated_Muon_pt = {};
    this->Isolated_Muon_eta = {};
    this->Isolated_Muon_phi = {};
    this->Isolated_Muon_mass = {};
    this->Isolated_Muon_charge = {};
    this->Isolated_Muon_PfRelIso03_all = {};
    this->JPsi_Muon_PfRelIso03_all = {};
}

MuonAndJPsiStatisticsAndKinematics::MuonAndJPsiStatisticsAndKinematics(RVec_f JPsi_Muon_pt, RVec_f JPsi_Muon_eta, RVec_f JPsi_Muon_phi, RVec_f JPsi_Muon_mass, RVec_i JPsi_Muon_charge, RVec_f Isolated_Muon_pt, RVec_f Isolated_Muon_eta, RVec_f Isolated_Muon_phi, RVec_f Isolated_Muon_mass, RVec_i Isolated_Muon_charge)
{
    this->JPsi_Muon_pt = JPsi_Muon_pt;
    this->JPsi_Muon_eta = JPsi_Muon_eta;
    this->JPsi_Muon_phi = JPsi_Muon_phi;
    this->JPsi_Muon_mass = JPsi_Muon_mass;
    this->JPsi_Muon_charge = JPsi_Muon_charge;
    this->Isolated_Muon_pt = Isolated_Muon_pt;
    this->Isolated_Muon_eta = Isolated_Muon_eta;
    this->Isolated_Muon_phi = Isolated_Muon_phi;
    this->Isolated_Muon_mass = Isolated_Muon_mass;
    this->Isolated_Muon_charge = Isolated_Muon_charge;
    this->Isolated_Muon_PfRelIso03_all = {};
    this->JPsi_Muon_PfRelIso03_all = {};
}

MuonAndJPsiStatisticsAndKinematics::MuonAndJPsiStatisticsAndKinematics(RVec_f JPsi_Muon_pt, RVec_f JPsi_Muon_eta, RVec_f JPsi_Muon_phi, RVec_f JPsi_Muon_mass, RVec_i JPsi_Muon_charge, RVec_f Isolated_Muon_pt, RVec_f Isolated_Muon_eta, RVec_f Isolated_Muon_phi, RVec_f Isolated_Muon_mass, RVec_i Isolated_Muon_charge, RVec_f Isolated_Muon_PfRelIso03_all, RVec_f JPsi_Muon_PfRelIso03_all)
{
    this->JPsi_Muon_pt = JPsi_Muon_pt;
    this->JPsi_Muon_eta = JPsi_Muon_eta;
    this->JPsi_Muon_phi = JPsi_Muon_phi;
    this->JPsi_Muon_mass = JPsi_Muon_mass;
    this->JPsi_Muon_charge = JPsi_Muon_charge;
    this->Isolated_Muon_pt = Isolated_Muon_pt;
    this->Isolated_Muon_eta = Isolated_Muon_eta;
    this->Isolated_Muon_phi = Isolated_Muon_phi;
    this->Isolated_Muon_mass = Isolated_Muon_mass;
    this->Isolated_Muon_charge = Isolated_Muon_charge;
    this->Isolated_Muon_PfRelIso03_all = Isolated_Muon_PfRelIso03_all;
    this->JPsi_Muon_PfRelIso03_all = JPsi_Muon_PfRelIso03_all;
}



RVec_FourVector MuonAndJPsiStatisticsAndKinematics::JPsiFourVectorCalculator()
{

    RVec_FourVector fvt = {};
    
    int FirstMuonCharge = 0;
    int SecondMuonCharge = 0;
    
    /* Loop over the set of muons to determine which muons have +1 charge, then match them with all the ones with -1 charge. */
    for(int i = 0; i < this->JPsi_Muon_charge.size(); i++)
    {    

        FirstMuonCharge = this->JPsi_Muon_charge[i];
        
        /* If charges are opposite, calculate the invariant mass of them */
        for(int j = i+1; j < this->JPsi_Muon_charge.size(); j++)
        {
            SecondMuonCharge = this->JPsi_Muon_charge[j];
            
            if(FirstMuonCharge * SecondMuonCharge == -1)
            {
                auto FourVector1 = ROOT::Math::PtEtaPhiMVector (this->JPsi_Muon_pt[i], this->JPsi_Muon_eta[i], this->JPsi_Muon_phi[i], this->JPsi_Muon_mass[i]);
                
                auto FourVector2 = ROOT::Math::PtEtaPhiMVector (this->JPsi_Muon_pt[j], this->JPsi_Muon_eta[j], this->JPsi_Muon_phi[j], this->JPsi_Muon_mass[j]);
                
                fvt.push_back(FourVector1 + FourVector2);
                    
                FourVector1 = {};
                FourVector2 = {};
            }
        }
    }
    
    return fvt;   
}

RVec_FourVector MuonAndJPsiStatisticsAndKinematics::JPsiForMuonsInPeakFourVectorCalculator()
{

    RVec_FourVector fvt = {};
    
    int FirstMuonCharge = 0;
    int SecondMuonCharge = 0;
    
    RVec_f pt {};
    RVec_f eta {};
    RVec_f phi {};
    RVec_f mass {};
    
    /* Loop over the set of muons to determine which muons have +1 charge, then match them with all the ones with -1 charge. */
    for(int i = 0; i < this->JPsi_Muon_charge.size(); i++)
    {    

        FirstMuonCharge = this->JPsi_Muon_charge[i];
        
        /* If charges are opposite, calculate the invariant mass of them */
        for(int j = i+1; j < this->JPsi_Muon_charge.size(); j++)
        {
            SecondMuonCharge = this->JPsi_Muon_charge[j];
            
            pt.push_back(this->JPsi_Muon_pt[i]);
            eta.push_back(this->JPsi_Muon_eta[i]);
            phi.push_back(this->JPsi_Muon_phi[i]);
            mass.push_back(this->JPsi_Muon_mass[i]);

            pt.push_back(this->JPsi_Muon_pt[j]);
            eta.push_back(this->JPsi_Muon_eta[j]);
            phi.push_back(this->JPsi_Muon_phi[j]);
            mass.push_back(this->JPsi_Muon_mass[j]);
            
            if((FirstMuonCharge * SecondMuonCharge == -1) && (ROOT::VecOps::InvariantMass(pt, eta, phi, mass) > 3.0) && (ROOT::VecOps::InvariantMass(pt, eta, phi, mass) < 3.2))
            {
                auto FourVector1 = ROOT::Math::PtEtaPhiMVector (this->JPsi_Muon_pt[i], this->JPsi_Muon_eta[i], this->JPsi_Muon_phi[i], this->JPsi_Muon_mass[i]);
                
                auto FourVector2 = ROOT::Math::PtEtaPhiMVector (this->JPsi_Muon_pt[j], this->JPsi_Muon_eta[j], this->JPsi_Muon_phi[j], this->JPsi_Muon_mass[j]);
                
                fvt.push_back(FourVector1 + FourVector2);
                    
                FourVector1 = {};
                FourVector2 = {};
            }
            
            pt.clear();
            eta.clear();
            phi.clear();
            mass.clear();
        }
    }
    
    return fvt;   
}

RVec_FourVector MuonAndJPsiStatisticsAndKinematics::TopQuarkInPeakFourVectorCalculator()
{

    RVec_FourVector fvt = {};
    
    int FirstMuonCharge = 0;
    int SecondMuonCharge = 0;
    
    RVec_f pt {};
    RVec_f eta {};
    RVec_f phi {};
    RVec_f mass {};
    
    for(int k = 0; k < this->Isolated_Muon_pt.size(); k++)
    {
        /* Loop over the set of muons to determine which muons have +1 charge, then match them with all the ones with -1 charge. */
        for(int i = 0; i < this->JPsi_Muon_charge.size(); i++)
        {    

            FirstMuonCharge = this->JPsi_Muon_charge[i];

            /* If charges are opposite, calculate the invariant mass of them */
            for(int j = i+1; j < this->JPsi_Muon_charge.size(); j++)
            {
                SecondMuonCharge = this->JPsi_Muon_charge[j];

                pt.push_back(this->JPsi_Muon_pt[i]);
                eta.push_back(this->JPsi_Muon_eta[i]);
                phi.push_back(this->JPsi_Muon_phi[i]);
                mass.push_back(this->JPsi_Muon_mass[i]);

                pt.push_back(this->JPsi_Muon_pt[j]);
                eta.push_back(this->JPsi_Muon_eta[j]);
                phi.push_back(this->JPsi_Muon_phi[j]);
                mass.push_back(this->JPsi_Muon_mass[j]);

                if((FirstMuonCharge * SecondMuonCharge == -1) && (ROOT::VecOps::InvariantMass(pt, eta, phi, mass) > 3.0) && (ROOT::VecOps::InvariantMass(pt, eta, phi, mass) < 3.2))
                {
                    auto FourVector1 = ROOT::Math::PtEtaPhiMVector (this->JPsi_Muon_pt[i], this->JPsi_Muon_eta[i], this->JPsi_Muon_phi[i], this->JPsi_Muon_mass[i]);

                    auto FourVector2 = ROOT::Math::PtEtaPhiMVector (this->JPsi_Muon_pt[j], this->JPsi_Muon_eta[j], this->JPsi_Muon_phi[j], this->JPsi_Muon_mass[j]);
                    
                    auto FourVector3 = ROOT::Math::PtEtaPhiMVector (this->Isolated_Muon_pt[k], this->Isolated_Muon_eta[k], this->Isolated_Muon_phi[k], this->Isolated_Muon_mass[k]);

                    fvt.push_back(FourVector1 + FourVector2 + FourVector3);

                    FourVector1 = {};
                    FourVector2 = {};
                    FourVector3 = {};
                }

                pt.clear();
                eta.clear();
                phi.clear();
                mass.clear();
            }
        }
    }
    
    return fvt;   
}



RVec_f MuonAndJPsiStatisticsAndKinematics::ReturnJPsiPt()
{
    RVec_f JPsiPt = {};
    
    RVec_FourVector FourVectorTotal = JPsiFourVectorCalculator();
    
    for(int i = 0; i < FourVectorTotal.size(); i++)
    {
        JPsiPt.push_back(FourVectorTotal.at(i).Pt());
    }
    
    return JPsiPt;
}

RVec_f MuonAndJPsiStatisticsAndKinematics::ReturnJPsiEta()
{
    RVec_f JPsiEta = {};
    
    RVec_FourVector FourVectorTotal = JPsiFourVectorCalculator();
    
    for(int i = 0; i < FourVectorTotal.size(); i++)
    {
        JPsiEta.push_back(FourVectorTotal.at(i).Eta());
    }
    
    return JPsiEta;
}

RVec_f MuonAndJPsiStatisticsAndKinematics::ReturnJPsiPhi()
{
    RVec_f JPsiPhi = {};
    
    RVec_FourVector FourVectorTotal = JPsiFourVectorCalculator();
    
    for(int i = 0; i < FourVectorTotal.size(); i++)
    {
        JPsiPhi.push_back(FourVectorTotal.at(i).Phi());
    }
    
    return JPsiPhi;
}

/* This function matches each muon with oppositely charged muons. */
RVec_f MuonAndJPsiStatisticsAndKinematics::JPsiMuonInvariantMassCalculator()
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
    for(int i = 0; i < this->JPsi_Muon_charge.size(); i++)
    {    

        FirstMuonCharge = this->JPsi_Muon_charge[i];
        
        /* If charges are opposite, calculate the invariant mass of them */
        for(int j = i+1; j < this->JPsi_Muon_charge.size(); j++)
        {
            SecondMuonCharge = this->JPsi_Muon_charge[j];
            
            if(FirstMuonCharge * SecondMuonCharge == -1)
            {
                pt.push_back(this->JPsi_Muon_pt[i]);
                eta.push_back(this->JPsi_Muon_eta[i]);
                phi.push_back(this->JPsi_Muon_phi[i]);
                mass.push_back(this->JPsi_Muon_mass[i]);
                    
                pt.push_back(this->JPsi_Muon_pt[j]);
                eta.push_back(this->JPsi_Muon_eta[j]);
                phi.push_back(this->JPsi_Muon_phi[j]);
                mass.push_back(this->JPsi_Muon_mass[j]);
                    
                im = ROOT::VecOps::InvariantMass(pt, eta, phi, mass);
                
                InvariantMasses.push_back(im);
                    
                pt.clear();
                eta.clear();
                phi.clear();
                mass.clear();
            }
        }
    }
        
    return InvariantMasses;
}


RVec_f MuonAndJPsiStatisticsAndKinematics::IsolatedAndJPsiMuonInvariantMassCalculator()
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
        for(int i = 0; i < this->JPsi_Muon_charge.size(); i++)
        {    

            FirstMuonCharge = this->JPsi_Muon_charge[i];

            /* If charges are opposite, calculate the invariant mass of them */
            for(int j = i+1; j < this->JPsi_Muon_charge.size(); j++)
            {
                SecondMuonCharge = this->JPsi_Muon_charge[j];

                if(FirstMuonCharge * SecondMuonCharge == -1)
                {
                    pt.push_back(this->JPsi_Muon_pt[i]);
                    eta.push_back(this->JPsi_Muon_eta[i]);
                    phi.push_back(this->JPsi_Muon_phi[i]);
                    mass.push_back(this->JPsi_Muon_mass[i]);

                    pt.push_back(this->JPsi_Muon_pt[j]);
                    eta.push_back(this->JPsi_Muon_eta[j]);
                    phi.push_back(this->JPsi_Muon_phi[j]);
                    mass.push_back(this->JPsi_Muon_mass[j]);
                    
                    pt.push_back(this->Isolated_Muon_pt[k]);
                    eta.push_back(this->Isolated_Muon_eta[k]);
                    phi.push_back(this->Isolated_Muon_phi[k]);
                    mass.push_back(this->Isolated_Muon_mass[k]);

                    im = ROOT::VecOps::InvariantMass(pt, eta, phi, mass);

                    InvariantMasses.push_back(im);

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





RVec_f MuonAndJPsiStatisticsAndKinematics::ReturnMuonPfRelIso03AllForIsolatedMuonInInvariantMass()
{
    int FirstMuonCharge = 0;
    int SecondMuonCharge = 0;
    
    float im = 0;
    
    RVec_f pt {};
    RVec_f eta {};
    RVec_f phi {};
    RVec_f mass {};
    
    RVec_f Isolations {};
    
    bool IsolatedMuonAlreadyAccountedFor = false;
 
    for(int k = 0; k < this->Isolated_Muon_pt.size(); k++)
    {
    
        /* Loop over the set of muons to determine which muons have +1 charge, then match them with all the ones with -1 charge. */
        for(int i = 0; i < this->JPsi_Muon_charge.size(); i++)
        {    

            FirstMuonCharge = this->JPsi_Muon_charge[i];

            /* If charges are opposite, calculate the invariant mass of them */
            for(int j = i+1; j < this->JPsi_Muon_charge.size(); j++)
            {
                SecondMuonCharge = this->JPsi_Muon_charge[j];

                if(FirstMuonCharge * SecondMuonCharge == -1)
                {
                
                    pt.push_back(this->JPsi_Muon_pt[i]);
                    eta.push_back(this->JPsi_Muon_eta[i]);
                    phi.push_back(this->JPsi_Muon_phi[i]);
                    mass.push_back(this->JPsi_Muon_mass[i]);

                    pt.push_back(this->JPsi_Muon_pt[j]);
                    eta.push_back(this->JPsi_Muon_eta[j]);
                    phi.push_back(this->JPsi_Muon_phi[j]);
                    mass.push_back(this->JPsi_Muon_mass[j]);
                    
                    pt.push_back(this->Isolated_Muon_pt[k]);
                    eta.push_back(this->Isolated_Muon_eta[k]);
                    phi.push_back(this->Isolated_Muon_phi[k]);
                    mass.push_back(this->Isolated_Muon_mass[k]);

                    im = ROOT::VecOps::InvariantMass(pt, eta, phi, mass);
                    
                    if((im >= 90) && (im <= 120))
                    {
                    
                        if(!IsolatedMuonAlreadyAccountedFor)
                        {
                            Isolations.push_back(this->Isolated_Muon_PfRelIso03_all[k]);
                            
                            IsolatedMuonAlreadyAccountedFor = true;
                        }
                    }
                    
                    pt.clear();
                    eta.clear();
                    phi.clear();
                    mass.clear();
                   
                }
            }
        }
        
        IsolatedMuonAlreadyAccountedFor = false;
    }
        
    return Isolations;
}






RVec_f MuonAndJPsiStatisticsAndKinematics::ReturnMuonPfRelIso03AllForJPsiMuonsInInvariantMass()
{
    int FirstMuonCharge = 0;
    int SecondMuonCharge = 0;
    
    float im = 0;
    
    RVec_f pt {};
    RVec_f eta {};
    RVec_f phi {};
    RVec_f mass {};
    
    RVec_f Isolations {};
 
    for(int k = 0; k < this->Isolated_Muon_pt.size(); k++)
    {
    
        /* Loop over the set of muons to determine which muons have +1 charge, then match them with all the ones with -1 charge. */
        for(int i = 0; i < this->JPsi_Muon_charge.size(); i++)
        {    

            FirstMuonCharge = this->JPsi_Muon_charge[i];

            /* If charges are opposite, calculate the invariant mass of them */
            for(int j = i+1; j < this->JPsi_Muon_charge.size(); j++)
            {
                SecondMuonCharge = this->JPsi_Muon_charge[j];

                if(FirstMuonCharge * SecondMuonCharge == -1)
                {
                
                    pt.push_back(this->JPsi_Muon_pt[i]);
                    eta.push_back(this->JPsi_Muon_eta[i]);
                    phi.push_back(this->JPsi_Muon_phi[i]);
                    mass.push_back(this->JPsi_Muon_mass[i]);

                    pt.push_back(this->JPsi_Muon_pt[j]);
                    eta.push_back(this->JPsi_Muon_eta[j]);
                    phi.push_back(this->JPsi_Muon_phi[j]);
                    mass.push_back(this->JPsi_Muon_mass[j]);
                    
                    pt.push_back(this->Isolated_Muon_pt[k]);
                    eta.push_back(this->Isolated_Muon_eta[k]);
                    phi.push_back(this->Isolated_Muon_phi[k]);
                    mass.push_back(this->Isolated_Muon_mass[k]);

                    im = ROOT::VecOps::InvariantMass(pt, eta, phi, mass);
                    
                    if((im >= 90) && (im <= 120))
                    {
                        Isolations.push_back(this->JPsi_Muon_PfRelIso03_all[i]);
                        Isolations.push_back(this->JPsi_Muon_PfRelIso03_all[j]);
                    }
                    
                    pt.clear();
                    eta.clear();
                    phi.clear();
                    mass.clear();
                   
                }
            }
        }
    }
        
    return Isolations;
}








RVec_f MuonAndJPsiStatisticsAndKinematics::ReturnJPsiMassforJPsiInInvariantMass()
{
    int FirstMuonCharge = 0;
    int SecondMuonCharge = 0;
    
    float imAll = 0;
    float imJPsi = 0;
    
    RVec_f ptAll {};
    RVec_f etaAll {};
    RVec_f phiAll {};
    RVec_f massAll {};
    
    RVec_f JPsiMasses {};
 
    for(int k = 0; k < this->Isolated_Muon_pt.size(); k++)
    {
    
        /* Loop over the set of muons to determine which muons have +1 charge, then match them with all the ones with -1 charge. */
        for(int i = 0; i < this->JPsi_Muon_charge.size(); i++)
        {    

            FirstMuonCharge = this->JPsi_Muon_charge[i];

            /* If charges are opposite, calculate the invariant mass of them */
            for(int j = i+1; j < this->JPsi_Muon_charge.size(); j++)
            {
                SecondMuonCharge = this->JPsi_Muon_charge[j];

                if(FirstMuonCharge * SecondMuonCharge == -1)
                {
                
                    ptAll.push_back(this->JPsi_Muon_pt[i]);
                    etaAll.push_back(this->JPsi_Muon_eta[i]);
                    phiAll.push_back(this->JPsi_Muon_phi[i]);
                    massAll.push_back(this->JPsi_Muon_mass[i]);

                    ptAll.push_back(this->JPsi_Muon_pt[j]);
                    etaAll.push_back(this->JPsi_Muon_eta[j]);
                    phiAll.push_back(this->JPsi_Muon_phi[j]);
                    massAll.push_back(this->JPsi_Muon_mass[j]);
                    
                    ptAll.push_back(this->Isolated_Muon_pt[k]);
                    etaAll.push_back(this->Isolated_Muon_eta[k]);
                    phiAll.push_back(this->Isolated_Muon_phi[k]);
                    massAll.push_back(this->Isolated_Muon_mass[k]);
                    
                    
                    imAll = ROOT::VecOps::InvariantMass(ptAll, etaAll, phiAll, massAll);
                    
                    if((imAll >= 90) && (imAll <= 120))
                    {
                    
                        auto FourVector1 = ROOT::Math::PtEtaPhiMVector (this->JPsi_Muon_pt[i], this->JPsi_Muon_eta[i], this->JPsi_Muon_phi[i], this->JPsi_Muon_mass[i]);

                        auto FourVector2 = ROOT::Math::PtEtaPhiMVector (this->JPsi_Muon_pt[j], this->JPsi_Muon_eta[j], this->JPsi_Muon_phi[j], this->JPsi_Muon_mass[j]);

                        auto FourVector3 = FourVector1 + FourVector2;
                        
                        JPsiMasses.push_back(FourVector3.M());

                        FourVector1 = {};
                        FourVector2 = {};
                    }
                    
                    ptAll.clear();
                    etaAll.clear();
                    phiAll.clear();
                    massAll.clear();
                   
                }
            }
        }
    }
        
    return JPsiMasses;
}











RVec_f MuonAndJPsiStatisticsAndKinematics::DeltaEtaBetweenIsolatedAndJPsiMuonCalculator()
{ 
    float DeltaEtaIndividual = 0;
    
    RVec_f DeltaEtaRVec {};
    
    for(int i = 0; i < this->Isolated_Muon_eta.size(); i++)
    {
    
        for(int j = 0; j < this->JPsi_Muon_eta.size(); j++)
        {    
                    DeltaEtaIndividual = this->Isolated_Muon_eta[i] - this->JPsi_Muon_eta[j];

                    DeltaEtaRVec.push_back(DeltaEtaIndividual);
        }
    }
        
    return DeltaEtaRVec;
}


RVec_f MuonAndJPsiStatisticsAndKinematics::DeltaPhiBetweenIsolatedAndJPsiMuonCalculator()
{ 
    float PhiIsolated = 0;
    float PhiJPsi = 0;
    
    float DeltaPhiIndividual = 0;
    
    RVec_f DeltaPhiRVec {};
    
    for(int i = 0; i < this->Isolated_Muon_phi.size(); i++)
    {
    
        for(int j = 0; j < this->JPsi_Muon_phi.size(); j++)
        {    

                    PhiIsolated = this->Isolated_Muon_phi[i];

                    PhiJPsi = this->JPsi_Muon_phi[j];

                    DeltaPhiIndividual = ROOT::VecOps::DeltaPhi(PhiIsolated, PhiJPsi);

                    DeltaPhiRVec.push_back(DeltaPhiIndividual);
        }
    }
        
    return DeltaPhiRVec;
}


RVec_f MuonAndJPsiStatisticsAndKinematics::DeltaRBetweenIsolatedAndJPsiMuonCalculator()
{ 
    float EtaIsolated = 0;
    float PhiIsolated = 0;
    float EtaJPsi = 0;
    float PhiJPsi = 0;
    
    float DeltaRIndividual = 0;
    
    RVec_f DeltaRRVec {};
    
    for(int i = 0; i < this->Isolated_Muon_phi.size(); i++)
    {
    
        for(int j = 0; j < this->JPsi_Muon_phi.size(); j++)
        {    
                    EtaIsolated = this->Isolated_Muon_eta[i];
                    PhiIsolated = this->Isolated_Muon_phi[i];
                    
                    EtaJPsi = this->JPsi_Muon_eta[j];
                    PhiJPsi = this->JPsi_Muon_phi[j];

                    DeltaRIndividual = ROOT::VecOps::DeltaR(EtaIsolated, EtaJPsi, PhiIsolated, PhiJPsi);

                    DeltaRRVec.push_back(DeltaRIndividual);
        }
    }
        
    return DeltaRRVec;
}


RVec_f MuonAndJPsiStatisticsAndKinematics::DeltaEtaBetweenJPsiMuonsCalculator()
{
    int FirstMuonCharge = 0;
    int SecondMuonCharge = 0;
    
    float DeltaEtaIndividual = 0;
    
    RVec_f DeltaEtaRVec {};
    
    /* Loop over the set of muons to determine which muons have +1 charge, then match them with all the ones with -1 charge. */
    for(int i = 0; i < this->JPsi_Muon_charge.size(); i++)
    {    
        FirstMuonCharge = this->JPsi_Muon_charge[i];
        
        /* If charges are opposite, calculate the invariant mass of them */
        for(int j = i+1; j < this->JPsi_Muon_charge.size(); j++)
        {
            SecondMuonCharge = this->JPsi_Muon_charge[j];
            
            if(FirstMuonCharge * SecondMuonCharge == -1)
            {                
                DeltaEtaIndividual = this->JPsi_Muon_eta[i] - this->JPsi_Muon_eta[j];

                DeltaEtaRVec.push_back(DeltaEtaIndividual);
            }
        }
    }
        
    return DeltaEtaRVec;
}


RVec_f MuonAndJPsiStatisticsAndKinematics::DeltaPhiBetweenJPsiMuonsCalculator()
{ 
    int FirstMuonCharge = 0;
    int SecondMuonCharge = 0;
    
    float PhiJPsi1 = 0;
    float PhiJPsi2 = 0;
    
    float DeltaPhiIndividual = 0;
    
    RVec_f DeltaPhiRVec {};
    
    for(int i = 0; i < this->JPsi_Muon_phi.size(); i++)
    {
        FirstMuonCharge = this->JPsi_Muon_charge[i];
        
        for(int j = i+1; j < this->JPsi_Muon_phi.size(); j++)
        {
            SecondMuonCharge = this->JPsi_Muon_charge[j];
            
            if(FirstMuonCharge * SecondMuonCharge == -1)
            {                
                PhiJPsi1 = this->JPsi_Muon_phi[i];

                PhiJPsi2 = this->JPsi_Muon_phi[j];

                DeltaPhiIndividual = ROOT::VecOps::DeltaPhi(PhiJPsi1, PhiJPsi2);

                DeltaPhiRVec.push_back(DeltaPhiIndividual);
            }
        }
    }
        
    return DeltaPhiRVec;
}


RVec_f MuonAndJPsiStatisticsAndKinematics::DeltaRBetweenJPsiMuonsCalculator()
{ 
    int FirstMuonCharge = 0;
    int SecondMuonCharge = 0;
    
    float EtaJPsi1 = 0;
    float PhiJPsi1 = 0;
    float EtaJPsi2 = 0;
    float PhiJPsi2 = 0;
    
    float DeltaRIndividual = 0;
    
    RVec_f DeltaRRVec {};
    
    for(int i = 0; i < this->JPsi_Muon_phi.size(); i++)
    {
        FirstMuonCharge = this->JPsi_Muon_charge[i];

            for(int j = i+1; j < this->JPsi_Muon_phi.size(); j++)
            {
                SecondMuonCharge = this->JPsi_Muon_charge[j];

                if(FirstMuonCharge * SecondMuonCharge == -1)
                {                
                    EtaJPsi1 = this->JPsi_Muon_eta[i];
                    PhiJPsi1 = this->JPsi_Muon_phi[i];
                    
                    EtaJPsi2 = this->JPsi_Muon_eta[j];
                    PhiJPsi2 = this->JPsi_Muon_phi[j];

                    DeltaRIndividual = ROOT::VecOps::DeltaR(EtaJPsi1, EtaJPsi2, PhiJPsi1, PhiJPsi2);

                    DeltaRRVec.push_back(DeltaRIndividual);
                }
            }
        }
        
    return DeltaRRVec;
}











RVec_f MuonAndJPsiStatisticsAndKinematics::ReturnDeltaEtaBetweenJPsiMuonsInPeak()
{
    int FirstMuonCharge = 0;
    int SecondMuonCharge = 0;
    
    float DeltaEtaIndividual = 0;
    
    RVec_f pt {};
    RVec_f eta {};
    RVec_f phi {};
    RVec_f mass {};
    
    RVec_f DeltaEtaRVec {};
    
    /* Loop over the set of muons to determine which muons have +1 charge, then match them with all the ones with -1 charge. */
    for(int i = 0; i < this->JPsi_Muon_charge.size(); i++)
    {    
        FirstMuonCharge = this->JPsi_Muon_charge[i];
        
        /* If charges are opposite, calculate the invariant mass of them */
        for(int j = i+1; j < this->JPsi_Muon_charge.size(); j++)
        {
            SecondMuonCharge = this->JPsi_Muon_charge[j];
            
            if(FirstMuonCharge * SecondMuonCharge == -1)
            {
                pt.push_back(this->JPsi_Muon_pt[i]);
                eta.push_back(this->JPsi_Muon_eta[i]);
                phi.push_back(this->JPsi_Muon_phi[i]);
                mass.push_back(this->JPsi_Muon_mass[i]);

                pt.push_back(this->JPsi_Muon_pt[j]);
                eta.push_back(this->JPsi_Muon_eta[j]);
                phi.push_back(this->JPsi_Muon_phi[j]);
                mass.push_back(this->JPsi_Muon_mass[j]);
                
                if((ROOT::VecOps::InvariantMass(pt, eta, phi, mass) > 3.0) && (ROOT::VecOps::InvariantMass(pt, eta, phi, mass) < 3.2))
                {            
                    DeltaEtaIndividual = this->JPsi_Muon_eta[i] - this->JPsi_Muon_eta[j];

                    DeltaEtaRVec.push_back(DeltaEtaIndividual);
                }
            }
            
            pt.clear();
            eta.clear();
            phi.clear();
            mass.clear();
        }
    }
        
    return DeltaEtaRVec;
}


RVec_f MuonAndJPsiStatisticsAndKinematics::ReturnDeltaPhiBetweenJPsiMuonsInPeak()
{ 
    int FirstMuonCharge = 0;
    int SecondMuonCharge = 0;
    
    float PhiJPsi1 = 0;
    float PhiJPsi2 = 0;
    
    float DeltaPhiIndividual = 0;
    
    RVec_f pt {};
    RVec_f eta {};
    RVec_f phi {};
    RVec_f mass {};
    
    RVec_f DeltaPhiRVec {};
    
    for(int i = 0; i < this->JPsi_Muon_phi.size(); i++)
    {
        FirstMuonCharge = this->JPsi_Muon_charge[i];
        
        for(int j = i+1; j < this->JPsi_Muon_phi.size(); j++)
        {
            SecondMuonCharge = this->JPsi_Muon_charge[j];
            
            if(FirstMuonCharge * SecondMuonCharge == -1)
            {
                pt.push_back(this->JPsi_Muon_pt[i]);
                eta.push_back(this->JPsi_Muon_eta[i]);
                phi.push_back(this->JPsi_Muon_phi[i]);
                mass.push_back(this->JPsi_Muon_mass[i]);

                pt.push_back(this->JPsi_Muon_pt[j]);
                eta.push_back(this->JPsi_Muon_eta[j]);
                phi.push_back(this->JPsi_Muon_phi[j]);
                mass.push_back(this->JPsi_Muon_mass[j]);
                
                if((ROOT::VecOps::InvariantMass(pt, eta, phi, mass) > 3.0) && (ROOT::VecOps::InvariantMass(pt, eta, phi, mass) < 3.2))
                {
                    PhiJPsi1 = this->JPsi_Muon_phi[i];

                    PhiJPsi2 = this->JPsi_Muon_phi[j];

                    DeltaPhiIndividual = ROOT::VecOps::DeltaPhi(PhiJPsi1, PhiJPsi2);

                    DeltaPhiRVec.push_back(DeltaPhiIndividual);
                }
            }
            
            pt.clear();
            eta.clear();
            phi.clear();
            mass.clear();
        }
    }
        
    return DeltaPhiRVec;
}


RVec_f MuonAndJPsiStatisticsAndKinematics::ReturnDeltaRBetweenJPsiMuonsInPeak()
{ 
    int FirstMuonCharge = 0;
    int SecondMuonCharge = 0;
    
    float EtaJPsi1 = 0;
    float PhiJPsi1 = 0;
    float EtaJPsi2 = 0;
    float PhiJPsi2 = 0;
    
    float DeltaRIndividual = 0;
    
    RVec_f pt {};
    RVec_f eta {};
    RVec_f phi {};
    RVec_f mass {};
    
    RVec_f DeltaRRVec {};
    
    for(int i = 0; i < this->JPsi_Muon_phi.size(); i++)
    {
        FirstMuonCharge = this->JPsi_Muon_charge[i];

            for(int j = i+1; j < this->JPsi_Muon_phi.size(); j++)
            {
                SecondMuonCharge = this->JPsi_Muon_charge[j];

                if(FirstMuonCharge * SecondMuonCharge == -1)
                {                
                    EtaJPsi1 = this->JPsi_Muon_eta[i];
                    PhiJPsi1 = this->JPsi_Muon_phi[i];
                    
                    EtaJPsi2 = this->JPsi_Muon_eta[j];
                    PhiJPsi2 = this->JPsi_Muon_phi[j];

                    DeltaRIndividual = ROOT::VecOps::DeltaR(EtaJPsi1, EtaJPsi2, PhiJPsi1, PhiJPsi2);
                    
                    pt.push_back(this->JPsi_Muon_pt[i]);
                    eta.push_back(this->JPsi_Muon_eta[i]);
                    phi.push_back(this->JPsi_Muon_phi[i]);
                    mass.push_back(this->JPsi_Muon_mass[i]);

                    pt.push_back(this->JPsi_Muon_pt[j]);
                    eta.push_back(this->JPsi_Muon_eta[j]);
                    phi.push_back(this->JPsi_Muon_phi[j]);
                    mass.push_back(this->JPsi_Muon_mass[j]);

                    if((ROOT::VecOps::InvariantMass(pt, eta, phi, mass) > 3.0) && (ROOT::VecOps::InvariantMass(pt, eta, phi, mass) < 3.2))
                    {
                        DeltaRRVec.push_back(DeltaRIndividual);
                    }
                }
                
                pt.clear();
                eta.clear();
                phi.clear();
                mass.clear();
            }
        }
        
    return DeltaRRVec;
}






RVec_f MuonAndJPsiStatisticsAndKinematics::DeltaEtaBetweenIsolatedMuonAndJPsiCalculator()
{
    RVec_FourVector FourVectorTotal = JPsiFourVectorCalculator();
    
    float DeltaEtaIndividual = 0;
    
    RVec_f DeltaEtaRVec {};
    
    for(int i = 0; i < this->Isolated_Muon_eta.size(); i++)
    {
    
        for(int j = 0; j < FourVectorTotal.size(); j++)
        {    
                    DeltaEtaIndividual = this->Isolated_Muon_eta[i] - FourVectorTotal.at(j).Eta();

                    DeltaEtaRVec.push_back(DeltaEtaIndividual);
        }
    }
        
    return DeltaEtaRVec;
}

RVec_f MuonAndJPsiStatisticsAndKinematics::DeltaPhiBetweenIsolatedMuonAndJPsiCalculator()
{
    RVec_FourVector FourVectorTotal = JPsiFourVectorCalculator();
    
    float PhiIsolated = 0;
    float PhiJPsi = 0;
    
    float DeltaPhiIndividual = 0;
    
    RVec_f DeltaPhiRVec {};
    
    for(int i = 0; i < this->Isolated_Muon_phi.size(); i++)
    {
    
        for(int j = 0; j < FourVectorTotal.size(); j++)
        {    
                    PhiIsolated = this->Isolated_Muon_phi[i];

                    PhiJPsi = FourVectorTotal.at(j).Phi();

                    DeltaPhiIndividual = ROOT::VecOps::DeltaPhi(PhiIsolated, PhiJPsi);

                    DeltaPhiRVec.push_back(DeltaPhiIndividual);
        }
    }
        
    return DeltaPhiRVec;   
}

RVec_f MuonAndJPsiStatisticsAndKinematics::DeltaRBetweenIsolatedMuonAndJPsiCalculator()
{
    RVec_FourVector FourVectorTotal = JPsiFourVectorCalculator();
    
    float EtaIsolated = 0;
    float PhiIsolated = 0;
    float EtaJPsi = 0;
    float PhiJPsi = 0;
    
    float DeltaRIndividual = 0;
    
    RVec_f DeltaRRVec {};
    
    for(int i = 0; i < this->Isolated_Muon_phi.size(); i++)
    {
    
        for(int j = 0; j < FourVectorTotal.size(); j++)
        {    
                    EtaIsolated = this->Isolated_Muon_eta[i];
                    PhiIsolated = this->Isolated_Muon_phi[i];
                    
                    EtaJPsi = FourVectorTotal.at(j).Eta();
                    PhiJPsi = FourVectorTotal.at(j).Phi();

                    DeltaRIndividual = ROOT::VecOps::DeltaR(EtaIsolated, EtaJPsi, PhiIsolated, PhiJPsi);
                    
                    DeltaRRVec.push_back(DeltaRIndividual);
        }
    }
        
    return DeltaRRVec;
}







RVec_f MuonAndJPsiStatisticsAndKinematics::ReturnLowDeltaRIsolatedMuonAndJPsiMuonPt(bool ReturnJPsiMuon)
{
    float EtaIsolated = 0;
    float PhiIsolated = 0;
    float EtaJPsiMuon = 0;
    float PhiJPsiMuon = 0;
    
    float DeltaRIndividual = 0;
    
    RVec_ff IsolatedMuonAndJPsiMuonPt {};
    
    RVec_f IsolatedMuonPt {};
    RVec_f JPsiMuonPt {};
    
    for(int i = 0; i < this->Isolated_Muon_pt.size(); i++)
    {    

        for(int j = 0; j < this->JPsi_Muon_pt.size(); j++)
        {
                    EtaIsolated = this->Isolated_Muon_eta[i];
                    PhiIsolated = this->Isolated_Muon_phi[i];
                    
                    EtaJPsiMuon = this->JPsi_Muon_eta[j];
                    PhiJPsiMuon = this->JPsi_Muon_phi[j];

                    DeltaRIndividual = ROOT::VecOps::DeltaR(EtaIsolated, EtaJPsiMuon, PhiIsolated, PhiJPsiMuon);

                    if(DeltaRIndividual < 0.3)
                    {
                        IsolatedMuonAndJPsiMuonPt.push_back(std::tuple(this->Isolated_Muon_pt[i], this->JPsi_Muon_pt[j]));
                    }
        }
    }
    
    if(!ReturnJPsiMuon)
    {
        for(int m = 0; m < IsolatedMuonAndJPsiMuonPt.size(); m++)
        {
            IsolatedMuonPt.push_back(std::get<0>(IsolatedMuonAndJPsiMuonPt[m]));
        }
        
        return IsolatedMuonPt;
    }
    
    else
    {
        for(int m = 0; m < IsolatedMuonAndJPsiMuonPt.size(); m++)
        {
            JPsiMuonPt.push_back(std::get<1>(IsolatedMuonAndJPsiMuonPt[m]));
        }
        
        return JPsiMuonPt;
    }
}





RVec_f MuonAndJPsiStatisticsAndKinematics::ReturnJPsiInPeakMass()
{
    RVec_f JPsiInPeakMass = {};
    
    RVec_FourVector FourVectorTotal = JPsiForMuonsInPeakFourVectorCalculator();
    
    for(int i = 0; i < FourVectorTotal.size(); i++)
    {
        JPsiInPeakMass.push_back(FourVectorTotal.at(i).M());
    }
    
    return JPsiInPeakMass;
}

RVec_f MuonAndJPsiStatisticsAndKinematics::ReturnJPsiMuonsInPeakCharge()
{
    int FirstMuonCharge = 0;
    int SecondMuonCharge = 0;
    
    RVec_f pt {};
    RVec_f eta {};
    RVec_f phi {};
    RVec_f mass {};
    
    RVec_f JPsiMuonChargeRVec {};
    
    for(int i = 0; i < this->JPsi_Muon_phi.size(); i++)
    {
        FirstMuonCharge = this->JPsi_Muon_charge[i];

            for(int j = i+1; j < this->JPsi_Muon_phi.size(); j++)
            {
                SecondMuonCharge = this->JPsi_Muon_charge[j];

                if(FirstMuonCharge * SecondMuonCharge == -1)
                {             
                    pt.push_back(this->JPsi_Muon_pt[i]);
                    eta.push_back(this->JPsi_Muon_eta[i]);
                    phi.push_back(this->JPsi_Muon_phi[i]);
                    mass.push_back(this->JPsi_Muon_mass[i]);

                    pt.push_back(this->JPsi_Muon_pt[j]);
                    eta.push_back(this->JPsi_Muon_eta[j]);
                    phi.push_back(this->JPsi_Muon_phi[j]);
                    mass.push_back(this->JPsi_Muon_mass[j]);

                    if((ROOT::VecOps::InvariantMass(pt, eta, phi, mass) > 3.0) && (ROOT::VecOps::InvariantMass(pt, eta, phi, mass) < 3.2))
                    {
                        JPsiMuonChargeRVec.push_back(FirstMuonCharge);
                        JPsiMuonChargeRVec.push_back(SecondMuonCharge);
                    }
                }
                
                pt.clear();
                eta.clear();
                phi.clear();
                mass.clear();
            }
    }
    
    return JPsiMuonChargeRVec;
}

RVec_f MuonAndJPsiStatisticsAndKinematics::ReturnJPsiInPeakPt()
{
    RVec_f JPsiPt = {};
    
    RVec_FourVector FourVectorTotal = JPsiForMuonsInPeakFourVectorCalculator();
    
    for(int i = 0; i < FourVectorTotal.size(); i++)
    {
        JPsiPt.push_back(FourVectorTotal.at(i).Pt());
    }
    
    return JPsiPt;
}





RVec_f MuonAndJPsiStatisticsAndKinematics::ReturnJPsiMuonPtInPeak(bool ReturnPositiveMuon)
{
    float EtaIsolated = 0;
    float PhiIsolated = 0;
    float EtaJPsiMuon = 0;
    float PhiJPsiMuon = 0;
    
    float DeltaRIndividual = 0;
    
    RVec_ff IsolatedMuonAndJPsiMuonPt {};
    
    RVec_f JPsiMuonPositivePt {};
    RVec_f JPsiMuonNegativePt {};
    
    RVec_f pt {};
    RVec_f eta {};
    RVec_f phi {};
    RVec_f mass {};
    
    int FirstMuonCharge = 0;
    int SecondMuonCharge = 0;
    
    for(int i = 0; i < this->JPsi_Muon_charge.size(); i++)
    {
        FirstMuonCharge = this->JPsi_Muon_charge[i];

        for(int j = i+1; j < this->JPsi_Muon_charge.size(); j++)
        {
                SecondMuonCharge = this->JPsi_Muon_charge[j];
                
                pt.push_back(this->JPsi_Muon_pt[i]);
                eta.push_back(this->JPsi_Muon_eta[i]);
                phi.push_back(this->JPsi_Muon_phi[i]);
                mass.push_back(this->JPsi_Muon_mass[i]);
                    
                pt.push_back(this->JPsi_Muon_pt[j]);
                eta.push_back(this->JPsi_Muon_eta[j]);
                phi.push_back(this->JPsi_Muon_phi[j]);
                mass.push_back(this->JPsi_Muon_mass[j]);

                if((ROOT::VecOps::InvariantMass(pt, eta, phi, mass) > 3.0) && (ROOT::VecOps::InvariantMass(pt, eta, phi, mass) < 3.2))
                {
                    if(FirstMuonCharge < 0)
                    {
                        JPsiMuonPositivePt.push_back(this->JPsi_Muon_pt[i]);
                        JPsiMuonNegativePt.push_back(this->JPsi_Muon_pt[j]);
                    }
                    
                    else
                    {
                        JPsiMuonPositivePt.push_back(this->JPsi_Muon_pt[j]);
                        JPsiMuonNegativePt.push_back(this->JPsi_Muon_pt[i]);
                    }
                }
                
                pt.clear();
                eta.clear();
                phi.clear();
                mass.clear();
        }
    }
    
    if(ReturnPositiveMuon)
    {        
        return JPsiMuonPositivePt;
    }
    
    else
    {        
        return JPsiMuonNegativePt;
    }
}









RVec_f MuonAndJPsiStatisticsAndKinematics::ReturnJPsiMuonInPeakPfRelIso03All()
{
    int FirstMuonCharge = 0;
    int SecondMuonCharge = 0;
    
    float im = 0;
    
    RVec_f pt {};
    RVec_f eta {};
    RVec_f phi {};
    RVec_f mass {};
    
    RVec_f Isolations {};
 
        /* Loop over the set of muons to determine which muons have +1 charge, then match them with all the ones with -1 charge. */
        for(int i = 0; i < this->JPsi_Muon_charge.size(); i++)
        {    

            FirstMuonCharge = this->JPsi_Muon_charge[i];

            /* If charges are opposite, calculate the invariant mass of them */
            for(int j = i+1; j < this->JPsi_Muon_charge.size(); j++)
            {
                SecondMuonCharge = this->JPsi_Muon_charge[j];

                if(FirstMuonCharge * SecondMuonCharge == -1)
                {
                
                    pt.push_back(this->JPsi_Muon_pt[i]);
                    eta.push_back(this->JPsi_Muon_eta[i]);
                    phi.push_back(this->JPsi_Muon_phi[i]);
                    mass.push_back(this->JPsi_Muon_mass[i]);

                    pt.push_back(this->JPsi_Muon_pt[j]);
                    eta.push_back(this->JPsi_Muon_eta[j]);
                    phi.push_back(this->JPsi_Muon_phi[j]);
                    mass.push_back(this->JPsi_Muon_mass[j]);

                    im = ROOT::VecOps::InvariantMass(pt, eta, phi, mass);
                    
                    if((im >= 3.0) && (im <= 3.2))
                    {
                        Isolations.push_back(this->JPsi_Muon_PfRelIso03_all[i]);
                        Isolations.push_back(this->JPsi_Muon_PfRelIso03_all[j]);
                    }
                    
                    pt.clear();
                    eta.clear();
                    phi.clear();
                    mass.clear();
                   
                }
            }
        }
        
    return Isolations;
}

RVec_f MuonAndJPsiStatisticsAndKinematics::ReturnTopQuarkInPeakMass()
{
    RVec_f TopQuarkInPeakMass = {};
    
    RVec_FourVector FourVectorTotal = TopQuarkInPeakFourVectorCalculator();
    
    for(int i = 0; i < FourVectorTotal.size(); i++)
    {
        TopQuarkInPeakMass.push_back(FourVectorTotal.at(i).M());
    }
    
    return TopQuarkInPeakMass;
}



RVec_f MuonAndJPsiStatisticsAndKinematics::DeltaEtaBetweenIsolatedMuonAndJPsiInPeakCalculator()
{
    RVec_FourVector FourVectorTotal = JPsiForMuonsInPeakFourVectorCalculator();
    
    float DeltaEtaIndividual = 0;
    
    RVec_f DeltaEtaRVec {};
    
    for(int i = 0; i < this->Isolated_Muon_eta.size(); i++)
    {
    
        for(int j = 0; j < FourVectorTotal.size(); j++)
        {    
                    DeltaEtaIndividual = this->Isolated_Muon_eta[i] - FourVectorTotal.at(j).Eta();

                    DeltaEtaRVec.push_back(DeltaEtaIndividual);
        }
    }
        
    return DeltaEtaRVec;
}

RVec_f MuonAndJPsiStatisticsAndKinematics::DeltaPhiBetweenIsolatedMuonAndJPsiInPeakCalculator()
{
    RVec_FourVector FourVectorTotal = JPsiForMuonsInPeakFourVectorCalculator();
    
    float PhiIsolated = 0;
    float PhiJPsi = 0;
    
    float DeltaPhiIndividual = 0;
    
    RVec_f DeltaPhiRVec {};
    
    for(int i = 0; i < this->Isolated_Muon_phi.size(); i++)
    {
    
        for(int j = 0; j < FourVectorTotal.size(); j++)
        {    
                    PhiIsolated = this->Isolated_Muon_phi[i];

                    PhiJPsi = FourVectorTotal.at(j).Phi();

                    DeltaPhiIndividual = ROOT::VecOps::DeltaPhi(PhiIsolated, PhiJPsi);

                    DeltaPhiRVec.push_back(DeltaPhiIndividual);
        }
    }
        
    return DeltaPhiRVec;   
}

RVec_f MuonAndJPsiStatisticsAndKinematics::DeltaRBetweenIsolatedMuonAndJPsiInPeakCalculator()
{
    RVec_FourVector FourVectorTotal = JPsiForMuonsInPeakFourVectorCalculator();
    
    float EtaIsolated = 0;
    float PhiIsolated = 0;
    float EtaJPsi = 0;
    float PhiJPsi = 0;
    
    float DeltaRIndividual = 0;
    
    RVec_f DeltaRRVec {};
    
    for(int i = 0; i < this->Isolated_Muon_phi.size(); i++)
    {
    
        for(int j = 0; j < FourVectorTotal.size(); j++)
        {    
                    EtaIsolated = this->Isolated_Muon_eta[i];
                    PhiIsolated = this->Isolated_Muon_phi[i];
                    
                    EtaJPsi = FourVectorTotal.at(j).Eta();
                    PhiJPsi = FourVectorTotal.at(j).Phi();

                    DeltaRIndividual = ROOT::VecOps::DeltaR(EtaIsolated, EtaJPsi, PhiIsolated, PhiJPsi);
                    
                    DeltaRRVec.push_back(DeltaRIndividual);
        }
    }
        
    return DeltaRRVec;
}

"""
ROOT.gInterpreter.Declare(cpp_code)


# In[7]:


def IsolatedLeptonSelection():

    for sample in dictOfListOfFiles:       

        rdfPassedIsolatedLeptonTrigger[sample] = rdf[sample].Filter("HLT_IsoMu24 == true ^ HLT_Ele32_WPTight_Gsf == true", "HLTLeptonTrigger")            .Define("LumiXS",wgtFormula[sample])            .Define("leading_isolated_muon_mask", LeadingIsolatedMuonMask)            .Define("leading_isolated_electron_mask", LeadingIsolatedElectronMask)

        if sample == 'Data':
            rdfIsolatedMuonNoHighWeights[sample] = rdfPassedIsolatedLeptonTrigger[sample].Filter("Sum(leading_isolated_muon_mask) == 1 && Sum(leading_isolated_electron_mask) == 0", "Exactly one isolated muon and exactly zero isolated electrons")                .Define("LeadingIsolatedMuon_pt", "Muon_pt[leading_isolated_muon_mask].at(0, -10)")                .Define("LeadingIsolatedMuon_eta", "Muon_eta[leading_isolated_muon_mask].at(0, -10)")                .Define("LeadingIsolatedMuon_phi", "Muon_phi[leading_isolated_muon_mask].at(0, -5)")                .Define("LeadingIsolatedMuon_mass", "Muon_mass[leading_isolated_muon_mask].at(0, -2.71)")                .Define("LeadingIsolatedMuon_charge", "Muon_charge[leading_isolated_muon_mask].at(0, -5)")                .Define("IsolatedMuon_pt", "Muon_pt[leading_isolated_muon_mask]")                .Define("IsolatedMuon_eta", "Muon_eta[leading_isolated_muon_mask]")                .Define("IsolatedMuon_phi", "Muon_phi[leading_isolated_muon_mask]")                .Define("IsolatedMuon_mass", "Muon_mass[leading_isolated_muon_mask]")                .Define("IsolatedMuon_charge", "Muon_charge[leading_isolated_muon_mask]")                .Define("IsolatedMuon_pdgId", "Muon_pdgId[leading_isolated_muon_mask]")                .Define("IsolatedElectron_pt", "Electron_pt[leading_isolated_electron_mask]")                .Define("IsolatedElectron_eta", "Electron_eta[leading_isolated_electron_mask]")                .Define("IsolatedElectron_phi", "Electron_phi[leading_isolated_electron_mask]")                .Define("IsolatedElectron_mass", "Electron_mass[leading_isolated_electron_mask]")                .Define("IsolatedElectron_charge", "Electron_charge[leading_isolated_electron_mask]")                .Define("IsolatedElectron_pdgId", "Electron_pdgId[leading_isolated_electron_mask]")                .Define("IsolatedLepton_pt", "Concatenate(IsolatedMuon_pt, IsolatedElectron_pt)")                .Define("IsolatedLepton_eta", "Concatenate(IsolatedMuon_eta, IsolatedElectron_eta)")                .Define("IsolatedLepton_phi", "Concatenate(IsolatedMuon_phi, IsolatedElectron_phi)")                .Define("IsolatedLepton_mass", "Concatenate(IsolatedMuon_mass, IsolatedElectron_mass)")                .Define("IsolatedLepton_charge", "Concatenate(IsolatedMuon_charge, IsolatedElectron_charge)")                .Define("IsolatedLepton_pdgid", "Concatenate(IsolatedMuon_pdgId, IsolatedElectron_pdgId)")                .Define("METBeforeMETCut", "MET_pt")
            
            rdfIsolatedElectronNoHighWeights[sample] = rdfPassedIsolatedLeptonTrigger[sample].Filter("Sum(leading_isolated_electron_mask) == 1 && Sum(leading_isolated_muon_mask) == 0", "Exactly one isolated electron and exactly zero isolated muon")                .Define("LeadingIsolatedElectron_pt", "Electron_pt[leading_isolated_electron_mask].at(0, -10)")                .Define("LeadingIsolatedElectron_eta", "Electron_eta[leading_isolated_electron_mask].at(0, -10)")                .Define("LeadingIsolatedElectron_phi", "Electron_phi[leading_isolated_electron_mask].at(0, -5)")                .Define("LeadingIsolatedElectron_mass", "Electron_mass[leading_isolated_electron_mask].at(0, -2.71)")                .Define("LeadingIsolatedElectron_charge", "Electron_charge[leading_isolated_electron_mask].at(0, -5)")                .Define("IsolatedMuon_pt", "Muon_pt[leading_isolated_muon_mask]")                .Define("IsolatedMuon_eta", "Muon_eta[leading_isolated_muon_mask]")                .Define("IsolatedMuon_phi", "Muon_phi[leading_isolated_muon_mask]")                .Define("IsolatedMuon_mass", "Muon_mass[leading_isolated_muon_mask]")                .Define("IsolatedMuon_charge", "Muon_charge[leading_isolated_muon_mask]")                .Define("IsolatedMuon_pdgId", "Muon_pdgId[leading_isolated_muon_mask]")                .Define("IsolatedElectron_pt", "Electron_pt[leading_isolated_electron_mask]")                .Define("IsolatedElectron_eta", "Electron_eta[leading_isolated_electron_mask]")                .Define("IsolatedElectron_phi", "Electron_phi[leading_isolated_electron_mask]")                .Define("IsolatedElectron_mass", "Electron_mass[leading_isolated_electron_mask]")                .Define("IsolatedElectron_charge", "Electron_charge[leading_isolated_electron_mask]")                .Define("IsolatedElectron_pdgId", "Electron_pdgId[leading_isolated_electron_mask]")                .Define("IsolatedLepton_pt", "Concatenate(IsolatedMuon_pt, IsolatedElectron_pt)")                .Define("IsolatedLepton_eta", "Concatenate(IsolatedMuon_eta, IsolatedElectron_eta)")                .Define("IsolatedLepton_phi", "Concatenate(IsolatedMuon_phi, IsolatedElectron_phi)")                .Define("IsolatedLepton_mass", "Concatenate(IsolatedMuon_mass, IsolatedElectron_mass)")                .Define("IsolatedLepton_charge", "Concatenate(IsolatedMuon_charge, IsolatedElectron_charge)")                .Define("IsolatedLepton_pdgid", "Concatenate(IsolatedMuon_pdgId, IsolatedElectron_pdgId)")                .Define("METBeforeMETCut", "MET_pt")


        # Change criteria to nothing below or above 4 standard deviations from the mean
        else:
            rdfIsolatedMuonNoHighWeights[sample] = rdfPassedIsolatedLeptonTrigger[sample].Filter("genWeight < 1000 && Sum(leading_isolated_muon_mask) == 1 && Sum(leading_isolated_electron_mask) == 0", "Exactly one isolated muon and exactly zero isolated electrons")                .Define("LeadingIsolatedMuon_pt", "Muon_pt[leading_isolated_muon_mask].at(0, -10)")                .Define("LeadingIsolatedMuon_eta", "Muon_eta[leading_isolated_muon_mask].at(0, -10)")                .Define("LeadingIsolatedMuon_phi", "Muon_phi[leading_isolated_muon_mask].at(0, -5)")                .Define("LeadingIsolatedMuon_mass", "Muon_mass[leading_isolated_muon_mask].at(0, -2.71)")                .Define("LeadingIsolatedMuon_charge", "Muon_charge[leading_isolated_muon_mask].at(0, -5)")                .Define("IsolatedMuon_pt", "Muon_pt[leading_isolated_muon_mask]")                .Define("IsolatedMuon_eta", "Muon_eta[leading_isolated_muon_mask]")                .Define("IsolatedMuon_phi", "Muon_phi[leading_isolated_muon_mask]")                .Define("IsolatedMuon_mass", "Muon_mass[leading_isolated_muon_mask]")                .Define("IsolatedMuon_charge", "Muon_charge[leading_isolated_muon_mask]")                .Define("IsolatedMuon_pdgId", "Muon_pdgId[leading_isolated_muon_mask]")                .Define("IsolatedElectron_pt", "Electron_pt[leading_isolated_electron_mask]")                .Define("IsolatedElectron_eta", "Electron_eta[leading_isolated_electron_mask]")                .Define("IsolatedElectron_phi", "Electron_phi[leading_isolated_electron_mask]")                .Define("IsolatedElectron_mass", "Electron_mass[leading_isolated_electron_mask]")                .Define("IsolatedElectron_charge", "Electron_charge[leading_isolated_electron_mask]")                .Define("IsolatedElectron_pdgId", "Electron_pdgId[leading_isolated_electron_mask]")                .Define("IsolatedLepton_pt", "Concatenate(IsolatedMuon_pt, IsolatedElectron_pt)")                .Define("IsolatedLepton_eta", "Concatenate(IsolatedMuon_eta, IsolatedElectron_eta)")                .Define("IsolatedLepton_phi", "Concatenate(IsolatedMuon_phi, IsolatedElectron_phi)")                .Define("IsolatedLepton_mass", "Concatenate(IsolatedMuon_mass, IsolatedElectron_mass)")                .Define("IsolatedLepton_charge", "Concatenate(IsolatedMuon_charge, IsolatedElectron_charge)")                .Define("IsolatedLepton_pdgid", "Concatenate(IsolatedMuon_pdgId, IsolatedElectron_pdgId)")                .Define("METBeforeMETCut", "MET_pt")

            rdfIsolatedElectronNoHighWeights[sample] = rdfPassedIsolatedLeptonTrigger[sample].Filter("genWeight < 1000 && Sum(leading_isolated_electron_mask) == 1 && Sum(leading_isolated_muon_mask) == 0", "Exactly one isolated electron and exactly zero isolated muon")                .Define("LeadingIsolatedElectron_pt", "Electron_pt[leading_isolated_electron_mask].at(0, -10)")                .Define("LeadingIsolatedElectron_eta", "Electron_eta[leading_isolated_electron_mask].at(0, -10)")                .Define("LeadingIsolatedElectron_phi", "Electron_phi[leading_isolated_electron_mask].at(0, -5)")                .Define("LeadingIsolatedElectron_mass", "Electron_mass[leading_isolated_electron_mask].at(0, -2.71)")                .Define("LeadingIsolatedElectron_charge", "Electron_charge[leading_isolated_electron_mask].at(0, -5)")                .Define("IsolatedMuon_pt", "Muon_pt[leading_isolated_muon_mask]")                .Define("IsolatedMuon_eta", "Muon_eta[leading_isolated_muon_mask]")                .Define("IsolatedMuon_phi", "Muon_phi[leading_isolated_muon_mask]")                .Define("IsolatedMuon_mass", "Muon_mass[leading_isolated_muon_mask]")                .Define("IsolatedMuon_charge", "Muon_charge[leading_isolated_muon_mask]")                .Define("IsolatedMuon_pdgId", "Muon_pdgId[leading_isolated_muon_mask]")                .Define("IsolatedElectron_pt", "Electron_pt[leading_isolated_electron_mask]")                .Define("IsolatedElectron_eta", "Electron_eta[leading_isolated_electron_mask]")                .Define("IsolatedElectron_phi", "Electron_phi[leading_isolated_electron_mask]")                .Define("IsolatedElectron_mass", "Electron_mass[leading_isolated_electron_mask]")                .Define("IsolatedElectron_charge", "Electron_charge[leading_isolated_electron_mask]")                .Define("IsolatedElectron_pdgId", "Electron_pdgId[leading_isolated_electron_mask]")                .Define("IsolatedLepton_pt", "Concatenate(IsolatedMuon_pt, IsolatedElectron_pt)")                .Define("IsolatedLepton_eta", "Concatenate(IsolatedMuon_eta, IsolatedElectron_eta)")                .Define("IsolatedLepton_phi", "Concatenate(IsolatedMuon_phi, IsolatedElectron_phi)")                .Define("IsolatedLepton_mass", "Concatenate(IsolatedMuon_mass, IsolatedElectron_mass)")                .Define("IsolatedLepton_charge", "Concatenate(IsolatedMuon_charge, IsolatedElectron_charge)")                .Define("IsolatedLepton_pdgid", "Concatenate(IsolatedMuon_pdgId, IsolatedElectron_pdgId)")                .Define("METBeforeMETCut", "MET_pt")

        rdfIsolatedMuonAfterMETCut[sample] = rdfIsolatedMuonNoHighWeights[sample].Filter("MET_pt > 30", "Muon MET Greater than 30 GeV")            .Define("METAfterMETCut", "MET_pt")            .Define("jet_mask", "ROOT::VecOps::RVec<Int_t> jmask = (Jet_pt >= 30 && abs(Jet_eta) <= 2.5 && Jet_jetId >= 2); "                        "for(int i=0; i < IsolatedLepton_pt.size(); ++i){"                            "ROOT::VecOps::RVec<Float_t> dr;"                            "for(int j=0; j < jmask.size(); ++j){"                                "dr.push_back(ROOT::VecOps::DeltaR(Jet_eta.at(j), IsolatedLepton_eta.at(i), Jet_phi.at(j), IsolatedLepton_phi.at(i)));}"                                "jmask = jmask && dr >= 0.4;"                                "dr.clear();}"                        "return jmask;")            .Define("jpsi_mu_candidate_mask", JPsiCandidateMask)            .Define("Num_Jets", "Jet_pt[jet_mask].size()")            .Define("Num_JPsi_Muons", "Muon_pt[jpsi_mu_candidate_mask].size()")
#             .Define("jpsi_muon_positive_candidate_mask", JPsiMuonPositiveCandidateMask)\
#             .Define("jpsi_muon_negative_candidate_mask", JPsiMuonNegativeCandidateMask)\
#             .Define("leading_isolated_muon_mask_pfisoid_medium_mask", LeadingIsolatedMuonMaskPfIsoMedium)\
#             .Define("leading_isolated_muon_mask_muon_id_medium_iso_tight", LeadingIsolatedMuonMaskMuonIdMediumIsoTight)\
#             .Define("leading_isolated_muon_mask_muon_id_loose_iso_tight", LeadingIsolatedMuonMaskMuonIdLooseIsoIdTight)\
#             .Define("leading_isolated_muon_mask_muon_id_medium_iso_medium", LeadingIsolatedMuonMaskMuonIdMediumIsoIdMedium)\
#             .Define("leading_isolated_muon_mask_muon_id_loose_iso_medium", LeadingIsolatedMuonMaskMuonIdLooseIsoIdMedium)\
#             .Define("jpsi_muon_positive_candidate_peak_only_mask", JPsiMuonPositiveCandidatePeakOnlyMask)\
#             .Define("jpsi_muon_negative_candidate_peak_only_mask", JPsiMuonNegativeCandidatePeakOnlyMask)\
#             .Define("jpsi_muon_positive_candidate_pfisoid_medium_mask", JPsiMuonPositiveCandidatePfIsoMediumMask)\
#             .Define("jpsi_muon_negative_candidate_pfisoid_medium_mask", JPsiMuonNegativeCandidatePfIsoMediumMask)\

            
        
        rdfIsolatedElectronAfterMETCut[sample] = rdfIsolatedElectronNoHighWeights[sample].Filter("MET_pt > 30", "Electron MET Greater than 30 GeV")            .Define("METAfterMETCut", "MET_pt")
        
    return rdfIsolatedMuonAfterMETCut, rdfIsolatedElectronAfterMETCut


# In[8]:


rdfIsolatedMuonAfterMETCut, rdfIsolatedElectronAfterMETCut = IsolatedLeptonSelection()


# In[9]:


def JetSelection():

    for sample in dictOfListOfFiles:
        
        rdfJetAndIsolatedLeptonFiltered[sample] = rdfIsolatedMuonAfterMETCut[sample].Filter("Num_Jets >= 2", "At Least Two Jets")            .Define("SJet1_pt", "Jet_pt[jet_mask].size() > 0 ? Jet_pt[jet_mask].at(0) : -500")            .Define("SJet2_pt", "Jet_pt[jet_mask].size() > 1 ? Jet_pt[jet_mask].at(1) : -500")            .Define("SJet1_eta", "Jet_eta[jet_mask].size() > 0 ? Jet_eta[jet_mask].at(0) : 500")            .Define("SJet2_eta", "Jet_eta[jet_mask].size() > 1 ? Jet_eta[jet_mask].at(1) : 500")            .Define("SJet1_phi", "Jet_phi[jet_mask].size() > 0 ? Jet_phi[jet_mask].at(0) : 500")            .Define("SJet2_phi", "Jet_phi[jet_mask].size() > 1 ? Jet_phi[jet_mask].at(1) : 500")            .Define("MTofMETandMu", "FTA::transverseMassMET(IsolatedMuon_pt, IsolatedMuon_phi, IsolatedMuon_mass, MET_pt, MET_phi)")            .Define("Ht", "Sum(Jet_pt[jet_mask])")
        
    return rdfJetAndIsolatedLeptonFiltered


# In[10]:


rdfJetAndIsolatedLeptonFiltered = JetSelection()


# In[11]:


def JPsiSelection():

    for sample in dictOfListOfFiles:

        rdfJPsiMuons[sample] = rdfJetAndIsolatedLeptonFiltered[sample].Filter("Sum(jpsi_mu_candidate_mask) >= 2", "At Least Two JPsi Candidates")            .Define("JPsiCandidate_pt", "Muon_pt[jpsi_mu_candidate_mask]")            .Define("JPsiCandidate_eta", "Muon_eta[jpsi_mu_candidate_mask]")            .Define("JPsiCandidate_phi", "Muon_phi[jpsi_mu_candidate_mask]")            .Define("JPsiCandidate_mass", "Muon_mass[jpsi_mu_candidate_mask]")            .Define("JPsiCandidate_charge", "Muon_charge[jpsi_mu_candidate_mask]")            .Define("InvariantMassJPsiMuons", "std::cout << rdfentry_ << std::endl; auto c = MuonAndJPsiStatisticsAndKinematics(JPsiCandidate_pt, JPsiCandidate_eta, JPsiCandidate_phi, JPsiCandidate_mass, JPsiCandidate_charge); return c.JPsiMuonInvariantMassCalculator();")            .Define("InvariantMassJPsiAndIsolatedMuons", "auto c = MuonAndJPsiStatisticsAndKinematics(JPsiCandidate_pt, JPsiCandidate_eta, JPsiCandidate_phi, JPsiCandidate_mass, JPsiCandidate_charge, IsolatedMuon_pt, IsolatedMuon_eta, IsolatedMuon_phi, IsolatedMuon_mass, IsolatedMuon_charge); return c.IsolatedAndJPsiMuonInvariantMassCalculator();")            .Define("DeltaEtaBetweenIsolatedAndJPsiMuon", "auto c = MuonAndJPsiStatisticsAndKinematics(JPsiCandidate_pt, JPsiCandidate_eta, JPsiCandidate_phi, JPsiCandidate_mass, JPsiCandidate_charge, IsolatedMuon_pt, IsolatedMuon_eta, IsolatedMuon_phi, IsolatedMuon_mass, IsolatedMuon_charge); return c.DeltaEtaBetweenIsolatedAndJPsiMuonCalculator();")            .Define("DeltaPhiBetweenIsolatedAndJPsiMuon", "auto c = MuonAndJPsiStatisticsAndKinematics(JPsiCandidate_pt, JPsiCandidate_eta, JPsiCandidate_phi, JPsiCandidate_mass, JPsiCandidate_charge, IsolatedMuon_pt, IsolatedMuon_eta, IsolatedMuon_phi, IsolatedMuon_mass, IsolatedMuon_charge); return c.DeltaPhiBetweenIsolatedAndJPsiMuonCalculator();")            .Define("DeltaRBetweenIsolatedAndJPsiMuon", "auto c = MuonAndJPsiStatisticsAndKinematics(JPsiCandidate_pt, JPsiCandidate_eta, JPsiCandidate_phi, JPsiCandidate_mass, JPsiCandidate_charge, IsolatedMuon_pt, IsolatedMuon_eta, IsolatedMuon_phi, IsolatedMuon_mass, IsolatedMuon_charge); return c.DeltaRBetweenIsolatedAndJPsiMuonCalculator();")            .Define("JPsi_pt", "auto c = MuonAndJPsiStatisticsAndKinematics(JPsiCandidate_pt, JPsiCandidate_eta, JPsiCandidate_phi, JPsiCandidate_mass, JPsiCandidate_charge); return c.ReturnJPsiPt();")            .Define("JPsi_eta", "auto c = MuonAndJPsiStatisticsAndKinematics(JPsiCandidate_pt, JPsiCandidate_eta, JPsiCandidate_phi, JPsiCandidate_mass, JPsiCandidate_charge); return c.ReturnJPsiEta();")            .Define("JPsi_phi", "auto c = MuonAndJPsiStatisticsAndKinematics(JPsiCandidate_pt, JPsiCandidate_eta, JPsiCandidate_phi, JPsiCandidate_mass, JPsiCandidate_charge); return c.ReturnJPsiPhi();")            .Define("DeltaEtaBetweenJPsiMuons", "auto c = MuonAndJPsiStatisticsAndKinematics(JPsiCandidate_pt, JPsiCandidate_eta, JPsiCandidate_phi, JPsiCandidate_mass, JPsiCandidate_charge); return c.DeltaEtaBetweenJPsiMuonsCalculator();")            .Define("DeltaPhiBetweenJPsiMuons", "auto c = MuonAndJPsiStatisticsAndKinematics(JPsiCandidate_pt, JPsiCandidate_eta, JPsiCandidate_phi, JPsiCandidate_mass, JPsiCandidate_charge); return c.DeltaPhiBetweenJPsiMuonsCalculator();")            .Define("DeltaRBetweenJPsiMuons", "auto c = MuonAndJPsiStatisticsAndKinematics(JPsiCandidate_pt, JPsiCandidate_eta, JPsiCandidate_phi, JPsiCandidate_mass, JPsiCandidate_charge); return c.DeltaRBetweenJPsiMuonsCalculator();")            .Define("DeltaRBetweenJPsiMuonsInPeak", "auto c = MuonAndJPsiStatisticsAndKinematics(JPsiCandidate_pt, JPsiCandidate_eta, JPsiCandidate_phi, JPsiCandidate_mass, JPsiCandidate_charge); return c.ReturnDeltaRBetweenJPsiMuonsInPeak();")            .Define("DeltaEtaBetweenIsolatedMuonAndJPsi", "auto c = MuonAndJPsiStatisticsAndKinematics(JPsiCandidate_pt, JPsiCandidate_eta, JPsiCandidate_phi, JPsiCandidate_mass, JPsiCandidate_charge, IsolatedMuon_pt, IsolatedMuon_eta, IsolatedMuon_phi, IsolatedMuon_mass, IsolatedMuon_charge); return c.DeltaEtaBetweenIsolatedMuonAndJPsiCalculator();")            .Define("DeltaPhiBetweenIsolatedMuonAndJPsi", "auto c = MuonAndJPsiStatisticsAndKinematics(JPsiCandidate_pt, JPsiCandidate_eta, JPsiCandidate_phi, JPsiCandidate_mass, JPsiCandidate_charge, IsolatedMuon_pt, IsolatedMuon_eta, IsolatedMuon_phi, IsolatedMuon_mass, IsolatedMuon_charge); return c.DeltaPhiBetweenIsolatedMuonAndJPsiCalculator();")            .Define("DeltaRBetweenIsolatedMuonAndJPsi", "auto c = MuonAndJPsiStatisticsAndKinematics(JPsiCandidate_pt, JPsiCandidate_eta, JPsiCandidate_phi, JPsiCandidate_mass, JPsiCandidate_charge, IsolatedMuon_pt, IsolatedMuon_eta, IsolatedMuon_phi, IsolatedMuon_mass, IsolatedMuon_charge); return c.DeltaRBetweenIsolatedMuonAndJPsiCalculator();")            .Define("JPsiMuon_pt", "Muon_pt[jpsi_mu_candidate_mask]")            .Define("JPsiMuon_pfRelIso03_all", "Muon_pfRelIso03_all[jpsi_mu_candidate_mask]")            .Define("JPsiMuon_pfIsoid", "Muon_pfIsoId[jpsi_mu_candidate_mask]")            .Define("IsolatedMuonsInRangeInvariantMassPlot_pfRelIso03_all", "auto c = MuonAndJPsiStatisticsAndKinematics(JPsiCandidate_pt, JPsiCandidate_eta, JPsiCandidate_phi, JPsiCandidate_mass, JPsiCandidate_charge, IsolatedMuon_pt, IsolatedMuon_eta, IsolatedMuon_phi, IsolatedMuon_mass, IsolatedMuon_charge, Muon_pfRelIso03_all[leading_isolated_muon_mask], Muon_pfRelIso03_all[jpsi_mu_candidate_mask]); return c.ReturnMuonPfRelIso03AllForIsolatedMuonInInvariantMass();")            .Define("JPsiMuonsInRangeInvariantMassPlot_pfRelIso03_all", "auto c = MuonAndJPsiStatisticsAndKinematics(JPsiCandidate_pt, JPsiCandidate_eta, JPsiCandidate_phi, JPsiCandidate_mass, JPsiCandidate_charge, IsolatedMuon_pt, IsolatedMuon_eta, IsolatedMuon_phi, IsolatedMuon_mass, IsolatedMuon_charge, Muon_pfRelIso03_all[leading_isolated_muon_mask], Muon_pfRelIso03_all[jpsi_mu_candidate_mask]); return c.ReturnMuonPfRelIso03AllForJPsiMuonsInInvariantMass();")            .Define("JPsiMassesInRangeInvariantMassPlot", "auto c = MuonAndJPsiStatisticsAndKinematics(JPsiCandidate_pt, JPsiCandidate_eta, JPsiCandidate_phi, JPsiCandidate_mass, JPsiCandidate_charge, IsolatedMuon_pt, IsolatedMuon_eta, IsolatedMuon_phi, IsolatedMuon_mass, IsolatedMuon_charge); return c.ReturnJPsiMassforJPsiInInvariantMass()")            .Define("IsolatedMuonPtAtLowDeltaR", "auto c = MuonAndJPsiStatisticsAndKinematics(JPsiCandidate_pt, JPsiCandidate_eta, JPsiCandidate_phi, JPsiCandidate_mass, JPsiCandidate_charge, IsolatedMuon_pt, IsolatedMuon_eta, IsolatedMuon_phi, IsolatedMuon_mass, IsolatedMuon_charge); return c.ReturnLowDeltaRIsolatedMuonAndJPsiMuonPt(false)")            .Define("JPsiMuonPtAtLowDeltaR", "auto c = MuonAndJPsiStatisticsAndKinematics(JPsiCandidate_pt, JPsiCandidate_eta, JPsiCandidate_phi, JPsiCandidate_mass, JPsiCandidate_charge, IsolatedMuon_pt, IsolatedMuon_eta, IsolatedMuon_phi, IsolatedMuon_mass, IsolatedMuon_charge); return c.ReturnLowDeltaRIsolatedMuonAndJPsiMuonPt(true)")            .Define("JPsiMassInPeak", "auto c = MuonAndJPsiStatisticsAndKinematics(JPsiCandidate_pt, JPsiCandidate_eta, JPsiCandidate_phi, JPsiCandidate_mass, JPsiCandidate_charge, IsolatedMuon_pt, IsolatedMuon_eta, IsolatedMuon_phi, IsolatedMuon_mass, IsolatedMuon_charge); return c.ReturnJPsiInPeakMass();")            .Define("DeltaEtaForJPsiMuonsInPeak", "auto c = MuonAndJPsiStatisticsAndKinematics(JPsiCandidate_pt, JPsiCandidate_eta, JPsiCandidate_phi, JPsiCandidate_mass, JPsiCandidate_charge, IsolatedMuon_pt, IsolatedMuon_eta, IsolatedMuon_phi, IsolatedMuon_mass, IsolatedMuon_charge); return c.ReturnDeltaEtaBetweenJPsiMuonsInPeak();")            .Define("DeltaPhiForJPsiMuonsInPeak", "auto c = MuonAndJPsiStatisticsAndKinematics(JPsiCandidate_pt, JPsiCandidate_eta, JPsiCandidate_phi, JPsiCandidate_mass, JPsiCandidate_charge, IsolatedMuon_pt, IsolatedMuon_eta, IsolatedMuon_phi, IsolatedMuon_mass, IsolatedMuon_charge); return c.ReturnDeltaPhiBetweenJPsiMuonsInPeak();")            .Define("DeltaRForJPsiMuonsInPeak", "auto c = MuonAndJPsiStatisticsAndKinematics(JPsiCandidate_pt, JPsiCandidate_eta, JPsiCandidate_phi, JPsiCandidate_mass, JPsiCandidate_charge, IsolatedMuon_pt, IsolatedMuon_eta, IsolatedMuon_phi, IsolatedMuon_mass, IsolatedMuon_charge); return c.ReturnDeltaRBetweenJPsiMuonsInPeak();")            .Define("JPsiMuonsChargeInPeak", "auto c = MuonAndJPsiStatisticsAndKinematics(JPsiCandidate_pt, JPsiCandidate_eta, JPsiCandidate_phi, JPsiCandidate_mass, JPsiCandidate_charge, IsolatedMuon_pt, IsolatedMuon_eta, IsolatedMuon_phi, IsolatedMuon_mass, IsolatedMuon_charge); return c.ReturnJPsiMuonsInPeakCharge();")            .Define("JPsiPtWithMuonInPeak", "auto c = MuonAndJPsiStatisticsAndKinematics(JPsiCandidate_pt, JPsiCandidate_eta, JPsiCandidate_phi, JPsiCandidate_mass, JPsiCandidate_charge, IsolatedMuon_pt, IsolatedMuon_eta, IsolatedMuon_phi, IsolatedMuon_mass, IsolatedMuon_charge); return c.ReturnJPsiInPeakPt();")            .Define("JPsiMuonPositivePt", "auto c = MuonAndJPsiStatisticsAndKinematics(JPsiCandidate_pt, JPsiCandidate_eta, JPsiCandidate_phi, JPsiCandidate_mass, JPsiCandidate_charge, IsolatedMuon_pt, IsolatedMuon_eta, IsolatedMuon_phi, IsolatedMuon_mass, IsolatedMuon_charge); return c.ReturnJPsiMuonPtInPeak(true);")            .Define("JPsiMuonNegativePt", "auto c = MuonAndJPsiStatisticsAndKinematics(JPsiCandidate_pt, JPsiCandidate_eta, JPsiCandidate_phi, JPsiCandidate_mass, JPsiCandidate_charge, IsolatedMuon_pt, IsolatedMuon_eta, IsolatedMuon_phi, IsolatedMuon_mass, IsolatedMuon_charge); return c.ReturnJPsiMuonPtInPeak(false);")            .Define("JPsiMuonIsolationInPeak", "auto c = MuonAndJPsiStatisticsAndKinematics(JPsiCandidate_pt, JPsiCandidate_eta, JPsiCandidate_phi, JPsiCandidate_mass, JPsiCandidate_charge, IsolatedMuon_pt, IsolatedMuon_eta, IsolatedMuon_phi, IsolatedMuon_mass, IsolatedMuon_charge); return c.ReturnJPsiMuonInPeakPfRelIso03All();")            .Define("TopMassInPeak", "auto c = MuonAndJPsiStatisticsAndKinematics(JPsiCandidate_pt, JPsiCandidate_eta, JPsiCandidate_phi, JPsiCandidate_mass, JPsiCandidate_charge, IsolatedMuon_pt, IsolatedMuon_eta, IsolatedMuon_phi, IsolatedMuon_mass, IsolatedMuon_charge); return c.ReturnTopQuarkInPeakMass();")            .Define("DeltaEtaBetweenIsolatedMuonAndJPsiInPeak", "auto c = MuonAndJPsiStatisticsAndKinematics(JPsiCandidate_pt, JPsiCandidate_eta, JPsiCandidate_phi, JPsiCandidate_mass, JPsiCandidate_charge, IsolatedMuon_pt, IsolatedMuon_eta, IsolatedMuon_phi, IsolatedMuon_mass, IsolatedMuon_charge); return c.DeltaEtaBetweenIsolatedMuonAndJPsiInPeakCalculator();")            .Define("DeltaPhiBetweenIsolatedMuonAndJPsiInPeak", "auto c = MuonAndJPsiStatisticsAndKinematics(JPsiCandidate_pt, JPsiCandidate_eta, JPsiCandidate_phi, JPsiCandidate_mass, JPsiCandidate_charge, IsolatedMuon_pt, IsolatedMuon_eta, IsolatedMuon_phi, IsolatedMuon_mass, IsolatedMuon_charge); return c.DeltaPhiBetweenIsolatedMuonAndJPsiInPeakCalculator();")            .Define("DeltaRBetweenIsolatedMuonAndJPsiInPeak", "auto c = MuonAndJPsiStatisticsAndKinematics(JPsiCandidate_pt, JPsiCandidate_eta, JPsiCandidate_phi, JPsiCandidate_mass, JPsiCandidate_charge, IsolatedMuon_pt, IsolatedMuon_eta, IsolatedMuon_phi, IsolatedMuon_mass, IsolatedMuon_charge); return c.DeltaEtaBetweenIsolatedMuonAndJPsiInPeakCalculator();")
#             .Define("IsolatedMuon_pfRelIso03_all", "Muon_pfRelIso03_all[leading_isolated_muon_mask]")\
#             .Define("IsolatedMuon_pfRelIso03_chg", "Muon_pfRelIso03_chg[leading_isolated_muon_mask]")\
#             .Define("JPsiMuon_pfRelIso03_all", "Muon_pfRelIso03_all[jpsi_mu_candidate_mask]")\
#             .Define("JPsiMuonPositive_pfRelIso03_all", "Muon_pfRelIso03_all[jpsi_muon_positive_candidate_mask]")\
#             .Define("JPsiMuonPositive_pfRelIso03_chg", "Muon_pfRelIso03_chg[jpsi_muon_positive_candidate_mask]")\
#             .Define("JPsiMuonNegative_pfRelIso03_all", "Muon_pfRelIso03_all[jpsi_muon_negative_candidate_mask]")\
#             .Define("JPsiMuonNegative_pfRelIso03_chg", "Muon_pfRelIso03_chg[jpsi_muon_negative_candidate_mask]")\
#             .Define("IsolatedMuon_pfRelIso03_all_mediumID", "Muon_pfRelIso03_all[leading_isolated_muon_mask_pfisoid_medium_mask]")\
#             .Define("IsolatedMuon_pfRelIso03_chg_mediumID", "Muon_pfRelIso03_chg[leading_isolated_muon_mask_pfisoid_medium_mask]")\
#             .Define("JPsiMuonPositive_pfRelIso03_all_mediumID", "Muon_pfRelIso03_all[jpsi_muon_positive_candidate_pfisoid_medium_mask]")\
#             .Define("JPsiMuonPositive_pfRelIso03_chg_mediumID", "Muon_pfRelIso03_chg[jpsi_muon_positive_candidate_pfisoid_medium_mask]")\
#             .Define("JPsiMuonNegative_pfRelIso03_all_mediumID", "Muon_pfRelIso03_all[jpsi_muon_negative_candidate_pfisoid_medium_mask]")\
#             .Define("JPsiMuonNegative_pfRelIso03_chg_mediumID", "Muon_pfRelIso03_chg[jpsi_muon_negative_candidate_pfisoid_medium_mask]")\
#             .Define("JPsiMuonPositive_peak_only_pfRelIso03_all", "Muon_pfRelIso03_all[jpsi_muon_positive_candidate_peak_only_mask]")\
#             .Define("JPsiMuonPositive_peak_only_pfRelIso03_chg", "Muon_pfRelIso03_chg[jpsi_muon_positive_candidate_peak_only_mask]")\
#             .Define("JPsiMuonNegative_peak_only_pfRelIso03_all", "Muon_pfRelIso03_all[jpsi_muon_negative_candidate_peak_only_mask]")\
#             .Define("JPsiMuonNegative_peak_only_pfRelIso03_chg", "Muon_pfRelIso03_chg[jpsi_muon_negative_candidate_peak_only_mask]")

        
    return rdfJPsiMuons


# In[12]:


rdfJPsiMuon = JPsiSelection()


# In[13]:


# def FurtherCuts():
    
#     for sample in dictOfListOfFiles:
        
#         rdfRemainingIsolatedMuonMuonIdMediumIsoIdTight[sample] = rdfJPsiMuons[sample].Filter("Sum(leading_isolated_muon_mask_muon_id_medium_iso_tight) >= 1", "At Least One Further Muon - Id Medium Iso Id Tight")\
#             .Define("FurtherMuonMediumTight_pt", "Muon_pt[leading_isolated_muon_mask_muon_id_medium_iso_tight]")\
#             .Define("FurtherMuonMediumTight_pfRelIso03_all", "Muon_pfRelIso03_all[leading_isolated_muon_mask_muon_id_medium_iso_tight]")\
#             .Define("FurtherMuonMediumTight_pfIsoid", "Muon_pfIsoId[leading_isolated_muon_mask_muon_id_medium_iso_tight]")

#         rdfRemainingIsolatedMuonMuonIdLooseIsoIdTight[sample] = rdfJPsiMuons[sample].Filter("Sum(leading_isolated_muon_mask_muon_id_loose_iso_tight) >= 1", "At Least One Further Muon - Id Loose Iso Id Tight")\
#             .Define("FurtherMuonLooseTight_pt", "Muon_pt[leading_isolated_muon_mask_muon_id_loose_iso_tight]")\
#             .Define("FurtherMuonLooseTight_pfRelIso03_all", "Muon_pfRelIso03_all[leading_isolated_muon_mask_muon_id_loose_iso_tight]")\
#             .Define("FurtherMuonLooseTight_pfIsoid", "Muon_pfIsoId[leading_isolated_muon_mask_muon_id_loose_iso_tight]")

#         rdfRemainingIsolatedMuonMuonIdMediumIsoIdMedium[sample] = rdfJPsiMuons[sample].Filter("Sum(leading_isolated_muon_mask_muon_id_medium_iso_medium) >= 1", "At Least One Further Muon - Id Medium Iso Id Medium")\
#             .Define("FurtherMuonMediumMedium_pt", "Muon_pt[leading_isolated_muon_mask_muon_id_medium_iso_medium]")\
#             .Define("FurtherMuonMediumMedium_pfRelIso03_all", "Muon_pfRelIso03_all[leading_isolated_muon_mask_muon_id_medium_iso_medium]")\
#             .Define("FurtherMuonMediumMedium_pfIsoid", "Muon_pfIsoId[leading_isolated_muon_mask_muon_id_medium_iso_medium]")

#         rdfRemainingIsolatedMuonMuonIdLooseIsoIdMedium[sample] = rdfJPsiMuons[sample].Filter("Sum(leading_isolated_muon_mask_muon_id_loose_iso_medium) >= 1", "At Least One Further Muon - Id Loose Iso Id Medium")\
#             .Define("FurtherMuonLooseMedium_pt", "Muon_pt[leading_isolated_muon_mask_muon_id_loose_iso_medium]")\
#             .Define("FurtherMuonLooseMedium_pfRelIso03_all", "Muon_pfRelIso03_all[leading_isolated_muon_mask_muon_id_loose_iso_medium]")\
#             .Define("FurtherMuonLooseMedium_pfIsoid", "Muon_pfIsoId[leading_isolated_muon_mask_muon_id_loose_iso_medium]")
    
#     return rdfRemainingIsolatedMuonMuonIdMediumIsoIdTight, rdfRemainingIsolatedMuonMuonIdLooseIsoIdTight, rdfRemainingIsolatedMuonMuonIdMediumIsoIdMedium, rdfRemainingIsolatedMuonMuonIdLooseIsoIdMedium


# In[14]:


# rdfRemainingIsolatedMuonMuonIdMediumIsoIdTight, rdfRemainingIsolatedMuonMuonIdLooseIsoIdTight, rdfRemainingIsolatedMuonMuonIdMediumIsoIdMedium, rdfRemainingIsolatedMuonMuonIdLooseIsoIdMedium = FurtherCuts()


# In[15]:


# Make this block a function with the inputs as the nodes we want to attach histograms to (2-17-22)
# rdfWithFourMomentum needs to be changed to whatever our final RDataFrame is (5-23-22)

for sample in dictOfListOfFiles:
    
    if sample not in hist.keys():
        hist[sample] = {}
        report[sample] = rdf[sample].Report()
        #nparraynode[sample] = rdfLeadingMuon[sample]
        
    if sample == 'Data':
        
        hist['Data']["leading_isolated_muon_pt_initial"] = rdfIsolatedMuonNoHighWeights['Data'].Histo1D(("leading_isolated_muon_pt_initial","Leading Isolated Muon Transverse Momentum (One Muon, No Electrons); Pt (GeV);Events",100,20,220),"LeadingIsolatedMuon_pt")
        hist['Data']["leading_isolated_muon_eta_initial"] = rdfIsolatedMuonNoHighWeights['Data'].Histo1D(("leading_isolated_muon_eta_initial", "Leading Isolated Muon Pseudorapidity (One Muon, No Electrons); Eta; Events",100,-3,3),"LeadingIsolatedMuon_eta")
        hist['Data']["leading_isolated_muon_phi_initial"] = rdfIsolatedMuonNoHighWeights['Data'].Histo1D(("leading_isolated_muon_phi_initial", "Leading Isolated Muon Angle (One Muon, No Electrons); Phi (Radians); Events",100,-3.5,3.5),"LeadingIsolatedMuon_phi")
        hist['Data']["leading_isolated_muon_mass_initial"] = rdfIsolatedMuonNoHighWeights['Data'].Histo1D(("leading_isolated_muon_mass_initial", "Leading Isolated Muon Mass (One Muon, No Electrons); Mass(Gev); Events",10,.1,.2),"LeadingIsolatedMuon_mass")
        hist['Data']["leading_isolated_muon_charge_initial"] = rdfIsolatedMuonNoHighWeights['Data'].Histo1D(("leading_isolated_muon_charge_initial", "Leading Isolated Muon Charge (One Muon, No Electrons); Charge; Events",5,-2,2),"LeadingIsolatedMuon_charge")
        hist['Data']["leading_isolated_electron_pt_initial"] = rdfIsolatedElectronNoHighWeights['Data'].Histo1D(("leading_isolated_electron_pt_initial","Leading Isolated Electron Transverse Momentum (One Electron, No Muons); Pt (GeV);Events",100,20,220),"LeadingIsolatedElectron_pt")
        hist['Data']["leading_isolated_electron_eta_initial"] = rdfIsolatedElectronNoHighWeights['Data'].Histo1D(("leading_isolated_electron_eta_initial", "Leading Isolated Electron Pseudorapidity (One Electron, No Muons); Eta; Events",100,-3,3),"LeadingIsolatedElectron_eta")
        hist['Data']["leading_isolated_electron_phi_initial"] = rdfIsolatedElectronNoHighWeights['Data'].Histo1D(("leading_isolated_electron_phi_initial", "Leading Isolated Electron Angle (One Electron, No Muons); Phi (Radians); Events",100,-3.5,3.5),"LeadingIsolatedElectron_phi")
        hist['Data']["leading_isolated_electron_mass_initial"] = rdfIsolatedElectronNoHighWeights['Data'].Histo1D(("leading_isolated_electron_mass_initial", "Leading Isolated Electron Mass (One Electron, No Muons); Mass(Gev); Events",10,.0001,.001),"LeadingIsolatedElectron_mass")
        hist['Data']["leading_isolated_electron_charge_initial"] = rdfIsolatedElectronNoHighWeights['Data'].Histo1D(("leading_isolated_electron_charge_initial", "Leading Isolated Electron Charge (One Electron, No Muons); Charge; Events",5,-2,2),"LeadingIsolatedElectron_charge")
        
        hist['Data']["met_before_met_cut"] = rdfIsolatedMuonNoHighWeights['Data'].Histo1D(("met_before_met_cut", "MET Before MET Cut; Pt (GeV); Events",100,0,250),"METBeforeMETCut")
        hist['Data']["met_after_met_cut"] = rdfIsolatedMuonAfterMETCut['Data'].Histo1D(("met_after_met_cut", "MET After MET Cut; Pt (GeV); Events",100,0,250),"METAfterMETCut")
        hist['Data']["leading_isolated_muon_pt_after_met_cut"] = rdfIsolatedMuonAfterMETCut['Data'].Histo1D(("leading_isolated_muon_pt_after_met_cut","Leading Isolated Muon Transverse Momentum After MET Cut; Pt (GeV);Events",100,20,220),"LeadingIsolatedMuon_pt")
        hist['Data']["leading_isolated_muon_eta_after_met_cut"] = rdfIsolatedMuonAfterMETCut['Data'].Histo1D(("leading_isolated_muon_eta_after_met_cut", "Leading Isolated Muon Pseudorapidity After MET Cut; Eta; Events",100,-3,3),"LeadingIsolatedMuon_eta")
        hist['Data']["leading_isolated_muon_phi_after_met_cut"] = rdfIsolatedMuonAfterMETCut['Data'].Histo1D(("leading_isolated_muon_phi_after_met_cut", "Leading Isolated Muon Angle After MET Cut; Phi (Radians); Events",100,-3.5,3.5),"LeadingIsolatedMuon_phi")
        hist['Data']["leading_isolated_muon_mass_after_met_cut"] = rdfIsolatedMuonAfterMETCut['Data'].Histo1D(("leading_isolated_muon_mass_after_met_cut", "Leading Isolated Muon Mass After MET Cut; Mass(Gev); Events",10,.1,.2),"LeadingIsolatedMuon_mass")
        hist['Data']["leading_isolated_muon_charge_after_met_cut"] = rdfIsolatedMuonAfterMETCut['Data'].Histo1D(("leading_isolated_muon_charge_after_met_cut", "Leading Isolated Muon Charge After MET Cut; Charge; Events",5,-2,2),"LeadingIsolatedMuon_charge")
        hist['Data']["leading_isolated_electron_pt_after_met_cut"] = rdfIsolatedElectronAfterMETCut['Data'].Histo1D(("leading_isolated_electron_pt_after_met_cut","Leading Isolated Electron Transverse Momentum After MET Cut; Pt (GeV);Events",100,20,220),"LeadingIsolatedElectron_pt")
        hist['Data']["leading_isolated_electron_eta_after_met_cut"] = rdfIsolatedElectronAfterMETCut['Data'].Histo1D(("leading_isolated_electron_eta_after_met_cut", "Leading Isolated Electron Pseudorapidity After MET Cut; Eta; Events",100,-3,3),"LeadingIsolatedElectron_eta")
        hist['Data']["leading_isolated_electron_phi_after_met_cut"] = rdfIsolatedElectronAfterMETCut['Data'].Histo1D(("leading_isolated_electron_phi_after_met_cut", "Leading Isolated Electron Angle After MET Cut; Phi (Radians); Events",100,-3.5,3.5),"LeadingIsolatedElectron_phi")
        hist['Data']["leading_isolated_electron_mass_after_met_cut"] = rdfIsolatedElectronAfterMETCut['Data'].Histo1D(("leading_isolated_electron_mass_after_met_cut", "Leading Isolated Electron Mass After MET Cut; Mass(Gev); Events",10,.0001,.001),"LeadingIsolatedElectron_mass")
        hist['Data']["leading_isolated_electron_charge_after_met_cut"] = rdfIsolatedElectronAfterMETCut['Data'].Histo1D(("leading_isolated_electron_charge_after_met_cut", "Leading Isolated Electron Charge After MET Cut; Charge; Events",5,-2,2),"LeadingIsolatedElectron_charge")
        
        hist['Data']["number_of_jets_initial"] = rdfIsolatedMuonAfterMETCut['Data'].Histo1D(("number_of_jets_initial", "Number of Jets Before Jet Cut; Number Of Jets; Events", 20, 0, 20), "Num_Jets")
        hist['Data']["number_of_jets"] = rdfJetAndIsolatedLeptonFiltered['Data'].Histo1D(("number_of_jets", "Number Of Jets; Number Of Jets; Events", 20, 0, 20), "Num_Jets")
        hist['Data']["leading_isolated_muon_pt_after_jet_cut"] = rdfJetAndIsolatedLeptonFiltered['Data'].Histo1D(("leading_isolated_muon_pt_after_jet_cut","Leading Isolated Muon Transverse Momentum After Jet Cut; Pt (GeV);Events",100,20,220),"LeadingIsolatedMuon_pt")
        hist['Data']["leading_isolated_muon_eta_after_jet_cut"] = rdfJetAndIsolatedLeptonFiltered['Data'].Histo1D(("leading_isolated_muon_eta_after_jet_cut", "Leading Isolated Muon Pseudorapidity After Jet Cut; Eta; Events",100,-3,3),"LeadingIsolatedMuon_eta")
        hist['Data']["leading_isolated_muon_phi_after_jet_cut"] = rdfJetAndIsolatedLeptonFiltered['Data'].Histo1D(("leading_isolated_muon_phi_after_jet_cut", "Leading Isolated Muon Angle After Jet Cut; Phi (Radians); Events",100,-3.5,3.5),"LeadingIsolatedMuon_phi")
        hist['Data']["leading_isolated_muon_mass_after_jet_cut"] = rdfJetAndIsolatedLeptonFiltered['Data'].Histo1D(("leading_isolated_muon_mass_after_jet_cut", "Leading Isolated Muon Mass After Jet Cut; Mass(Gev); Events",10,.1,.2),"LeadingIsolatedMuon_mass")
        hist['Data']["leading_isolated_muon_charge_after_jet_cut"] = rdfJetAndIsolatedLeptonFiltered['Data'].Histo1D(("leading_isolated_muon_charge_after_jet_cut", "Leading Isolated Muon Charge After Jet Cut; Charge; Events",5,-2,2),"LeadingIsolatedMuon_charge")
        
        hist['Data']["number_of_jpsi_muons_initial"] = rdfJetAndIsolatedLeptonFiltered['Data'].Histo1D(("number_of_jpsi_muons_initial", "Number of JPsi Muons Before JPsi Muon Cut; Number of Muons; Events",10,0,9), "Num_JPsi_Muons")
        hist['Data']["number_of_jpsi_muons"] = rdfJPsiMuons['Data'].Histo1D(("number_of_jpsi_muons", "Number of JPsi Muons; Number of Muons; Events",10,0,9), "Num_JPsi_Muons")

        
        # Final Results
        hist['Data']["leading_isolated_muon_pt"] = rdfJPsiMuons['Data'].Histo1D(("leading_isolated_muon_pt","Leading Isolated Muon Transverse Momentum; Pt (GeV);Events",100,20,220),"LeadingIsolatedMuon_pt")
        hist['Data']["leading_isolated_muon_eta"] = rdfJPsiMuons['Data'].Histo1D(("leading_isolated_muon_eta", "Leading Isolated Muon Pseudorapidity; Eta; Events",100,-3,3),"LeadingIsolatedMuon_eta")
        hist['Data']["leading_isolated_muon_phi"] = rdfJPsiMuons['Data'].Histo1D(("leading_isolated_muon_phi", "Leading Isolated Muon Angle; Phi (Radians); Events",100,-3.5,3.5),"LeadingIsolatedMuon_phi")
        hist['Data']["leading_isolated_muon_mass"] = rdfJPsiMuons['Data'].Histo1D(("leading_isolated_muon_mass", "Leading Isolated Muon Mass; Mass(Gev); Events",10,.1,.2),"LeadingIsolatedMuon_mass")
        hist['Data']["leading_isolated_muon_charge"] = rdfJPsiMuons['Data'].Histo1D(("leading_isolated_muon_charge", "Leading Isolated Muon Charge; Charge; Events",5,-2,2),"LeadingIsolatedMuon_charge")
        
        hist['Data']["jet1_pt"] = rdfJPsiMuons['Data'].Histo1D(("jet1_pt", "Jet Transverse Momentum for Leading Jet; Pt (GeV); Events", 100, 20, 200), "SJet1_pt")         
        hist['Data']["jet2_pt"] = rdfJPsiMuons['Data'].Histo1D(("jet2_pt", "Jet Transverse Momentum for Subleading Jet; Pt (GeV); Events", 100, 20, 200), "SJet2_pt")
        hist['Data']["jet1_eta"] = rdfJPsiMuons['Data'].Histo1D(("jet1_eta", "Jet Pseudorapidity for Leading Jet; Eta; Events", 100, -3, 3), "SJet1_eta")
        hist['Data']["jet2_eta"] = rdfJPsiMuons['Data'].Histo1D(("jet2_eta", "Jet Pseudorapidity for Subleading Jet; Eta; Events", 100, -3, 3), "SJet2_eta")
        hist['Data']["jet1_phi"] = rdfJPsiMuons['Data'].Histo1D(("jet1_phi", "Jet Angle for Leading Jet; Phi (Radians); Events", 100, -3.5, 3.5), "SJet1_phi")
        hist['Data']["jet2_phi"] = rdfJPsiMuons['Data'].Histo1D(("jet2_phi", "Jet Angle for Subleading Jet; Phi (Radians); Events", 100, -3.5, 3.5), "SJet2_phi")
        hist['Data']["transverse_mass"] = rdfJPsiMuons['Data'].Histo1D(("transverse_mass", "Transverse Mass; Transverse Mass (GeV); Events", 150, 0, 150), "MTofMETandMu")  
        hist['Data']["ht"] = rdfJPsiMuons['Data'].Histo1D(("ht", "Ht; Ht; Events", 300, 0, 1500), "Ht")
               
        hist['Data']["jpsi_muons_pt"] = rdfJPsiMuons['Data'].Histo1D(("jpsi_muons_pt", "Transverse Momentum for JPsi Muons; Pt; Events", 150, 0, 50), "JPsiCandidate_pt")
        hist['Data']["jpsi_muons_eta"] = rdfJPsiMuons['Data'].Histo1D(("jpsi_muons_eta", "Pseudorapidity for JPsi Muons; Eta; Events", 50, -3, 3), "JPsiCandidate_eta")
        hist['Data']["jpsi_muons_phi"] = rdfJPsiMuons['Data'].Histo1D(("jpsi_muons_phi", "Angle for JPsi Muons; Phi; Events", 50, -3.5, 3.5), "JPsiCandidate_phi")
        hist['Data']["jpsi_muons_charge"] = rdfJPsiMuons['Data'].Histo1D(("jpsi_muons_charge", "Charge of JPsi Muons; Charge; Events", 5, -2, 2), "JPsiCandidate_charge")
        hist['Data']["invariant_mass_jpsi_muons"] = rdfJPsiMuons['Data'].Histo1D(("invariant_mass_jpsi_muons", "Invariant Masses for J/Psi Candidate Muons (Oppositely Charged); Invariant Masses; Events", 100, .5, 12), "InvariantMassJPsiMuons")
        hist['Data']["invariant_masses_zoomed"] = rdfJPsiMuons['Data'].Histo1D(("invariant_masses_zoomed", "Invariant Masses for J/Psi Candidate Muons (Oppositely Charged); Invariant Masses; Events", 50, 2.8, 3.4), "InvariantMassJPsiMuons")
        hist['Data']["invariant_mass_jpsi_and_isolated_muons"] = rdfJPsiMuons['Data'].Histo1D(("invariant_mass_jpsi_and_isolated_muons", "Invariant Masses for J/Psi Candidate And Isolated Muons; Invariant Masses; Events", 100, 0, 200), "InvariantMassJPsiAndIsolatedMuons")
        hist['Data']["delta_eta_between_isolated_and_jpsi_muons"] = rdfJPsiMuons['Data'].Histo1D(("delta_eta_between_isolated_and_jpsi_muons", "Delta Eta for Isolated Muon - JPsi Muons; Delta Eta; Events", 50, 0, 6), "DeltaEtaBetweenIsolatedAndJPsiMuon")
        hist['Data']["delta_phi_between_isolated_and_jpsi_muons"] = rdfJPsiMuons['Data'].Histo1D(("delta_phi_between_isolated_and_jpsi_muons", "Delta Phi for Isolated Muon - JPsi Muons; Delta Phi; Events", 50, -3.5, 3.5), "DeltaPhiBetweenIsolatedAndJPsiMuon")
        hist['Data']["delta_r_between_isolated_and_jpsi_muons"] = rdfJPsiMuons['Data'].Histo1D(("delta_r_between_isolated_and_jpsi_muons", "Delta R for Isolated and JPsi Muons; Delta R; Events", 50, 0, 6), "DeltaRBetweenIsolatedAndJPsiMuon")
        hist['Data']["jpsi_pt"] = rdfJPsiMuons['Data'].Histo1D(("jpsi_pt", "Transverse Momentum for JPsi; Pt; Events", 240, 0, 120), "JPsi_pt")
        hist['Data']["jpsi_eta"] = rdfJPsiMuons['Data'].Histo1D(("jpsi_eta", "Pseudorapidity for JPsi; Eta; Events", 50, -3, 3), "JPsi_eta")
        hist['Data']["jpsi_phi"] = rdfJPsiMuons['Data'].Histo1D(("jpsi_phi", "Angle for JPsi; Phi; Events", 50, -3.5, 3.5), "JPsi_phi")
        hist['Data']["delta_eta_between_jpsi_muons"] = rdfJPsiMuons['Data'].Histo1D(("delta_eta_between_jpsi_muons", "Delta Eta for JPsi Muons; Delta Eta; Events", 50, 0, 6), "DeltaEtaBetweenJPsiMuons")
        hist['Data']["delta_phi_between_jpsi_muons"] = rdfJPsiMuons['Data'].Histo1D(("delta_phi_between_jpsi_muons", "Delta Phi for JPsi Muons; Delta Phi; Events", 50, -3.5, 3.5), "DeltaPhiBetweenJPsiMuons")
        hist['Data']["delta_r_between_jpsi_muons"] = rdfJPsiMuons['Data'].Histo1D(("delta_r_between_jpsi_muons", "Delta R for JPsi Muons; Delta R; Events", 50, 0, 6), "DeltaRBetweenJPsiMuons")
        hist['Data']["delta_r_between_jpsi_muons_in_peak"] = rdfJPsiMuons['Data'].Histo1D(("delta_r_between_jpsi_muons_in_peak", "Delta R for JPsi Muons In Peak; Delta R; Events", 50, 0, 1), "DeltaRBetweenJPsiMuonsInPeak")
        hist['Data']["delta_eta_between_isolated_muon_and_jpsi"] = rdfJPsiMuons['Data'].Histo1D(("delta_eta_between_isolated_muon_and_jpsi", "Delta Eta for Isolated Muon And JPsi; Delta Eta; Events", 50, 0, 6), "DeltaEtaBetweenIsolatedMuonAndJPsi")
        hist['Data']["delta_phi_between_isolated_muon_and_jpsi"] = rdfJPsiMuons['Data'].Histo1D(("delta_phi_between_isolated_muon_and_jpsi", "Delta Phi for Isolated Muon And JPsi; Delta Phi; Events", 50, -3.5, 3.5), "DeltaPhiBetweenIsolatedMuonAndJPsi")
        hist['Data']["delta_r_between_isolated_muon_and_jpsi"] = rdfJPsiMuons['Data'].Histo1D(("delta_r_between_isolated_muon_and_jpsi", "Delta R for Isolated Muon And JPsi; Delta R; Events", 50, 0, 6), "DeltaRBetweenIsolatedMuonAndJPsi")
        

        hist['Data']["jpsi_muon_pt"] = rdfJPsiMuons['Data'].Histo1D(("jpsi_muon_pt", "JPsi Muon Pt - Muon Id Medium, Iso Id Tight; JPsi Muon Pt; Events", 100, 0, 100), "JPsiMuon_pt")
        hist['Data']["jpsi_muon_pf_rel_iso_03_all"] = rdfJPsiMuons['Data'].Histo1D(("jpsi_muon_pf_rel_iso_03_all", "JPsi Muon Pf Rel Iso 03 All - Muon Id Medium, Iso Id Tight; JPsi Muon Pf Rel Iso 03 All; Events", 50, 0, .5), "JPsiMuon_pfRelIso03_all")
        hist['Data']["jpsi_muon_pf_iso_id"] = rdfJPsiMuons['Data'].Histo1D(("jpsi_muon_pf_iso_id", "JPsi Muon Pf Iso Id - Muon Id Medium, Iso Id Tight; JPsi Muon Pf Iso Id; Events", 6, .5, 6.5), "JPsiMuon_pfIsoid")
        hist['Data']["isolated_muons_in_range_pfRelIso03_all"] = rdfJPsiMuons['Data'].Histo1D(("isolated_muons_in_range_pfRelIso03_all", "PfRelIso03_All For Isolated Muons In 90-120 GeV Range; PfRelIso03_All; Events", 50, 0, .5), "IsolatedMuonsInRangeInvariantMassPlot_pfRelIso03_all")
        hist['Data']["jpsi_muons_in_range_pfRelIso03_all"] = rdfJPsiMuons['Data'].Histo1D(("jpsi_muons_in_range_pfRelIso03_all", "PfRelIso03_All For JPsi Muons In 90-120 GeV Range; PfRelIso03_All; Events", 50, 0, .5), "JPsiMuonsInRangeInvariantMassPlot_pfRelIso03_all")
        hist['Data']["jpsi_mass_muons_in_range"] = rdfJPsiMuons['Data'].Histo1D(("jpsi_mass_muons_in_range", "JPsi Mass For JPsis In 90-120 GeV Range; JPsi Mass; Events", 200, 0, 100), "JPsiMassesInRangeInvariantMassPlot")
        hist['Data']["pt_at_low_delta_r"] = rdfJPsiMuons['Data'].Histo2D(("pt_at_low_delta_r", "Isolated Muon and JPsi Muon Pt at Delta R < 0.3; Isolated Muon Pt; JPsi Muon Pt; Events", 100, 20, 220, 100, 0, 100), "IsolatedMuonPtAtLowDeltaR", "JPsiMuonPtAtLowDeltaR")
        
        
        hist['Data']["jpsi_mass_in_peak"] = rdfJPsiMuons['Data'].Histo1D(("jpsi_mass_in_peak", "J/Psi Mass Made from Muons In Peak; J/Psi Mass; Events", 80, 2.9, 3.3), "JPsiMassInPeak")
        hist['Data']["delta_eta_for_jpsi_muons_in_peak"] = rdfJPsiMuons['Data'].Histo1D(("delta_eta_for_jpsi_muons_in_peak", "Delta Eta For J/Psi Muons In Peak; Delta Eta; Events", 50, 0, 6), "DeltaEtaForJPsiMuonsInPeak")
        hist['Data']["delta_phi_for_jpsi_muons_in_peak"] = rdfJPsiMuons['Data'].Histo1D(("delta_phi_for_jpsi_muons_in_peak", "Delta Phi For J/Psi Muons In Peak; Delta Phi; Events", 50, -3.5, 3.5), "DeltaPhiForJPsiMuonsInPeak")
        hist['Data']["delta_r_for_jpsi_muons_in_peak"] = rdfJPsiMuons['Data'].Histo1D(("delta_r_for_jpsi_muons_in_peak", "Delta R For J/Psi Muons In Peak; Delta R; Events", 50, 0, 6), "DeltaRForJPsiMuonsInPeak")
        hist['Data']["jpsi_muons_charge_in_peak"] = rdfJPsiMuons['Data'].Histo1D(("jpsi_muons_charge_in_peak", "J/Psi Muons In Peak Charge; Charge; Events", 5, -2, 2), "JPsiMuonsChargeInPeak")
        hist['Data']["jpsi_pt_with_muon_in_peak"] = rdfJPsiMuons['Data'].Histo1D(("jpsi_pt_with_muon_in_peak", "J/Psi Pt With Muons In Peak; Pt; Events", 200, 0, 200), "JPsiPtWithMuonInPeak")
        hist['Data']["jpsi_muon_pts_in_peak"] = rdfJPsiMuons['Data'].Histo2D(("jpsi_muon_pts_in_peak", "J/Psi Muon Pts In Peak; Positive Muon Pt; Negative Muon Pt; Events", 100, 0, 100, 100, 0, 100), "JPsiMuonPositivePt", "JPsiMuonNegativePt")
        hist['Data']["jpsi_muon_isolation_in_peak"] = rdfJPsiMuons['Data'].Histo1D(("jpsi_muon_isolation_in_peak", "J/Psi Pf Rel Iso 03 All Made From Muons In Peak; Pf Rel Iso 03 All; Events", 50, 0, .5), "JPsiMuonIsolationInPeak")
        hist['Data']["top_mass_in_peak"] = rdfJPsiMuons['Data'].Histo1D(("top_mass_in_peak", "Top Quark Mass Made From Muons In Peak; Mass; Events", 300, 0, 300), "TopMassInPeak")
        hist['Data']["delta_eta_between_isolated_muon_and_jpsi_in_peak"] = rdfJPsiMuons['Data'].Histo1D(("delta_eta_between_isolated_muon_and_jpsi_in_peak", "Delta Eta Between Isolated Muon And J/Psi; Delta Eta; Events", 50, 0, 6), "DeltaEtaBetweenIsolatedMuonAndJPsiInPeak")
        hist['Data']["delta_phi_between_isolated_muon_and_jpsi_in_peak"] = rdfJPsiMuons['Data'].Histo1D(("delta_phi_between_isolated_muon_and_jpsi_in_peak", "Delta Phi Between Isolated Muon And J/Psi; Delta Phi; Events", 50, -3.5, 3.5), "DeltaPhiBetweenIsolatedMuonAndJPsiInPeak")
        hist['Data']["delta_r_between_isolated_muon_and_jpsi_in_peak"] = rdfJPsiMuons['Data'].Histo1D(("delta_r_between_isolated_and_jpsi_in_peak", "Delta R Between Isolated Muon And J/Psi; Delta R; Events", 50, 0, 6), "DeltaRBetweenIsolatedMuonAndJPsiInPeak")
#         hist['Data']["isolated_muon_pf_rel_iso_03_all"] = rdfJPsiMuons['Data'].Histo1D(("isolated_muon_pf_rel_iso_03_all", "Isolated Muon PfRelIso03_All; Isolated Muon PfRelIso03_All; Events", 50, 0, .5), "IsolatedMuon_pfRelIso03_all")
#         hist['Data']["isolated_muon_pf_rel_iso_03_all_zoomed"] = rdfJPsiMuons['Data'].Histo1D(("isolated_muon_pf_rel_iso_03_all_zoomed", "Isolated Muon PfRelIso03_All Zoomed; Isolated Muon PfRelIso03_All; Events", 30, 0, .15), "IsolatedMuon_pfRelIso03_all")
#         hist['Data']["isolated_muon_pf_rel_iso_03_chg"] = rdfJPsiMuons['Data'].Histo1D(("isolated_muon_pf_rel_iso_03_chg", "Isolated Muon PfRelIso03_Chg; Isolated Muon PfRelIso03_Chg; Events", 50, 0, .5), "IsolatedMuon_pfRelIso03_chg")
#         hist['Data']["jpsi_muon_positive_pf_rel_iso_03_all"] = rdfJPsiMuons['Data'].Histo1D(("jpsi_muon_positive_pf_rel_iso_03_all", "JPsi Muon Positive PfRelIso03_All; JPsi Muon Positive PfRelIso03_All; Events", 50, 0, .5), "JPsiMuonPositive_pfRelIso03_all")
#         hist['Data']["jpsi_muon_positive_pf_rel_iso_03_all_zoomed"] = rdfJPsiMuons['Data'].Histo1D(("jpsi_muon_positive_pf_rel_iso_03_all_zoomed", "JPsi Muon Positive PfRelIso03_All Zoomed; JPsi Muon Positive PfRelIso03_All; Events", 30, 0, .15), "JPsiMuonPositive_pfRelIso03_all")
#         hist['Data']["jpsi_muon_positive_pf_rel_iso_03_chg"] = rdfJPsiMuons['Data'].Histo1D(("jpsi_muon_positive_pf_rel_iso_03_chg", "JPsi Muon Positive PfRelIso03_Chg; JPsi Muon Positive PfRelIso03_Chg; Events", 50, 0, .5), "JPsiMuonPositive_pfRelIso03_chg")
#         hist['Data']["jpsi_muon_negative_pf_rel_iso_03_all"] = rdfJPsiMuons['Data'].Histo1D(("jpsi_muon_negative_pf_rel_iso_03_all", "JPsi Muon Negative PfRelIso03_All; JPsi Muon Negative PfRelIso03_All; Events", 50, 0, .5), "JPsiMuonNegative_pfRelIso03_all")
#         hist['Data']["jpsi_muon_negative_pf_rel_iso_03_all_zoomed"] = rdfJPsiMuons['Data'].Histo1D(("jpsi_muon_negative_pf_rel_iso_03_all_zoomed", "JPsi Muon Negative PfRelIso03_All Zoomed; JPsi Muon Negative PfRelIso03_All; Events", 30, 0, .15), "JPsiMuonNegative_pfRelIso03_all")
#         hist['Data']["jpsi_muon_negative_pf_rel_iso_03_chg"] = rdfJPsiMuons['Data'].Histo1D(("jpsi_muon_negative_pf_rel_iso_03_chg", "JPsi Muon Negative PfRelIso03_Chg; JPsi Muon Negative PfRelIso03_Chg; Events", 50, 0, .5), "JPsiMuonNegative_pfRelIso03_chg")
        
#         # Medium ID means Medium Muon ID, but still tight or higher Iso ID
#         hist['Data']["isolated_muon_pf_rel_iso_03_all_medium_id"] = rdfJPsiMuons['Data'].Histo1D(("isolated_muon_pf_rel_iso_03_all_medium_id", "Isolated Muon PfRelIso03_All Medium ID; Isolated Muon PfRelIso03_All Medium ID; Events", 50, 0, .5), "IsolatedMuon_pfRelIso03_all_mediumID")
#         hist['Data']["isolated_muon_pf_rel_iso_03_all_medium_id_zoomed"] = rdfJPsiMuons['Data'].Histo1D(("isolated_muon_pf_rel_iso_03_all_medium_id_zoomed", "Isolated Muon PfRelIso03_All Medium ID Zoomed; Isolated Muon PfRelIso03_All Medium ID; Events", 30, 0, .15), "IsolatedMuon_pfRelIso03_all_mediumID")
#         hist['Data']["isolated_muon_pf_rel_iso_03_chg_medium_id"] = rdfJPsiMuons['Data'].Histo1D(("isolated_muon_pf_rel_iso_03_chg_medium_id", "Isolated Muon PfRelIso03_Chg Medium ID; Isolated Muon PfRelIso03_Chg Medium ID; Events", 50, 0, .5), "IsolatedMuon_pfRelIso03_chg_mediumID")
#         hist['Data']["jpsi_muon_positive_pf_rel_iso_03_all_medium_id"] = rdfJPsiMuons['Data'].Histo1D(("jpsi_muon_positive_pf_rel_iso_03_all_medium_id", "JPsi Muon Positive PfRelIso03_All Medium ID; JPsi Muon Positive PfRelIso03_All Medium ID; Events", 50, 0, .5), "JPsiMuonPositive_pfRelIso03_all_mediumID")
#         hist['Data']["jpsi_muon_positive_pf_rel_iso_03_all_medium_id_zoomed"] = rdfJPsiMuons['Data'].Histo1D(("jpsi_muon_positive_pf_rel_iso_03_all_medium_id_zoomed", "JPsi Muon Positive PfRelIso03_All Zoomed Medium ID; JPsi Muon Positive PfRelIso03_All Medium ID; Events", 30, 0, .15), "JPsiMuonPositive_pfRelIso03_all_mediumID")
#         hist['Data']["jpsi_muon_positive_pf_rel_iso_03_chg_medium_id"] = rdfJPsiMuons['Data'].Histo1D(("jpsi_muon_positive_pf_rel_iso_03_chg_medium_id", "JPsi Muon Positive PfRelIso03_Chg Medium ID; JPsi Muon Positive PfRelIso03_Chg Medium ID; Events", 50, 0, .5), "JPsiMuonPositive_pfRelIso03_chg_mediumID")
#         hist['Data']["jpsi_muon_negative_pf_rel_iso_03_all_medium_id"] = rdfJPsiMuons['Data'].Histo1D(("jpsi_muon_negative_pf_rel_iso_03_all_medium_id", "JPsi Muon Negative PfRelIso03_All Medium ID; JPsi Muon Negative PfRelIso03_All Medium ID; Events", 50, 0, .5), "JPsiMuonNegative_pfRelIso03_all_mediumID")
#         hist['Data']["jpsi_muon_negative_pf_rel_iso_03_all_medium_id_zoomed"] = rdfJPsiMuons['Data'].Histo1D(("jpsi_muon_negative_pf_rel_iso_03_all_medium_id_zoomed", "JPsi Muon Negative PfRelIso03_All Zoomed Medium ID; JPsi Muon Negative PfRelIso03_All Medium ID; Events", 30, 0, .15), "JPsiMuonNegative_pfRelIso03_all_mediumID")
#         hist['Data']["jpsi_muon_negative_pf_rel_iso_03_chg_medium_id"] = rdfJPsiMuons['Data'].Histo1D(("jpsi_muon_negative_pf_rel_iso_03_chg_medium_id", "JPsi Muon Negative PfRelIso03_Chg Medium ID; JPsi Muon Negative PfRelIso03_Chg Medium ID; Events", 50, 0, .5), "JPsiMuonNegative_pfRelIso03_chg_mediumID")
    
    
#         hist['Data']["jpsi_muon_positive_peak_only_pf_rel_iso_03_all"] = rdfJPsiMuons['Data'].Histo1D(("jpsi_muon_positive_peak_only_pf_rel_iso_03_all", "JPsi Muon Positive Peak Only PfRelIso03_All; JPsi Muon Positive Peak Only PfRelIso03_All; Events", 50, 0, .5), "JPsiMuonPositive_peak_only_pfRelIso03_all")
#         hist['Data']["jpsi_muon_positive_peak_only_pf_rel_iso_03_chg"] = rdfJPsiMuons['Data'].Histo1D(("jpsi_muon_positive_peak_only_pf_rel_iso_03_chg", "JPsi Muon Positive Peak Only PfRelIso03_Chg; JPsi Muon Positive Peak Only PfRelIso03_Chg; Events", 50, 0, .5), "JPsiMuonPositive_peak_only_pfRelIso03_chg")
#         hist['Data']["jpsi_muon_negative_peak_only_pf_rel_iso_03_all"] = rdfJPsiMuons['Data'].Histo1D(("jpsi_muon_negative_peak_only_pf_rel_iso_03_all", "JPsi Muon Negative Peak Only PfRelIso03_All; JPsi Muon Negative Peak Only PfRelIso03_All; Events", 50, 0, .5), "JPsiMuonNegative_peak_only_pfRelIso03_all")
#         hist['Data']["jpsi_muon_negative_peak_only_pf_rel_iso_03_chg"] = rdfJPsiMuons['Data'].Histo1D(("jpsi_muon_negative_peak_only_pf_rel_iso_03_chg", "JPsi Muon Negative Peak Only PfRelIso03_Chg; JPsi Muon Negative Peak Only PfRelIso03_Chg; Events", 50, 0, .5), "JPsiMuonNegative_peak_only_pfRelIso03_chg")
        
        
        
        

#         hist['Data']["further_muon_medium_tight_pt"] = rdfRemainingIsolatedMuonMuonIdMediumIsoIdTight['Data'].Histo1D(("further_muon_medium_tight_pt", "Further Muon Pt - Muon Id Medium, Iso Id Tight; Further Muon Pt - Muon Id Medium, Iso Id Tight; Events", 100, 20, 220), "FurtherMuonMediumTight_pt")
#         hist['Data']["further_muon_medium_tight_pf_rel_iso_03_all"] = rdfRemainingIsolatedMuonMuonIdMediumIsoIdTight['Data'].Histo1D(("further_muon_medium_tight_pf_rel_iso_03_all", "Further Muon Pf Rel Iso 03 All - Muon Id Medium, Iso Id Tight; Further Muon Pf Rel Iso 03 All - Muon Id Medium, Iso Id Tight; Events", 50, 0, .5), "FurtherMuonMediumTight_pfRelIso03_all")
#         hist['Data']["further_muon_medium_tight_pf_iso_id"] = rdfRemainingIsolatedMuonMuonIdMediumIsoIdTight['Data'].Histo1D(("further_muon_medium_tight_pf_iso_id", "Further Muon Pf Iso Id - Muon Id Medium, Iso Id Tight; Further Muon Pf Iso Id - Muon Id Medium, Iso Id Tight; Events", 6, 1, 6), "FurtherMuonMediumTight_pfIsoid")
        
#         hist['Data']["further_muon_loose_tight_pt"] = rdfRemainingIsolatedMuonMuonIdLooseIsoIdTight['Data'].Histo1D(("further_muon_loose_tight_pt", "Further Muon Pt - Muon Id Loose Iso Id Tight; Further Muon Pt - Muon Id Loose, Iso Id Tight; Events", 100, 20, 220), "FurtherMuonLooseTight_pt")
#         hist['Data']["further_muon_loose_tight_pf_rel_iso_03_all"] = rdfRemainingIsolatedMuonMuonIdLooseIsoIdTight['Data'].Histo1D(("further_muon_loose_tight_pf_rel_iso_03_all", "Further Muon Pf Rel Iso 03 All - Muon Id Loose, Iso Id Tight; Further Muon Pf Rel Iso 03 All - Muon Id Loose, Iso Id Tight; Events", 50, 0, .5), "FurtherMuonLooseTight_pfRelIso03_all")
#         hist['Data']["further_muon_loose_tight_pf_iso_id"] = rdfRemainingIsolatedMuonMuonIdLooseIsoIdTight['Data'].Histo1D(("further_muon_loose_tight_pf_iso_id", "Further Muon Pf Iso Id - Muon Id Loose, Iso Id Tight; Further Muon Pf Iso Id - Muon Id Loose, Iso Id Tight; Events", 6, 1, 6), "FurtherMuonLooseTight_pfIsoid")
        
#         hist['Data']["further_muon_medium_medium_pt"] = rdfRemainingIsolatedMuonMuonIdMediumIsoIdMedium['Data'].Histo1D(("further_muon_medium_medium_pt", "Further Muon Pt - Muon Id Medium, Iso Id Medium; Further Muon Pt - Muon Id Medium, Iso Id Medium; Events", 100, 20, 220), "FurtherMuonMediumMedium_pt")
#         hist['Data']["further_muon_medium_medium_pf_rel_iso_03_all"] = rdfRemainingIsolatedMuonMuonIdMediumIsoIdMedium['Data'].Histo1D(("further_muon_medium_medium_pf_rel_iso_03_all", "Further Muon Pf Rel Iso 03 All - Muon Id Medium, Iso Id Medium; Further Muon Pf Rel Iso 03 All - Muon Id Medium, Iso Id Medium; Events", 50, 0, .5), "FurtherMuonMediumMedium_pfRelIso03_all")
#         hist['Data']["further_muon_medium_medium_pf_iso_id"] = rdfRemainingIsolatedMuonMuonIdMediumIsoIdMedium['Data'].Histo1D(("further_muon_medium_medium_pf_iso_id", "Further Muon Pf Iso Id - Muon Id Medium, Iso Id Medium; Further Muon Pf Iso Id - Muon Id Medium, Iso Id Medium; Events", 6, 1, 6), "FurtherMuonMediumMedium_pfIsoid")
        
#         hist['Data']["further_muon_loose_medium_pt"] = rdfRemainingIsolatedMuonMuonIdLooseIsoIdMedium['Data'].Histo1D(("further_muon_loose_medium_pt", "Further Muon Pt - Muon Id Loose, Iso Id Medium; Further Muon Pt - Muon Id Loose, Iso Id Medium; Events", 100, 20, 220), "FurtherMuonLooseMedium_pt")
#         hist['Data']["further_muon_loose_medium_pf_rel_iso_03_all"] = rdfRemainingIsolatedMuonMuonIdLooseIsoIdMedium['Data'].Histo1D(("further_muon_loose_medium_pf_rel_iso_03_all", "Further Muon Pf Rel Iso 03 All - Muon Id Loose, Iso Id Medium; Further Muon Pf Rel Iso 03 All - Muon Id Loose, Iso Id Medium; Events", 50, 0, .5), "FurtherMuonLooseMedium_pfRelIso03_all")
#         hist['Data']["further_muon_loose_medium_pf_iso_id"] = rdfRemainingIsolatedMuonMuonIdLooseIsoIdMedium['Data'].Histo1D(("further_muon_loose_medium_pf_iso_id", "Further Muon Pf Iso Id - Muon Id Loose, Iso Id Medium; Further Muon Pf Iso Id - Muon Id Loose, Iso Id Medium; Events", 6, 1, 6), "FurtherMuonLooseMedium_pfIsoid")
        
        
        
        
        
        
        

    else:
        
        hist[sample]["leading_isolated_muon_pt_initial"] = rdfIsolatedMuonNoHighWeights[sample].Histo1D((sample + "_" + "leading_isolated_muon_pt_initial", "Monte Carlo " + sample + ";Leading Isolated Muon Transverse Momentum (One Muon, No Electrons); Pt (GeV);Events",100,20,220),"LeadingIsolatedMuon_pt","LumiXS")
        hist[sample]["leading_isolated_muon_eta_initial"] = rdfIsolatedMuonNoHighWeights[sample].Histo1D((sample + "_" + "leading_isolated_muon_eta_initial", "Monte Carlo " + sample + ";Leading Isolated Muon Pseudorapidity (One Muon, No Electrons); Eta; Events",100,-3,3),"LeadingIsolatedMuon_eta","LumiXS")
        hist[sample]["leading_isolated_muon_phi_initial"] = rdfIsolatedMuonNoHighWeights[sample].Histo1D((sample + "_" + "leading_isolated_muon_phi_initial", "Monte Carlo " + sample + ";Leading Isolated Muon Angle (One Muon, No Electrons); Phi (Radians); Events",100,-3.5,3.5),"LeadingIsolatedMuon_phi","LumiXS")
        hist[sample]["leading_isolated_muon_mass_initial"] = rdfIsolatedMuonNoHighWeights[sample].Histo1D((sample + "_" + "leading_isolated_muon_mass_initial", "Monte Carlo " + sample + ";Leading Isolated Muon Mass (One Muon, No Electrons); Mass(Gev); Events",10,.1,.2),"LeadingIsolatedMuon_mass","LumiXS")
        hist[sample]["leading_isolated_muon_charge_initial"] = rdfIsolatedMuonNoHighWeights[sample].Histo1D((sample + "_" + "leading_isolated_muon_charge_initial", "Monte Carlo " + sample + ";Leading Isolated Muon Charge (One Muon, No Electrons); Charge; Events",5,-2,2),"LeadingIsolatedMuon_charge","LumiXS")
        hist[sample]["leading_isolated_electron_pt_initial"] = rdfIsolatedElectronNoHighWeights[sample].Histo1D((sample + "_" + "leading_isolated_electron_pt_initial", "Monte Carlo " + sample + ";Leading Isolated Electron Transverse Momentum (One Electron, No Muons); Pt (GeV);Events",100,20,220),"LeadingIsolatedElectron_pt","LumiXS")
        hist[sample]["leading_isolated_electron_eta_initial"] = rdfIsolatedElectronNoHighWeights[sample].Histo1D((sample + "_" + "leading_isolated_electron_eta_initial", "Monte Carlo " + sample + ";Leading Isolated Electron Pseudorapidity (One Electron, No Muons); Eta; Events",100,-3,3),"LeadingIsolatedElectron_eta","LumiXS")
        hist[sample]["leading_isolated_electron_phi_initial"] = rdfIsolatedElectronNoHighWeights[sample].Histo1D((sample + "_" + "leading_isolated_electron_phi_initial", "Monte Carlo " + sample + ";Leading Isolated Electron Angle (One Electron, No Muons); Phi (Radians); Events",100,-3.5,3.5),"LeadingIsolatedElectron_phi","LumiXS")
        hist[sample]["leading_isolated_electron_mass_initial"] = rdfIsolatedElectronNoHighWeights[sample].Histo1D((sample + "_" + "leading_isolated_electron_mass_initial", "Monte Carlo " + sample + ";Leading Isolated Electron Mass (One Electron, No Muons); Mass(Gev); Events",10,.0001,.001),"LeadingIsolatedElectron_mass","LumiXS")
        hist[sample]["leading_isolated_electron_charge_initial"] = rdfIsolatedElectronNoHighWeights[sample].Histo1D((sample + "_" + "leading_isolated_electron_charge_initial", "Monte Carlo " + sample + ";Leading Isolated Electron Charge (One Electron, No Muons); Charge; Events",5,-2,2),"LeadingIsolatedElectron_charge","LumiXS")
        
        hist[sample]["met_before_met_cut"] = rdfIsolatedMuonNoHighWeights[sample].Histo1D((sample + "_" + "met_before_met_cut", "Monte Carlo " + sample + "; MET Before MET Cut; Pt (GeV); Events",100,0,250), "METBeforeMETCut", "LumiXS")
        hist[sample]["met_after_met_cut"] = rdfIsolatedMuonAfterMETCut[sample].Histo1D((sample + "_" + "met_after_met_cut", "Monte Carlo " + sample + "; MET After MET Cut; Pt (GeV); Events",100,0,250), "METAfterMETCut", "LumiXS")
        hist[sample]["leading_isolated_muon_pt_after_met_cut"] = rdfIsolatedMuonAfterMETCut[sample].Histo1D((sample + "_" + "leading_isolated_muon_pt_after_met_cut", "Monte Carlo " + sample + ";Leading Isolated Muon Transverse Momentum After MET Cut; Pt (GeV);Events",100,20,220),"LeadingIsolatedMuon_pt","LumiXS")
        hist[sample]["leading_isolated_muon_eta_after_met_cut"] = rdfIsolatedMuonAfterMETCut[sample].Histo1D((sample + "_" + "leading_isolated_muon_eta_after_met_cut", "Monte Carlo " + sample + ";Leading Isolated Muon Pseudorapidity After MET Cut; Eta; Events",100,-3,3),"LeadingIsolatedMuon_eta","LumiXS")
        hist[sample]["leading_isolated_muon_phi_after_met_cut"] = rdfIsolatedMuonAfterMETCut[sample].Histo1D((sample + "_" + "leading_isolated_muon_phi_after_met_cut", "Monte Carlo " + sample + ";Leading Isolated Muon Angle After MET Cut; Phi (Radians); Events",100,-3.5,3.5),"LeadingIsolatedMuon_phi","LumiXS")
        hist[sample]["leading_isolated_muon_mass_after_met_cut"] = rdfIsolatedMuonAfterMETCut[sample].Histo1D((sample + "_" + "leading_isolated_muon_mass_after_met_cut", "Monte Carlo " + sample + ";Leading Isolated Muon Mass After MET Cut; Mass(Gev); Events",10,.1,.2),"LeadingIsolatedMuon_mass","LumiXS")
        hist[sample]["leading_isolated_muon_charge_after_met_cut"] = rdfIsolatedMuonAfterMETCut[sample].Histo1D((sample + "_" + "leading_isolated_muon_charge_after_met_cut", "Monte Carlo " + sample + ";Leading Isolated Muon Charge After MET Cut; Charge; Events",5,-2,2),"LeadingIsolatedMuon_charge","LumiXS")
        hist[sample]["leading_isolated_electron_pt_after_met_cut"] = rdfIsolatedElectronAfterMETCut[sample].Histo1D((sample + "_" + "leading_isolated_electron_pt_after_met_cut", "Monte Carlo " + sample + ";Leading Isolated Electron Transverse Momentum After MET Cut; Pt (GeV);Events",100,20,220),"LeadingIsolatedElectron_pt","LumiXS")
        hist[sample]["leading_isolated_electron_eta_after_met_cut"] = rdfIsolatedElectronAfterMETCut[sample].Histo1D((sample + "_" + "leading_isolated_electron_eta_after_met_cut", "Monte Carlo " + sample + ";Leading Isolated Electron Pseudorapidity After MET Cut; Eta; Events",100,-3,3),"LeadingIsolatedElectron_eta","LumiXS")
        hist[sample]["leading_isolated_electron_phi_after_met_cut"] = rdfIsolatedElectronAfterMETCut[sample].Histo1D((sample + "_" + "leading_isolated_electron_phi_after_met_cut", "Monte Carlo " + sample + ";Leading Isolated Electron Angle After MET Cut; Phi (Radians); Events",100,-3.5,3.5),"LeadingIsolatedElectron_phi","LumiXS")
        hist[sample]["leading_isolated_electron_mass_after_met_cut"] = rdfIsolatedElectronAfterMETCut[sample].Histo1D((sample + "_" + "leading_isolated_electron_mass_after_met_cut", "Monte Carlo " + sample + ";Leading Isolated Electron Mass After MET Cut; Mass(Gev); Events",10,.0001,.001),"LeadingIsolatedElectron_mass","LumiXS")
        hist[sample]["leading_isolated_electron_charge_after_met_cut"] = rdfIsolatedElectronAfterMETCut[sample].Histo1D((sample + "_" + "leading_isolated_electron_charge_after_met_cut", "Monte Carlo " + sample + ";Leading Isolated Electron Charge After MET Cut; Charge; Events",5,-2,2),"LeadingIsolatedElectron_charge","LumiXS")
        
        hist[sample]["number_of_jets_initial"] = rdfIsolatedMuonAfterMETCut[sample].Histo1D((sample + "_" + "number_of_jets_initial", "Monte Carlo " + sample +"; Number Of Jets Before Jet Cut; Number Of Jets; Events", 20, 0, 20), "Num_Jets", "LumiXS")
        hist[sample]["number_of_jets"] = rdfJetAndIsolatedLeptonFiltered[sample].Histo1D((sample + "_" + "number_of_jets", "Monte Carlo " + sample +"; Number Of Jets; Number Of Jets; Events", 20, 0, 20), "Num_Jets", "LumiXS")
        hist[sample]["leading_isolated_muon_pt_after_jet_cut"] = rdfJetAndIsolatedLeptonFiltered[sample].Histo1D((sample + "_" + "leading_isolated_muon_pt_after_jet_cut", "Monte Carlo " + sample + ";Leading Isolated Muon Transverse Momentum After Jet Cut; Pt (GeV);Events",100,20,220),"LeadingIsolatedMuon_pt","LumiXS")
        hist[sample]["leading_isolated_muon_eta_after_jet_cut"] = rdfJetAndIsolatedLeptonFiltered[sample].Histo1D((sample + "_" + "leading_isolated_muon_eta_after_jet_cut", "Monte Carlo " + sample + ";Leading Isolated Muon Pseudorapidity After Jet Cut; Eta; Events",100,-3,3),"LeadingIsolatedMuon_eta","LumiXS")
        hist[sample]["leading_isolated_muon_phi_after_jet_cut"] = rdfJetAndIsolatedLeptonFiltered[sample].Histo1D((sample + "_" + "leading_isolated_muon_phi_after_jet_cut", "Monte Carlo " + sample + ";Leading Isolated Muon Angle After Jet Cut; Phi (Radians); Events",100,-3.5,3.5),"LeadingIsolatedMuon_phi","LumiXS")
        hist[sample]["leading_isolated_muon_mass_after_jet_cut"] = rdfJetAndIsolatedLeptonFiltered[sample].Histo1D((sample + "_" + "leading_isolated_muon_mass_after_jet_cut", "Monte Carlo " + sample + ";Leading Isolated Muon Mass After Jet Cut; Mass(Gev); Events",10,.1,.2),"LeadingIsolatedMuon_mass","LumiXS")
        hist[sample]["leading_isolated_muon_charge_after_jet_cut"] = rdfJetAndIsolatedLeptonFiltered[sample].Histo1D((sample + "_" + "leading_isolated_muon_charge_after_jet_cut", "Monte Carlo " + sample + ";Leading Isolated Muon Charge After Jet Cut; Charge; Events",5,-2,2),"LeadingIsolatedMuon_charge","LumiXS")
        
        hist[sample]["number_of_jpsi_muons_initial"] = rdfJetAndIsolatedLeptonFiltered[sample].Histo1D((sample + "_" + "number_of_jpsi_muons_initial", "Monte Carlo " + sample +"; Number Of JPsi Muons Before JPsi Muon Cut; Number Of Muons; Events", 10, 0, 9), "Num_JPsi_Muons", "LumiXS")
        hist[sample]["number_of_jpsi_muons"] = rdfJPsiMuons[sample].Histo1D((sample + "_" + "number_of_jpsi_muons", "Monte Carlo " + sample +"; Number Of JPsi Muons; Number Of Muons; Events", 10, 0, 9), "Num_JPsi_Muons", "LumiXS")
        
        
        
        # Final Results
        hist[sample]["leading_isolated_muon_pt"] = rdfJPsiMuons[sample].Histo1D((sample + "_" + "leading_isolated_muon_pt", "Monte Carlo " + sample + ";Leading Isolated Muon Transverse Momentum (One Muon, No Electrons); Pt (GeV);Events",100,20,220),"LeadingIsolatedMuon_pt","LumiXS")
        hist[sample]["leading_isolated_muon_eta"] = rdfJPsiMuons[sample].Histo1D((sample + "_" + "leading_isolated_muon_eta", "Monte Carlo " + sample + ";Leading Isolated Muon Pseudorapidity (One Muon, No Electrons); Eta; Events",100,-3,3),"LeadingIsolatedMuon_eta","LumiXS")
        hist[sample]["leading_isolated_muon_phi"] = rdfJPsiMuons[sample].Histo1D((sample + "_" + "leading_isolated_muon_phi", "Monte Carlo " + sample + ";Leading Isolated Muon Angle (One Muon, No Electrons); Phi (Radians); Events",100,-3.5,3.5),"LeadingIsolatedMuon_phi","LumiXS")
        hist[sample]["leading_isolated_muon_mass"] = rdfJPsiMuons[sample].Histo1D((sample + "_" + "leading_isolated_muon_mass", "Monte Carlo " + sample + ";Leading Isolated Muon Mass (One Muon, No Electrons); Mass(Gev); Events",10,.1,.2),"LeadingIsolatedMuon_mass","LumiXS")
        hist[sample]["leading_isolated_muon_charge"] = rdfJPsiMuons[sample].Histo1D((sample + "_" + "leading_isolated_muon_charge", "Monte Carlo " + sample + ";Leading Isolated Muon Charge (One Muon, No Electrons); Charge; Events",5,-2,2),"LeadingIsolatedMuon_charge","LumiXS")
               
        hist[sample]["jet1_pt"] = rdfJPsiMuons[sample].Histo1D((sample + "_" + "jet1_pt", "Monte Carlo " + sample + "; Jet Transverse Momentum for Leading Jet; Pt (GeV); Events", 100, 20, 200), "SJet1_pt", "LumiXS")
        hist[sample]["jet2_pt"] = rdfJPsiMuons[sample].Histo1D((sample + "_" + "jet2_pt", "Monte Carlo " + sample + "; Jet Transverse Momentum for Subeading Jet; Pt (GeV); Events", 100, 20, 200), "SJet2_pt", "LumiXS")
        hist[sample]["jet1_eta"] = rdfJPsiMuons[sample].Histo1D((sample + "_" + "jet1_eta", "Monte Carlo " + sample + "; Jet Pseudorapidity for Leading Jet; Eta; Events", 100, -3, 3), "SJet1_eta", "LumiXS")
        hist[sample]["jet2_eta"] = rdfJPsiMuons[sample].Histo1D((sample + "_" + "jet2_eta", "Monte Carlo " + sample +"; Jet Pseudorapidity for Subeading Jet; Eta; Events", 100, -3, 3), "SJet2_eta", "LumiXS")
        hist[sample]["jet1_phi"] = rdfJPsiMuons[sample].Histo1D((sample + "_" + "jet1_phi", "Monte Carlo " + sample + "; Jet Angle for Leading Jet; Phi (Radians); Events", 100, -3.5, 3.5), "SJet1_phi", "LumiXS")
        hist[sample]["jet2_phi"] = rdfJPsiMuons[sample].Histo1D((sample + "_" + "jet2_phi", "Monte Carlo " + sample +"; Jet Angle for Subleading Jet; Phi (Radians); Events", 100, -3.5, 3.5), "SJet2_phi", "LumiXS")
        hist[sample]["transverse_mass"] = rdfJPsiMuons[sample].Histo1D((sample + "_" + "transverse_mass", "Monte Carlo " + sample +"; Transverse Mass; Transverse Mass (GeV); Events", 150, 0, 150), "MTofMETandMu", "LumiXS")
        hist[sample]["ht"] = rdfJPsiMuons[sample].Histo1D((sample + "_" + "ht", "Monte Carlo " + sample + "; Ht; Ht; Events", 300, 0, 1500), "Ht", "LumiXS")
                
        hist[sample]["jpsi_muons_pt"] = rdfJPsiMuons[sample].Histo1D((sample + "_" + "jpsi_muons_pt", "Monte Carlo " + sample + "; Transverse Momentum for JPsi Muons; Pt; Events", 150, 0, 50), "JPsiCandidate_pt", "LumiXS")
        hist[sample]["jpsi_muons_eta"] = rdfJPsiMuons[sample].Histo1D((sample + "_" + "jpsi_muons_eta", "Monte Carlo " + sample + "; Pseudorapidity for JPsi Muons; Eta; Events", 50, -3, 3), "JPsiCandidate_eta", "LumiXS")
        hist[sample]["jpsi_muons_phi"] = rdfJPsiMuons[sample].Histo1D((sample + "_" + "jpsi_muons_phi", "Monte Carlo " + sample + "; Angle for JPsi Muons; Phi; Events", 50, -3.5, 3.5), "JPsiCandidate_phi", "LumiXS")
        hist[sample]["jpsi_muons_charge"] = rdfJPsiMuons[sample].Histo1D((sample + "_" + "jpsi_muon_charge", "Monte Carlo" + sample + "; Charge of JPsi Muons; Charge; Events", 5, -2, 2), "JPsiCandidate_charge", "LumiXS")
        hist[sample]["invariant_mass_jpsi_muons"] = rdfJPsiMuons[sample].Histo1D((sample + "_" + "invariant_mass_jpsi_muons", "Monte Carlo " + sample + "; Invariant Masses for J/Psi Candidate Muons (Oppositely Charged); Invariant Masses; Events", 100, .5, 12), "InvariantMassJPsiMuons", "LumiXS")
        hist[sample]["invariant_masses_zoomed"] = rdfJPsiMuons[sample].Histo1D((sample + "_" + "invariant_masses_zoomed", "Monte Carlo " + sample + "; Invariant Masses for J/Psi Candidate Muons (Oppositely Charged); Invariant Masses; Events", 50, 2.8, 3.4), "InvariantMassJPsiMuons", "LumiXS")
        hist[sample]["invariant_mass_jpsi_and_isolated_muons"] = rdfJPsiMuons[sample].Histo1D((sample + "_" + "invariant_mass_jpsi_and_isolated_muons", "Monte Carlo " + sample + "; Invariant Masses for J/Psi Candidate And Isolated Muons; Invariant Masses; Events", 100, 0, 200), "InvariantMassJPsiAndIsolatedMuons", "LumiXS")
        hist[sample]["delta_eta_between_isolated_and_jpsi_muons"] = rdfJPsiMuons[sample].Histo1D((sample + "_" + "delta_eta_betweeen_isolated_and_jpsi_muons", "Monte Carlo " + sample + "; Delta Eta for Isolated Muon - JPsi Muons; Delta Eta; Events", 50, 0, 6), "DeltaEtaBetweenIsolatedAndJPsiMuon", "LumiXS")
        hist[sample]["delta_phi_between_isolated_and_jpsi_muons"] = rdfJPsiMuons[sample].Histo1D((sample + "_" + "delta_phi_betweeen_isolated_and_jpsi_muons", "Monte Carlo " + sample + "; Delta Phi for Isolated Muon - JPsi Muons; Delta Phi; Events", 50, -3.5, 3.5), "DeltaPhiBetweenIsolatedAndJPsiMuon", "LumiXS")
        hist[sample]["delta_r_between_isolated_and_jpsi_muons"] = rdfJPsiMuons[sample].Histo1D((sample + "_" + "delta_r_betweeen_isolated_and_jpsi_muons", "Monte Carlo " + sample + "; Delta R for Isolated and JPsi Muons; Delta R; Events", 50, 0, 6), "DeltaRBetweenIsolatedAndJPsiMuon", "LumiXS")
        hist[sample]["jpsi_pt"] = rdfJPsiMuons[sample].Histo1D((sample + "_" + "jpsi_pt", "Monte Carlo " + sample + "; Transverse Momentum for JPsi; Pt; Events", 240, 0, 120), "JPsi_pt", "LumiXS")
        hist[sample]["jpsi_eta"] = rdfJPsiMuons[sample].Histo1D((sample + "_" + "jpsi_eta", "Monte Carlo " + sample + "; Pseudorapidity for JPsi; Eta; Events", 50, -3, 3), "JPsi_eta", "LumiXS")
        hist[sample]["jpsi_phi"] = rdfJPsiMuons[sample].Histo1D((sample + "_" + "jpsi_phi", "Monte Carlo " + sample + "; Angle for JPsi; Phi; Events", 50, -3.5, 3.5), "JPsi_phi", "LumiXS")
        hist[sample]["delta_eta_between_jpsi_muons"] = rdfJPsiMuons[sample].Histo1D((sample + "_" + "delta_eta_between_jpsi_muons", "Monte Carlo " + sample + "; Delta Eta for JPsi Muons; Delta Eta; Events", 50, 0, 6), "DeltaEtaBetweenJPsiMuons", "LumiXS")
        hist[sample]["delta_phi_between_jpsi_muons"] = rdfJPsiMuons[sample].Histo1D((sample + "_" + "delta_phi_between_jpsi_muons", "Monte Carlo " + sample + "; Delta Phi for JPsi Muons; Delta Phi; Events", 50, -3.5, 3.5), "DeltaPhiBetweenJPsiMuons", "LumiXS")
        hist[sample]["delta_r_between_jpsi_muons"] = rdfJPsiMuons[sample].Histo1D((sample + "_" + "delta_r_between_jpsi_muons", "Monte Carlo " + sample + "; Delta R for JPsi Muons; Delta R; Events", 50, 0, 6), "DeltaRBetweenJPsiMuons", "LumiXS")
        hist[sample]["delta_r_between_jpsi_muons_in_peak"] = rdfJPsiMuons[sample].Histo1D((sample + "_" + "delta_r_between_jpsi_muons_in_peak", "Monte Carlo " + sample + "; Delta R for JPsi Muons In Peak; Delta R; Events", 50, 0, 1), "DeltaRBetweenJPsiMuonsInPeak", "LumiXS")
        hist[sample]["delta_eta_between_isolated_muon_and_jpsi"] = rdfJPsiMuons[sample].Histo1D((sample + "_" + "delta_eta_between_isolated_muon_and_jpsi", "Monte Carlo " + sample + "; Delta Eta for Isolated Muon And JPsi; Delta Eta; Events", 50, 0, 6), "DeltaEtaBetweenIsolatedMuonAndJPsi", "LumiXS")
        hist[sample]["delta_phi_between_isolated_muon_and_jpsi"] = rdfJPsiMuons[sample].Histo1D((sample + "_" + "delta_phi_between_isolated_muon_and_jpsi", "Monte Carlo " + sample + "; Delta Phi for Isolated Muon And JPsi; Delta Phi; Events", 50, -3.5, 3.5), "DeltaPhiBetweenIsolatedMuonAndJPsi", "LumiXS")
        hist[sample]["delta_r_between_isolated_muon_and_jpsi"] = rdfJPsiMuons[sample].Histo1D((sample + "_" + "delta_r_between_isolated_muon_and_jpsi", "Monte Carlo " + sample + "; Delta R for Isolated Muon And JPsi; Delta R; Events", 50, 0, 6), "DeltaRBetweenIsolatedMuonAndJPsi", "LumiXS")
        
        
        hist[sample]["jpsi_muon_pt"] = rdfJPsiMuons[sample].Histo1D((sample + "_" + "jpsi_muon_pt", "Monte Carlo" + sample + "; JPsi Muon Pt - Muon Id Medium, Iso Id Tight; JPsi Muon Pt; Events", 100, 00, 100), "JPsiMuon_pt", "LumiXS")
        hist[sample]["jpsi_muon_pf_rel_iso_03_all"] = rdfJPsiMuons[sample].Histo1D((sample + "_" + "jpsi_muon_pf_rel_iso_03_all", "Monte Carlo" + sample + "; JPsi Muon Pf Rel Iso 03 All - Muon Id Medium, Iso Id Tight; JPsi Muon Pf Rel Iso 03 All; Events", 50, 0, .5), "JPsiMuon_pfRelIso03_all", "LumiXS")
        hist[sample]["jpsi_muon_pf_iso_id"] = rdfJPsiMuons[sample].Histo1D((sample + "_" + "jpsi_muon_pf_iso_id", "Monte Carlo" + sample + "; JPsi Muon Pf Iso Id - Muon Id Medium, Iso Id Tight; JPsi Muon Pf Iso Id; Events", 6, .5, 6.5), "JPsiMuon_pfIsoid", "LumiXS")
        hist[sample]["isolated_muons_in_range_pfRelIso03_all"] = rdfJPsiMuons[sample].Histo1D((sample + "_" + "isolated_muons_in_range_pfRelIso03_all", "Monte Carlo" + sample + "; PfRelIso03_All For Isolated Muons In 90-120 GeV Range; PfRelIso03_All; Events", 50, 0, .5), "IsolatedMuonsInRangeInvariantMassPlot_pfRelIso03_all", "LumiXS")
        hist[sample]["jpsi_muons_in_range_pfRelIso03_all"] = rdfJPsiMuons[sample].Histo1D((sample + "_" + "jpsi_muons_in_range_pfRelIso03_all", "Monte Carlo" + sample + "; PfRelIso03_All For JPsi Muons In 90-120 GeV Range; PfRelIso03_All; Events", 50, 0, .5), "JPsiMuonsInRangeInvariantMassPlot_pfRelIso03_all", "LumiXS")
        hist[sample]["jpsi_mass_muons_in_range"] = rdfJPsiMuons[sample].Histo1D((sample + "_" + "jpsi_mass_muons_in_range", "Monte Carlo" + sample + "; JPsi Mass For JPsis In 90-120 GeV Range; JPsi Mass; Events", 200, 0, 100), "JPsiMassesInRangeInvariantMassPlot")
        hist[sample]["pt_at_low_delta_r"] = rdfJPsiMuons[sample].Histo2D((sample + "_" + "pt_at_low_delta_r", "Monte Carlo" + sample + "; Isolated Muon and JPsi Muon Pt at Delta R < 0.3; Isolated Muon Pt; JPsi Muon Pt; Events", 100, 20, 220, 100, 0, 100), "IsolatedMuonPtAtLowDeltaR", "JPsiMuonPtAtLowDeltaR", "LumiXS")
        
        
        hist[sample]["jpsi_mass_in_peak"] = rdfJPsiMuons[sample].Histo1D((sample + "_" + "jpsi_mass_in_peak", "Monte Carlo" + sample + "; J/Psi Mass Made from Muons In Peak; J/Psi Mass; Events", 80, 2.9, 3.3), "JPsiMassInPeak", "LumiXS")
        hist[sample]["delta_eta_for_jpsi_muons_in_peak"] = rdfJPsiMuons[sample].Histo1D((sample + "_" + "delta_eta_for_jpsi_muons_in_peak", "Monte Carlo" + sample + "; Delta Eta For J/Psi Muons In Peak; Delta Eta; Events", 50, 0, 6), "DeltaEtaForJPsiMuonsInPeak", "LumiXS")
        hist[sample]["delta_phi_for_jpsi_muons_in_peak"] = rdfJPsiMuons[sample].Histo1D((sample + "_" + "delta_phi_for_jpsi_muons_in_peak", "Monte Carlo" + sample + "; Delta Phi For J/Psi Muons In Peak; Delta Phi; Events", 50, -3.5, 3.5), "DeltaPhiForJPsiMuonsInPeak", "LumiXS")
        hist[sample]["delta_r_for_jpsi_muons_in_peak"] = rdfJPsiMuons[sample].Histo1D((sample + "_" + "delta_r_for_jpsi_muons_in_peak", "Monte Carlo" + sample + "; Delta R For J/Psi Muons In Peak; Delta R; Events", 50, 0, 6), "DeltaRForJPsiMuonsInPeak", "LumiXS")
        hist[sample]["jpsi_muons_charge_in_peak"] = rdfJPsiMuons[sample].Histo1D((sample + "_" + "jpsi_muons_charge_in_peak", "Monte Carlo" + sample + "; J/Psi Muons In Peak Charge; Charge; Events", 5, -2, 2), "JPsiMuonsChargeInPeak", "LumiXS")
        hist[sample]["jpsi_pt_with_muon_in_peak"] = rdfJPsiMuons[sample].Histo1D((sample + "_" + "jpsi_pt_with_muon_in_peak", "Monte Carlo" + sample + "; J/Psi Pt With Muons In Peak; Pt; Events", 200, 0, 200), "JPsiPtWithMuonInPeak", "LumiXS")
        hist[sample]["jpsi_muon_pts_in_peak"] = rdfJPsiMuons[sample].Histo2D((sample + "_" + "jpsi_muon_pts_in_peak", "Monte Carlo" + sample + "; J/Psi Muon Pts In Peak; Positive Muon Pts; Negative Muon Pts; Events", 200, 0, 200, 100, 0, 100), "JPsiMuonPositivePt", "JPsiMuonNegativePt", "LumiXS")
        hist[sample]["jpsi_muon_isolation_in_peak"] = rdfJPsiMuons[sample].Histo1D((sample + "_" + "jpsi_muon_isolation_in_peak", "Monte Carlo" + sample + "; J/Psi Pf Rel Iso 03 All Made From Muons In Peak; Pf Rel Iso 03 All; Events", 50, 0, .5), "JPsiMuonIsolationInPeak", "LumiXS")
        hist[sample]["top_mass_in_peak"] = rdfJPsiMuons[sample].Histo1D((sample + "_" + "top_mass_in_peak", "Monte Carlo" + sample + "; Top Quark Mass Made From Muons In Peak; Mass; Events", 300, 0, 300), "TopMassInPeak", "LumiXS")
        hist[sample]["delta_eta_between_isolated_muon_and_jpsi_in_peak"] = rdfJPsiMuons[sample].Histo1D((sample + "_" + "delta_eta_between_isolated_muon_and_jpsi_in_peak", "Monte Carlo" + sample + "; Delta Eta Between Isolated Muon And J/Psi; Delta Eta; Events", 50, 0, 6), "DeltaEtaBetweenIsolatedMuonAndJPsiInPeak", "LumiXS")
        hist[sample]["delta_phi_between_isolated_muon_and_jpsi_in_peak"] = rdfJPsiMuons[sample].Histo1D((sample + "_" + "delta_phi_between_isolated_muon_and_jpsi_in_peak", "Monte Carlo" + sample + "; Delta Phi Between Isolated Muon And J/Psi; Delta Phi; Events", 50, -3.5, 3.5), "DeltaPhiBetweenIsolatedMuonAndJPsiInPeak", "LumiXS")
        hist[sample]["delta_r_between_isolated_muon_and_jpsi_in_peak"] = rdfJPsiMuons[sample].Histo1D((sample + "_" + "delta_r_between_isolated_muon_and_jpsi_in_peak", "Motne Carlo" + sample + "Delta R Between Isolated Muon And J/Psi; Delta R; Events", 50, 0, 6), "DeltaRBetweenIsolatedMuonAndJPsiInPeak", "LumiXS")
#         hist[sample]["isolated_muon_pf_rel_iso_03_all"] = rdfJPsiMuons[sample].Histo1D((sample + "_" + "isolated_muon_pf_rel_iso_03_all", "Monte Carlo" + sample + "; Isolated Muon PfRelIso03_All; Isolated Muon PfRelIso03_All; Events", 50, 0, .5), "IsolatedMuon_pfRelIso03_all", "LumiXS")
#         hist[sample]["isolated_muon_pf_rel_iso_03_all_zoomed"] = rdfJPsiMuons[sample].Histo1D((sample + "_" + "isolated_muon_pf_rel_iso_03_all_zoomed", "Monte Carlo" + sample + "; Isolated Muon PfRelIso03_All Zoomed; Isolated Muon PfRelIso03_All Zoomed; Events", 30, 0, .15), "IsolatedMuon_pfRelIso03_all", "LumiXS")
#         hist[sample]["isolated_muon_pf_rel_iso_03_chg"] = rdfJPsiMuons[sample].Histo1D((sample + "_" + "isolated_muon_pf_rel_iso_03_chg", "Monte Carlo" + sample + "; Isolated Muon PfRelIso03_Chg; Isolated Muon PfRelIso03_Chg; Events", 50, 0, .5), "IsolatedMuon_pfRelIso03_chg", "LumiXS")
#         hist[sample]["jpsi_muon_positive_pf_rel_iso_03_all"] = rdfJPsiMuons[sample].Histo1D((sample + "_" + "jpsi_muon_positive_pf_rel_iso_03_all", "Monte Carlo" + sample + "; JPsi Muon Positive PfRelIso03_All; JPsi Muon Positive PfRelIso03_All; Events", 50, 0, .5), "JPsiMuonPositive_pfRelIso03_all", "LumiXS")
#         hist[sample]["jpsi_muon_positive_pf_rel_iso_03_all_zoomed"] = rdfJPsiMuons[sample].Histo1D((sample + "_" + "jpsi_muon_positive_pf_rel_iso_03_all_zoomed", "Monte Carlo" + sample + "; JPsi Muon Positive PfRelIso03_All Zoomed; JPsi Muon Positive PfRelIso03_All Zoomed; Events", 30, 0, .15), "JPsiMuonPositive_pfRelIso03_all", "LumiXS")
#         hist[sample]["jpsi_muon_positive_pf_rel_iso_03_chg"] = rdfJPsiMuons[sample].Histo1D((sample + "_" + "jpsi_muon_positive_pf_rel_iso_03_chg", "Monte Carlo" + sample + "; JPsi Muon Positive PfRelIso03_Chg; JPsi Muon Positive PfRelIso03_Chg; Events", 50, 0, .5), "JPsiMuonPositive_pfRelIso03_chg", "LumiXS")
#         hist[sample]["jpsi_muon_negative_pf_rel_iso_03_all"] = rdfJPsiMuons[sample].Histo1D((sample + "_" + "jpsi_muon_negative_pf_rel_iso_03_all", "Monte Carlo" + sample + "; JPsi Muon Negative PfRelIso03_All; JPsi Muon Negative PfRelIso03_All; Events", 50, 0, .5), "JPsiMuonNegative_pfRelIso03_all", "LumiXS")
#         hist[sample]["jpsi_muon_negative_pf_rel_iso_03_all_zoomed"] = rdfJPsiMuons[sample].Histo1D((sample + "_" + "jpsi_muon_negative_pf_rel_iso_03_all_zoomed", "Monte Carlo" + sample + "; JPsi Muon Negative PfRelIso03_All Zoomed; JPsi Muon Negative PfRelIso03_All Zoomed; Events", 30, 0, .15), "JPsiMuonNegative_pfRelIso03_all", "LumiXS")
#         hist[sample]["jpsi_muon_negative_pf_rel_iso_03_chg"] = rdfJPsiMuons[sample].Histo1D((sample + "_" + "jpsi_muon_negative_pf_rel_iso_03_chg", "Monte Carlo" + sample + "; JPsi Muon Negative PfRelIso03_Chg; JPsi Muon Negative PfRelIso03_Chg; Events", 50, 0, .5), "JPsiMuonNegative_pfRelIso03_chg", "LumiXS")
        
        
#         hist[sample]["isolated_muon_pf_rel_iso_03_all_medium_id"] = rdfJPsiMuons[sample].Histo1D((sample + "_" + "isolated_muon_pf_rel_iso_03_all_medium_id", "Monte Carlo" + sample + "; Isolated Muon PfRelIso03_All Medium ID; Isolated Muon PfRelIso03_All Medium ID; Events", 50, 0, .5), "IsolatedMuon_pfRelIso03_all_mediumID", "LumiXS")
#         hist[sample]["isolated_muon_pf_rel_iso_03_all_medium_id_zoomed"] = rdfJPsiMuons[sample].Histo1D((sample + "_" + "isolated_muon_pf_rel_iso_03_all_medium_id_zoomed", "Monte Carlo" + sample +"; Isolated Muon PfRelIso03_All Medium ID Zoomed; Isolated Muon PfRelIso03_All Medium ID; Events", 30, 0, .15), "IsolatedMuon_pfRelIso03_all_mediumID", "LumiXS")
#         hist[sample]["isolated_muon_pf_rel_iso_03_chg_medium_id"] = rdfJPsiMuons[sample].Histo1D((sample + "_" + "isolated_muon_pf_rel_iso_03_chg_medium_id", "Monte Carlo" + sample + "; Isolated Muon PfRelIso03_Chg Medium ID; Isolated Muon PfRelIso03_Chg Medium ID; Events", 50, 0, .5), "IsolatedMuon_pfRelIso03_chg_mediumID", "LumiXS")
#         hist[sample]["jpsi_muon_positive_pf_rel_iso_03_all_medium_id"] = rdfJPsiMuons[sample].Histo1D((sample + "_" + "jpsi_muon_positive_pf_rel_iso_03_all_medium_id", "Monte Carlo" + sample + "; JPsi Muon Positive PfRelIso03_All Medium ID; JPsi Muon Positive PfRelIso03_All Medium ID; Events", 50, 0, .5), "JPsiMuonPositive_pfRelIso03_all_mediumID", "LumiXS")
#         hist[sample]["jpsi_muon_positive_pf_rel_iso_03_all_medium_id_zoomed"] = rdfJPsiMuons[sample].Histo1D((sample + "_" + "jpsi_muon_positive_pf_rel_iso_03_all_medium_id_zoomed", "Monte Carlo" + sample + "; JPsi Muon Positive PfRelIso03_All Zoomed Medium ID; JPsi Muon Positive PfRelIso03_All Medium ID; Events", 30, 0, .15), "JPsiMuonPositive_pfRelIso03_all_mediumID", "LumiXS")
#         hist[sample]["jpsi_muon_positive_pf_rel_iso_03_chg_medium_id"] = rdfJPsiMuons[sample].Histo1D((sample + "_" + "jpsi_muon_positive_pf_rel_iso_03_chg_medium_id", "Monte Carlo" + sample + "; JPsi Muon Positive PfRelIso03_Chg Medium ID; JPsi Muon Positive PfRelIso03_Chg Medium ID; Events", 50, 0, .5), "JPsiMuonPositive_pfRelIso03_chg_mediumID", "LumiXS")
#         hist[sample]["jpsi_muon_negative_pf_rel_iso_03_all_medium_id"] = rdfJPsiMuons[sample].Histo1D((sample + "_" + "jpsi_muon_negative_pf_rel_iso_03_all_medium_id", "Monte Carlo" + sample + "; JPsi Muon Negative PfRelIso03_All Medium ID; JPsi Muon Negative PfRelIso03_All Medium ID; Events", 50, 0, .5), "JPsiMuonNegative_pfRelIso03_all_mediumID", "LumiXS")
#         hist[sample]["jpsi_muon_negative_pf_rel_iso_03_all_medium_id_zoomed"] = rdfJPsiMuons[sample].Histo1D((sample + "_" + "jpsi_muon_negative_pf_rel_iso_03_all_medium_id_zoomed", "JPsi Muon Negative PfRelIso03_All Zoomed Medium ID; JPsi Muon Negative PfRelIso03_All Medium ID; Events", 30, 0, .15), "JPsiMuonNegative_pfRelIso03_all_mediumID", "LumiXS")
#         hist[sample]["jpsi_muon_negative_pf_rel_iso_03_chg_medium_id"] = rdfJPsiMuons[sample].Histo1D((sample + "_" + "jpsi_muon_negative_pf_rel_iso_03_chg_medium_id", "Monte Carlo" + sample + ";JPsi Muon Negative PfRelIso03_Chg Medium ID; JPsi Muon Negative PfRelIso03_Chg Medium ID; Events", 50, 0, .5), "JPsiMuonNegative_pfRelIso03_chg_mediumID", "LumiXS")
    
    
#         hist[sample]["jpsi_muon_positive_peak_only_pf_rel_iso_03_all"] = rdfJPsiMuons[sample].Histo1D((sample + "_" + "jpsi_muon_positive_peak_only_pf_rel_iso_03_all", "Monte Carlo" + sample + "; JPsi Muon Positive Peak Only PfRelIso03_All; JPsi Muon Positive Peak Only PfRelIso03_All; Events", 50, 0, .5), "JPsiMuonPositive_peak_only_pfRelIso03_all", "LumiXS")
#         hist[sample]["jpsi_muon_positive_peak_only_pf_rel_iso_03_chg"] = rdfJPsiMuons[sample].Histo1D((sample + "_" + "jpsi_muon_positive_peak_only_pf_rel_iso_03_chg", "Monte Carlo" + sample + "; JPsi Muon Positive Peak Only PfRelIso03_Chg; JPsi Muon Positive Peak Only PfRelIso03_Chg; Events", 50, 0, .5), "JPsiMuonPositive_peak_only_pfRelIso03_chg", "LumiXS")
#         hist[sample]["jpsi_muon_negative_peak_only_pf_rel_iso_03_all"] = rdfJPsiMuons[sample].Histo1D((sample + "_" + "jpsi_muon_negative_peak_only_pf_rel_iso_03_all", "Monte Carlo" + sample + "; JPsi Muon Negative Peak Only PfRelIso03_All; JPsi Muon Negative Peak Only PfRelIso03_All; Events", 50, 0, .5), "JPsiMuonNegative_peak_only_pfRelIso03_all", "LumiXS")
#         hist[sample]["jpsi_muon_negative_peak_only_pf_rel_iso_03_chg"] = rdfJPsiMuons[sample].Histo1D((sample + "_" + "jpsi_muon_negative_peak_only_pf_rel_iso_03_chg", "Monte Carlo" + sample + "; JPsi Muon Negative Peak Only PfRelIso03_Chg; JPsi Muon Negative Peak Only PfRelIso03_Chg; Events", 50, 0, .5), "JPsiMuonNegative_peak_only_pfRelIso03_chg", "LumiXS")
        
        
        
        
        
        
        
        
#         hist[sample]["further_muon_medium_tight_pt"] = rdfRemainingIsolatedMuonMuonIdMediumIsoIdTight[sample].Histo1D((sample + "_" + "further_muon_medium_tight_pt", "Monte Carlo" + sample + "; Further Muon Pt - Muon Id Medium, Iso Id Tight; Further Muon Pt - Muon Id Medium, Iso Id Tight; Events", 100, 20, 220), "FurtherMuonMediumTight_pt", "LumiXS")
#         hist[sample]["further_muon_medium_tight_pf_rel_iso_03_all"] = rdfRemainingIsolatedMuonMuonIdMediumIsoIdTight[sample].Histo1D((sample + "_" + "further_muon_medium_tight_pf_rel_iso_03_all", "Monte Carlo" + sample + "; Further Muon Pf Rel Iso 03 All - Muon Id Medium, Iso Id Tight; Further Muon Pf Rel Iso 03 All - Muon Id Medium, Iso Id Tight; Events", 50, 0, .5), "FurtherMuonMediumTight_pfRelIso03_all", "LumiXS")
#         hist[sample]["further_muon_medium_tight_pf_iso_id"] = rdfRemainingIsolatedMuonMuonIdMediumIsoIdTight[sample].Histo1D((sample + "_" + "further_muon_medium_tight_pf_iso_id", "Monte Carlo" + sample + "; Further Muon Pf Iso Id - Muon Id Medium, Iso Id Tight; Further Muon Pf Iso Id - Muon Id Medium, Iso Id Tight; Events", 6, 1, 6), "FurtherMuonMediumTight_pfIsoid", "LumiXS")
        
#         hist[sample]["further_muon_loose_tight_pt"] = rdfRemainingIsolatedMuonMuonIdLooseIsoIdTight[sample].Histo1D((sample + "_" + "further_muon_loose_tight_pt", "Monte Carlo" + sample + "; Further Muon Pt - Muon Id Loose Iso Id Tight; Further Muon Pt - Muon Id Loose, Iso Id Tight; Events", 100, 20, 220), "FurtherMuonLooseTight_pt", "LumiXS")
#         hist[sample]["further_muon_loose_tight_pf_rel_iso_03_all"] = rdfRemainingIsolatedMuonMuonIdLooseIsoIdTight[sample].Histo1D((sample + "_" + "further_muon_loose_tight_pf_rel_iso_03_all", "Monte Carlo" + sample + "; Further Muon Pf Rel Iso 03 All - Muon Id Loose, Iso Id Tight; Further Muon Pf Rel Iso 03 All - Muon Id Loose, Iso Id Tight; Events", 50, 0, .5), "FurtherMuonLooseTight_pfRelIso03_all", "LumiXS")
#         hist[sample]["further_muon_loose_tight_pf_iso_id"] = rdfRemainingIsolatedMuonMuonIdLooseIsoIdTight[sample].Histo1D((sample + "_" + "further_muon_loose_tight_pf_iso_id", "Monte Carlo" + sample + "; Further Muon Pf Iso Id - Muon Id Loose, Iso Id Tight; Further Muon Pf Iso Id - Muon Id Loose, Iso Id Tight; Events", 6, 1, 6), "FurtherMuonLooseTight_pfIsoid", "LumiXS")
        
#         hist[sample]["further_muon_medium_medium_pt"] = rdfRemainingIsolatedMuonMuonIdMediumIsoIdMedium[sample].Histo1D((sample + "_" + "further_muon_medium_medium_pt", "Monte Carlo" + sample + "; Further Muon Pt - Muon Id Medium, Iso Id Medium; Further Muon Pt - Muon Id Medium, Iso Id Medium; Events", 100, 20, 220), "FurtherMuonMediumMedium_pt", "LumiXS")
#         hist[sample]["further_muon_medium_medium_pf_rel_iso_03_all"] = rdfRemainingIsolatedMuonMuonIdMediumIsoIdMedium[sample].Histo1D((sample + "_" + "further_muon_medium_medium_pf_rel_iso_03_all", "Monte Carlo" + sample + "; Further Muon Pf Rel Iso 03 All - Muon Id Medium, Iso Id Medium; Further Muon Pf Rel Iso 03 All - Muon Id Medium, Iso Id Medium; Events", 50, 0, .5), "FurtherMuonMediumMedium_pfRelIso03_all", "LumiXS")
#         hist[sample]["further_muon_medium_medium_pf_iso_id"] = rdfRemainingIsolatedMuonMuonIdMediumIsoIdMedium[sample].Histo1D((sample + "_" + "further_muon_medium_medium_pf_iso_id", "Monte Carlo" + sample + "; Further Muon Pf Iso Id - Muon Id Medium, Iso Id Medium; Further Muon Pf Iso Id - Muon Id Medium, Iso Id Medium; Events", 6, 1, 6), "FurtherMuonMediumMedium_pfIsoid", "LumiXS")
        
#         hist[sample]["further_muon_loose_medium_pt"] = rdfRemainingIsolatedMuonMuonIdLooseIsoIdMedium[sample].Histo1D((sample + "_" + "further_muon_loose_medium_pt", "Monte Carlo" + sample + "; Further Muon Pt - Muon Id Loose, Iso Id Medium; Further Muon Pt - Muon Id Loose, Iso Id Medium; Events", 100, 20, 220), "FurtherMuonLooseMedium_pt", "LumiXS")
#         hist[sample]["further_muon_loose_medium_pf_rel_iso_03_all"] = rdfRemainingIsolatedMuonMuonIdLooseIsoIdMedium[sample].Histo1D((sample + "_" + "further_muon_loose_medium_pf_rel_iso_03_all", "Monte Carlo" + sample + "; Further Muon Pf Rel Iso 03 All - Muon Id Loose, Iso Id Medium; Further Muon Pf Rel Iso 03 All - Muon Id Loose, Iso Id Medium; Events", 50, 0, .5), "FurtherMuonLooseMedium_pfRelIso03_all", "LumiXS")
#         hist[sample]["further_muon_loose_medium_pf_iso_id"] = rdfRemainingIsolatedMuonMuonIdLooseIsoIdMedium[sample].Histo1D((sample + "_" + "further_muon_loose_medium_pf_iso_id", "Monte Carlo" + sample + "; Further Muon Pf Iso Id - Muon Id Loose, Iso Id Medium; Further Muon Pf Iso Id - Muon Id Loose, Iso Id Medium; Events", 6, 1, 6), "FurtherMuonLooseMedium_pfIsoid", "LumiXS")


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



