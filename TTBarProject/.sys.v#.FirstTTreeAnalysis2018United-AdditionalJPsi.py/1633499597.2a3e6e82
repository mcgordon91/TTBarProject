#!/usr/bin/env python
# coding: utf-8

# In[1]:


import ROOT
import glob
import numpy as np
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
rdfTriggersAndMasks = {}
rdfOnePlusMuons = {}
rdfJetVariables = {}
rdfTwoPlusJets = {}
rdfNoMuTrigger = {}
rdfmuOneBTaggedJet = {}
rdfmu = {}
rdfOtherMuons = {}
rdfAllMuons = {}
hist = {}
report = {}
stats = {}

npArrays = {}

for sample, fileList in dictOfListOfFiles.items():
    
    if sample == 'Data':
        vecList = ROOT.std.vector(str)()

        for element in dictOfListOfFiles['Data']:
            vecList.push_back(element)

        rdf['Data'] = ROOT.ROOT.RDataFrame("Events", vecList)    
        #mureport['Data'] = rdf['Data'].Report()
        print(rdf['Data'].Report())
        mrdf['Data'] = ROOT.ROOT.RDataFrame("Runs", vecList)
        
        if useRange:
            rdf['Data'] = rdf['Data'].Range(5000)
            nrange = 5000
            printcode = ' if(rdfentry_ % 5000 == 0) { std::cout << "Processed entry " << rdfentry_ << "/' + str(nrange) + '" << std::endl; } return rdfentry_;'
            print("Data run")
            rdf['Data'] = rdf['Data'].Define("my_rdfentry", printcode)
            
        #s = rdf['Data'].Min('my_rdfentry')
                         
        #c = rdf['Data'].Count()
            

        
    else:
        chain[sample] = ROOT.TChain("Events")
        meta[sample] = ROOT.TChain("Runs")

        for file in fileList:
            #print(file, end=" ")
            
            chain[sample].Add(file)
            meta[sample].Add(file)

        rdf[sample] = ROOT.ROOT.RDataFrame(chain[sample])
        mrdf[sample] = ROOT.ROOT.RDataFrame(meta[sample])
        
        if useRange:
            rdf[sample] = rdf[sample].Range(5000)
            
            printcode = ' if(rdfentry_ % 5000 == 0) { std::cout << "Processed entry " << rdfentry_ << " " << rdfslot_ << std::endl; } return rdfentry_;'
            print("MC run")
            rdf[sample] = rdf[sample].Define("my_rdfentry", printcode)
        
        nevents[sample] = mrdf[sample].Sum("genEventCount")
        sumweight[sample] = mrdf[sample].Sum("genEventSumw")


# In[4]:


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
wgtFormula['MonteCarloST_tW_top'] = "{XS:f} * {lumi:f} * 1000 * genWeight / {sW:f}".format(XS=61526.7, lumi=lumiDict["2018"], sW=float(sumweight['MonteCarloST_tW_top'].GetValue()))
wgtFormula['MonteCarloST_tchannel_top'] = "{XS:f} * {lumi:f} * 1000 * genWeight / {sW:f}".format(XS=61526.7, lumi=lumiDict["2018"], sW=float(sumweight['MonteCarloST_tchannel_top'].GetValue()))


# In[5]:


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

int MatchOppositelyChargedMuons::ReturnFlagPair()
{
    return this->FlagPair;
}
"""

ROOT.gInterpreter.Declare(cpp_code)

#rvf = ROOT.VecOps.RVec(float)([2.71,3.14])
#rvi = ROOT.VecOps.RVec(int)([1,-1])

#x = ROOT.MatchOppositelyChargedMuons(rvf, rvf, rvf, rvf, rvi, 123456789)
#print(x.InvariantMassCalculator())


# In[6]:


# Object selection (masks) and cuts (filters)
# LumiXS is lumi * xs
#.Define("Junk", 'if(rdfentry_ % 200000 == 0){ std::cout << "Processed " << rdfentry_ << " entries on slot " << rdfslot_ << std::endl;} return rdfentry_;')

# Initial cuts for finding muons, isolated or not.
for sample in dictOfListOfFiles:
    rdfTriggersAndMasks[sample] = rdf[sample].Filter("HLT_IsoMu27 == true", "HLTTriggerSingleMuonTrigger").Define("LumiXS",wgtFormula[sample])        .Define("mu_mask", "Muon_pt > 30 && abs(Muon_eta) < 2.4 && Muon_tightId == true && Muon_pfIsoId >= 4")        .Define("mu_veto", "Muon_pt > 10 && abs(Muon_eta) < 2.4 && Muon_looseId == true && Muon_pfIsoId >= 4 && mu_mask == false")        .Define("jpsi_mu_candidate_mask", "Muon_pt > 3 && Muon_mediumId == true && mu_mask == false")
    rdfOnePlusMuons[sample] = rdfTriggersAndMasks[sample].Filter("Sum(mu_mask) == 1", "Exactly One Good Muon")        .Define("SMuon_pt", "Muon_pt[mu_mask]")        .Define("SMuon_eta", "Muon_eta[mu_mask]")        .Define("SMuon_phi", "Muon_phi[mu_mask]")        .Define("SMuon_mass", "Muon_mass[mu_mask]")        .Define("SMuon_charge", "Muon_charge[mu_mask]")        .Define("SMuon_pfRelIso03_all", "Muon_pfRelIso03_all[mu_mask]")        .Define("SMuon_pfRelIso03_chg", "Muon_pfRelIso03_chg[mu_mask]")        .Define("SMuon_pfRelIso04_all", "Muon_pfRelIso04_all[mu_mask]")
    rdfJetVariables[sample] = rdfOnePlusMuons[sample]        .Define("jet_mask", "ROOT::VecOps::RVec<Int_t> jmask = (Jet_pt >= 30 && abs(Jet_eta) <= 2.5 && Jet_jetId >= 2); "                          "for(int i=0; i < SMuon_pt.size(); ++i){"                              "ROOT::VecOps::RVec<Float_t> dr;"                              "for(int j=0; j < jmask.size(); ++j){"                                  "dr.push_back(ROOT::VecOps::DeltaR(Jet_eta.at(j), SMuon_eta.at(i), Jet_phi.at(j), SMuon_phi.at(i)));}"                                  "jmask = jmask && dr >= 0.4;"                                  "dr.clear();}"                          "return jmask;")        .Define("MediumBJetMask", "Jet_btagDeepFlavB > 0.3033 && jet_mask" )        .Define("MTofMETandMu", "FTA::transverseMassMET(SMuon_pt, SMuon_phi, SMuon_mass, MET_pt, MET_phi)")        .Define("Num_Jets", "Jet_pt[jet_mask].size()")        .Define("SJet1_pt", "Jet_pt[jet_mask].size() > 0 ? Jet_pt[jet_mask].at(0) : -500")        .Define("SJet2_pt", "Jet_pt[jet_mask].size() > 1 ? Jet_pt[jet_mask].at(1) : -500")        .Define("SJet1_eta", "Jet_eta[jet_mask].size() > 0 ? Jet_eta[jet_mask].at(0) : 500")        .Define("SJet2_eta", "Jet_eta[jet_mask].size() > 1 ? Jet_eta[jet_mask].at(1) : 500")        .Define("SJet1_phi", "Jet_phi[jet_mask].size() > 0 ? Jet_phi[jet_mask].at(0) : 500")        .Define("SJet2_phi", "Jet_phi[jet_mask].size() > 1 ? Jet_phi[jet_mask].at(1) : 500")        .Define("SJet_btagDeepFlavB", "Jet_btagDeepFlavB[jet_mask]")        .Define("Num_BTaggedJets", "Sum(MediumBJetMask)")        .Define("Ht", "Sum(Jet_pt[jet_mask])")
    rdfTwoPlusJets[sample] = rdfJetVariables[sample].Filter("Num_Jets >= 4", "At Least Four Jets")        .Define("DeepJetB", "Jet_pt[jet_mask].size() > 0 ? Jet_btagDeepFlavB[jet_mask].at(0) : 0")        .Define("Num_Muons", "Muon_pt[mu_mask].size()")
    rdfmuOneBTaggedJet[sample] = rdfTwoPlusJets[sample].Filter("Num_BTaggedJets >= 1", "At Least One B-Tagged Jet")
    rdfmu[sample] = rdfmuOneBTaggedJet[sample].Filter("Sum(mu_veto) == 0", "No Vetoed Muons")
    
    rdfOtherMuons[sample] = rdfmu[sample].Filter("Sum(jpsi_mu_candidate_mask) >= 1", "JPsi Candidate")        .Define("InvariantMasses", "std::cout << rdfentry_ << std::endl; int FlagPair = -1; auto c = MatchOppositelyChargedMuons(Muon_pt, Muon_eta, Muon_phi, Muon_mass, Muon_charge, event, FlagPair); return c.OppositelyChargedMuonInvariantMassCalculator()")        .Define("InvariantMasses_0", "InvariantMasses.at(0, -1.0)")        .Define("InvariantMasses_1", "InvariantMasses.at(1, -1.0)")#if (FlagPair == 1) {print("winning");}; if (FlagPair == 1) {return c.InvariantMassCalculator();}")
    
    rdfAllMuons[sample] = rdfmu[sample].Filter("Sum(jpsi_mu_candidate_mask) >= 1", "JPsi Candidate")        .Define("InvariantMassesAllMuons", "std::cout << rdfentry_ << std::endl; int FlagPair = -1; auto c = MatchOppositelyChargedMuons(Muon_pt, Muon_eta, Muon_phi, Muon_mass, Muon_charge, event, FlagPair, SMuon_pt, SMuon_eta, SMuon_phi, SMuon_mass, SMuon_charge); return c.ThreeMuonInvariantMassCalculator()")
        #.Define("InvariantMasses_0", "InvariantMasses.at(0, -1.0)")\
        #.Define("InvariantMasses_1", "InvariantMasses.at(1, -1.0)")#if (FlagPair == 1) {print("winning");}; if (FlagPair == 1) {return c.InvariantMassCalculator();}")


# In[7]:


for sample in dictOfListOfFiles:
    
    if sample not in hist.keys():
        hist[sample] = {}
        stats[sample] = {}
        report[sample] = rdf[sample].Report()
        
    #stats[sample]["countercode"] = rdf[sample].Min('my_rdfentry')
        
    if sample == 'Data':
        
        hist['Data']["mu_pt"] = rdfmu['Data'].Histo1D(("mu_pt","Muon Transverse Momentum; Pt (GeV);Events",100,20,220),"SMuon_pt")
        hist['Data']["mu_eta"] = rdfmu['Data'].Histo1D(("mu_eta", "Muon Pseudorapidity; Eta; Events",100,-3,3),"SMuon_eta")
        hist['Data']["mu_phi"] = rdfmu['Data'].Histo1D(("mu_phi", "Muon Angle; Phi (Radians); Events",100,-3.5,3.5),"SMuon_phi")
        hist['Data']["jet1_pt"] = rdfmu['Data'].Histo1D(("jet1_pt", "Jet Transverse Momentum for Leading Jet; Pt (GeV); Events", 100, 20, 200), "SJet1_pt")
        hist['Data']["jet2_pt"] = rdfmu['Data'].Histo1D(("jet2_pt", "Jet Transverse Momentum for Subleading Jet; Pt (GeV); Events", 100, 20, 200), "SJet2_pt")
        hist['Data']["jet1_eta"] = rdfmu['Data'].Histo1D(("jet1_eta", "Jet Pseudorapidity for Leading Jet; Eta; Events", 100, -3, 3), "SJet1_eta")
        hist['Data']["jet2_eta"] = rdfmu['Data'].Histo1D(("jet2_eta", "Jet Pseudorapidity for Subleading Jet; Eta; Events", 100, -3, 3), "SJet2_eta")
        hist['Data']["jet1_phi"] = rdfmu['Data'].Histo1D(("jet1_phi", "Jet Angle for Leading Jet; Phi (Radians); Events", 100, -3.5, 3.5), "SJet1_phi")
        hist['Data']["jet2_phi"] = rdfmu['Data'].Histo1D(("jet2_phi", "Jet Angle for Subleading Jet; Phi (Radians); Events", 100, -3.5, 3.5), "SJet2_phi")
        hist['Data']["jet_deep"] = rdfmu['Data'].Histo1D(("jet_deep", "Deep Jet B Discriminator; Discriminant Value; Events", 100, 0, 1), "DeepJetB")
        hist['Data']["number_of_jets"] = rdfmu['Data'].Histo1D(("number_of_jets", "Number of Jets; Number Of Jets; Events", 20, 0, 20), "Num_Jets") 
        hist['Data']["number_of_muons"] = rdfmu['Data'].Histo1D(("number_of_muons", "Number of Muons; Number of Muons; Events", 5, 0, 5), "Num_Muons")
        hist['Data']["transverse_mass"] = rdfmu['Data'].Histo1D(("transverse_mass", "Transverse Mass; Transverse Mass (GeV); Events", 150, 0, 150), "MTofMETandMu")
        hist['Data']["missing_transverse_momentum"] = rdfmu['Data'].Histo1D(("missing_transverse_momentum", "Missing Transverse Momentum; Missing Transverse Momentum(GeV); Events", 150, 0, 300), "MET_pt")
        hist['Data']["ht"] = rdfmu['Data'].Histo1D(("ht", "Ht; Ht; Events", 300, 0, 1500), "Ht")
        hist['Data']["muon_pfRelIso03_all"] = rdfmu['Data'].Histo1D(("muon_pfRelIso03_all", "Muon Pf Rel Iso 03 (All); Muon Pf Rel Iso 03 (All); Events", 60, 0, .3), "SMuon_pfRelIso03_all")
        hist['Data']["muon_pfRelIso03_chg"] = rdfmu['Data'].Histo1D(("muon_pfRelIso03_chg", "Muon Pf Rel Iso 03 (Chg); Muon Pf Rel Iso 03 (Chg); Events", 60, 0, .3), "SMuon_pfRelIso03_chg")
        hist['Data']["muon_pfRelIso04_all"] = rdfmu['Data'].Histo1D(("muon_pfRelIso04_all", "Muon Pf Rel Iso 04 (All); Muon Pf Rel Iso 04 (All); Events", 60, 0, .3), "SMuon_pfRelIso04_all")
        hist['Data']["invariant_masses"] = rdfOtherMuons['Data'].Histo1D(("invariant_masses", "Invariant Masses for Two Muons (Oppositely Charged); Invariant Masses; Events", 50, .5, 12), "InvariantMasses")        
        hist['Data']["invariant_masses_zoomed"] = rdfOtherMuons['Data'].Histo1D(("invariant_masses_zoomed", "Invariant Masses for Two Muons (Oppositely Charged); Invariant Masses; Events", 50, 2.8, 3.4), "InvariantMasses")
        hist['Data']["invariant_masses_all_muons"] = rdfAllMuons['Data'].Histo1D(("invariant_masses_all_muons", "Invariant Masses for Three Muons (Isolated and Paired, Oppositely Charged); Invariant Masses; Events", 50, .5, 200), "InvariantMassesAllMuons")

    else:
        
        hist[sample]["mu_pt"] = rdfmu[sample].Histo1D((sample + "_" + "mu_pt","Monte Carlo " + sample + ";Pt (GeV);Events",100,20,220),"SMuon_pt","LumiXS")
        hist[sample]["mu_eta"] = rdfmu[sample].Histo1D((sample + "_" + "mu_eta", "Monte Carlo " + sample + "; Eta; Events",100,-3,3),"SMuon_eta","LumiXS")
        hist[sample]["mu_phi"] = rdfmu[sample].Histo1D((sample + "_" + "mu_phi", "Monte Carlo " + sample + "; Phi (Radians); Events",100,-3.5,3.5),"SMuon_phi","LumiXS")
        hist[sample]["jet1_pt"] = rdfmu[sample].Histo1D((sample + "_" + "jet1_pt", "Monte Carlo " + sample + "; Pt (GeV); Events", 100, 20, 200), "SJet1_pt", "LumiXS")
        hist[sample]["jet2_pt"] = rdfmu[sample].Histo1D((sample + "_" + "jet2_pt", "Monte Carlo " + sample + "; Pt (GeV); Events", 100, 20, 200), "SJet2_pt", "LumiXS")
        hist[sample]["jet1_eta"] = rdfmu[sample].Histo1D((sample + "_" + "jet1_eta", "Monte Carlo " + sample + "; Eta; Events", 100, -3, 3), "SJet1_eta", "LumiXS")
        hist[sample]["jet2_eta"] = rdfmu[sample].Histo1D((sample + "_" + "jet2_eta", "Monte Carlo " + sample +"; Eta; Events", 100, -3, 3), "SJet2_eta", "LumiXS")
        hist[sample]["jet1_phi"] = rdfmu[sample].Histo1D((sample + "_" + "jet1_phi", "Monte Carlo " + sample + "; Phi (Radians); Events", 100, -3.5, 3.5), "SJet1_phi", "LumiXS")
        hist[sample]["jet2_phi"] = rdfmu[sample].Histo1D((sample + "_" + "jet2_phi", "Monte Carlo " + sample +"; Phi (Radians); Events", 100, -3.5, 3.5), "SJet2_phi", "LumiXS")
        hist[sample]["jet_deep"] = rdfmu[sample].Histo1D((sample + "_" + "jet_deep", "Monte Carlo " + sample +"; Discriminant Value; Events", 100, 0, 1), "DeepJetB", "LumiXS")
        hist[sample]["number_of_jets"] = rdfmu[sample].Histo1D((sample + "_" + "number_of_jets", "Monte Carlo " + sample +"; Number Of Jets; Events", 20, 0, 20), "Num_Jets", "LumiXS") 
        hist[sample]["number_of_muons"] = rdfmu[sample].Histo1D((sample + "_" + "number_of_muons", "Monte Carlo " + sample +"; Number of Muons; Events", 5, 0, 5), "Num_Muons", "LumiXS")
        hist[sample]["transverse_mass"] = rdfmu[sample].Histo1D((sample + "_" + "transverse_mass", "Monte Carlo " + sample +"; Transverse Mass (GeV); Events", 150, 0, 150), "MTofMETandMu", "LumiXS")
        hist[sample]["missing_transverse_momentum"] = rdfmu[sample].Histo1D((sample + "_" + "missing_transverse_momentum", "Monte Carlo" + sample + "; Missing Transverse Momentum(GeV); Events", 150, 0, 300), "MET_pt", "LumiXS")
        hist[sample]["ht"] = rdfmu[sample].Histo1D((sample + "_" + "ht", "Monte Carlo " + sample + "; Ht; Events", 300, 0, 1500), "Ht", "LumiXS")
        hist[sample]["muon_pfRelIso03_all"] = rdfmu[sample].Histo1D((sample + "_" + "muon_pfRelIso03_all", "Monte Carlo " + sample + "; Muon Pf Rel Iso 03 (All); Events", 60, 0, .3), "SMuon_pfRelIso03_all", "LumiXS")
        hist[sample]["muon_pfRelIso03_chg"] = rdfmu[sample].Histo1D((sample + "_" + "muon_pfRelIso03_chg", "Monte Carlo " + sample + "; Muon Pf Rel Iso 03 (Chg); Events", 60, 0, .3), "SMuon_pfRelIso03_chg", "LumiXS")
        hist[sample]["muon_pfRelIso04_all"] = rdfmu[sample].Histo1D((sample + "_" + "muon_pfRelIso04_all", "Monte Carlo " + sample + "; Muon Pf Rel Iso 04 (All); Events", 60, 0, .3), "SMuon_pfRelIso04_all", "LumiXS")
        hist[sample]["invariant_masses"] = rdfOtherMuons[sample].Histo1D((sample + "_" + "invariant_masses", "Monte Carlo " + sample + "; Invariant Masses for Two Muons (Oppositely Charged); Invariant Masses; Events", 50, .5, 12), "InvariantMasses", "LumiXS")
        hist[sample]["invariant_masses_zoomed"] = rdfOtherMuons[sample].Histo1D((sample + "_" + "invariant_masses_zoomed", "Monte Carlo " + sample + "; Invariant Masses for Two Muons (Oppositely Charged); Events", 50, 2.8, 3.4), "InvariantMasses", "LumiXS")
        hist[sample]["invariant_masses_all_muons"] = rdfAllMuons[sample].Histo1D((sample + "_" + "invariant_masses_all_muons", "Monte Carlo" + sample + "; Invariant Masses for Three Muons (Isolated and Paired, Oppositely Charged); Invariant Masses; Events", 50, .5, 200), "InvariantMassesAllMuons", "LumiXS")


# In[8]:


for sample in dictOfListOfFiles:

    ROOT.RDF.SaveGraph(rdf[sample], str(sample) + ".dot")


# In[ ]:


for sample in dictOfListOfFiles:
    #npArrays[sample] = rdfOtherMuons[sample].AsNumpy(["InvariantMasses_0", "InvariantMasses_1"])
            
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


# In[12]:


npArrays['Data']['InvariantMasses_0']


# In[ ]:




