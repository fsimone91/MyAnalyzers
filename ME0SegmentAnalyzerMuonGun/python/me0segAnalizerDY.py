import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.load('Configuration.StandardSequences.Services_cff')
process.load('Configuration.Geometry.GeometryExtended2023HGCalMuonReco_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_cff')
process.load('Configuration.Geometry.GeometryExtended2023HGCalMuon_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.load('RecoMuon.MuonIdentification.me0MuonReco_cff')



from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'PH2_1K_FB_V6::All', '')
#process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:upgrade2019', '')


process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source = cms.Source("PoolSource",
                            # replace 'myfile.root',',',',',' with the source file you want to use
                            fileNames = cms.untracked.vstring(
                                                              
    'file:/lustre/cms/store/user/rosma/DYToMuMu_M_20_TuneZ2star_14TeV_pythia6_GEN-SIM_2023HGCalME0GeomV2_SP5_v4/crab_DYToMuMu_M20_14TeV_2023HGCalFastTime_PU0_1ns_500um_1cm_NeutrBkg_7ns_SP5_costPhiSmear_RECO_v2/160321_161927/0000/step3_v4_34.root',
    'file:/lustre/cms/store/user/rosma/DYToMuMu_M_20_TuneZ2star_14TeV_pythia6_GEN-SIM_2023HGCalME0GeomV2_SP5_v4/crab_DYToMuMu_M20_14TeV_2023HGCalFastTime_PU0_1ns_500um_1cm_NeutrBkg_7ns_SP5_costPhiSmear_RECO_v2/160321_161927/0000/step3_v4_340.root',
#    'file:/lustre/cms/store/user/rosma/DYToMuMu_M_20_TuneZ2star_14TeV_pythia6_GEN-SIM_2023HGCalME0GeomV2_SP5_v4/crab_DYToMuMu_M20_14TeV_2023HGCalFastTime_PU0_1ns_500um_1cm_NeutrBkg_7ns_SP5_costPhiSmear_RECO_v2/160321_161927/0000/step3_v4_341.root',
    'file:/lustre/cms/store/user/rosma/DYToMuMu_M_20_TuneZ2star_14TeV_pythia6_GEN-SIM_2023HGCalME0GeomV2_SP5_v4/crab_DYToMuMu_M20_14TeV_2023HGCalFastTime_PU0_1ns_500um_1cm_NeutrBkg_7ns_SP5_costPhiSmear_RECO_v2/160321_161927/0000/step3_v4_342.root',
    'file:/lustre/cms/store/user/rosma/DYToMuMu_M_20_TuneZ2star_14TeV_pythia6_GEN-SIM_2023HGCalME0GeomV2_SP5_v4/crab_DYToMuMu_M20_14TeV_2023HGCalFastTime_PU0_1ns_500um_1cm_NeutrBkg_7ns_SP5_costPhiSmear_RECO_v2/160321_161927/0000/step3_v4_343.root',
     'file:/lustre/cms/store/user/rosma/DYToMuMu_M_20_TuneZ2star_14TeV_pythia6_GEN-SIM_2023HGCalME0GeomV2_SP5_v4/crab_DYToMuMu_M20_14TeV_2023HGCalFastTime_PU0_1ns_500um_1cm_NeutrBkg_7ns_SP5_costPhiSmear_RECO_v2/160321_161927/0000/step3_v4_344.root',
    'file:/lustre/cms/store/user/rosma/DYToMuMu_M_20_TuneZ2star_14TeV_pythia6_GEN-SIM_2023HGCalME0GeomV2_SP5_v4/crab_DYToMuMu_M20_14TeV_2023HGCalFastTime_PU0_1ns_500um_1cm_NeutrBkg_7ns_SP5_costPhiSmear_RECO_v2/160321_161927/0000/step3_v4_345.root',
     'file:/lustre/cms/store/user/rosma/DYToMuMu_M_20_TuneZ2star_14TeV_pythia6_GEN-SIM_2023HGCalME0GeomV2_SP5_v4/crab_DYToMuMu_M20_14TeV_2023HGCalFastTime_PU0_1ns_500um_1cm_NeutrBkg_7ns_SP5_costPhiSmear_RECO_v2/160321_161927/0000/step3_v4_346.root',
    'file:/lustre/cms/store/user/rosma/DYToMuMu_M_20_TuneZ2star_14TeV_pythia6_GEN-SIM_2023HGCalME0GeomV2_SP5_v4/crab_DYToMuMu_M20_14TeV_2023HGCalFastTime_PU0_1ns_500um_1cm_NeutrBkg_7ns_SP5_costPhiSmear_RECO_v2/160321_161927/0000/step3_v4_347.root',
    'file:/lustre/cms/store/user/rosma/DYToMuMu_M_20_TuneZ2star_14TeV_pythia6_GEN-SIM_2023HGCalME0GeomV2_SP5_v4/crab_DYToMuMu_M20_14TeV_2023HGCalFastTime_PU0_1ns_500um_1cm_NeutrBkg_7ns_SP5_costPhiSmear_RECO_v2/160321_161927/0000/step3_v4_348.root',
     'file:/lustre/cms/store/user/rosma/DYToMuMu_M_20_TuneZ2star_14TeV_pythia6_GEN-SIM_2023HGCalME0GeomV2_SP5_v4/crab_DYToMuMu_M20_14TeV_2023HGCalFastTime_PU0_1ns_500um_1cm_NeutrBkg_7ns_SP5_costPhiSmear_RECO_v2/160321_161927/0000/step3_v4_349.root',
      'file:/lustre/cms/store/user/rosma/DYToMuMu_M_20_TuneZ2star_14TeV_pythia6_GEN-SIM_2023HGCalME0GeomV2_SP5_v4/crab_DYToMuMu_M20_14TeV_2023HGCalFastTime_PU0_1ns_500um_1cm_NeutrBkg_7ns_SP5_costPhiSmear_RECO_v2/160321_161927/0000/step3_v4_350.root',
    'file:/lustre/cms/store/user/rosma/DYToMuMu_M_20_TuneZ2star_14TeV_pythia6_GEN-SIM_2023HGCalME0GeomV2_SP5_v4/crab_DYToMuMu_M20_14TeV_2023HGCalFastTime_PU0_1ns_500um_1cm_NeutrBkg_7ns_SP5_costPhiSmear_RECO_v2/160321_161927/0000/step3_v4_351.root',
                                                              
                                                              
        ))
                                                              
process.me0SegAna = cms.EDAnalyzer('ME0SegmentAnalyzerMuonGun',
                                  wp = cms.string("tight"),
                                  minEta = cms.double(2.0),
                                  maxEta = cms.double(2.8),
                                  #    etaMin = cms.double(2.0),
                                  #    etaMax = cms.double(2.8),
                                  #    dr = cms.double(0.25),
                                  #    ptMin = cms.double(5.0),
                                  timeMin = cms.double(5.5),
                                  timeMax = cms.double(30.5)
)


process.TFileService = cms.Service("TFileService",
fileName = cms.string("histoME0SegDY.root")
)

process.p = cms.Path(process.me0SegAna)


