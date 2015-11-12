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
    # replace 'myfile.root' with the source file you want to use
                            fileNames = cms.untracked.vstring( 'file:/lustre/cms/store/user/piet/ME0Segment_Time/DYToMuMu_M-20_TuneZ2star_14TeV-pythia6-tauola/DYToMuMu_M-20_TuneZ2star_14TeV-pythia6-tauola_HGCALGS_PU140_100ps_1p5ns_v2_RECO/151021_081029/0000/out_reco_1.root',
                                                              'file:/lustre/cms/store/user/piet/ME0Segment_Time/DYToMuMu_M-20_TuneZ2star_14TeV-pythia6-tauola/DYToMuMu_M-20_TuneZ2star_14TeV-pythia6-tauola_HGCALGS_PU140_100ps_1p5ns_v2_RECO/151021_081029/0000/out_reco_10.root',
                                                              'file:/lustre/cms/store/user/piet/ME0Segment_Time/DYToMuMu_M-20_TuneZ2star_14TeV-pythia6-tauola/DYToMuMu_M-20_TuneZ2star_14TeV-pythia6-tauola_HGCALGS_PU140_100ps_1p5ns_v2_RECO/151021_081029/0000/out_reco_101.root',
                                                              'file:/lustre/cms/store/user/piet/ME0Segment_Time/DYToMuMu_M-20_TuneZ2star_14TeV-pythia6-tauola/DYToMuMu_M-20_TuneZ2star_14TeV-pythia6-tauola_HGCALGS_PU140_100ps_1p5ns_v2_RECO/151021_081029/0000/out_reco_102.root',
                                                              'file:/lustre/cms/store/user/piet/ME0Segment_Time/DYToMuMu_M-20_TuneZ2star_14TeV-pythia6-tauola/DYToMuMu_M-20_TuneZ2star_14TeV-pythia6-tauola_HGCALGS_PU140_100ps_1p5ns_v2_RECO/151021_081029/0000/out_reco_103.root',
                                                              'file:/lustre/cms/store/user/piet/ME0Segment_Time/DYToMuMu_M-20_TuneZ2star_14TeV-pythia6-tauola/DYToMuMu_M-20_TuneZ2star_14TeV-pythia6-tauola_HGCALGS_PU140_100ps_1p5ns_v2_RECO/151021_081029/0000/out_reco_104.root',
                                                              'file:/lustre/cms/store/user/piet/ME0Segment_Time/DYToMuMu_M-20_TuneZ2star_14TeV-pythia6-tauola/DYToMuMu_M-20_TuneZ2star_14TeV-pythia6-tauola_HGCALGS_PU140_100ps_1p5ns_v2_RECO/151021_081029/0000/out_reco_105.root',
                                                              'file:/lustre/cms/store/user/piet/ME0Segment_Time/DYToMuMu_M-20_TuneZ2star_14TeV-pythia6-tauola/DYToMuMu_M-20_TuneZ2star_14TeV-pythia6-tauola_HGCALGS_PU140_100ps_1p5ns_v2_RECO/151021_081029/0000/out_reco_107.root',
                                                              'file:/lustre/cms/store/user/piet/ME0Segment_Time/DYToMuMu_M-20_TuneZ2star_14TeV-pythia6-tauola/DYToMuMu_M-20_TuneZ2star_14TeV-pythia6-tauola_HGCALGS_PU140_100ps_1p5ns_v2_RECO/151021_081029/0000/out_reco_108.root',
                                                              'file:/lustre/cms/store/user/piet/ME0Segment_Time/DYToMuMu_M-20_TuneZ2star_14TeV-pythia6-tauola/DYToMuMu_M-20_TuneZ2star_14TeV-pythia6-tauola_HGCALGS_PU140_100ps_1p5ns_v2_RECO/151021_081029/0000/out_reco_11.root',
                                                              'file:/lustre/cms/store/user/piet/ME0Segment_Time/DYToMuMu_M-20_TuneZ2star_14TeV-pythia6-tauola/DYToMuMu_M-20_TuneZ2star_14TeV-pythia6-tauola_HGCALGS_PU140_100ps_1p5ns_v2_RECO/151021_081029/0000/out_reco_110.root',
                                                              'file:/lustre/cms/store/user/piet/ME0Segment_Time/DYToMuMu_M-20_TuneZ2star_14TeV-pythia6-tauola/DYToMuMu_M-20_TuneZ2star_14TeV-pythia6-tauola_HGCALGS_PU140_100ps_1p5ns_v2_RECO/151021_081029/0000/out_reco_111.root',
                                                              'file:/lustre/cms/store/user/piet/ME0Segment_Time/DYToMuMu_M-20_TuneZ2star_14TeV-pythia6-tauola/DYToMuMu_M-20_TuneZ2star_14TeV-pythia6-tauola_HGCALGS_PU140_100ps_1p5ns_v2_RECO/151021_081029/0000/out_reco_112.root',
                                                              'file:/lustre/cms/store/user/piet/ME0Segment_Time/DYToMuMu_M-20_TuneZ2star_14TeV-pythia6-tauola/DYToMuMu_M-20_TuneZ2star_14TeV-pythia6-tauola_HGCALGS_PU140_100ps_1p5ns_v2_RECO/151021_081029/0000/out_reco_113.root',
                                                              'file:/lustre/cms/store/user/piet/ME0Segment_Time/DYToMuMu_M-20_TuneZ2star_14TeV-pythia6-tauola/DYToMuMu_M-20_TuneZ2star_14TeV-pythia6-tauola_HGCALGS_PU140_100ps_1p5ns_v2_RECO/151021_081029/0000/out_reco_114.root',
                                                              'file:/lustre/cms/store/user/piet/ME0Segment_Time/DYToMuMu_M-20_TuneZ2star_14TeV-pythia6-tauola/DYToMuMu_M-20_TuneZ2star_14TeV-pythia6-tauola_HGCALGS_PU140_100ps_1p5ns_v2_RECO/151021_081029/0000/out_reco_115.root',
                                                              'file:/lustre/cms/store/user/piet/ME0Segment_Time/DYToMuMu_M-20_TuneZ2star_14TeV-pythia6-tauola/DYToMuMu_M-20_TuneZ2star_14TeV-pythia6-tauola_HGCALGS_PU140_100ps_1p5ns_v2_RECO/151021_081029/0000/out_reco_116.root',
                                                              'file:/lustre/cms/store/user/piet/ME0Segment_Time/DYToMuMu_M-20_TuneZ2star_14TeV-pythia6-tauola/DYToMuMu_M-20_TuneZ2star_14TeV-pythia6-tauola_HGCALGS_PU140_100ps_1p5ns_v2_RECO/151021_081029/0000/out_reco_117.root',
                                                              'file:/lustre/cms/store/user/piet/ME0Segment_Time/DYToMuMu_M-20_TuneZ2star_14TeV-pythia6-tauola/DYToMuMu_M-20_TuneZ2star_14TeV-pythia6-tauola_HGCALGS_PU140_100ps_1p5ns_v2_RECO/151021_081029/0000/out_reco_118.root',
                                                              'file:/lustre/cms/store/user/piet/ME0Segment_Time/DYToMuMu_M-20_TuneZ2star_14TeV-pythia6-tauola/DYToMuMu_M-20_TuneZ2star_14TeV-pythia6-tauola_HGCALGS_PU140_100ps_1p5ns_v2_RECO/151021_081029/0000/out_reco_119.root',
                                                              'file:/lustre/cms/store/user/piet/ME0Segment_Time/DYToMuMu_M-20_TuneZ2star_14TeV-pythia6-tauola/DYToMuMu_M-20_TuneZ2star_14TeV-pythia6-tauola_HGCALGS_PU140_100ps_1p5ns_v2_RECO/151021_081029/0000/out_reco_120.root',
                                                              'file:/lustre/cms/store/user/piet/ME0Segment_Time/DYToMuMu_M-20_TuneZ2star_14TeV-pythia6-tauola/DYToMuMu_M-20_TuneZ2star_14TeV-pythia6-tauola_HGCALGS_PU140_100ps_1p5ns_v2_RECO/151021_081029/0000/out_reco_121.root',
                                                              'file:/lustre/cms/store/user/piet/ME0Segment_Time/DYToMuMu_M-20_TuneZ2star_14TeV-pythia6-tauola/DYToMuMu_M-20_TuneZ2star_14TeV-pythia6-tauola_HGCALGS_PU140_100ps_1p5ns_v2_RECO/151021_081029/0000/out_reco_122.root',
                                                              'file:/lustre/cms/store/user/piet/ME0Segment_Time/DYToMuMu_M-20_TuneZ2star_14TeV-pythia6-tauola/DYToMuMu_M-20_TuneZ2star_14TeV-pythia6-tauola_HGCALGS_PU140_100ps_1p5ns_v2_RECO/151021_081029/0000/out_reco_123.root',
                                                              'file:/lustre/cms/store/user/piet/ME0Segment_Time/DYToMuMu_M-20_TuneZ2star_14TeV-pythia6-tauola/DYToMuMu_M-20_TuneZ2star_14TeV-pythia6-tauola_HGCALGS_PU140_100ps_1p5ns_v2_RECO/151021_081029/0000/out_reco_124.root',
                                                              'file:/lustre/cms/store/user/piet/ME0Segment_Time/DYToMuMu_M-20_TuneZ2star_14TeV-pythia6-tauola/DYToMuMu_M-20_TuneZ2star_14TeV-pythia6-tauola_HGCALGS_PU140_100ps_1p5ns_v2_RECO/151021_081029/0000/out_reco_125.root',
                                                              'file:/lustre/cms/store/user/piet/ME0Segment_Time/DYToMuMu_M-20_TuneZ2star_14TeV-pythia6-tauola/DYToMuMu_M-20_TuneZ2star_14TeV-pythia6-tauola_HGCALGS_PU140_100ps_1p5ns_v2_RECO/151021_081029/0000/out_reco_126.root',
                                                              'file:/lustre/cms/store/user/piet/ME0Segment_Time/DYToMuMu_M-20_TuneZ2star_14TeV-pythia6-tauola/DYToMuMu_M-20_TuneZ2star_14TeV-pythia6-tauola_HGCALGS_PU140_100ps_1p5ns_v2_RECO/151021_081029/0000/out_reco_127.root',
                                                              'file:/lustre/cms/store/user/piet/ME0Segment_Time/DYToMuMu_M-20_TuneZ2star_14TeV-pythia6-tauola/DYToMuMu_M-20_TuneZ2star_14TeV-pythia6-tauola_HGCALGS_PU140_100ps_1p5ns_v2_RECO/151021_081029/0000/out_reco_128.root',
                                                              'file:/lustre/cms/store/user/piet/ME0Segment_Time/DYToMuMu_M-20_TuneZ2star_14TeV-pythia6-tauola/DYToMuMu_M-20_TuneZ2star_14TeV-pythia6-tauola_HGCALGS_PU140_100ps_1p5ns_v2_RECO/151021_081029/0000/out_reco_129.root',
                                                              'file:/lustre/cms/store/user/piet/ME0Segment_Time/DYToMuMu_M-20_TuneZ2star_14TeV-pythia6-tauola/DYToMuMu_M-20_TuneZ2star_14TeV-pythia6-tauola_HGCALGS_PU140_100ps_1p5ns_v2_RECO/151021_081029/0000/out_reco_130.root',
                                                              'file:/lustre/cms/store/user/piet/ME0Segment_Time/DYToMuMu_M-20_TuneZ2star_14TeV-pythia6-tauola/DYToMuMu_M-20_TuneZ2star_14TeV-pythia6-tauola_HGCALGS_PU140_100ps_1p5ns_v2_RECO/151021_081029/0000/out_reco_131.root',
                                                              'file:/lustre/cms/store/user/piet/ME0Segment_Time/DYToMuMu_M-20_TuneZ2star_14TeV-pythia6-tauola/DYToMuMu_M-20_TuneZ2star_14TeV-pythia6-tauola_HGCALGS_PU140_100ps_1p5ns_v2_RECO/151021_081029/0000/out_reco_132.root',
                                                              'file:/lustre/cms/store/user/piet/ME0Segment_Time/DYToMuMu_M-20_TuneZ2star_14TeV-pythia6-tauola/DYToMuMu_M-20_TuneZ2star_14TeV-pythia6-tauola_HGCALGS_PU140_100ps_1p5ns_v2_RECO/151021_081029/0000/out_reco_133.root',
                                                              
                                                              
                                                              
                )
                            )
process.looseSeq = cms.EDAnalyzer('ME0PhiAnalyzer',
    wp = cms.string("loose"),
#    etaMin = cms.double(2.0),
#    etaMax = cms.double(2.8),
#    dr = cms.double(0.25),
#    ptMin = cms.double(5.0),
    timeMin = cms.double(5.5),
    timeMax = cms.double(30.5)
)

process.tightSeq = cms.EDAnalyzer('ME0PhiAnalyzer',
    wp = cms.string("tight"),
#    etaMin = cms.double(2.0),
#    etaMax = cms.double(2.8),
#    dr = cms.double(0.25),
#    ptMin = cms.double(5.0),
    timeMin = cms.double(5.5),
    timeMax = cms.double(30.5)
)


process.TFileService = cms.Service("TFileService",
fileName = cms.string("histoME0.root")
)

#process.p = cms.Path(process.looseSeq )
process.p = cms.Path(process.looseSeq + process.tightSeq)


