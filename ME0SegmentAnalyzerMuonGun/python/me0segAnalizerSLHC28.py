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

#'file:/lustre/cms/store/user/rosma/SingleMuPlusPt30_14TeV_GEN-SIM_2023HGCalME0Geomv2_AllEtaPhi/crab_SingleMuPt30_14TeV_2023HGCalME0Geomv2_RECO_500um_500um_5ns_AllEtaPhi_PU200/160509_110657/0000/step3_1um_10.root',
#'file:/lustre/cms/store/user/rosma/SingleMuPlusPt30_14TeV_GEN-SIM_2023HGCalME0Geomv2_AllEtaPhi/crab_SingleMuPt30_14TeV_2023HGCalME0Geomv2_RECO_500um_500um_5ns_AllEtaPhi_PU200/160509_110657/0000/step3_1um_100.root',
#'file:/lustre/cms/store/user/rosma/SingleMuPlusPt30_14TeV_GEN-SIM_2023HGCalME0Geomv2_AllEtaPhi/crab_SingleMuPt30_14TeV_2023HGCalME0Geomv2_RECO_500um_500um_5ns_AllEtaPhi_PU200/160509_110657/0000/step3_1um_101.root',
#'file:/lustre/cms/store/user/rosma/SingleMuPlusPt30_14TeV_GEN-SIM_2023HGCalME0Geomv2_AllEtaPhi/crab_SingleMuPt30_14TeV_2023HGCalME0Geomv2_RECO_500um_500um_5ns_AllEtaPhi_PU200/160509_110657/0000/step3_1um_102.root',
#'file:/lustre/cms/store/user/rosma/SingleMuPlusPt30_14TeV_GEN-SIM_2023HGCalME0Geomv2_AllEtaPhi/crab_SingleMuPt30_14TeV_2023HGCalME0Geomv2_RECO_500um_500um_5ns_AllEtaPhi_PU200/160509_110657/0000/step3_1um_103.root',
#'file:/lustre/cms/store/user/rosma/SingleMuPlusPt30_14TeV_GEN-SIM_2023HGCalME0Geomv2_AllEtaPhi/crab_SingleMuPt30_14TeV_2023HGCalME0Geomv2_RECO_500um_500um_5ns_AllEtaPhi_PU200/160509_110657/0000/step3_1um_104.root',
#'file:/lustre/cms/store/user/rosma/SingleMuPlusPt30_14TeV_GEN-SIM_2023HGCalME0Geomv2_AllEtaPhi/crab_SingleMuPt30_14TeV_2023HGCalME0Geomv2_RECO_500um_500um_5ns_AllEtaPhi_PU200/160509_110657/0000/step3_1um_105.root',
#'file:/lustre/cms/store/user/rosma/SingleMuPlusPt30_14TeV_GEN-SIM_2023HGCalME0Geomv2_AllEtaPhi/crab_SingleMuPt30_14TeV_2023HGCalME0Geomv2_RECO_500um_500um_5ns_AllEtaPhi_PU200/160509_110657/0000/step3_1um_106.root',
#'file:/lustre/cms/store/user/rosma/SingleMuPlusPt30_14TeV_GEN-SIM_2023HGCalME0Geomv2_AllEtaPhi/crab_SingleMuPt30_14TeV_2023HGCalME0Geomv2_RECO_500um_500um_5ns_AllEtaPhi_PU200/160509_110657/0000/step3_1um_107.root',
#'file:/lustre/cms/store/user/rosma/SingleMuPlusPt30_14TeV_GEN-SIM_2023HGCalME0Geomv2_AllEtaPhi/crab_SingleMuPt30_14TeV_2023HGCalME0Geomv2_RECO_500um_500um_5ns_AllEtaPhi_PU200/160509_110657/0000/step3_1um_108.root',
                                                              'file:/lustre/cms/store/user/rosma/SingleMuPlusPt30_14TeV_GEN-SIM_2023HGCalME0Geomv2_AllEtaPhi/crab_SingleMuPt30_14TeV_2023HGCalME0Geomv2_RECO_500um_500um_100ps_AllEtaPhi_PU200/160506_083146/0000/step3_10.root',
                                                              'file:/lustre/cms/store/user/rosma/SingleMuPlusPt30_14TeV_GEN-SIM_2023HGCalME0Geomv2_AllEtaPhi/crab_SingleMuPt30_14TeV_2023HGCalME0Geomv2_RECO_500um_500um_100ps_AllEtaPhi_PU200/160506_083146/0000/step3_100.root',
                                                              'file:/lustre/cms/store/user/rosma/SingleMuPlusPt30_14TeV_GEN-SIM_2023HGCalME0Geomv2_AllEtaPhi/crab_SingleMuPt30_14TeV_2023HGCalME0Geomv2_RECO_500um_500um_100ps_AllEtaPhi_PU200/160506_083146/0000/step3_101.root',
                                                              'file:/lustre/cms/store/user/rosma/SingleMuPlusPt30_14TeV_GEN-SIM_2023HGCalME0Geomv2_AllEtaPhi/crab_SingleMuPt30_14TeV_2023HGCalME0Geomv2_RECO_500um_500um_100ps_AllEtaPhi_PU200/160506_083146/0000/step3_102.root',
                                                              'file:/lustre/cms/store/user/rosma/SingleMuPlusPt30_14TeV_GEN-SIM_2023HGCalME0Geomv2_AllEtaPhi/crab_SingleMuPt30_14TeV_2023HGCalME0Geomv2_RECO_500um_500um_100ps_AllEtaPhi_PU200/160506_083146/0000/step3_103.root',
                                                              'file:/lustre/cms/store/user/rosma/SingleMuPlusPt30_14TeV_GEN-SIM_2023HGCalME0Geomv2_AllEtaPhi/crab_SingleMuPt30_14TeV_2023HGCalME0Geomv2_RECO_500um_500um_100ps_AllEtaPhi_PU200/160506_083146/0000/step3_104.root',
                                                              'file:/lustre/cms/store/user/rosma/SingleMuPlusPt30_14TeV_GEN-SIM_2023HGCalME0Geomv2_AllEtaPhi/crab_SingleMuPt30_14TeV_2023HGCalME0Geomv2_RECO_500um_500um_100ps_AllEtaPhi_PU200/160506_083146/0000/step3_105.root',
                                                              'file:/lustre/cms/store/user/rosma/SingleMuPlusPt30_14TeV_GEN-SIM_2023HGCalME0Geomv2_AllEtaPhi/crab_SingleMuPt30_14TeV_2023HGCalME0Geomv2_RECO_500um_500um_100ps_AllEtaPhi_PU200/160506_083146/0000/step3_106.root',
                                                              'file:/lustre/cms/store/user/rosma/SingleMuPlusPt30_14TeV_GEN-SIM_2023HGCalME0Geomv2_AllEtaPhi/crab_SingleMuPt30_14TeV_2023HGCalME0Geomv2_RECO_500um_500um_100ps_AllEtaPhi_PU200/160506_083146/0000/step3_107.root',
                                                              
                                                              
        ))
                                                              
process.me0SegAna = cms.EDAnalyzer('ME0SegmentAnalyzerMuonGun',
                                  wp = cms.string("tight"),
                                  minEta = cms.double(2.0),
                                  maxEta = cms.double(2.8),
                                  #    etaMin = cms.double(2.0),
                                  #    etaMax = cms.double(2.8),
                                  #    dr = cms.double(0.25),
                                  #    ptMin = cms.double(5.0),
                                  timeMin = cms.double(18.3 - 3*0.1),
                                  timeMax = cms.double(18.3 + 3*0.1)
)


process.TFileService = cms.Service("TFileService",
fileName = cms.string("histoME0SegMuGun.root")
)

process.p = cms.Path(process.me0SegAna)


