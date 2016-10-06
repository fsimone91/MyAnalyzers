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
                            # replace 'myfile.root',',',',',',' with the source file you want to use
                            fileNames = cms.untracked.vstring(
                                                              #'file:./out_sim.root','
                                                              'file:/lustre/cms/store/user/rosma/TauTo3mu/crab_Tau3Mu_14TeV_GEN-SIM_2023HGCalME0Geomv2_RECO_500um_500um_10ps_PU200/160719_103707/0000/step3_10.root',
                                                              'file:/lustre/cms/store/user/rosma/TauTo3mu/crab_Tau3Mu_14TeV_GEN-SIM_2023HGCalME0Geomv2_RECO_500um_500um_10ps_PU200/160719_103707/0000/step3_100.root',
                                                              'file:/lustre/cms/store/user/rosma/TauTo3mu/crab_Tau3Mu_14TeV_GEN-SIM_2023HGCalME0Geomv2_RECO_500um_500um_10ps_PU200/160719_103707/0000/step3_101.root',
                                                              'file:/lustre/cms/store/user/rosma/TauTo3mu/crab_Tau3Mu_14TeV_GEN-SIM_2023HGCalME0Geomv2_RECO_500um_500um_10ps_PU200/160719_103707/0000/step3_103.root',
                                                              'file:/lustre/cms/store/user/rosma/TauTo3mu/crab_Tau3Mu_14TeV_GEN-SIM_2023HGCalME0Geomv2_RECO_500um_500um_10ps_PU200/160719_103707/0000/step3_104.root',
                                                              'file:/lustre/cms/store/user/rosma/TauTo3mu/crab_Tau3Mu_14TeV_GEN-SIM_2023HGCalME0Geomv2_RECO_500um_500um_10ps_PU200/160719_103707/0000/step3_105.root',
                                                              'file:/lustre/cms/store/user/rosma/TauTo3mu/crab_Tau3Mu_14TeV_GEN-SIM_2023HGCalME0Geomv2_RECO_500um_500um_10ps_PU200/160719_103707/0000/step3_106.root',
                                                              'file:/lustre/cms/store/user/rosma/TauTo3mu/crab_Tau3Mu_14TeV_GEN-SIM_2023HGCalME0Geomv2_RECO_500um_500um_10ps_PU200/160719_103707/0000/step3_107.root',
                                                              'file:/lustre/cms/store/user/rosma/TauTo3mu/crab_Tau3Mu_14TeV_GEN-SIM_2023HGCalME0Geomv2_RECO_500um_500um_10ps_PU200/160719_103707/0000/step3_108.root',
                                                              'file:/lustre/cms/store/user/rosma/TauTo3mu/crab_Tau3Mu_14TeV_GEN-SIM_2023HGCalME0Geomv2_RECO_500um_500um_10ps_PU200/160719_103707/0000/step3_109.root',
                                                              'file:/lustre/cms/store/user/rosma/TauTo3mu/crab_Tau3Mu_14TeV_GEN-SIM_2023HGCalME0Geomv2_RECO_500um_500um_10ps_PU200/160719_103707/0000/step3_11.root',
                                                              'file:/lustre/cms/store/user/rosma/TauTo3mu/crab_Tau3Mu_14TeV_GEN-SIM_2023HGCalME0Geomv2_RECO_500um_500um_10ps_PU200/160719_103707/0000/step3_110.root',
                                                              'file:/lustre/cms/store/user/rosma/TauTo3mu/crab_Tau3Mu_14TeV_GEN-SIM_2023HGCalME0Geomv2_RECO_500um_500um_10ps_PU200/160719_103707/0000/step3_111.root',
                                                              'file:/lustre/cms/store/user/rosma/TauTo3mu/crab_Tau3Mu_14TeV_GEN-SIM_2023HGCalME0Geomv2_RECO_500um_500um_10ps_PU200/160719_103707/0000/step3_112.root',
                                                              'file:/lustre/cms/store/user/rosma/TauTo3mu/crab_Tau3Mu_14TeV_GEN-SIM_2023HGCalME0Geomv2_RECO_500um_500um_10ps_PU200/160719_103707/0000/step3_113.root',
                                                              'file:/lustre/cms/store/user/rosma/TauTo3mu/crab_Tau3Mu_14TeV_GEN-SIM_2023HGCalME0Geomv2_RECO_500um_500um_10ps_PU200/160719_103707/0000/step3_114.root',
                                                              'file:/lustre/cms/store/user/rosma/TauTo3mu/crab_Tau3Mu_14TeV_GEN-SIM_2023HGCalME0Geomv2_RECO_500um_500um_10ps_PU200/160719_103707/0000/step3_115.root',
                                                              'file:/lustre/cms/store/user/rosma/TauTo3mu/crab_Tau3Mu_14TeV_GEN-SIM_2023HGCalME0Geomv2_RECO_500um_500um_10ps_PU200/160719_103707/0000/step3_116.root',
                                                              'file:/lustre/cms/store/user/rosma/TauTo3mu/crab_Tau3Mu_14TeV_GEN-SIM_2023HGCalME0Geomv2_RECO_500um_500um_10ps_PU200/160719_103707/0000/step3_117.root',
                                                              'file:/lustre/cms/store/user/rosma/TauTo3mu/crab_Tau3Mu_14TeV_GEN-SIM_2023HGCalME0Geomv2_RECO_500um_500um_10ps_PU200/160719_103707/0000/step3_118.root',
                                                              'file:/lustre/cms/store/user/rosma/TauTo3mu/crab_Tau3Mu_14TeV_GEN-SIM_2023HGCalME0Geomv2_RECO_500um_500um_10ps_PU200/160719_103707/0000/step3_119.root',
                                                              'file:/lustre/cms/store/user/rosma/TauTo3mu/crab_Tau3Mu_14TeV_GEN-SIM_2023HGCalME0Geomv2_RECO_500um_500um_10ps_PU200/160719_103707/0000/step3_12.root',
                                                              'file:/lustre/cms/store/user/rosma/TauTo3mu/crab_Tau3Mu_14TeV_GEN-SIM_2023HGCalME0Geomv2_RECO_500um_500um_10ps_PU200/160719_103707/0000/step3_120.root',
                                                              'file:/lustre/cms/store/user/rosma/TauTo3mu/crab_Tau3Mu_14TeV_GEN-SIM_2023HGCalME0Geomv2_RECO_500um_500um_10ps_PU200/160719_103707/0000/step3_121.root',
                                                              'file:/lustre/cms/store/user/rosma/TauTo3mu/crab_Tau3Mu_14TeV_GEN-SIM_2023HGCalME0Geomv2_RECO_500um_500um_10ps_PU200/160719_103707/0000/step3_122.root',
                                                              'file:/lustre/cms/store/user/rosma/TauTo3mu/crab_Tau3Mu_14TeV_GEN-SIM_2023HGCalME0Geomv2_RECO_500um_500um_10ps_PU200/160719_103707/0000/step3_123.root',
                                                              'file:/lustre/cms/store/user/rosma/TauTo3mu/crab_Tau3Mu_14TeV_GEN-SIM_2023HGCalME0Geomv2_RECO_500um_500um_10ps_PU200/160719_103707/0000/step3_124.root',
                                                              'file:/lustre/cms/store/user/rosma/TauTo3mu/crab_Tau3Mu_14TeV_GEN-SIM_2023HGCalME0Geomv2_RECO_500um_500um_10ps_PU200/160719_103707/0000/step3_125.root',
                                                              'file:/lustre/cms/store/user/rosma/TauTo3mu/crab_Tau3Mu_14TeV_GEN-SIM_2023HGCalME0Geomv2_RECO_500um_500um_10ps_PU200/160719_103707/0000/step3_126.root',
                                                              'file:/lustre/cms/store/user/rosma/TauTo3mu/crab_Tau3Mu_14TeV_GEN-SIM_2023HGCalME0Geomv2_RECO_500um_500um_10ps_PU200/160719_103707/0000/step3_127.root',
                                                              'file:/lustre/cms/store/user/rosma/TauTo3mu/crab_Tau3Mu_14TeV_GEN-SIM_2023HGCalME0Geomv2_RECO_500um_500um_10ps_PU200/160719_103707/0000/step3_128.root',
                                                              'file:/lustre/cms/store/user/rosma/TauTo3mu/crab_Tau3Mu_14TeV_GEN-SIM_2023HGCalME0Geomv2_RECO_500um_500um_10ps_PU200/160719_103707/0000/step3_129.root',
                                                              'file:/lustre/cms/store/user/rosma/TauTo3mu/crab_Tau3Mu_14TeV_GEN-SIM_2023HGCalME0Geomv2_RECO_500um_500um_10ps_PU200/160719_103707/0000/step3_13.root',
                                                              'file:/lustre/cms/store/user/rosma/TauTo3mu/crab_Tau3Mu_14TeV_GEN-SIM_2023HGCalME0Geomv2_RECO_500um_500um_10ps_PU200/160719_103707/0000/step3_130.root',
                                                              'file:/lustre/cms/store/user/rosma/TauTo3mu/crab_Tau3Mu_14TeV_GEN-SIM_2023HGCalME0Geomv2_RECO_500um_500um_10ps_PU200/160719_103707/0000/step3_131.root',
                                                              'file:/lustre/cms/store/user/rosma/TauTo3mu/crab_Tau3Mu_14TeV_GEN-SIM_2023HGCalME0Geomv2_RECO_500um_500um_10ps_PU200/160719_103707/0000/step3_132.root',
                                                              'file:/lustre/cms/store/user/rosma/TauTo3mu/crab_Tau3Mu_14TeV_GEN-SIM_2023HGCalME0Geomv2_RECO_500um_500um_10ps_PU200/160719_103707/0000/step3_133.root',
                                                              'file:/lustre/cms/store/user/rosma/TauTo3mu/crab_Tau3Mu_14TeV_GEN-SIM_2023HGCalME0Geomv2_RECO_500um_500um_10ps_PU200/160719_103707/0000/step3_134.root',
                                                              'file:/lustre/cms/store/user/rosma/TauTo3mu/crab_Tau3Mu_14TeV_GEN-SIM_2023HGCalME0Geomv2_RECO_500um_500um_10ps_PU200/160719_103707/0000/step3_135.root',
                                                              'file:/lustre/cms/store/user/rosma/TauTo3mu/crab_Tau3Mu_14TeV_GEN-SIM_2023HGCalME0Geomv2_RECO_500um_500um_10ps_PU200/160719_103707/0000/step3_136.root',
                                                              'file:/lustre/cms/store/user/rosma/TauTo3mu/crab_Tau3Mu_14TeV_GEN-SIM_2023HGCalME0Geomv2_RECO_500um_500um_10ps_PU200/160719_103707/0000/step3_137.root',
                                                              'file:/lustre/cms/store/user/rosma/TauTo3mu/crab_Tau3Mu_14TeV_GEN-SIM_2023HGCalME0Geomv2_RECO_500um_500um_10ps_PU200/160719_103707/0000/step3_138.root',
                                                              'file:/lustre/cms/store/user/rosma/TauTo3mu/crab_Tau3Mu_14TeV_GEN-SIM_2023HGCalME0Geomv2_RECO_500um_500um_10ps_PU200/160719_103707/0000/step3_139.root',
                                                              'file:/lustre/cms/store/user/rosma/TauTo3mu/crab_Tau3Mu_14TeV_GEN-SIM_2023HGCalME0Geomv2_RECO_500um_500um_10ps_PU200/160719_103707/0000/step3_14.root',
                                                              'file:/lustre/cms/store/user/rosma/TauTo3mu/crab_Tau3Mu_14TeV_GEN-SIM_2023HGCalME0Geomv2_RECO_500um_500um_10ps_PU200/160719_103707/0000/step3_140.root',
                                                              'file:/lustre/cms/store/user/rosma/TauTo3mu/crab_Tau3Mu_14TeV_GEN-SIM_2023HGCalME0Geomv2_RECO_500um_500um_10ps_PU200/160719_103707/0000/step3_141.root',
                                                              'file:/lustre/cms/store/user/rosma/TauTo3mu/crab_Tau3Mu_14TeV_GEN-SIM_2023HGCalME0Geomv2_RECO_500um_500um_10ps_PU200/160719_103707/0000/step3_142.root',
                                                              'file:/lustre/cms/store/user/rosma/TauTo3mu/crab_Tau3Mu_14TeV_GEN-SIM_2023HGCalME0Geomv2_RECO_500um_500um_10ps_PU200/160719_103707/0000/step3_143.root',
                                                              'file:/lustre/cms/store/user/rosma/TauTo3mu/crab_Tau3Mu_14TeV_GEN-SIM_2023HGCalME0Geomv2_RECO_500um_500um_10ps_PU200/160719_103707/0000/step3_144.root',
                                                              'file:/lustre/cms/store/user/rosma/TauTo3mu/crab_Tau3Mu_14TeV_GEN-SIM_2023HGCalME0Geomv2_RECO_500um_500um_10ps_PU200/160719_103707/0000/step3_145.root',
                                                              'file:/lustre/cms/store/user/rosma/TauTo3mu/crab_Tau3Mu_14TeV_GEN-SIM_2023HGCalME0Geomv2_RECO_500um_500um_10ps_PU200/160719_103707/0000/step3_146.root',
                                                              'file:/lustre/cms/store/user/rosma/TauTo3mu/crab_Tau3Mu_14TeV_GEN-SIM_2023HGCalME0Geomv2_RECO_500um_500um_10ps_PU200/160719_103707/0000/step3_147.root',
                                                              'file:/lustre/cms/store/user/rosma/TauTo3mu/crab_Tau3Mu_14TeV_GEN-SIM_2023HGCalME0Geomv2_RECO_500um_500um_10ps_PU200/160719_103707/0000/step3_148.root',
                                                              'file:/lustre/cms/store/user/rosma/TauTo3mu/crab_Tau3Mu_14TeV_GEN-SIM_2023HGCalME0Geomv2_RECO_500um_500um_10ps_PU200/160719_103707/0000/step3_149.root',
                                                              'file:/lustre/cms/store/user/rosma/TauTo3mu/crab_Tau3Mu_14TeV_GEN-SIM_2023HGCalME0Geomv2_RECO_500um_500um_10ps_PU200/160719_103707/0000/step3_15.root',
                                                              'file:/lustre/cms/store/user/rosma/TauTo3mu/crab_Tau3Mu_14TeV_GEN-SIM_2023HGCalME0Geomv2_RECO_500um_500um_10ps_PU200/160719_103707/0000/step3_150.root',
                                                              'file:/lustre/cms/store/user/rosma/TauTo3mu/crab_Tau3Mu_14TeV_GEN-SIM_2023HGCalME0Geomv2_RECO_500um_500um_10ps_PU200/160719_103707/0000/step3_151.root',
                                                              'file:/lustre/cms/store/user/rosma/TauTo3mu/crab_Tau3Mu_14TeV_GEN-SIM_2023HGCalME0Geomv2_RECO_500um_500um_10ps_PU200/160719_103707/0000/step3_152.root',
                                                              'file:/lustre/cms/store/user/rosma/TauTo3mu/crab_Tau3Mu_14TeV_GEN-SIM_2023HGCalME0Geomv2_RECO_500um_500um_10ps_PU200/160719_103707/0000/step3_153.root',
                                                              'file:/lustre/cms/store/user/rosma/TauTo3mu/crab_Tau3Mu_14TeV_GEN-SIM_2023HGCalME0Geomv2_RECO_500um_500um_10ps_PU200/160719_103707/0000/step3_154.root',
                                                              'file:/lustre/cms/store/user/rosma/TauTo3mu/crab_Tau3Mu_14TeV_GEN-SIM_2023HGCalME0Geomv2_RECO_500um_500um_10ps_PU200/160719_103707/0000/step3_155.root',
                                                              'file:/lustre/cms/store/user/rosma/TauTo3mu/crab_Tau3Mu_14TeV_GEN-SIM_2023HGCalME0Geomv2_RECO_500um_500um_10ps_PU200/160719_103707/0000/step3_156.root',
                                                              'file:/lustre/cms/store/user/rosma/TauTo3mu/crab_Tau3Mu_14TeV_GEN-SIM_2023HGCalME0Geomv2_RECO_500um_500um_10ps_PU200/160719_103707/0000/step3_157.root',
                                                              'file:/lustre/cms/store/user/rosma/TauTo3mu/crab_Tau3Mu_14TeV_GEN-SIM_2023HGCalME0Geomv2_RECO_500um_500um_10ps_PU200/160719_103707/0000/step3_158.root',
                                                              'file:/lustre/cms/store/user/rosma/TauTo3mu/crab_Tau3Mu_14TeV_GEN-SIM_2023HGCalME0Geomv2_RECO_500um_500um_10ps_PU200/160719_103707/0000/step3_159.root',
                                                              'file:/lustre/cms/store/user/rosma/TauTo3mu/crab_Tau3Mu_14TeV_GEN-SIM_2023HGCalME0Geomv2_RECO_500um_500um_10ps_PU200/160719_103707/0000/step3_16.root',
                                                              'file:/lustre/cms/store/user/rosma/TauTo3mu/crab_Tau3Mu_14TeV_GEN-SIM_2023HGCalME0Geomv2_RECO_500um_500um_10ps_PU200/160719_103707/0000/step3_160.root',
                                                              'file:/lustre/cms/store/user/rosma/TauTo3mu/crab_Tau3Mu_14TeV_GEN-SIM_2023HGCalME0Geomv2_RECO_500um_500um_10ps_PU200/160719_103707/0000/step3_161.root',
                                                              'file:/lustre/cms/store/user/rosma/TauTo3mu/crab_Tau3Mu_14TeV_GEN-SIM_2023HGCalME0Geomv2_RECO_500um_500um_10ps_PU200/160719_103707/0000/step3_162.root',
                                                              'file:/lustre/cms/store/user/rosma/TauTo3mu/crab_Tau3Mu_14TeV_GEN-SIM_2023HGCalME0Geomv2_RECO_500um_500um_10ps_PU200/160719_103707/0000/step3_164.root',
                                                              'file:/lustre/cms/store/user/rosma/TauTo3mu/crab_Tau3Mu_14TeV_GEN-SIM_2023HGCalME0Geomv2_RECO_500um_500um_10ps_PU200/160719_103707/0000/step3_165.root',
                                                              'file:/lustre/cms/store/user/rosma/TauTo3mu/crab_Tau3Mu_14TeV_GEN-SIM_2023HGCalME0Geomv2_RECO_500um_500um_10ps_PU200/160719_103707/0000/step3_166.root',
                                                              'file:/lustre/cms/store/user/rosma/TauTo3mu/crab_Tau3Mu_14TeV_GEN-SIM_2023HGCalME0Geomv2_RECO_500um_500um_10ps_PU200/160719_103707/0000/step3_167.root',
                                                              'file:/lustre/cms/store/user/rosma/TauTo3mu/crab_Tau3Mu_14TeV_GEN-SIM_2023HGCalME0Geomv2_RECO_500um_500um_10ps_PU200/160719_103707/0000/step3_168.root',
                                                              'file:/lustre/cms/store/user/rosma/TauTo3mu/crab_Tau3Mu_14TeV_GEN-SIM_2023HGCalME0Geomv2_RECO_500um_500um_10ps_PU200/160719_103707/0000/step3_169.root',
                                                              'file:/lustre/cms/store/user/rosma/TauTo3mu/crab_Tau3Mu_14TeV_GEN-SIM_2023HGCalME0Geomv2_RECO_500um_500um_10ps_PU200/160719_103707/0000/step3_17.root',
                                                              'file:/lustre/cms/store/user/rosma/TauTo3mu/crab_Tau3Mu_14TeV_GEN-SIM_2023HGCalME0Geomv2_RECO_500um_500um_10ps_PU200/160719_103707/0000/step3_170.root',
                                                              'file:/lustre/cms/store/user/rosma/TauTo3mu/crab_Tau3Mu_14TeV_GEN-SIM_2023HGCalME0Geomv2_RECO_500um_500um_10ps_PU200/160719_103707/0000/step3_171.root',
                                                              'file:/lustre/cms/store/user/rosma/TauTo3mu/crab_Tau3Mu_14TeV_GEN-SIM_2023HGCalME0Geomv2_RECO_500um_500um_10ps_PU200/160719_103707/0000/step3_172.root',
                                                              'file:/lustre/cms/store/user/rosma/TauTo3mu/crab_Tau3Mu_14TeV_GEN-SIM_2023HGCalME0Geomv2_RECO_500um_500um_10ps_PU200/160719_103707/0000/step3_173.root',


                                                              
        
        ))
                                                              
process.me0SegAna = cms.EDAnalyzer('ME0SegmentAnalyzerSim',
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
fileName = cms.string("histoME0SimTau3Mu.root")
)

process.p = cms.Path(process.me0SegAna)


