import FWCore.ParameterSet.Config as cms

from Configuration.StandardSequences.Eras import eras
process = cms.Process("Demo",eras.Phase2C1)
process.load("FWCore.MessageService.MessageLogger_cfi")
process.load('Configuration.StandardSequences.Services_cff')
#process.load('Configuration.Geometry.GeometryExtended2023D1Reco_cff')
process.load('Configuration.Geometry.GeometryExtended2023D6Reco_cff')
#process.load('Configuration.StandardSequences.MagneticField_38T_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
#process.load('RecoMuon.MuonIdentification.me0MuonReco_cff')



from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:run2_mc', '')


process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(200) )

process.source = cms.Source("PoolSource",
                            # replace 'myfile.root',',',',',',',' with the source file you want to use
                            fileNames = cms.untracked.vstring(
                                                              'root://cmsxrootd.fnal.gov//store/user/lpcgem/ME0TDRStudies/MuGun_0p5_30/90X_DIGI_PU0_WithoutNoise/job_1.root',
                                                              'root://cmsxrootd.fnal.gov//store/user/lpcgem/ME0TDRStudies/MuGun_0p5_30/90X_DIGI_PU0_WithoutNoise/job_2.root',
                                                              'root://cmsxrootd.fnal.gov//store/user/lpcgem/ME0TDRStudies/MuGun_0p5_30/90X_DIGI_PU0_WithoutNoise/job_3.root',
                                                              'root://cmsxrootd.fnal.gov//store/user/lpcgem/ME0TDRStudies/MuGun_0p5_30/90X_DIGI_PU0_WithoutNoise/job_4.root',
                                                              'root://cmsxrootd.fnal.gov//store/user/lpcgem/ME0TDRStudies/MuGun_0p5_30/90X_DIGI_PU0_WithoutNoise/job_5.root',
                                                              'root://cmsxrootd.fnal.gov//store/user/lpcgem/ME0TDRStudies/MuGun_0p5_30/90X_DIGI_PU0_WithoutNoise/job_6.root',
                                                              'root://cmsxrootd.fnal.gov//store/user/lpcgem/ME0TDRStudies/MuGun_0p5_30/90X_DIGI_PU0_WithoutNoise/job_7.root',
                                                              'root://cmsxrootd.fnal.gov//store/user/lpcgem/ME0TDRStudies/MuGun_0p5_30/90X_DIGI_PU0_WithoutNoise/job_8.root',
                                                              'root://cmsxrootd.fnal.gov//store/user/lpcgem/ME0TDRStudies/MuGun_0p5_30/90X_DIGI_PU0_WithoutNoise/job_9.root',
                                                              'root://cmsxrootd.fnal.gov//store/user/lpcgem/ME0TDRStudies/MuGun_0p5_30/90X_DIGI_PU0_WithoutNoise/job_10.root',
                                                              'root://cmsxrootd.fnal.gov//store/user/lpcgem/ME0TDRStudies/MuGun_0p5_30/90X_DIGI_PU0_WithoutNoise/job_11.root',
                                                              'root://cmsxrootd.fnal.gov//store/user/lpcgem/ME0TDRStudies/MuGun_0p5_30/90X_DIGI_PU0_WithoutNoise/job_12.root',
                                                              'root://cmsxrootd.fnal.gov//store/user/lpcgem/ME0TDRStudies/MuGun_0p5_30/90X_DIGI_PU0_WithoutNoise/job_13.root',
                                                              'root://cmsxrootd.fnal.gov//store/user/lpcgem/ME0TDRStudies/MuGun_0p5_30/90X_DIGI_PU0_WithoutNoise/job_14.root',
                                                              'root://cmsxrootd.fnal.gov//store/user/lpcgem/ME0TDRStudies/MuGun_0p5_30/90X_DIGI_PU0_WithoutNoise/job_15.root',
                                                              'root://cmsxrootd.fnal.gov//store/user/lpcgem/ME0TDRStudies/MuGun_0p5_30/90X_DIGI_PU0_WithoutNoise/job_16.root',
                                                              'root://cmsxrootd.fnal.gov//store/user/lpcgem/ME0TDRStudies/MuGun_0p5_30/90X_DIGI_PU0_WithoutNoise/job_17.root',
                                                              'root://cmsxrootd.fnal.gov//store/user/lpcgem/ME0TDRStudies/MuGun_0p5_30/90X_DIGI_PU0_WithoutNoise/job_18.root',
                                                              'root://cmsxrootd.fnal.gov//store/user/lpcgem/ME0TDRStudies/MuGun_0p5_30/90X_DIGI_PU0_WithoutNoise/job_19.root',
                                                              'root://cmsxrootd.fnal.gov//store/user/lpcgem/ME0TDRStudies/MuGun_0p5_30/90X_DIGI_PU0_WithoutNoise/job_20.root',
                                                              'root://cmsxrootd.fnal.gov//store/user/lpcgem/ME0TDRStudies/MuGun_0p5_30/90X_DIGI_PU0_WithoutNoise/job_21.root',
                                                              'root://cmsxrootd.fnal.gov//store/user/lpcgem/ME0TDRStudies/MuGun_0p5_30/90X_DIGI_PU0_WithoutNoise/job_22.root',
                                                              'root://cmsxrootd.fnal.gov//store/user/lpcgem/ME0TDRStudies/MuGun_0p5_30/90X_DIGI_PU0_WithoutNoise/job_23.root',
                                                              'root://cmsxrootd.fnal.gov//store/user/lpcgem/ME0TDRStudies/MuGun_0p5_30/90X_DIGI_PU0_WithoutNoise/job_24.root',
                                                              'root://cmsxrootd.fnal.gov//store/user/lpcgem/ME0TDRStudies/MuGun_0p5_30/90X_DIGI_PU0_WithoutNoise/job_25.root',
                                                              'root://cmsxrootd.fnal.gov//store/user/lpcgem/ME0TDRStudies/MuGun_0p5_30/90X_DIGI_PU0_WithoutNoise/job_26.root',
                                                              'root://cmsxrootd.fnal.gov//store/user/lpcgem/ME0TDRStudies/MuGun_0p5_30/90X_DIGI_PU0_WithoutNoise/job_27.root',
                                                              'root://cmsxrootd.fnal.gov//store/user/lpcgem/ME0TDRStudies/MuGun_0p5_30/90X_DIGI_PU0_WithoutNoise/job_28.root',
                                                              'root://cmsxrootd.fnal.gov//store/user/lpcgem/ME0TDRStudies/MuGun_0p5_30/90X_DIGI_PU0_WithoutNoise/job_29.root',
                                                              'root://cmsxrootd.fnal.gov//store/user/lpcgem/ME0TDRStudies/MuGun_0p5_30/90X_DIGI_PU0_WithoutNoise/job_30.root',
                                                              'root://cmsxrootd.fnal.gov//store/user/lpcgem/ME0TDRStudies/MuGun_0p5_30/90X_DIGI_PU0_WithoutNoise/job_31.root',
                                                              'root://cmsxrootd.fnal.gov//store/user/lpcgem/ME0TDRStudies/MuGun_0p5_30/90X_DIGI_PU0_WithoutNoise/job_32.root',
                                                              'root://cmsxrootd.fnal.gov//store/user/lpcgem/ME0TDRStudies/MuGun_0p5_30/90X_DIGI_PU0_WithoutNoise/job_33.root',
                                                              'root://cmsxrootd.fnal.gov//store/user/lpcgem/ME0TDRStudies/MuGun_0p5_30/90X_DIGI_PU0_WithoutNoise/job_34.root',
                                                              'root://cmsxrootd.fnal.gov//store/user/lpcgem/ME0TDRStudies/MuGun_0p5_30/90X_DIGI_PU0_WithoutNoise/job_35.root',
                                                              'root://cmsxrootd.fnal.gov//store/user/lpcgem/ME0TDRStudies/MuGun_0p5_30/90X_DIGI_PU0_WithoutNoise/job_36.root',
                                                              'root://cmsxrootd.fnal.gov//store/user/lpcgem/ME0TDRStudies/MuGun_0p5_30/90X_DIGI_PU0_WithoutNoise/job_37.root',
                                                              'root://cmsxrootd.fnal.gov//store/user/lpcgem/ME0TDRStudies/MuGun_0p5_30/90X_DIGI_PU0_WithoutNoise/job_38.root',
                                                              'root://cmsxrootd.fnal.gov//store/user/lpcgem/ME0TDRStudies/MuGun_0p5_30/90X_DIGI_PU0_WithoutNoise/job_39.root',
                                                              'root://cmsxrootd.fnal.gov//store/user/lpcgem/ME0TDRStudies/MuGun_0p5_30/90X_DIGI_PU0_WithoutNoise/job_40.root',
                                                              'root://cmsxrootd.fnal.gov//store/user/lpcgem/ME0TDRStudies/MuGun_0p5_30/90X_DIGI_PU0_WithoutNoise/job_41.root',
                                                              'root://cmsxrootd.fnal.gov//store/user/lpcgem/ME0TDRStudies/MuGun_0p5_30/90X_DIGI_PU0_WithoutNoise/job_42.root',
                                                              'root://cmsxrootd.fnal.gov//store/user/lpcgem/ME0TDRStudies/MuGun_0p5_30/90X_DIGI_PU0_WithoutNoise/job_43.root',
                                                              'root://cmsxrootd.fnal.gov//store/user/lpcgem/ME0TDRStudies/MuGun_0p5_30/90X_DIGI_PU0_WithoutNoise/job_44.root',
                                                              'root://cmsxrootd.fnal.gov//store/user/lpcgem/ME0TDRStudies/MuGun_0p5_30/90X_DIGI_PU0_WithoutNoise/job_45.root',
                                                              'root://cmsxrootd.fnal.gov//store/user/lpcgem/ME0TDRStudies/MuGun_0p5_30/90X_DIGI_PU0_WithoutNoise/job_46.root',
                                                              'root://cmsxrootd.fnal.gov//store/user/lpcgem/ME0TDRStudies/MuGun_0p5_30/90X_DIGI_PU0_WithoutNoise/job_47.root',
                                                              'root://cmsxrootd.fnal.gov//store/user/lpcgem/ME0TDRStudies/MuGun_0p5_30/90X_DIGI_PU0_WithoutNoise/job_48.root',
                                                              'root://cmsxrootd.fnal.gov//store/user/lpcgem/ME0TDRStudies/MuGun_0p5_30/90X_DIGI_PU0_WithoutNoise/job_49.root',
                                                              'root://cmsxrootd.fnal.gov//store/user/lpcgem/ME0TDRStudies/MuGun_0p5_30/90X_DIGI_PU0_WithoutNoise/job_50.root',
                                                              'root://cmsxrootd.fnal.gov//store/user/lpcgem/ME0TDRStudies/MuGun_0p5_30/90X_DIGI_PU0_WithoutNoise/job_51.root',
                                                              'root://cmsxrootd.fnal.gov//store/user/lpcgem/ME0TDRStudies/MuGun_0p5_30/90X_DIGI_PU0_WithoutNoise/job_52.root',
                                                              'root://cmsxrootd.fnal.gov//store/user/lpcgem/ME0TDRStudies/MuGun_0p5_30/90X_DIGI_PU0_WithoutNoise/job_53.root',
                                                              'root://cmsxrootd.fnal.gov//store/user/lpcgem/ME0TDRStudies/MuGun_0p5_30/90X_DIGI_PU0_WithoutNoise/job_54.root',
                                                              'root://cmsxrootd.fnal.gov//store/user/lpcgem/ME0TDRStudies/MuGun_0p5_30/90X_DIGI_PU0_WithoutNoise/job_55.root',
                                                              'root://cmsxrootd.fnal.gov//store/user/lpcgem/ME0TDRStudies/MuGun_0p5_30/90X_DIGI_PU0_WithoutNoise/job_56.root',
                                                              'root://cmsxrootd.fnal.gov//store/user/lpcgem/ME0TDRStudies/MuGun_0p5_30/90X_DIGI_PU0_WithoutNoise/job_57.root',
                                                              'root://cmsxrootd.fnal.gov//store/user/lpcgem/ME0TDRStudies/MuGun_0p5_30/90X_DIGI_PU0_WithoutNoise/job_58.root',
                                                              'root://cmsxrootd.fnal.gov//store/user/lpcgem/ME0TDRStudies/MuGun_0p5_30/90X_DIGI_PU0_WithoutNoise/job_59.root',
                                                              'root://cmsxrootd.fnal.gov//store/user/lpcgem/ME0TDRStudies/MuGun_0p5_30/90X_DIGI_PU0_WithoutNoise/job_60.root',
                                                              'root://cmsxrootd.fnal.gov//store/user/lpcgem/ME0TDRStudies/MuGun_0p5_30/90X_DIGI_PU0_WithoutNoise/job_61.root',
                                                              'root://cmsxrootd.fnal.gov//store/user/lpcgem/ME0TDRStudies/MuGun_0p5_30/90X_DIGI_PU0_WithoutNoise/job_62.root',
                                                              'root://cmsxrootd.fnal.gov//store/user/lpcgem/ME0TDRStudies/MuGun_0p5_30/90X_DIGI_PU0_WithoutNoise/job_63.root',
                                                              'root://cmsxrootd.fnal.gov//store/user/lpcgem/ME0TDRStudies/MuGun_0p5_30/90X_DIGI_PU0_WithoutNoise/job_64.root',
                                                              'root://cmsxrootd.fnal.gov//store/user/lpcgem/ME0TDRStudies/MuGun_0p5_30/90X_DIGI_PU0_WithoutNoise/job_65.root',
                                                              'root://cmsxrootd.fnal.gov//store/user/lpcgem/ME0TDRStudies/MuGun_0p5_30/90X_DIGI_PU0_WithoutNoise/job_66.root',
                                                              'root://cmsxrootd.fnal.gov//store/user/lpcgem/ME0TDRStudies/MuGun_0p5_30/90X_DIGI_PU0_WithoutNoise/job_67.root',
                                                              'root://cmsxrootd.fnal.gov//store/user/lpcgem/ME0TDRStudies/MuGun_0p5_30/90X_DIGI_PU0_WithoutNoise/job_68.root',
                                                              'root://cmsxrootd.fnal.gov//store/user/lpcgem/ME0TDRStudies/MuGun_0p5_30/90X_DIGI_PU0_WithoutNoise/job_69.root',
                                                              'root://cmsxrootd.fnal.gov//store/user/lpcgem/ME0TDRStudies/MuGun_0p5_30/90X_DIGI_PU0_WithoutNoise/job_70.root',
                                                              'root://cmsxrootd.fnal.gov//store/user/lpcgem/ME0TDRStudies/MuGun_0p5_30/90X_DIGI_PU0_WithoutNoise/job_71.root',
                                                              'root://cmsxrootd.fnal.gov//store/user/lpcgem/ME0TDRStudies/MuGun_0p5_30/90X_DIGI_PU0_WithoutNoise/job_72.root',
                                                              'root://cmsxrootd.fnal.gov//store/user/lpcgem/ME0TDRStudies/MuGun_0p5_30/90X_DIGI_PU0_WithoutNoise/job_73.root',
                                                              'root://cmsxrootd.fnal.gov//store/user/lpcgem/ME0TDRStudies/MuGun_0p5_30/90X_DIGI_PU0_WithoutNoise/job_74.root',
                                                              'root://cmsxrootd.fnal.gov//store/user/lpcgem/ME0TDRStudies/MuGun_0p5_30/90X_DIGI_PU0_WithoutNoise/job_75.root',
                                                              'root://cmsxrootd.fnal.gov//store/user/lpcgem/ME0TDRStudies/MuGun_0p5_30/90X_DIGI_PU0_WithoutNoise/job_76.root',
                                                              'root://cmsxrootd.fnal.gov//store/user/lpcgem/ME0TDRStudies/MuGun_0p5_30/90X_DIGI_PU0_WithoutNoise/job_77.root',
                                                              'root://cmsxrootd.fnal.gov//store/user/lpcgem/ME0TDRStudies/MuGun_0p5_30/90X_DIGI_PU0_WithoutNoise/job_78.root',
                                                              'root://cmsxrootd.fnal.gov//store/user/lpcgem/ME0TDRStudies/MuGun_0p5_30/90X_DIGI_PU0_WithoutNoise/job_79.root',
                                                              'root://cmsxrootd.fnal.gov//store/user/lpcgem/ME0TDRStudies/MuGun_0p5_30/90X_DIGI_PU0_WithoutNoise/job_80.root',
                                                              'root://cmsxrootd.fnal.gov//store/user/lpcgem/ME0TDRStudies/MuGun_0p5_30/90X_DIGI_PU0_WithoutNoise/job_81.root',
                                                              'root://cmsxrootd.fnal.gov//store/user/lpcgem/ME0TDRStudies/MuGun_0p5_30/90X_DIGI_PU0_WithoutNoise/job_82.root',
                                                              'root://cmsxrootd.fnal.gov//store/user/lpcgem/ME0TDRStudies/MuGun_0p5_30/90X_DIGI_PU0_WithoutNoise/job_83.root',
                                                              'root://cmsxrootd.fnal.gov//store/user/lpcgem/ME0TDRStudies/MuGun_0p5_30/90X_DIGI_PU0_WithoutNoise/job_84.root',
                                                              'root://cmsxrootd.fnal.gov//store/user/lpcgem/ME0TDRStudies/MuGun_0p5_30/90X_DIGI_PU0_WithoutNoise/job_85.root',
                                                              'root://cmsxrootd.fnal.gov//store/user/lpcgem/ME0TDRStudies/MuGun_0p5_30/90X_DIGI_PU0_WithoutNoise/job_86.root',
                                                              'root://cmsxrootd.fnal.gov//store/user/lpcgem/ME0TDRStudies/MuGun_0p5_30/90X_DIGI_PU0_WithoutNoise/job_87.root',
                                                              'root://cmsxrootd.fnal.gov//store/user/lpcgem/ME0TDRStudies/MuGun_0p5_30/90X_DIGI_PU0_WithoutNoise/job_88.root',
                                                              'root://cmsxrootd.fnal.gov//store/user/lpcgem/ME0TDRStudies/MuGun_0p5_30/90X_DIGI_PU0_WithoutNoise/job_89.root',
                                                              'root://cmsxrootd.fnal.gov//store/user/lpcgem/ME0TDRStudies/MuGun_0p5_30/90X_DIGI_PU0_WithoutNoise/job_90.root',
                            
               
        ),
        duplicateCheckMode = cms.untracked.string("noDuplicateCheck"),
                            #eventsToSkip = cms.untracked.VEventRange('1:233-1:235'),
                            )
                            
                                                              
process.me0SegAna = cms.EDAnalyzer('ME0SegmentAnalyzerMuonGun',
                                  me0SegmentInputTag = cms.InputTag("me0SegmentsSmeared","", "RECO"),
                                  me0DigiInputTag = cms.InputTag("simMuonME0ReDigis"),
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
fileName = cms.string("histoSingleMu_PU0_10k_NoNoise_RERECO_STDSEGMENT_REDIGI.root")
)

process.p = cms.Path(process.me0SegAna)


