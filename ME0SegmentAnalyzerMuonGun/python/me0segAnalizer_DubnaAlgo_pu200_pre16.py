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


process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(3000) )

process.source = cms.Source("PoolSource",
                            # replace 'myfile.root',',',',',',',' with the source file you want to use
                            fileNames = cms.untracked.vstring(
        #'root://cmsxrootd.fnal.gov//store/user/lpcgem/ME0TDRStudies/MuGun_0p5_30/90X_DIGI_PU0_WithoutNoise/job_1.root',
        'root://cmsxrootd.fnal.gov//store/user/lpcgem/ME0TDRStudies/MuGun_0p5_30/90X_RECO_PU200_WithNoise_BXFilter/job_0.root',
        'root://cmsxrootd.fnal.gov//store/user/lpcgem/ME0TDRStudies/MuGun_0p5_30/90X_RECO_PU200_WithNoise_BXFilter/job_1.root',
        'root://cmsxrootd.fnal.gov//store/user/lpcgem/ME0TDRStudies/MuGun_0p5_30/90X_RECO_PU200_WithNoise_BXFilter/job_2.root',
        'root://cmsxrootd.fnal.gov//store/user/lpcgem/ME0TDRStudies/MuGun_0p5_30/90X_RECO_PU200_WithNoise_BXFilter/job_3.root',
        'root://cmsxrootd.fnal.gov//store/user/lpcgem/ME0TDRStudies/MuGun_0p5_30/90X_RECO_PU200_WithNoise_BXFilter/job_4.root',
        'root://cmsxrootd.fnal.gov//store/user/lpcgem/ME0TDRStudies/MuGun_0p5_30/90X_RECO_PU200_WithNoise_BXFilter/job_5.root',
        'root://cmsxrootd.fnal.gov//store/user/lpcgem/ME0TDRStudies/MuGun_0p5_30/90X_RECO_PU200_WithNoise_BXFilter/job_6.root',
        'root://cmsxrootd.fnal.gov//store/user/lpcgem/ME0TDRStudies/MuGun_0p5_30/90X_RECO_PU200_WithNoise_BXFilter/job_7.root',
        'root://cmsxrootd.fnal.gov//store/user/lpcgem/ME0TDRStudies/MuGun_0p5_30/90X_RECO_PU200_WithNoise_BXFilter/job_8.root',
        'root://cmsxrootd.fnal.gov//store/user/lpcgem/ME0TDRStudies/MuGun_0p5_30/90X_RECO_PU200_WithNoise_BXFilter/job_9.root',
        'root://cmsxrootd.fnal.gov//store/user/lpcgem/ME0TDRStudies/MuGun_0p5_30/90X_RECO_PU200_WithNoise_BXFilter/job_10.root',
      'root://cmsxrootd.fnal.gov//store/user/lpcgem/ME0TDRStudies/MuGun_0p5_30/90X_RECO_PU200_WithNoise_BXFilter/job_11.root',
      'root://cmsxrootd.fnal.gov//store/user/lpcgem/ME0TDRStudies/MuGun_0p5_30/90X_RECO_PU200_WithNoise_BXFilter/job_12.root',
      'root://cmsxrootd.fnal.gov//store/user/lpcgem/ME0TDRStudies/MuGun_0p5_30/90X_RECO_PU200_WithNoise_BXFilter/job_13.root',
      'root://cmsxrootd.fnal.gov//store/user/lpcgem/ME0TDRStudies/MuGun_0p5_30/90X_RECO_PU200_WithNoise_BXFilter/job_14.root',
      'root://cmsxrootd.fnal.gov//store/user/lpcgem/ME0TDRStudies/MuGun_0p5_30/90X_RECO_PU200_WithNoise_BXFilter/job_15.root',
      'root://cmsxrootd.fnal.gov//store/user/lpcgem/ME0TDRStudies/MuGun_0p5_30/90X_RECO_PU200_WithNoise_BXFilter/job_16.root',
      'root://cmsxrootd.fnal.gov//store/user/lpcgem/ME0TDRStudies/MuGun_0p5_30/90X_RECO_PU200_WithNoise_BXFilter/job_17.root',
      'root://cmsxrootd.fnal.gov//store/user/lpcgem/ME0TDRStudies/MuGun_0p5_30/90X_RECO_PU200_WithNoise_BXFilter/job_18.root',
      'root://cmsxrootd.fnal.gov//store/user/lpcgem/ME0TDRStudies/MuGun_0p5_30/90X_RECO_PU200_WithNoise_BXFilter/job_19.root',
      'root://cmsxrootd.fnal.gov//store/user/lpcgem/ME0TDRStudies/MuGun_0p5_30/90X_RECO_PU200_WithNoise_BXFilter/job_20.root',
      'root://cmsxrootd.fnal.gov//store/user/lpcgem/ME0TDRStudies/MuGun_0p5_30/90X_RECO_PU200_WithNoise_BXFilter/job_21.root',
      'root://cmsxrootd.fnal.gov//store/user/lpcgem/ME0TDRStudies/MuGun_0p5_30/90X_RECO_PU200_WithNoise_BXFilter/job_22.root',
      'root://cmsxrootd.fnal.gov//store/user/lpcgem/ME0TDRStudies/MuGun_0p5_30/90X_RECO_PU200_WithNoise_BXFilter/job_23.root',
      'root://cmsxrootd.fnal.gov//store/user/lpcgem/ME0TDRStudies/MuGun_0p5_30/90X_RECO_PU200_WithNoise_BXFilter/job_24.root',
      'root://cmsxrootd.fnal.gov//store/user/lpcgem/ME0TDRStudies/MuGun_0p5_30/90X_RECO_PU200_WithNoise_BXFilter/job_25.root',
      'root://cmsxrootd.fnal.gov//store/user/lpcgem/ME0TDRStudies/MuGun_0p5_30/90X_RECO_PU200_WithNoise_BXFilter/job_26.root',
      'root://cmsxrootd.fnal.gov//store/user/lpcgem/ME0TDRStudies/MuGun_0p5_30/90X_RECO_PU200_WithNoise_BXFilter/job_27.root',
      'root://cmsxrootd.fnal.gov//store/user/lpcgem/ME0TDRStudies/MuGun_0p5_30/90X_RECO_PU200_WithNoise_BXFilter/job_28.root',
      'root://cmsxrootd.fnal.gov//store/user/lpcgem/ME0TDRStudies/MuGun_0p5_30/90X_RECO_PU200_WithNoise_BXFilter/job_29.root',
      'root://cmsxrootd.fnal.gov//store/user/lpcgem/ME0TDRStudies/MuGun_0p5_30/90X_RECO_PU200_WithNoise_BXFilter/job_30.root',
      'root://cmsxrootd.fnal.gov//store/user/lpcgem/ME0TDRStudies/MuGun_0p5_30/90X_RECO_PU200_WithNoise_BXFilter/job_31.root',
      'root://cmsxrootd.fnal.gov//store/user/lpcgem/ME0TDRStudies/MuGun_0p5_30/90X_RECO_PU200_WithNoise_BXFilter/job_32.root',
      'root://cmsxrootd.fnal.gov//store/user/lpcgem/ME0TDRStudies/MuGun_0p5_30/90X_RECO_PU200_WithNoise_BXFilter/job_33.root',
      'root://cmsxrootd.fnal.gov//store/user/lpcgem/ME0TDRStudies/MuGun_0p5_30/90X_RECO_PU200_WithNoise_BXFilter/job_34.root',
      'root://cmsxrootd.fnal.gov//store/user/lpcgem/ME0TDRStudies/MuGun_0p5_30/90X_RECO_PU200_WithNoise_BXFilter/job_35.root',
      'root://cmsxrootd.fnal.gov//store/user/lpcgem/ME0TDRStudies/MuGun_0p5_30/90X_RECO_PU200_WithNoise_BXFilter/job_36.root',
      'root://cmsxrootd.fnal.gov//store/user/lpcgem/ME0TDRStudies/MuGun_0p5_30/90X_RECO_PU200_WithNoise_BXFilter/job_37.root',
      'root://cmsxrootd.fnal.gov//store/user/lpcgem/ME0TDRStudies/MuGun_0p5_30/90X_RECO_PU200_WithNoise_BXFilter/job_38.root',
      'root://cmsxrootd.fnal.gov//store/user/lpcgem/ME0TDRStudies/MuGun_0p5_30/90X_RECO_PU200_WithNoise_BXFilter/job_39.root',
      'root://cmsxrootd.fnal.gov//store/user/lpcgem/ME0TDRStudies/MuGun_0p5_30/90X_RECO_PU200_WithNoise_BXFilter/job_40.root',
      'root://cmsxrootd.fnal.gov//store/user/lpcgem/ME0TDRStudies/MuGun_0p5_30/90X_RECO_PU200_WithNoise_BXFilter/job_41.root',
      'root://cmsxrootd.fnal.gov//store/user/lpcgem/ME0TDRStudies/MuGun_0p5_30/90X_RECO_PU200_WithNoise_BXFilter/job_42.root',
      'root://cmsxrootd.fnal.gov//store/user/lpcgem/ME0TDRStudies/MuGun_0p5_30/90X_RECO_PU200_WithNoise_BXFilter/job_43.root',
      'root://cmsxrootd.fnal.gov//store/user/lpcgem/ME0TDRStudies/MuGun_0p5_30/90X_RECO_PU200_WithNoise_BXFilter/job_44.root',
      'root://cmsxrootd.fnal.gov//store/user/lpcgem/ME0TDRStudies/MuGun_0p5_30/90X_RECO_PU200_WithNoise_BXFilter/job_45.root',
      'root://cmsxrootd.fnal.gov//store/user/lpcgem/ME0TDRStudies/MuGun_0p5_30/90X_RECO_PU200_WithNoise_BXFilter/job_46.root',
      'root://cmsxrootd.fnal.gov//store/user/lpcgem/ME0TDRStudies/MuGun_0p5_30/90X_RECO_PU200_WithNoise_BXFilter/job_47.root',
      'root://cmsxrootd.fnal.gov//store/user/lpcgem/ME0TDRStudies/MuGun_0p5_30/90X_RECO_PU200_WithNoise_BXFilter/job_48.root',
      'root://cmsxrootd.fnal.gov//store/user/lpcgem/ME0TDRStudies/MuGun_0p5_30/90X_RECO_PU200_WithNoise_BXFilter/job_49.root',
      'root://cmsxrootd.fnal.gov//store/user/lpcgem/ME0TDRStudies/MuGun_0p5_30/90X_RECO_PU200_WithNoise_BXFilter/job_50.root',
      'root://cmsxrootd.fnal.gov//store/user/lpcgem/ME0TDRStudies/MuGun_0p5_30/90X_RECO_PU200_WithNoise_BXFilter/job_51.root',
      'root://cmsxrootd.fnal.gov//store/user/lpcgem/ME0TDRStudies/MuGun_0p5_30/90X_RECO_PU200_WithNoise_BXFilter/job_52.root',
      'root://cmsxrootd.fnal.gov//store/user/lpcgem/ME0TDRStudies/MuGun_0p5_30/90X_RECO_PU200_WithNoise_BXFilter/job_53.root',
      'root://cmsxrootd.fnal.gov//store/user/lpcgem/ME0TDRStudies/MuGun_0p5_30/90X_RECO_PU200_WithNoise_BXFilter/job_54.root',
      'root://cmsxrootd.fnal.gov//store/user/lpcgem/ME0TDRStudies/MuGun_0p5_30/90X_RECO_PU200_WithNoise_BXFilter/job_55.root',
      'root://cmsxrootd.fnal.gov//store/user/lpcgem/ME0TDRStudies/MuGun_0p5_30/90X_RECO_PU200_WithNoise_BXFilter/job_56.root',
      'root://cmsxrootd.fnal.gov//store/user/lpcgem/ME0TDRStudies/MuGun_0p5_30/90X_RECO_PU200_WithNoise_BXFilter/job_57.root',
      'root://cmsxrootd.fnal.gov//store/user/lpcgem/ME0TDRStudies/MuGun_0p5_30/90X_RECO_PU200_WithNoise_BXFilter/job_58.root',
      'root://cmsxrootd.fnal.gov//store/user/lpcgem/ME0TDRStudies/MuGun_0p5_30/90X_RECO_PU200_WithNoise_BXFilter/job_59.root',
      'root://cmsxrootd.fnal.gov//store/user/lpcgem/ME0TDRStudies/MuGun_0p5_30/90X_RECO_PU200_WithNoise_BXFilter/job_60.root',
      'root://cmsxrootd.fnal.gov//store/user/lpcgem/ME0TDRStudies/MuGun_0p5_30/90X_RECO_PU200_WithNoise_BXFilter/job_61.root',
      'root://cmsxrootd.fnal.gov//store/user/lpcgem/ME0TDRStudies/MuGun_0p5_30/90X_RECO_PU200_WithNoise_BXFilter/job_62.root',
      'root://cmsxrootd.fnal.gov//store/user/lpcgem/ME0TDRStudies/MuGun_0p5_30/90X_RECO_PU200_WithNoise_BXFilter/job_63.root',
      'root://cmsxrootd.fnal.gov//store/user/lpcgem/ME0TDRStudies/MuGun_0p5_30/90X_RECO_PU200_WithNoise_BXFilter/job_64.root',
      'root://cmsxrootd.fnal.gov//store/user/lpcgem/ME0TDRStudies/MuGun_0p5_30/90X_RECO_PU200_WithNoise_BXFilter/job_65.root',
      'root://cmsxrootd.fnal.gov//store/user/lpcgem/ME0TDRStudies/MuGun_0p5_30/90X_RECO_PU200_WithNoise_BXFilter/job_66.root',
      'root://cmsxrootd.fnal.gov//store/user/lpcgem/ME0TDRStudies/MuGun_0p5_30/90X_RECO_PU200_WithNoise_BXFilter/job_67.root',
      'root://cmsxrootd.fnal.gov//store/user/lpcgem/ME0TDRStudies/MuGun_0p5_30/90X_RECO_PU200_WithNoise_BXFilter/job_68.root',
      'root://cmsxrootd.fnal.gov//store/user/lpcgem/ME0TDRStudies/MuGun_0p5_30/90X_RECO_PU200_WithNoise_BXFilter/job_69.root',
      'root://cmsxrootd.fnal.gov//store/user/lpcgem/ME0TDRStudies/MuGun_0p5_30/90X_RECO_PU200_WithNoise_BXFilter/job_70.root',
      'root://cmsxrootd.fnal.gov//store/user/lpcgem/ME0TDRStudies/MuGun_0p5_30/90X_RECO_PU200_WithNoise_BXFilter/job_71.root',
      'root://cmsxrootd.fnal.gov//store/user/lpcgem/ME0TDRStudies/MuGun_0p5_30/90X_RECO_PU200_WithNoise_BXFilter/job_72.root',
      'root://cmsxrootd.fnal.gov//store/user/lpcgem/ME0TDRStudies/MuGun_0p5_30/90X_RECO_PU200_WithNoise_BXFilter/job_73.root',
      'root://cmsxrootd.fnal.gov//store/user/lpcgem/ME0TDRStudies/MuGun_0p5_30/90X_RECO_PU200_WithNoise_BXFilter/job_74.root',
      'root://cmsxrootd.fnal.gov//store/user/lpcgem/ME0TDRStudies/MuGun_0p5_30/90X_RECO_PU200_WithNoise_BXFilter/job_75.root',
      'root://cmsxrootd.fnal.gov//store/user/lpcgem/ME0TDRStudies/MuGun_0p5_30/90X_RECO_PU200_WithNoise_BXFilter/job_76.root',
      'root://cmsxrootd.fnal.gov//store/user/lpcgem/ME0TDRStudies/MuGun_0p5_30/90X_RECO_PU200_WithNoise_BXFilter/job_77.root',
      'root://cmsxrootd.fnal.gov//store/user/lpcgem/ME0TDRStudies/MuGun_0p5_30/90X_RECO_PU200_WithNoise_BXFilter/job_78.root',
      'root://cmsxrootd.fnal.gov//store/user/lpcgem/ME0TDRStudies/MuGun_0p5_30/90X_RECO_PU200_WithNoise_BXFilter/job_79.root',
      'root://cmsxrootd.fnal.gov//store/user/lpcgem/ME0TDRStudies/MuGun_0p5_30/90X_RECO_PU200_WithNoise_BXFilter/job_80.root',
      'root://cmsxrootd.fnal.gov//store/user/lpcgem/ME0TDRStudies/MuGun_0p5_30/90X_RECO_PU200_WithNoise_BXFilter/job_81.root',
      'root://cmsxrootd.fnal.gov//store/user/lpcgem/ME0TDRStudies/MuGun_0p5_30/90X_RECO_PU200_WithNoise_BXFilter/job_82.root',
      'root://cmsxrootd.fnal.gov//store/user/lpcgem/ME0TDRStudies/MuGun_0p5_30/90X_RECO_PU200_WithNoise_BXFilter/job_83.root',
      'root://cmsxrootd.fnal.gov//store/user/lpcgem/ME0TDRStudies/MuGun_0p5_30/90X_RECO_PU200_WithNoise_BXFilter/job_84.root',
      'root://cmsxrootd.fnal.gov//store/user/lpcgem/ME0TDRStudies/MuGun_0p5_30/90X_RECO_PU200_WithNoise_BXFilter/job_85.root',
      'root://cmsxrootd.fnal.gov//store/user/lpcgem/ME0TDRStudies/MuGun_0p5_30/90X_RECO_PU200_WithNoise_BXFilter/job_86.root',
      'root://cmsxrootd.fnal.gov//store/user/lpcgem/ME0TDRStudies/MuGun_0p5_30/90X_RECO_PU200_WithNoise_BXFilter/job_87.root',
      'root://cmsxrootd.fnal.gov//store/user/lpcgem/ME0TDRStudies/MuGun_0p5_30/90X_RECO_PU200_WithNoise_BXFilter/job_88.root',
      'root://cmsxrootd.fnal.gov//store/user/lpcgem/ME0TDRStudies/MuGun_0p5_30/90X_RECO_PU200_WithNoise_BXFilter/job_89.root',
      'root://cmsxrootd.fnal.gov//store/user/lpcgem/ME0TDRStudies/MuGun_0p5_30/90X_RECO_PU200_WithNoise_BXFilter/job_90.root',
      'root://cmsxrootd.fnal.gov//store/user/lpcgem/ME0TDRStudies/MuGun_0p5_30/90X_RECO_PU200_WithNoise_BXFilter/job_91.root',
      'root://cmsxrootd.fnal.gov//store/user/lpcgem/ME0TDRStudies/MuGun_0p5_30/90X_RECO_PU200_WithNoise_BXFilter/job_92.root',
      'root://cmsxrootd.fnal.gov//store/user/lpcgem/ME0TDRStudies/MuGun_0p5_30/90X_RECO_PU200_WithNoise_BXFilter/job_93.root',
      'root://cmsxrootd.fnal.gov//store/user/lpcgem/ME0TDRStudies/MuGun_0p5_30/90X_RECO_PU200_WithNoise_BXFilter/job_94.root',
      'root://cmsxrootd.fnal.gov//store/user/lpcgem/ME0TDRStudies/MuGun_0p5_30/90X_RECO_PU200_WithNoise_BXFilter/job_95.root',
      'root://cmsxrootd.fnal.gov//store/user/lpcgem/ME0TDRStudies/MuGun_0p5_30/90X_RECO_PU200_WithNoise_BXFilter/job_96.root',
      'root://cmsxrootd.fnal.gov//store/user/lpcgem/ME0TDRStudies/MuGun_0p5_30/90X_RECO_PU200_WithNoise_BXFilter/job_97.root',
      'root://cmsxrootd.fnal.gov//store/user/lpcgem/ME0TDRStudies/MuGun_0p5_30/90X_RECO_PU200_WithNoise_BXFilter/job_98.root',
      'root://cmsxrootd.fnal.gov//store/user/lpcgem/ME0TDRStudies/MuGun_0p5_30/90X_RECO_PU200_WithNoise_BXFilter/job_99.root',
      'root://cmsxrootd.fnal.gov//store/user/lpcgem/ME0TDRStudies/MuGun_0p5_30/90X_RECO_PU200_WithNoise_BXFilter/job_100.root',
                                                 
                   
        ),
                            duplicateCheckMode = cms.untracked.string("noDuplicateCheck"),
                            # eventsToSkip = cms.untracked.VEventRange('1:233-1:235'),
)
                            

process.me0SegAna = cms.EDAnalyzer('ME0SegmentAnalyzerMuonGun',
                                  me0SegmentInputTag = cms.InputTag("me0SegmentsPerfectRU","", "RECO"),
                                  me0DigiInputTag = cms.InputTag("simMuonME0Digis"),
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
fileName = cms.string("histoSingleMu_PU200_10k_RERECO_RUSEGMENT.root")
)

process.p = cms.Path(process.me0SegAna)


