import FWCore.ParameterSet.Config as cms

from Configuration.StandardSequences.Eras import eras
process = cms.Process("Demo",eras.Phase2C1)
process.load("FWCore.MessageService.MessageLogger_cfi")
process.load('Configuration.StandardSequences.Services_cff')
process.load('Configuration.Geometry.GeometryExtended2023D1Reco_cff')
#process.load('Configuration.StandardSequences.MagneticField_38T_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
#process.load('RecoMuon.MuonIdentification.me0MuonReco_cff')



from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:run2_mc', '')


process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(1000) )

process.source = cms.Source("PoolSource",
                            # replace 'myfile.root',',',',',' with the source file you want to use
                            fileNames = cms.untracked.vstring(
                                                              'file:/lustre/cms/store/user/rosma/step3.root'
                                                              #'root://xrootd-cms.infn.it:1194//store/relval/CMSSW_8_1_0_pre11/RelValTenMuExtendedE_0_200/GEN-SIM-RECO/81X_mcRun2_asymptotic_v5_2023D1-v1/00000/044DA3BD-7770-E611-B8C0-0CC47A78A3D8.root',
                                                              #'root://xrootd-cms.infn.it:1194//store/relval/CMSSW_8_1_0_pre11/RelValTenMuExtendedE_0_200/GEN-SIM-RECO/81X_mcRun2_asymptotic_v5_2023D1-v1/00000/4624D6C7-5370-E611-B006-0025905A607E.root',
                                                              #'root://xrootd-cms.infn.it:1194//store/relval/CMSSW_8_1_0_pre11/RelValTenMuExtendedE_0_200/GEN-SIM-RECO/81X_mcRun2_asymptotic_v5_2023D1-v1/00000/B0FB198C-7870-E611-915A-0CC47A4C8E0E.root'
                                                              
                                                              
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
fileName = cms.string("histoME0SegMuGun_81X_1k.root")
)

process.p = cms.Path(process.me0SegAna)


