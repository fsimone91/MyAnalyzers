import FWCore.ParameterSet.Config as cms

# Give the process a name
process = cms.Process("PickEvent")

# Tell the process which files to use as the source
process.source = cms.Source ("PoolSource",
          fileNames = cms.untracked.vstring(
        "file:./SingleMu_PU0_10k_NoNoise_RERECO_RUSEGMENT.root",
        
	),
          eventsToProcess = cms.untracked.VEventRange( '1:1:18',
                                                      '1:1:8',
                                                      '1:1:24',
                                                      '1:1:185',
                                                      '1:1:236',
                                                      '1:1:253',
                                                      '1:1:39',
                                                      '1:1:46',
                                                      '1:1:166',
                                                      '1:1:183',
                                                      '1:1:304',
                                                      '1:1:311',
                                                      '1:1:27',
                                                      '1:1:108',
                                                      '1:1:150',
                                                      '1:1:218',
                                                      '1:1:251',
                                                      '1:1:256',
                                                      '1:1:274',
                                                      '1:1:16',
                                                      '1:1:76',
                                                      '1:1:144',
                                                      '1:1:45',
                                                      '1:1:272',
                                                      '1:1:131',
                                                       '1:1:103',
                                                      '1:1:34',
                                                      '1:1:312',

                                                      
						)

)

# tell the process to only run over 100 events (-1 would mean run over
#  everything
process.maxEvents = cms.untracked.PSet(
            input = cms.untracked.int32(-1)

)

# Tell the process what filename to use to save the output
process.Out = cms.OutputModule("PoolOutputModule",
         fileName = cms.untracked.string("LostMuon_events.root")
)

# make sure everything is hooked up
process.end = cms.EndPath(process.Out)

