import FWCore.ParameterSet.Config as cms

# ---------------------------------------------------------
from FWCore.ParameterSet.VarParsing import VarParsing
options = VarParsing('analysis')


#options.outputFile = '/eos/home-m/mmalucch/dnntuple_output/output_tot_cutidx_cutpca_pt3.root'
options.outputFile = '/scratchnvme/malucchi/output_tot_cutidx_cutpca_pt4.root'


#options.inputFiles = '/store/mc/RunIISummer19UL17MiniAOD/BulkGravitonToHHTo4Q_MX-600to6000_MH-15to250_part2_TuneCP5_13TeV-madgraph_pythia8/MINIAODSIM/multigridpack_106X_mc2017_realistic_v6-v1/50000/FB46C2C2-73A4-A64C-A3D7-FC47C6A48871.root'
# options.inputFiles = '/store/mc/RunIISummer19UL17MiniAOD/QCD_Pt-15to7000_TuneCP5_Flat_13TeV_pythia8/MINIAODSIM/106X_mc2017_realistic_v6_ext2-v2/60000/E0502E95-AA7E-B441-9277-498113BA458C.root'
# options.inputFiles = '/store/mc/RunIISummer19UL17MiniAOD/TTToSemiLeptonic_mtop171p5_TuneCP5_13TeV-powheg-pythia8/MINIAODSIM/106X_mc2017_realistic_v6-v2/60000/F06BD23E-48AB-1F4D-A0FB-80DD370F4868.root'

'''options.inputFiles = ['/store/mc/RunIISummer20UL18MiniAODv2/TTJets_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/106X_upgrade2018_realistic_v16_L1v1-v2/2430000/0984A792-8B13-1543-AA86-063CC14B1678.root',
                      '/store/mc/RunIISummer20UL18MiniAODv2/TTJets_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/106X_upgrade2018_realistic_v16_L1v1-v2/2430000/1069D98E-0622-5642-BE60-EC3C8C9CF87D.root',
                      '/store/mc/RunIISummer20UL18MiniAODv2/TTJets_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/106X_upgrade2018_realistic_v16_L1v1-v2/2430000/125C21CC-F228-6547-9A38-FC135A8D9764.root',
                      '/store/mc/RunIISummer20UL18MiniAODv2/TTJets_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/106X_upgrade2018_realistic_v16_L1v1-v2/2430000/195947E9-1F3A-D84A-8ED8-34D7F38467D4.root',
                      '/store/mc/RunIISummer20UL18MiniAODv2/TTJets_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/106X_upgrade2018_realistic_v16_L1v1-v2/2430000/1C772400-8C89-CB4E-93DC-463AE40E4EC1.root']

options.inputFiles = ['/store/mc/RunIISummer20UL18MiniAODv2/TTJets_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/106X_upgrade2018_realistic_v16_L1v1-v2/2430000/1CEEDB6C-987C-FE49-B96E-13FC35BF1558.root',
                      '/store/mc/RunIISummer20UL18MiniAODv2/TTJets_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/106X_upgrade2018_realistic_v16_L1v1-v2/2430000/1FA712C2-D07D-5447-951D-D424992820C4.root',
                      '/store/mc/RunIISummer20UL18MiniAODv2/TTJets_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/106X_upgrade2018_realistic_v16_L1v1-v2/2430000/29667DEB-8CD9-814F-B8EB-1C7A049A0E7A.root',
                      '/store/mc/RunIISummer20UL18MiniAODv2/TTJets_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/106X_upgrade2018_realistic_v16_L1v1-v2/2430000/2ACB3396-DA76-444C-B65F-F79BE6732F58.root',
                      '/store/mc/RunIISummer20UL18MiniAODv2/TTJets_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/106X_upgrade2018_realistic_v16_L1v1-v2/2430000/359158CA-4820-354F-BEF3-F8B1EA63A1B4.root',
                      '/store/mc/RunIISummer20UL18MiniAODv2/TTJets_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/106X_upgrade2018_realistic_v16_L1v1-v2/2430000/3638C320-DA2B-A74F-A417-981BB01E852E.root',
                      '/store/mc/RunIISummer20UL18MiniAODv2/TTJets_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/106X_upgrade2018_realistic_v16_L1v1-v2/2430000/36B941BB-9B5A-EA43-A8A2-77E3858B1B13.root']
'''

"""options.inputFiles = ['/store/mc/RunIISummer20UL18MiniAODv2/TTJets_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/106X_upgrade2018_realistic_v16_L1v1-v2/2430000/49973BC0-3D16-EB4E-81E1-5CE4B1FFFA0F.root',
                    '/store/mc/RunIISummer20UL18MiniAODv2/TTJets_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/106X_upgrade2018_realistic_v16_L1v1-v2/2430000/4A170998-CF71-1941-AFF0-8DD32AED6B46.root',
                    '/store/mc/RunIISummer20UL18MiniAODv2/TTJets_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/106X_upgrade2018_realistic_v16_L1v1-v2/2430000/506DCAD2-574B-F24A-80A1-63E96711850E.root',
                    '/store/mc/RunIISummer20UL18MiniAODv2/TTJets_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/106X_upgrade2018_realistic_v16_L1v1-v2/2430000/5E52CEEA-3D4B-2342-8465-592D2BC0AD16.root',
                    '/store/mc/RunIISummer20UL18MiniAODv2/TTJets_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/106X_upgrade2018_realistic_v16_L1v1-v2/2430000/6930178C-2478-F743-B538-DD1BCF14592D.root',
                    '/store/mc/RunIISummer20UL18MiniAODv2/TTJets_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/106X_upgrade2018_realistic_v16_L1v1-v2/2430000/6B153BB1-B12D-D741-9D2B-6B48A8D69467.root',
                    '/store/mc/RunIISummer20UL18MiniAODv2/TTJets_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/106X_upgrade2018_realistic_v16_L1v1-v2/2430000/7FB4B911-32D3-6648-9B0E-B112B092492C.root',
                    '/store/mc/RunIISummer20UL18MiniAODv2/TTJets_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/106X_upgrade2018_realistic_v16_L1v1-v2/2430000/9739BF6A-D21B-A546-B8E9-216458C68831.root',
                    '/store/mc/RunIISummer20UL18MiniAODv2/TTJets_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/106X_upgrade2018_realistic_v16_L1v1-v2/2430000/9CB6344C-40BA-3146-B75A-2124D6C039BE.root',
                    '/store/mc/RunIISummer20UL18MiniAODv2/TTJets_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/106X_upgrade2018_realistic_v16_L1v1-v2/2430000/9D30A9AA-CF6E-EA45-81AF-F4571C09D431.root',
                    '/store/mc/RunIISummer20UL18MiniAODv2/TTJets_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/106X_upgrade2018_realistic_v16_L1v1-v2/2430000/9F4DD37F-2F34-3A49-B552-F0100162A7C2.root',
                    '/store/mc/RunIISummer20UL18MiniAODv2/TTJets_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/106X_upgrade2018_realistic_v16_L1v1-v2/2430000/AFF1A7DD-6B79-9A45-956C-D3385BC6EB63.root',
                    '/store/mc/RunIISummer20UL18MiniAODv2/TTJets_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/106X_upgrade2018_realistic_v16_L1v1-v2/2430000/B0FFB691-1562-224F-B8E8-1F98C8B5C78D.root',
                    '/store/mc/RunIISummer20UL18MiniAODv2/TTJets_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/106X_upgrade2018_realistic_v16_L1v1-v2/2430000/CF3A46EE-61A2-454B-A4CE-350A8920AAC5.root',
                    '/store/mc/RunIISummer20UL18MiniAODv2/TTJets_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/106X_upgrade2018_realistic_v16_L1v1-v2/2430000/F68CE2D2-F078-DB42-BC30-4D73566DB39D.root',
                    '/store/mc/RunIISummer20UL18MiniAODv2/TTJets_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/106X_upgrade2018_realistic_v16_L1v1-v2/2430000/FBB49546-F1ED-EA4B-8F12-0466361B8404.root',
]"""

with open("files.txt", 'r') as f:
    options.inputFiles = [line.strip() for line in f][-70:-60]
    print(options.inputFiles)

options.maxEvents = -1

options.register('skipEvents', 0, VarParsing.multiplicity.singleton, VarParsing.varType.int, "skip N events")
options.register('inputDataset',
                 '',
                 VarParsing.multiplicity.singleton,
                 VarParsing.varType.string,
                 "Input dataset")
options.register('isTrainSample', True, VarParsing.multiplicity.singleton, VarParsing.varType.bool, "if the sample is used for training")

options.parseArguments()

globalTagMap = {
    'Summer19UL17': '106X_mc2017_realistic_v6',
    'Summer19UL18': '106X_upgrade2018_realistic_v11_L1v1',
    'Summer19UL16': '',
}

era = None if options.inputDataset else 'Summer19UL17'
for k in globalTagMap:
    if k in options.inputDataset:
        era = k
# ---------------------------------------------------------
process = cms.Process("DNNFiller")

process.options.numberOfThreads=cms.untracked.uint32(32)
process.options.numberOfStreams=cms.untracked.uint32(32)

process.load('FWCore.MessageService.MessageLogger_cfi')
process.MessageLogger.cerr.FwkReport.reportEvery = 1000

process.options = cms.untracked.PSet(
    allowUnscheduled = cms.untracked.bool(True),
    wantSummary=cms.untracked.bool(False)
)

print ('Using output file ' + options.outputFile)

process.TFileService = cms.Service("TFileService",
                                   fileName=cms.string(options.outputFile))

process.maxEvents = cms.untracked.PSet(input=cms.untracked.int32(options.maxEvents))

process.source = cms.Source('PoolSource',
    fileNames=cms.untracked.vstring(options.inputFiles),
    skipEvents=cms.untracked.uint32(options.skipEvents)
)
# ---------------------------------------------------------
process.load("FWCore.MessageService.MessageLogger_cfi")
process.load("Configuration.EventContent.EventContent_cff")
process.load('Configuration.StandardSequences.Services_cff')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, globalTagMap[era], '')
print('Using global tag', process.GlobalTag.globaltag)
process.TransientTrackBuilderESProducer = cms.ESProducer("TransientTrackBuilderESProducer",
    ComponentName=cms.string('TransientTrackBuilder')
)
# ---------------------------------------------------------
# read JEC from sqlite
if era == 'Summer19UL17':
    import os
    jecTag = 'Summer19UL17_V5_MC'
    jecFile = '%s.db' % jecTag
    if not os.path.exists(jecFile):
        os.symlink('../data/'+jecFile, jecFile)
    from CondCore.CondDB.CondDB_cfi import CondDB
    CondDBJECFile = CondDB.clone(connect = cms.string( 'sqlite:%s'%jecFile ) )
    process.jec = cms.ESSource('PoolDBESSource',
        CondDBJECFile,
        toGet = cms.VPSet(
            cms.PSet(
                record = cms.string('JetCorrectionsRecord'),
                tag    = cms.string('JetCorrectorParametersCollection_%s_AK4PFchs' % jecTag),
                label  = cms.untracked.string('AK4PFchs')
            ),
            cms.PSet(
                record = cms.string('JetCorrectionsRecord'),
                tag    = cms.string('JetCorrectorParametersCollection_%s_AK4PFPuppi' % jecTag),
                label  = cms.untracked.string('AK4PFPuppi')
            ),
            # ...and so on for all jet types you need
        )
    )
    print(jecTag, process.jec.toGet)
    # Add an ESPrefer to override JEC that might be available from the global tag
    process.es_prefer_jec = cms.ESPrefer('PoolDBESSource', 'jec')
# ---------------------------------------------------------
from PhysicsTools.PatAlgos.tools.jetTools import updateJetCollection

isPuppiJets = True
jetR = 0.4

bTagDiscriminators = [
    'pfDeepFlavourJetTags:probb',
    'pfDeepFlavourJetTags:probbb',
    'pfDeepFlavourJetTags:problepb',
    'pfDeepFlavourJetTags:probc',
    'pfDeepFlavourJetTags:probuds',
    'pfDeepFlavourJetTags:probg',
]

JETCorrLevels = ['L2Relative', 'L3Absolute']

from DeepNTuples.Ntupler.jetToolbox_cff import jetToolbox
jetToolbox(process, 'ak4', 'dummySeq', 'noOutput',
           PUMethod='Puppi', JETCorrPayload='AK4PFPuppi', JETCorrLevels=JETCorrLevels,
           Cut='pt > 10',
           runOnMC=True,
           bTagDiscriminators=['None'], subjetBTagDiscriminators=['None'])

updateJetCollection(
    process,
    jetSource=cms.InputTag('selectedPatJetsAK4PFPuppi'),
    jetCorrections=('AK4PFPuppi', cms.vstring(JETCorrLevels), 'None'),
    btagDiscriminators=bTagDiscriminators,
    postfix='AK4Puppi',
)
srcJets = cms.InputTag('selectedUpdatedPatJetsAK4Puppi')
# ---------------------------------------------------------
from RecoJets.JetProducers.QGTagger_cfi import QGTagger
process.qgtagger = QGTagger.clone(srcJets=srcJets, srcVertexCollection="offlineSlimmedPrimaryVertices")
process.updatedJetsAK4PuppiWithUserData = cms.EDProducer("PATJetUserDataEmbedder",
    src=srcJets,
    userFloats=cms.PSet(qgl=cms.InputTag('qgtagger:qgLikelihood')),
    )
srcJets = cms.InputTag('updatedJetsAK4PuppiWithUserData')
process.qgTask = cms.Task(
    process.qgtagger,
    process.updatedJetsAK4PuppiWithUserData,
)
# ---------------------------------------------------------
from PhysicsTools.PatAlgos.tools.helpers import getPatAlgosToolsTask
patTask = getPatAlgosToolsTask(process)

from RecoJets.JetProducers.ak4GenJets_cfi import ak4GenJets
process.ak4GenJetsWithNu = ak4GenJets.clone(
    src='packedGenParticles'
    )
process.ak4GenJetsWithNuMatch = cms.EDProducer("GenJetMatcher",  # cut on deltaR; pick best by deltaR
    src=srcJets,  # RECO jets (any View<Jet> is ok)
    matched=cms.InputTag("ak4GenJetsWithNu"),  # GEN jets  (must be GenJetCollection)
    mcPdgId=cms.vint32(),  # n/a
    mcStatus=cms.vint32(),  # n/a
    checkCharge=cms.bool(False),  # n/a
    maxDeltaR=cms.double(jetR),  # Minimum deltaR for the match
    # maxDPtRel   = cms.double(3.0),                  # Minimum deltaPt/Pt for the match (not used in GenJetMatcher)
    resolveAmbiguities=cms.bool(True),  # Forbid two RECO objects to match to the same GEN object
    resolveByMatchQuality=cms.bool(False),  # False = just match input in order; True = pick lowest deltaR pair first
)
process.genJetTask = cms.Task(
    process.ak4GenJetsWithNu,
    process.ak4GenJetsWithNuMatch,
)

# DeepNtuplizer
process.load("DeepNTuples.Ntupler.DeepNtuplizer_cfi")
process.deepntuplizer.jets = srcJets
process.deepntuplizer.isPuppiJets = isPuppiJets
process.deepntuplizer.bDiscriminators = bTagDiscriminators

process.deepntuplizer.genJetsMatch = 'ak4GenJetsWithNuMatch'

process.deepntuplizer.isQCDSample = '/QCD_' in options.inputDataset
process.deepntuplizer.isPythia = 'pythia' in options.inputDataset.lower()
process.deepntuplizer.isHerwig = 'herwig' in options.inputDataset.lower()
process.deepntuplizer.isMadGraph = 'madgraph' in options.inputDataset.lower()  # note: MG can be interfaced w/ either pythia or herwig

process.deepntuplizer.isTrainSample = options.isTrainSample
if not options.inputDataset:
    # interactive running
    process.deepntuplizer.isTrainSample = False
#==============================================================================================================================#
process.p = cms.Path(process.deepntuplizer)
process.p.associate(patTask)
process.p.associate(process.genJetTask)
process.p.associate(process.qgTask)
