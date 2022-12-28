/*
 * PFCompleteFiller.h
 *
 *  Created on: Sep 25, 2017
 *      Author: hqu
 */

#ifndef NTUPLER_INTERFACE_PFCOMPLETEFILLER_H_
#define NTUPLER_INTERFACE_PFCOMPLETEFILLER_H_

#include <memory>

#include "DeepNTuples/BTagHelpers/interface/TrackInfoBuilder.h"
#include "DataFormats/Candidate/interface/VertexCompositePtrCandidate.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/Common/interface/ValueMap.h"

#include "DataFormats/PatCandidates/interface/PackedGenParticle.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

#include "DeepNTuples/NtupleCommons/interface/NtupleBase.h"

namespace deepntuples {

class PFCompleteFiller: public NtupleBase {
public:
  PFCompleteFiller() : PFCompleteFiller("", 0.8) {}
  PFCompleteFiller(std::string branchName, double jetR=0.8) : NtupleBase(branchName, jetR) {}
  virtual ~PFCompleteFiller() {}

  // get input parameters from the cfg file
  virtual void readConfig(const edm::ParameterSet& iConfig, edm::ConsumesCollector && cc) override;

  // read event content or event setup for each event
  virtual void readEvent(const edm::Event& iEvent, const edm::EventSetup& iSetup) override;

  // check if particle comes from a given parton
  bool containParton(const reco::Candidate * pruned_part, int pdgid);

protected:
  // declare the data branches (name, type, default values)
  virtual void book() override;
  // fill the branches
  virtual bool fill(const pat::Jet &jet, size_t jetidx, const JetHelper &jet_helper) override;

private:
  edm::EDGetTokenT<reco::VertexCollection> vtxToken_;
  edm::Handle<reco::VertexCollection> vertices;

  edm::EDGetTokenT<reco::VertexCompositePtrCandidateCollection> svToken_;
  edm::Handle<reco::VertexCompositePtrCandidateCollection> SVs;

  edm::EDGetTokenT<std::vector<pat::PackedGenParticle>> packedToken_;
  edm::Handle<std::vector<pat::PackedGenParticle>> packed;

  edm::EDGetTokenT<std::vector<reco::GenParticle>> prunedToken_;
  edm::Handle<std::vector<reco::GenParticle>> pruned;

  edm::ESHandle<TransientTrackBuilder> builder_;
};

} /* namespace deepntuples */

#endif /* NTUPLER_INTERFACE_PFCOMPLETEFILLER_H_ */
