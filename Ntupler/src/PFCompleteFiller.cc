/*
 * PFCompleteFiller.cc
 *
 *  Created on: Sep 25, 2017
 *      Author: hqu
 */

#include <unordered_map>
#include "TrackingTools/Records/interface/TransientTrackRecord.h"

#include "DeepNTuples/Ntupler/interface/PFCompleteFiller.h"

namespace deepntuples {

void PFCompleteFiller::readConfig(const edm::ParameterSet& iConfig, edm::ConsumesCollector&& cc) {
  vtxToken_ = cc.consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("vertices"));
  svToken_ = cc.consumes<reco::VertexCompositePtrCandidateCollection>(iConfig.getParameter<edm::InputTag>("SVs"));
  packedToken_ = cc.consumes<std::vector<pat::PackedGenParticle>>(iConfig.getParameter<edm::InputTag>("packed"));
  prunedToken_ = cc.consumes<std::vector<reco::GenParticle>>(iConfig.getParameter<edm::InputTag>("pruned"));
}

void PFCompleteFiller::readEvent(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
  iEvent.getByToken(vtxToken_, vertices);
  iEvent.getByToken(svToken_, SVs);
  iEvent.getByToken(packedToken_, packed);
  iEvent.getByToken(prunedToken_, pruned);
  iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder", builder_);
}

bool PFCompleteFiller::containParton(const reco::Candidate * pruned_part, int pdgid) {
  if (abs(pruned_part->pdgId())==pdgid) return true;
  for(size_t i=0;i< pruned_part->numberOfMothers();i++){
    if (containParton(pruned_part->mother(i), pdgid)) return true;
  }
  return false;
}


void PFCompleteFiller::book() {

  data.add<int>("n_pfcands", 0);
  data.add<float>("npfcands", 0);

  // no puppi scaled
  data.addMulti<float>("pfcand_pt_nopuppi");
  data.addMulti<float>("pfcand_pt_log_nopuppi");
  data.addMulti<float>("pfcand_e_log_nopuppi");

  data.addMulti<float>("pfcand_phirel");
  data.addMulti<float>("pfcand_etarel");
  // data.addMulti<float>("pfcand_deltaR");
  data.addMulti<float>("pfcand_puppiw");
  data.addMulti<float>("pfcand_abseta");

  data.addMulti<float>("pfcand_drminsvin"); // restricted to within the jet cone

  data.addMulti<float>("pfcand_charge");
  data.addMulti<float>("pfcand_isMu");
  data.addMulti<float>("pfcand_isEl");
  data.addMulti<float>("pfcand_isChargedHad");
  data.addMulti<float>("pfcand_isGamma");
  data.addMulti<float>("pfcand_isNeutralHad");

  // for neutral
  data.addMulti<float>("pfcand_hcalFrac");
  data.addMulti<float>("pfcand_hcalFracCalib");

  // for charged
  data.addMulti<float>("pfcand_VTX_ass");
  data.addMulti<float>("pfcand_fromPV");
  data.addMulti<float>("pfcand_nValidHits");
  data.addMulti<float>("pfcand_nValidPixelHits");
  data.addMulti<float>("pfcand_lostInnerHits");
  data.addMulti<float>("pfcand_trackHighPurity");

  // impact parameters
  data.addMulti<float>("pfcand_dz");
  data.addMulti<float>("pfcand_dzsig");
  data.addMulti<float>("pfcand_dxy");
  data.addMulti<float>("pfcand_dxysig");

  // track quality
  data.addMulti<float>("pfcand_normchi2");
  data.addMulti<float>("pfcand_quality");

  // track covariance
  data.addMulti<float>("pfcand_dptdpt");
  data.addMulti<float>("pfcand_detadeta");
  data.addMulti<float>("pfcand_dphidphi");
  data.addMulti<float>("pfcand_dxydxy");
  data.addMulti<float>("pfcand_dzdz");
  data.addMulti<float>("pfcand_dxydz");
  data.addMulti<float>("pfcand_dphidxy");
  data.addMulti<float>("pfcand_dlambdadz");

  // track btag info
  // data.addMulti<float>("pfcand_btagMomentum");
  // data.addMulti<float>("pfcand_btagEta");
  data.addMulti<float>("pfcand_btagEtaRel");
  data.addMulti<float>("pfcand_btagPtRel");
  // data.addMulti<float>("pfcand_btagPPar");
  // data.addMulti<float>("pfcand_btagDeltaR");
  data.addMulti<float>("pfcand_btagPtRatio");
  data.addMulti<float>("pfcand_btagPParRatio");
  data.addMulti<float>("pfcand_btagSip2dVal");
  data.addMulti<float>("pfcand_btagSip2dSig");
  data.addMulti<float>("pfcand_btagSip3dVal");
  data.addMulti<float>("pfcand_btagSip3dSig");
  data.addMulti<float>("pfcand_btagJetDistVal");
  data.addMulti<float>("pfcand_btagDecayLengthVal");
  data.addMulti<float>("pfcand_btagDecayLengthSig");

  data.addMulti<int>("pfcand_from_b");
  data.addMulti<int>("pfcand_from_c");
  data.addMulti<int>("pfcand_from_g");

  data.addMulti<float>("pfcand_vtx_x");
  data.addMulti<float>("pfcand_vtx_y");
  data.addMulti<float>("pfcand_vtx_z");

  data.addMulti<float>("pfcand_dist_from_pv");

}

bool PFCompleteFiller::fill(const pat::Jet& jet, size_t jetidx, const JetHelper& jet_helper) {

  const auto& pfCands = jet_helper.getJetConstituents();

  data.fill<int>("n_pfcands", pfCands.size());
  data.fill<float>("npfcands", pfCands.size());

  float etasign = jet.eta()>0 ? 1 : -1;

  //float pv_x=-1,pv_y=-1,pv_z=-1;
  for(const auto &pruned_part : *pruned){
    if(pruned_part.pdgId()!=2212) {
      const auto pv = pruned_part.vertex();
      //pv_x= pv.x();
      //pv_y= pv.y();
      //pv_z= pv.z();

      for (const auto& cand : pfCands){

        const auto *packed_cand = dynamic_cast<const pat::PackedCandidate *>(&(*cand));

        // basic kinematics, valid for both charged and neutral
        // not puppi weighted
        data.fillMulti<float>("pfcand_pt_nopuppi", packed_cand->pt());
        data.fillMulti<float>("pfcand_pt_log_nopuppi", catchInfs(std::log(packed_cand->pt()), -99));
        data.fillMulti<float>("pfcand_e_log_nopuppi", catchInfs(std::log(packed_cand->energy()), -99));

        data.fillMulti<float>("pfcand_phirel", reco::deltaPhi(*packed_cand, jet));
        data.fillMulti<float>("pfcand_etarel", etasign * (packed_cand->eta() - jet.eta()));
        // data.fillMulti<float>("pfcand_deltaR", reco::deltaR(*packed_cand, jet));
        data.fillMulti<float>("pfcand_abseta", std::abs(packed_cand->eta()));

        data.fillMulti<float>("pfcand_puppiw", jet_helper.getPuppiWeight(cand));

        double minDRin = 2.*jetR_;
        for (const auto &sv : *SVs){
          double dr = reco::deltaR(*packed_cand, sv);
          if (dr < minDRin && reco::deltaR(jet, sv) < jetR_) minDRin = dr;
        }
        data.fillMulti<float>("pfcand_drminsvin", minDRin);

        data.fillMulti<float>("pfcand_charge", packed_cand->charge());
        data.fillMulti<float>("pfcand_isEl", std::abs(packed_cand->pdgId())==11);
        data.fillMulti<float>("pfcand_isMu", std::abs(packed_cand->pdgId())==13);
        data.fillMulti<float>("pfcand_isChargedHad", std::abs(packed_cand->pdgId())==211);
        data.fillMulti<float>("pfcand_isGamma", std::abs(packed_cand->pdgId())==22);
        data.fillMulti<float>("pfcand_isNeutralHad", std::abs(packed_cand->pdgId())==130);

        // for neutral
        float hcal_fraction = 0.;
        if (packed_cand->pdgId() == 1 || packed_cand->pdgId() == 130) {
          hcal_fraction = packed_cand->hcalFraction();
        } else if (packed_cand->isIsolatedChargedHadron()) {
          hcal_fraction = packed_cand->rawHcalFraction();
        }
        data.fillMulti<float>("pfcand_hcalFrac", hcal_fraction);
        data.fillMulti<float>("pfcand_hcalFracCalib", packed_cand->hcalFraction());

        // for charged
        data.fillMulti<float>("pfcand_VTX_ass", packed_cand->pvAssociationQuality());
        data.fillMulti<float>("pfcand_fromPV", packed_cand->fromPV());
        data.fillMulti<float>("pfcand_lostInnerHits", packed_cand->lostInnerHits());
        data.fillMulti<float>("pfcand_trackHighPurity", packed_cand->trackHighPurity());

        // impact parameters
        data.fillMulti<float>("pfcand_dz", catchInfs(packed_cand->dz()));
        data.fillMulti<float>("pfcand_dzsig", packed_cand->bestTrack() ? catchInfs(packed_cand->dz()/packed_cand->dzError()) : 0);
        data.fillMulti<float>("pfcand_dxy", catchInfs(packed_cand->dxy()));
        data.fillMulti<float>("pfcand_dxysig", packed_cand->bestTrack() ? catchInfs(packed_cand->dxy()/packed_cand->dxyError()) : 0);

        if (packed_cand->bestTrack()){
          const auto *trk = packed_cand->bestTrack();
          data.fillMulti<float>("pfcand_normchi2", catchInfs(trk->normalizedChi2()));
          data.fillMulti<float>("pfcand_quality", trk->qualityMask());
          data.fillMulti<float>("pfcand_nValidHits", trk->hitPattern().numberOfValidHits());
          data.fillMulti<float>("pfcand_nValidPixelHits", trk->hitPattern().numberOfValidPixelHits());

          // track covariance
          auto cov = [&](unsigned i, unsigned j) {
            return catchInfs(trk->covariance(i, j));
          };
          data.fillMulti<float>("pfcand_dptdpt", cov(0,0));
          data.fillMulti<float>("pfcand_detadeta", cov(1,1));
          data.fillMulti<float>("pfcand_dphidphi", cov(2,2));
          data.fillMulti<float>("pfcand_dxydxy", cov(3,3));
          data.fillMulti<float>("pfcand_dzdz", cov(4,4));
          data.fillMulti<float>("pfcand_dxydz", cov(3,4));
          data.fillMulti<float>("pfcand_dphidxy", cov(2,3));
          data.fillMulti<float>("pfcand_dlambdadz", cov(1,4));
        }else{
          data.fillMulti<float>("pfcand_normchi2", 999);
          data.fillMulti<float>("pfcand_quality", 0);
          data.fillMulti<float>("pfcand_nValidHits", 0);
          data.fillMulti<float>("pfcand_nValidPixelHits", 0);

          data.fillMulti<float>("pfcand_dptdpt", 0);
          data.fillMulti<float>("pfcand_detadeta", 0);
          data.fillMulti<float>("pfcand_dphidphi", 0);
          data.fillMulti<float>("pfcand_dxydxy", 0);
          data.fillMulti<float>("pfcand_dzdz", 0);
          data.fillMulti<float>("pfcand_dxydz", 0);
          data.fillMulti<float>("pfcand_dphidxy", 0);
          data.fillMulti<float>("pfcand_dlambdadz", 0);
        }

        // build track info map
        TrackInfoBuilder trkinfo;
        trkinfo.buildTrackInfo(builder_, *packed_cand, jet, vertices->at(0));

        // data.fillMulti<float>("pfcand_btagMomentum", catchInfs(trkinfo.getTrackMomentum()));
        // data.fillMulti<float>("pfcand_btagEta", catchInfs(trkinfo.getTrackEta()));
        data.fillMulti<float>("pfcand_btagEtaRel", catchInfs(trkinfo.getTrackEtaRel()));
        data.fillMulti<float>("pfcand_btagPtRel", catchInfs(trkinfo.getTrackPtRel()));
        // data.fillMulti<float>("pfcand_btagPPar", catchInfs(trkinfo.getTrackPPar()));
        // data.fillMulti<float>("pfcand_btagDeltaR", catchInfs(trkinfo.getTrackDeltaR()));
        data.fillMulti<float>("pfcand_btagPtRatio", catchInfs(trkinfo.getTrackPtRatio()));
        data.fillMulti<float>("pfcand_btagPParRatio", catchInfs(trkinfo.getTrackPParRatio()));
        data.fillMulti<float>("pfcand_btagSip2dVal", catchInfs(trkinfo.getTrackSip2dVal()));
        data.fillMulti<float>("pfcand_btagSip2dSig", catchInfs(trkinfo.getTrackSip2dSig()));
        data.fillMulti<float>("pfcand_btagSip3dVal", catchInfs(trkinfo.getTrackSip3dVal()));
        data.fillMulti<float>("pfcand_btagSip3dSig", catchInfs(trkinfo.getTrackSip3dSig()));
        data.fillMulti<float>("pfcand_btagJetDistVal", catchInfs(trkinfo.getTrackJetDistVal()));
        data.fillMulti<float>("pfcand_btagDecayLengthVal", catchInfs(trkinfo.getTrackDecayLengthVal()));
        data.fillMulti<float>("pfcand_btagDecayLengthSig", catchInfs(trkinfo.getTrackDecayLengthSig()));

        int b_tag=-1, c_tag=-1, g_tag=-1;
        float dist_from_pv=-1;
        float vtx_x=0, vtx_y=0, vtx_z=0;

        double dR_min=pow10(6);

        const reco::Candidate * pruned_part_match=nullptr;
        for (const auto &packed_part : *packed){
          double dR = reco::deltaR(*packed_cand, packed_part);
          double dpt = std::abs((packed_cand->pt()- packed_part.pt())/packed_cand->pt());

          if(dR<0.01 && dpt<0.1 && packed_cand->charge()==packed_part.charge()){
            if (dR<dR_min) {
              pruned_part_match=packed_part.lastPrunedRef().get();
              //pruned_part_match=packed_part.mother(0);
            }
          }
        }

        if (pruned_part_match != nullptr){
          c_tag=containParton(pruned_part_match, 4)? 1 : 0;
          b_tag=containParton(pruned_part_match, 5)? 1 : 0;
          g_tag=containParton(pruned_part_match, 21)? 1 : 0;

          dist_from_pv= sqrt((pv- pruned_part_match->vertex()).mag2());

          vtx_x=pruned_part_match->vertex().x();
          vtx_y=pruned_part_match->vertex().y();
          vtx_z=pruned_part_match->vertex().z();
        }

        data.fillMulti<int>("pfcand_from_b", b_tag);
        data.fillMulti<int>("pfcand_from_c", c_tag);
        data.fillMulti<int>("pfcand_from_g", g_tag);

        data.fillMulti<float>("pfcand_vtx_x", vtx_x);
        data.fillMulti<float>("pfcand_vtx_y", vtx_y);
        data.fillMulti<float>("pfcand_vtx_z", vtx_z);

        data.fillMulti<float>("pfcand_dist_from_pv", dist_from_pv);
      }
      break;
    }
  }
  return true;
}

} /* namespace deepntuples */
