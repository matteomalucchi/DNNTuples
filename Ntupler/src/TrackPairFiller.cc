#include <unordered_map>
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"

#include "DeepNTuples/Ntupler/interface/TrackPairFiller.h"

namespace deepntuples {

void TrackPairFiller::readConfig(const edm::ParameterSet& iConfig, edm::ConsumesCollector&& cc) {
  vtxToken_ = cc.consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("vertices"));
  svToken_ = cc.consumes<reco::VertexCompositePtrCandidateCollection>(iConfig.getParameter<edm::InputTag>("SVs"));
  packedToken_ = cc.consumes<std::vector<pat::PackedGenParticle>>(iConfig.getParameter<edm::InputTag>("packed"));
}

void TrackPairFiller::readEvent(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
  iEvent.getByToken(vtxToken_, vertices);
  iEvent.getByToken(svToken_, SVs);
  iEvent.getByToken(packedToken_, packed);
  iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder", builder_);
}

void TrackPairFiller::book() {

  //test
  data.add<int>("n_pfcandidates", 0);

  data.addMulti<int>("track1_index");
  data.addMulti<int>("track2_index");


  data.addMulti<float>("pt_1");
  data.addMulti<float>("pt_2");


  // pca (1 to 2)
  data.addMulti<float>("pca_distance");
  data.addMulti<float>("pca_significance");

  //   pcas (1 and 2) poistions
  data.addMulti<float>("pcaSeed_x1");
  data.addMulti<float>("pcaSeed_y1");
  data.addMulti<float>("pcaSeed_z1");

  data.addMulti<float>("pcaSeed_x2");
  data.addMulti<float>("pcaSeed_y2");
  data.addMulti<float>("pcaSeed_z2");

  data.addMulti<float>("pcaSeed_xerr1");
  data.addMulti<float>("pcaSeed_yerr1");
  data.addMulti<float>("pcaSeed_zerr1");

  data.addMulti<float>("pcaSeed_xerr2");
  data.addMulti<float>("pcaSeed_yerr2");
  data.addMulti<float>("pcaSeed_zerr2");

  // dot prod betweeen track and pca direction
  data.addMulti<float>("dotprod1");
  data.addMulti<float>("dotprod2");

  //pca distance form PV on both tracks
  data.addMulti<float>("pca_dist1");
  data.addMulti<float>("pca_dist2");

  //track track or dir dir dotprod
  data.addMulti<float>("dotprod12_2D");
  data.addMulti<float>("dotprod12_2DV");
  data.addMulti<float>("dotprod12_3D");
  data.addMulti<float>("dotprod12_3DV");

  //jet pca relative
  data.addMulti<float>("pca_jetAxis_dist");
  data.addMulti<float>("pca_jetAxis_dotprod");
  data.addMulti<float>("pca_jetAxis_dEta");
  data.addMulti<float>("pca_jetAxis_dPhi_");

  data.addMulti<int>("index_pf1");
  data.addMulti<int>("index_pf2");
  data.addMulti<float>("dist_vtx_12");

}

bool TrackPairFiller::fill(const pat::Jet& jet, size_t jetidx, const JetHelper& jet_helper) {

  const auto& pfCands = jet_helper.getJetConstituents();

  data.fill<int>("n_pfcandidates", pfCands.size());

  std::vector<reco::TransientTrack> selectedTracks;
  std::vector<int> selectedTracks_pfidx;
  int counter = 0;


  for (const auto& cand : pfCands){

    //const auto *packed_cand = dynamic_cast<const pat::PackedCandidate *>(&(*cand));

    if (cand->bestTrack()){
        selectedTracks.push_back( builder_->build(cand) );
        selectedTracks_pfidx.push_back(counter);
    }
    counter++;

  }

  // for each possible pfCand pair save the distance between the vertices of the
  // last pruned ancestors of each particle
  for(std::vector<reco::CandidatePtr>::const_iterator it1 = pfCands.begin(); it1 != pfCands.end(); it1++){
    int index1 = it1 - pfCands.begin();
    const auto *packed_cand1 = dynamic_cast<const pat::PackedCandidate *>(&(*(*it1)));

    for(std::vector<reco::CandidatePtr>::const_iterator it2 = pfCands.begin(); it2 != pfCands.end(); it2++){
      int index2 = it2 - pfCands.begin();
      const auto *packed_cand2 = dynamic_cast<const pat::PackedCandidate *>(&(*(*it2)));
      data.fillMulti<int>("index_pf1", index1);
      data.fillMulti<int>("index_pf2", index2);

      float dist_vtx_12 = -1;
      for (const auto &packed_part1 : *packed){
        double dR1 = reco::deltaR(*packed_cand1, packed_part1);
        double dpt1 = std::abs((packed_cand1->pt()- packed_part1.pt())/packed_cand1->pt());
        if(dR1<0.01 && dpt1<0.1 && packed_cand1->charge()==packed_part1.charge()){

          for (const auto &packed_part2 : *packed){
            double dR2 = reco::deltaR(*packed_cand2, packed_part2);
            double dpt2 = std::abs((packed_cand2->pt()- packed_part2.pt())/packed_cand2->pt());
            if(dR2<0.01 && dpt2<0.1 && packed_cand2->charge()==packed_part2.charge()){

              const reco::Candidate * pruned_part1=packed_part1.lastPrunedRef().get();
              const reco::Candidate * pruned_part2=packed_part2.lastPrunedRef().get();

              dist_vtx_12= sqrt((pruned_part1->vertex()- pruned_part2->vertex()).mag2());

              break;
            }
          }
          break;
        }
      }
      data.fillMulti<float>("dist_vtx_12", dist_vtx_12);
    }
  }



  for(std::vector<reco::TransientTrack>::const_iterator it = selectedTracks.begin(); it != selectedTracks.end(); it++){

      int index1 = it - selectedTracks.begin();

      for(std::vector<reco::TransientTrack>::const_iterator tt = selectedTracks.begin(); tt != selectedTracks.end(); tt++){

          int index2 = tt - selectedTracks.begin();

          if (index1!=index2 ){
            TrackPairInfoBuilder trkpairinfo;
            trkpairinfo.buildTrackPairInfo(&(*it),&(*tt),vertices->at(0),jet);

            if (trkpairinfo.pca_distance()/trkpairinfo.pcaSeed_dist()<20. && selectedTracks_pfidx[index1]<50 && selectedTracks_pfidx[index2]<50){

            data.fillMulti<int>("track1_index", selectedTracks_pfidx[index1]);
            data.fillMulti<int>("track2_index", selectedTracks_pfidx[index2]);

            data.fillMulti<float>("pt_1", trkpairinfo.track_i_pt());
            data.fillMulti<float>("pt_2", trkpairinfo.track_t_pt());


            data.fillMulti<float>("pca_distance", trkpairinfo.pca_distance());
            data.fillMulti<float>("pca_significance", trkpairinfo.pca_significance());

            //   pcas (1 and 2) poistions
            data.fillMulti<float>("pcaSeed_x1", trkpairinfo.pcaSeed_x());
            data.fillMulti<float>("pcaSeed_y1", trkpairinfo.pcaSeed_y());
            data.fillMulti<float>("pcaSeed_z1", trkpairinfo.pcaSeed_z());

            data.fillMulti<float>("pcaSeed_x2", trkpairinfo.pcaTrack_x());
            data.fillMulti<float>("pcaSeed_y2", trkpairinfo.pcaTrack_y());
            data.fillMulti<float>("pcaSeed_z2", trkpairinfo.pcaTrack_z());

            data.fillMulti<float>("pcaSeed_xerr1", trkpairinfo.pcaSeed_xerr());
            data.fillMulti<float>("pcaSeed_yerr1", trkpairinfo.pcaSeed_yerr());
            data.fillMulti<float>("pcaSeed_zerr1", trkpairinfo.pcaSeed_zerr());

            data.fillMulti<float>("pcaSeed_xerr2", trkpairinfo.pcaTrack_xerr());
            data.fillMulti<float>("pcaSeed_yerr2", trkpairinfo.pcaTrack_yerr());
            data.fillMulti<float>("pcaSeed_zerr2", trkpairinfo.pcaTrack_zerr());

            // dot prod betweeen track and pca direction
            data.fillMulti<float>("dotprod1", trkpairinfo.dotprodTrack());
            data.fillMulti<float>("dotprod2", trkpairinfo.dotprodSeed());

            //pca distance form PV on both tracks
            data.fillMulti<float>("pca_dist1", trkpairinfo.pcaSeed_dist());
            data.fillMulti<float>("pca_dist2", trkpairinfo.pcaTrack_dist());

            //track track or dir dir dotprod
            data.fillMulti<float>("dotprod12_2D", trkpairinfo.dotprodTrackSeed2D());
            data.fillMulti<float>("dotprod12_2DV", trkpairinfo.dotprodTrackSeed2DV());
            data.fillMulti<float>("dotprod12_3D", trkpairinfo.dotprodTrackSeed3D());
            data.fillMulti<float>("dotprod12_3DV", trkpairinfo.dotprodTrackSeed3DV());

            //jet pca relative
            data.fillMulti<float>("pca_jetAxis_dist", trkpairinfo.pca_jetAxis_dist());
            data.fillMulti<float>("pca_jetAxis_dotprod", trkpairinfo.pca_jetAxis_dotprod());
            data.fillMulti<float>("pca_jetAxis_dEta", trkpairinfo.pca_jetAxis_dEta());
            data.fillMulti<float>("pca_jetAxis_dPhi_", trkpairinfo.pca_jetAxis_dPhi());


            /*for (const auto &packed_part : *packed){
              double dR = reco::deltaR(*packed_cand, packed_part);
              double dpt = std::abs((packed_cand->pt()- packed_part.pt())/packed_cand->pt());


            }*/
            }
          }


      }


  }


  return true;
}

} /* namespace deepntuples */
