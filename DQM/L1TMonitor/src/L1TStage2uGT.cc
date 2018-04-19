/**
 * \class L1TStage2uGT
 *
 * Description: DQM for L1 Micro Global Trigger.
 *
 * \author Mateusz Zarucki 2016
 * \author J. Berryhill, I. Mikulec
 * \author Vasile Mihai Ghete - HEPHY Vienna
 *
 */

#include "DQM/L1TMonitor/interface/L1TStage2uGT.h"

// Constructor
L1TStage2uGT::L1TStage2uGT(const edm::ParameterSet& params):
   l1tStage2uGtSource_(consumes<GlobalAlgBlkBxCollection>(params.getParameter<edm::InputTag>("l1tStage2uGtSource"))),
   monitorDir_(params.getUntrackedParameter<std::string> ("monitorDir", "")),
   verbose_(params.getUntrackedParameter<bool>("verbose", false)),
   gtUtil_(new l1t::L1TGlobalUtil(params, consumesCollector(), *this, params.getParameter<edm::InputTag>("l1tStage2uGtSource"), params.getParameter<edm::InputTag>("l1tStage2uGtSource"))),
   numAlgs_(0)
{
}

// Destructor
L1TStage2uGT::~L1TStage2uGT() {}

void L1TStage2uGT::dqmBeginRun(edm::Run const& iRun, edm::EventSetup const& evtSetup) {
   int algoBitPre_=-1; 
   int algoBitUnpre_=-1; 
   // Get the trigger menu information
   gtUtil_->retrieveL1Setup(evtSetup);
   // Find the number of algos defined
   numAlgs_ = static_cast<int>(gtUtil_->decisionsInitial().size());
}


void L1TStage2uGT::bookHistograms(DQMStore::IBooker &ibooker, edm::Run const&, edm::EventSetup const& evtSetup) {
   
   // Book histograms
   const int numLS = 2000;
   const double numLS_d = static_cast<double>(numLS);
   const double numAlgs_d = static_cast<double>(numAlgs_);
   const int numBx = 3564; 
   const double numBx_d = static_cast<double>(numBx);

   ibooker.setCurrentFolder(monitorDir_);
   
   // Algorithm bits 
   algoBits_before_bxmask_ = ibooker.book1D("algoBits_before_bxmask", "uGT: Algorithm Trigger Bits (before AlgoBX mask)", numAlgs_, -0.5, numAlgs_d-0.5);
   algoBits_before_bxmask_->setAxisTitle("Algorithm Trigger Bits (before AlgoBX mask)", 1);
   
   algoBits_before_prescale_ = ibooker.book1D("algoBits_before_prescale", "uGT: Algorithm Trigger Bits (before prescale)", numAlgs_, -0.5, numAlgs_d-0.5);
   algoBits_before_prescale_->setAxisTitle("Algorithm Trigger Bits (before prescale)", 1);
   
   algoBits_after_prescale_ = ibooker.book1D("algoBits_after_prescale", "uGT: Algorithm Trigger Bits (after prescale)", numAlgs_, -0.5, numAlgs_d-0.5);
   algoBits_after_prescale_->setAxisTitle("Algorithm Trigger Bits (after prescale)", 1);
  
   // Algorithm bits correlation 
   algoBits_before_bxmask_corr_ = ibooker.book2D("algoBits_before_bxmask_corr","uGT: Algorithm Trigger Bit Correlation (before AlgoBX mask)", numAlgs_, -0.5, numAlgs_d-0.5, numAlgs_, -0.5, numAlgs_d-0.5);
   algoBits_before_bxmask_corr_->setAxisTitle("Algorithm Trigger Bits (before AlgoBX mask)", 1);
   algoBits_before_bxmask_corr_->setAxisTitle("Algorithm Trigger Bits (before AlgoBX mask)", 2);
   
   algoBits_before_prescale_corr_ = ibooker.book2D("algoBits_before_prescale_corr","uGT: Algorithm Trigger Bit Correlation (before prescale)", numAlgs_, -0.5, numAlgs_d-0.5, numAlgs_, -0.5, numAlgs_d-0.5);
   algoBits_before_prescale_corr_->setAxisTitle("Algorithm Trigger Bits (before prescale)", 1);
   algoBits_before_prescale_corr_->setAxisTitle("Algorithm Trigger Bits (before prescale)", 2);
   
   algoBits_after_prescale_corr_ = ibooker.book2D("algoBits_after_prescale_corr","uGT: Algorithm Trigger Bit Correlation (after prescale)", numAlgs_, -0.5, numAlgs_d-0.5, numAlgs_, -0.5, numAlgs_d-0.5);
   algoBits_after_prescale_corr_->setAxisTitle("Algorithm Trigger Bits (after prescale)", 1);
   algoBits_after_prescale_corr_->setAxisTitle("Algorithm Trigger Bits (after prescale)", 2);
  
   // Algorithm bits vs global BX number
   algoBits_before_bxmask_bx_global_ = ibooker.book2D("algoBits_before_bxmask_bx_global", "uGT: Algorithm Trigger Bits (before AlgoBX mask) vs. Global BX Number", numBx, 0.5, numBx_d + 0.5, numAlgs_, -0.5, numAlgs_d-0.5);
   algoBits_before_bxmask_bx_global_->setAxisTitle("Global Bunch Crossing Number", 1); 
   algoBits_before_bxmask_bx_global_->setAxisTitle("Algorithm Trigger Bits (before AlgoBX mask)", 2);
   
   algoBits_before_prescale_bx_global_ = ibooker.book2D("algoBits_before_prescale_bx_global", "uGT: Algorithm Trigger Bits (before prescale) vs. Global BX Number", numBx, 0.5, numBx_d + 0.5, numAlgs_, -0.5, numAlgs_d-0.5);
   algoBits_before_prescale_bx_global_->setAxisTitle("Global Bunch Crossing Number", 1); 
   algoBits_before_prescale_bx_global_->setAxisTitle("Algorithm Trigger Bits (before prescale)", 2);
   
   algoBits_after_prescale_bx_global_ = ibooker.book2D("algoBits_after_prescale_bx_global", "uGT: Algorithm Trigger Bits (after prescale) vs. Global BX Number", numBx, 0.5, numBx_d + 0.5, numAlgs_, -0.5, numAlgs_d-0.5);
   algoBits_after_prescale_bx_global_->setAxisTitle("Global Bunch Crossing Number", 1); 
   algoBits_after_prescale_bx_global_->setAxisTitle("Algorithm Trigger Bits (after prescale)", 2);
  
   // Algorithm bits vs BX number in event
   algoBits_before_bxmask_bx_inEvt_ = ibooker.book2D("algoBits_before_bxmask_bx_inEvt", "uGT: Algorithm Trigger Bits (before AlgoBX mask) vs. BX Number in Event", 5, -2.5, 2.5, numAlgs_, -0.5, numAlgs_d-0.5);
   algoBits_before_bxmask_bx_inEvt_->setAxisTitle("Bunch Crossing Number in Event", 1);
   algoBits_before_bxmask_bx_inEvt_->setAxisTitle("Algorithm Trigger Bits (before AlgoBX mask)", 2);
   
   algoBits_before_prescale_bx_inEvt_ = ibooker.book2D("algoBits_before_prescale_bx_inEvt", "uGT: Algorithm Trigger Bits (before prescale) vs. BX Number in Event", 5, -2.5, 2.5, numAlgs_, -0.5, numAlgs_d-0.5);
   algoBits_before_prescale_bx_inEvt_->setAxisTitle("Bunch Crossing Number in Event", 1);
   algoBits_before_prescale_bx_inEvt_->setAxisTitle("Algorithm Trigger Bits (before prescale)", 2);
   
   algoBits_after_prescale_bx_inEvt_ = ibooker.book2D("algoBits_after_prescale_bx_inEvt", "uGT: Algorithm Trigger Bits (after prescale) vs. BX Number in Event", 5, -2.5, 2.5, numAlgs_, -0.5, numAlgs_d-0.5);
   algoBits_after_prescale_bx_inEvt_->setAxisTitle("Bunch Crossing Number in Event", 1);
   algoBits_after_prescale_bx_inEvt_->setAxisTitle("Algorithm Trigger Bits (after prescale)", 2);
  
   // Algorithm bits vs LS
   algoBits_before_bxmask_lumi_ = ibooker.book2D("algoBits_before_bxmask_lumi","uGT: Algorithm Trigger Bits (before AlgoBX mask) vs. LS", numLS, 0., numLS_d, numAlgs_, -0.5, numAlgs_d-0.5);
   algoBits_before_bxmask_lumi_->setAxisTitle("Luminosity Segment", 1);
   algoBits_before_bxmask_lumi_->setAxisTitle("Algorithm Trigger Bits (before AlgoBX mask)", 2);
   
   algoBits_before_prescale_lumi_ = ibooker.book2D("algoBits_before_prescale_lumi","uGT: Algorithm Trigger Bits (before prescale) vs. LS", numLS, 0., numLS_d, numAlgs_, -0.5, numAlgs_d-0.5);
   algoBits_before_prescale_lumi_->setAxisTitle("Luminosity Segment", 1);
   algoBits_before_prescale_lumi_->setAxisTitle("Algorithm Trigger Bits (before prescale)", 2);
  
   algoBits_after_prescale_lumi_ = ibooker.book2D("algoBits_after_prescale_lumi","uGT: Algorithm Trigger Bits (after prescale) vs. LS", numLS, 0., numLS_d, numAlgs_, -0.5, numAlgs_d-0.5);
   algoBits_after_prescale_lumi_->setAxisTitle("Luminosity Segment", 1);
   algoBits_after_prescale_lumi_->setAxisTitle("Algorithm Trigger Bits (after prescale)", 2);

   // Prescale factor index 
   prescaleFactorSet_ = ibooker.book2D("prescaleFactorSet", "uGT: Index of Prescale Factor Set vs. LS", numLS, 0., numLS_d, 25, 0., 25.);
   prescaleFactorSet_->setAxisTitle("Luminosity Segment", 1);
   prescaleFactorSet_->setAxisTitle("Prescale Factor Set Index", 2);

  // Prescaled and Unprescaled Algo Trigger Bits
  prescaled_first_collision_run_ = ibooker.book2D("prescaled_first_collision_run_", "uGT: Prescaled Algorithm Trigger Bits  vs. BX Number In Event", 5, -2.5, 2.5, prescaledAlgoBitName_.size(), -0.5, preAlgs_d-0.5);
  prescaled_first_collision_run_->setAxisTitle("Bunch Crossing Number In Event", 1);
  prescaled_first_collision_run_->setAxisTitle("Algorithm Trigger Bits + Names ", 2);
  for(unsigned int algo=1; algo<prescaledAlgoBitName_.size(); algo++) {
    prescaled_first_collision_run_->setBinLabel(algo,prescaledAlgoBitName_.at(algo).first+"("+std::to_string(prescaledAlgoBitName_.at(algo).second)+")",2);
  }

  den_prescaled_first_collision_run_ = ibooker.book2D("den_prescaled_first_collision_run_", "uGT: Prescaled Algorithm Trigger Bits Deno vs. BX Number In Event", 5, -2.5, 2.5, numAlgs, -0.5, numAlgs_d-0.5);
  den_prescaled_first_collision_run_->setAxisTitle("Bunch Crossing Number In Event", 1);
  den_prescaled_first_collision_run_->setAxisTitle("Algorithm Trigger Bits + Names", 2);
  for(unsigned int algo=1; algo<prescaledAlgoBitName_.size(); algo++) {
    den_prescaled_first_collision_run_->setBinLabel(algo,prescaledAlgoBitName_.at(algo).first+"("+std::to_string(prescaledAlgoBitName_.at(algo).second)+")",2);
  }

  unprescaled_algo_first_collision_in_train_ = ibooker.book2D("unprescaled_algo_first_collision_in_train_", "uGT: Unprescaled Algorithm Trigger Bits  vs. BX Number In Event", 5, -2.5, 2.5, unprescaledAlgoBitName_.size(), -0.5, unpreAlgs_d-0.5);
  unprescaled_algo_first_collision_in_train_->setAxisTitle("Bunch Crossing Number In Event", 1);
  unprescaled_algo_first_collision_in_train_->setAxisTitle("Algorithm Trigger Bits + Names", 2);
  for(unsigned int algo=0; algo<unprescaledAlgoBitName_.size(); algo++) {
    unprescaled_first_collision_run_->setBinLabel(algo,unprescaledAlgoBitName_.at(algo).first+"("+std::to_string(unprescaledAlgoBitName_.at(algo).second)+")",2);
  }

  den_unprescaled_first_collision_run_ = ibooker.book2D("den_unprescaled_first_collision_run_", "uGT: Unprescaled Algorithm Trigger Bits Deno vs. BX Number In Event", 5, -2.5, 2.5, numAlgs, -0.5, numAlgs_d-0.5);
  den_unprescaled_first_collision_run_->setAxisTitle("Bunch Crossing Number In Event", 1);
  den_unprescaled_first_collision_run_->setAxisTitle("Algorithm Trigger Bits + Names", 2);
  for(unsigned int algo=1; algo<unprescaledAlgoBitName_.size(); algo++) {
    den_unprescaled_first_collision_run_->setBinLabel(algo,unprescaledAlgoBitName_.at(algo).first+"("+std::to_string(unprescaledAlgoBitName_.at(algo).second)+")",2);
  }

  prescaled_algo_isolated_collision_in_train_ = ibooker.book2D("prescaled_algo_isolated_collision_in_train_", "uGT: Prescaled Algorithm Trigger Bits vs. BX Number In Event", 5, -2.5, 2.5, prescaledAlgoBitName_.size(), -0.5, preAlgs_d-0.5);
  prescaled_algo_isolated_collision_in_train_->setAxisTitle("Bunch Crossing Number In Event", 1);
  prescaled_algo_isolated_collision_in_train_->setAxisTitle("Algorithm Trigger Bits + Names", 2);
  for(unsigned int algo=0; algo<prescaledAlgoBitName_.size(); algo++) {
    prescaled_isolated_collision_run_->setBinLabel(algo,prescaledAlgoBitName_.at(algo).first+"("+std::to_string(prescaledAlgoBitName_.at(algo).second)+")",2);
  }

  den_prescaled_isolated_collision_run_ = ibooker.book2D("den_prescaled_isolated_collision_run_", "uGT: Prescaled Algorithm Trigger Bits Deno vs. BX Number In Event", 5, -2.5, 2.5, numAlgs, -0.5, numAlgs_d-0.5);
  den_prescaled_isolated_collision_run_->setAxisTitle("Bunch Crossing Number In Event", 1);
  den_prescaled_isolated_collision_run_->setAxisTitle("Algorithm Trigger Bits + Names", 2);
  for(unsigned int algo=1; algo<prescaledAlgoBitName_.size(); algo++) {
    den_prescaled_isolated_collision_run_->setBinLabel(algo,prescaledAlgoBitName_.at(algo).first+"("+std::to_string(prescaledAlgoBitName_.at(algo).second)+")",2);
  }

  unprescaled_algo_isolated_collision_in_train_ = ibooker.book2D("unprescaled_algo_isolated_collision_in_train_", "uGT: Unprescaled Algorithm Trigger Bits vs. BX Number In Event", 5, -2.5, 2.5, unprescaledAlgoBitName_.size(), -0.5, unpreAlgs_d-0.5);
  unprescaled_algo_isolated_collision_in_train_->setAxisTitle("Bunch Crossing Number In Event", 1);
  unprescaled_algo_isolated_collision_in_train_->setAxisTitle("Algorithm Trigger Bits + Names", 2);
  for(unsigned int algo=0; algo<unprescaledAlgoBitName_.size(); algo++) {
    unprescaled_isolated_collision_run_->setBinLabel(algo,unprescaledAlgoBitName_.at(algo).first+"("+std::to_string(unprescaledAlgoBitName_.at(algo).second)+")",2);
  }

  unprescaled_algo_isolated_collision_in_train_den_ = ibooker.book2D("unprescaled_algo_isolated_collision_in_train_den_", "uGT: Unprescaled Algorithm Trigger Bits Deno vs. BX Number In Event", 5, -2.5, 2.5, unprescaledAlgoBitName_.size(), -0.5, unpreAlgs_d-0.5);
  unprescaled_algo_isolated_collision_in_train_den_->setAxisTitle("Bunch Crossing Number In Event", 1);
  unprescaled_algo_isolated_collision_in_train_den_->setAxisTitle("Algorithm Trigger Bits + Names", 2);
  for(unsigned int algo=0; algo<unprescaledAlgoBitName_.size(); algo++) {
    den_unprescaled_isolated_collision_run_->setBinLabel(algo,unprescaledAlgoBitName_.at(algo).first+"("+std::to_string(unprescaledAlgoBitName_.at(algo).second)+")",2);
  }

  prescaled_algo_last_collision_in_train_ = ibooker.book2D("prescaled_algo_last_collision_in_train_", "uGT: Prescaled Algorithm Trigger Bits vs. BX Number In Event", 5, -2.5, 2.5, prescaledAlgoBitName_.size(), -0.5, preAlgs_d-0.5);
  prescaled_algo_last_collision_in_train_->setAxisTitle("Bunch Crossing Number In Event", 1);
  prescaled_algo_last_collision_in_train_->setAxisTitle("Algorithm Trigger Bits + Names", 2);
  for(unsigned int algo=0; algo<prescaledAlgoBitName_.size(); algo++) {
    prescaled_last_collision_run_->setBinLabel(algo,prescaledAlgoBitName_.at(algo).first+"("+std::to_string(prescaledAlgoBitName_.at(algo).second)+")",2);
  }

  prescaled_algo_last_collision_in_train_den_ = ibooker.book2D("prescaled_algo_last_collision_in_train_den_", "uGT: Prescaled Algorithm Trigger Bits Deno vs. BX Number In Event", 5, -2.5, 2.5, prescaledAlgoBitName_.size(), -0.5, preAlgs_d-0.5);
  prescaled_algo_last_collision_in_train_den_->setAxisTitle("Bunch Crossing Number In Event", 1);
  prescaled_algo_last_collision_in_train_den_->setAxisTitle("Algorithm Trigger Bits + Names", 2);
  for(unsigned int algo=0; algo<prescaledAlgoBitName_.size(); algo++) {
    den_prescaled_last_collision_run_->setBinLabel(algo,prescaledAlgoBitName_.at(algo).first+"("+std::to_string(prescaledAlgoBitName_.at(algo).second)+")",2);
  }

  unprescaled_algo_last_collision_in_train_ = ibooker.book2D("unprescaled_algo_last_collision_in_train_", "uGT: Unprescaled Algorithm Trigger Bits vs. BX Number In Event", 5, -2.5, 2.5, unprescaledAlgoBitName_.size(), -0.5, unpreAlgs_d-0.5);
  unprescaled_algo_last_collision_in_train_->setAxisTitle("Bunch Crossing Number In Event", 1);
  unprescaled_algo_last_collision_in_train_->setAxisTitle("Algorithm Trigger Bits + Names", 2);
  for(unsigned int algo=0; algo<unprescaledAlgoBitName_.size(); algo++) {
    unprescaled_last_collision_run_->setBinLabel(algo,unprescaledAlgoBitName_.at(algo).first+"("+std::to_string(unprescaledAlgoBitName_.at(algo).second)+")",2);
  }

  unprescaled_algo_last_collision_in_train_den_ = ibooker.book2D("unprescaled_algo_last_collision_in_train_den_", "uGT: Unprescaled Algorithm Trigger Bits Deno vs. BX Number In Event", 5, -2.5, 2.5, unprescaledAlgoBitName_.size(), -0.5, unpreAlgs_d-0.5);
  unprescaled_algo_last_collision_in_train_den_->setAxisTitle("Bunch Crossing Number In Event", 1);
  unprescaled_algo_last_collision_in_train_den_->setAxisTitle("Algorithm Trigger Bits + Names", 2);
  for(unsigned int algo=0; algo<unprescaledAlgoBitName_.size(); algo++) {
    den_unprescaled_last_collision_run_->setBinLabel(algo,unprescaledAlgoBitName_.at(algo).first+"("+std::to_string(unprescaledAlgoBitName_.at(algo).second)+")",2);
  }

}

void L1TStage2uGT::analyze(const edm::Event& evt, const edm::EventSetup& evtSetup) {
   if (verbose_) {
      edm::LogInfo("L1TStage2uGT") << "L1TStage2uGT DQM: Analyzing.." << std::endl;
   }
   
   // Get standard event parameters 
   int lumi = evt.luminosityBlock();
   int bx = evt.bunchCrossing();
      
   // Open uGT readout record
   edm::Handle<GlobalAlgBlkBxCollection> uGtAlgs;
   evt.getByToken(l1tStage2uGtSource_, uGtAlgs);
   
   if (!uGtAlgs.isValid()) {
      edm::LogInfo("L1TStage2uGT") << "Cannot find uGT readout record.";
      return;
   }
   
   //algoBits_->Fill(-1.); // fill underflow to normalize // FIXME: needed?
   for (int ibx = uGtAlgs->getFirstBX(); ibx <= uGtAlgs->getLastBX(); ++ibx) {
      for (auto itr = uGtAlgs->begin(ibx); itr != uGtAlgs->end(ibx); ++itr) { // FIXME: redundant loop?

         // Fills prescale factor set histogram
         prescaleFactorSet_->Fill(lumi, itr->getPreScColumn());

         // Fills algorithm bits histograms
         for(int algoBit = 0; algoBit < numAlgs_; ++algoBit) {

            // Algorithm bits before AlgoBX mask
            if(itr->getAlgoDecisionInitial(algoBit)) {
               algoBits_before_bxmask_->Fill(algoBit);
               algoBits_before_bxmask_lumi_->Fill(lumi, algoBit);
               algoBits_before_bxmask_bx_inEvt_->Fill(ibx, algoBit); // FIXME: or itr->getbxInEventNr()/getbxNr()?
               algoBits_before_bxmask_bx_global_->Fill(bx + ibx, algoBit);

               for(int algoBit2 = 0; algoBit2 < numAlgs_; ++algoBit2) {
                  if(itr->getAlgoDecisionInitial(algoBit2)) {
                     algoBits_before_bxmask_corr_->Fill(algoBit, algoBit2);
                  }
               }
            }

            // Algorithm bits before prescale
            if(itr->getAlgoDecisionInterm(algoBit)) {
               algoBits_before_prescale_->Fill(algoBit);
               algoBits_before_prescale_lumi_->Fill(lumi, algoBit);
               algoBits_before_prescale_bx_inEvt_->Fill(ibx, algoBit);
               algoBits_before_prescale_bx_global_->Fill(bx + ibx, algoBit);

               for(int algoBit2 = 0; algoBit2 < numAlgs_; ++algoBit2) {
                  if(itr->getAlgoDecisionInterm(algoBit2)) {
                     algoBits_before_prescale_corr_->Fill(algoBit, algoBit2);
                  }
               }
            }

            // Algorithm bits after prescale
            if(itr->getAlgoDecisionFinal(algoBit)) {
               algoBits_after_prescale_->Fill(algoBit);
               algoBits_after_prescale_lumi_->Fill(lumi, algoBit);
               algoBits_after_prescale_bx_inEvt_->Fill(ibx, algoBit);
               algoBits_after_prescale_bx_global_->Fill(bx + ibx, algoBit);

               for(int algoBit2 = 0; algoBit2 < numAlgs_; ++algoBit2) {
                  if(itr->getAlgoDecisionFinal(algoBit2)) {
                     algoBits_after_prescale_corr_->Fill(algoBit, algoBit2);
                  }
               }
            }
         }
      }
   }
}
      
// End section
