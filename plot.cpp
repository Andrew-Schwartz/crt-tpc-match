#include <TFile.h>
#include <TTree.h>
#include <TPolyLine3D.h>
#include <TPolyMarker3D.h>
#include <TMarker3DBox.h>
#include <TCanvas.h>
#include <iostream>

void plot(const char *hitDumperTreeDir = ".", bool strips = true) {
  std::string dir(hitDumperTreeDir);
  auto *c1 = new TCanvas("c1", "c1");
  auto *file = TFile::Open(dir.append("/hitdumper_tree.root").c_str());
  auto *hitdumper = (TDirectoryFile *) file->Get("hitdumper");
  hitdumper->cd();
  auto *tree = (TTree *) hitdumper->Get("hitdumpertree");

  // Tracks
  double ct_x1[50], ct_x2[50],
      ct_y1[50], ct_y2[50],
      ct_z1[50], ct_z2[50];
  int ncts;
  double ctrk_x1[50], ctrk_x2[50],
      ctrk_y1[50], ctrk_y2[50],
      ctrk_z1[50], ctrk_z2[50];
  int nctrks;

  // Hits
  double x[500], y[500], z[500], sipm_dist[500],
      cx[500], cy[500], cz[500], hw[500], hh[500], hl[500];
  int nhits;

  // SiPMs
  int nstrips;
  double sx1[3000], sy1[3000], sz1[3000], sx2[3000], sy2[3000], sz2[3000];

  tree->SetBranchAddress("ct_x1", &ct_x1);
  tree->SetBranchAddress("ct_y1", &ct_y1);
  tree->SetBranchAddress("ct_z1", &ct_z1);
  tree->SetBranchAddress("ct_x2", &ct_x2);
  tree->SetBranchAddress("ct_y2", &ct_y2);
  tree->SetBranchAddress("ct_z2", &ct_z2);
  tree->SetBranchAddress("ncts", &ncts);

  tree->SetBranchAddress("ctrk_x1", &ctrk_x1);
  tree->SetBranchAddress("ctrk_y1", &ctrk_y1);
  tree->SetBranchAddress("ctrk_z1", &ctrk_z1);
  tree->SetBranchAddress("ctrk_x2", &ctrk_x2);
  tree->SetBranchAddress("ctrk_y2", &ctrk_y2);
  tree->SetBranchAddress("ctrk_z2", &ctrk_z2);
  tree->SetBranchAddress("nctrks", &nctrks);

  tree->SetBranchAddress("chit_x", &x);
  tree->SetBranchAddress("chit_y", &y);
  tree->SetBranchAddress("chit_z", &z);
  int set_sipm = tree->SetBranchAddress("chit_sipm_dist", &sipm_dist);
  bool has_sipm = set_sipm != TTree::ESetBranchAddressStatus::kMissingBranch;
  /*tree->SetBranchAddress("chit_vol_x", &cx);
  tree->SetBranchAddress("chit_vol_y", &cy);
  tree->SetBranchAddress("chit_vol_z", &cz);
  tree->SetBranchAddress("chit_vol_hwidth", &hw);
  tree->SetBranchAddress("chit_vol_hheight", &hh);
  tree->SetBranchAddress("chit_vol_hlength", &hl);
  */
  tree->SetBranchAddress("nchits", &nhits);

  if (strips) {
    tree->SetBranchAddress("tot_strips", &nstrips);
    tree->SetBranchAddress("crt_min_x", &sx1);
    tree->SetBranchAddress("crt_max_x", &sx2);
    tree->SetBranchAddress("crt_min_y", &sy1);
    tree->SetBranchAddress("crt_max_y", &sy2);
    tree->SetBranchAddress("crt_min_z", &sz1);
    tree->SetBranchAddress("crt_max_z", &sz2);
  }

  auto *hits = new TPolyMarker3D(100);
  auto *inAds = new TPolyMarker3D(100);

  long len = tree->GetEntries();
//  for (long i = 0; i < 1; ++i) {
  for (long i = 0; i < len; ++i) {
    tree->GetEntry(i);

    for (int j = 0; j < ncts; ++j) {
      double xs[2] {ct_x1[j], ct_x2[j]},
          ys[2] {ct_y1[j], ct_y2[j]},
          zs[2] {ct_z1[j], ct_z2[j]};
      auto *track = new TPolyLine3D(2, xs, ys, zs);
      track->SetLineColor(2);
      if (strips) track->SetLineWidth(2);
      track->Draw();
    }
    for (int j = 0; j < nctrks; ++j) {
      double xs[] {ctrk_x1[j], ctrk_x2[j]},
          ys[] {ctrk_y1[j], ctrk_y2[j]},
          zs[] {ctrk_z1[j], ctrk_z2[j]};
      auto *track = new TPolyLine3D(2, xs, ys, zs);
      track->SetLineColor(kRed - 7);
      track->Draw();
    }
    for (int j = 0; j < nhits; ++j) {
      if (has_sipm && sipm_dist[i] != -9999.9) {
        inAds->SetNextPoint(x[j], y[j], z[j]);
      } else {
        hits->SetNextPoint(x[j], y[j], z[j]);
      }
      auto vol = new TMarker3DBox(cx[j], cy[j], cz[j], hw[j], hh[j], hl[j], 0.0, 0.0);
      vol->Draw();
    }
    if (strips) {
      for (int j = 0; j < nstrips; ++j) {
        double cx = (sx2[j] + sx1[j]) / 2,
            cy = (sy2[j] + sy1[j]) / 2,
            cz = (sz2[j] + sz1[j]) / 2;
        double dx = std::abs(sx2[j] - cx),
            dy = std::abs(sy2[j] - cy),
            dz = std::abs(sz2[j] - cz);
        auto *box = new TMarker3DBox(
            cx, cy, cz,
            dx, dy, dz,
            0.0, 0.0
        );
        box->Draw();
        //auto *strip = new TPolyLine3D(2, xs, ys, zs);
        //strip->Draw();
      }
    }
  }
  hits->SetMarkerStyle(kFullDotMedium);
  hits->SetMarkerColor(kGreen + 3);
  hits->Draw();
  inAds->SetMarkerStyle(kFullDotSmall);
  inAds->SetMarkerColor(kBlue);
  inAds->Draw();
  // Axes
  //auto axes = new TAxis3D();
  //axes->Draw();
}
