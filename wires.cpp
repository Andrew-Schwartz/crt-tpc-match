#pragma clang diagnostic push
#pragma ide diagnostic ignored "Simplify"

#include <vector>
#include <fstream>
#include <iostream>
#include <cmath>
#include <TTree.h>
#include <TPolyLine3D.h>
#include <TFile.h>
#include <TPolyMarker3D.h>
#include <TCanvas.h>
#include <TVector3.h>
#include <TMarker3DBox.h>
#include <unordered_set>
#include <numeric>

#pragma clang diagnostic push
#pragma ide diagnostic ignored "OCSimplifyInspection"

template<class T>
using hashset = std::unordered_set<T>;

constexpr bool FLIP = true;
constexpr size_t CHITS = 100;
constexpr size_t WHITS = 30000;
constexpr double ytop = 866.42;
constexpr double ymid = 662.42;
constexpr double ybot = -358.9;

bool near_top(double y, double eps = 0.1) {
  return std::abs(y - ytop) < eps;
}

bool near_mid(double y, double eps = 0.1) {
  return std::abs(y - ymid) < eps;
}

bool near_bot(double y, double eps = 0.1) {
  return std::abs(y - ybot) < eps;
}

void flip(double *ys, double *zs, int len) {
  for (int i = 0; i < len; ++i) {
    double temp = ys[i];
    ys[i] = -zs[i];
    zs[i] = temp;
  }
}

struct Wire {
  int channel;
  int cryo, tpc, plane, wire;
  double x1, y1, z1;
  double x2, y2, z2;

  [[nodiscard]] double len() const {
    return std::hypot(x2 - x1, y2 - y1, z2 - z1);
  }

  void draw() const {
    double xs[] {x1, x2},
        ys[] {y1, y2},
        zs[] {z1, z2};
    if constexpr (FLIP) {
      flip(ys, zs, 2);
    }
    auto *line = new TPolyLine3D(2, xs, ys, zs);
    line->Draw();
  }

  friend std::ostream &operator<<(std::ostream &os, const Wire &wire) {
    os << "channel: " << wire.channel << ", cryo: " << wire.cryo << ", tpc: " << wire.tpc << ", plane: "
       << wire.plane << ", wire: " << wire.wire << ", x1: " << wire.x1 << ", y1: " << wire.y1 << ", z1: "
       << wire.z1 << ", x2: " << wire.x2 << ", y2: " << wire.y2 << ", z2: " << wire.z2;
    return os;
  }
};

struct WireHit {
  int channel;
  int cryo, tpc, plane, wire;
  double peakTick;

  WireHit(
      int channel,
      int cryo,
      int tpc,
      int plane,
      int wire,
      double peakTick
  ) : channel(channel), cryo(cryo), tpc(tpc), plane(plane), wire(wire), peakTick(peakTick) {}

  [[nodiscard]] double peakNs() const {
    return peakTick * 500;
  }

  friend std::ostream &operator<<(std::ostream &os, const WireHit &hit) {
    os << "channel: " << hit.channel << ", cryo: " << hit.cryo << ", tpc: " << hit.tpc << ", plane: " << hit.plane
       << ", wire: " << hit.wire << ", peakTick: " << hit.peakTick/* << ", peakNs: " << hit.peakNs*/;
    return os;
  }
};

template<class T>
inline void hash_combine(std::size_t &s, const T &v) {
  std::hash<T> h;
  s ^= h(v) + 0x9e3779b9 + (s << 6) + (s >> 2);
}

namespace std {
  template<>
  struct hash<WireHit> {
    std::size_t operator()(WireHit const &hit) const noexcept {
      std::size_t ret = 0;
      hash_combine(ret, hit.channel);
      hash_combine(ret, hit.cryo);
      hash_combine(ret, hit.tpc);
      hash_combine(ret, hit.plane);
      hash_combine(ret, hit.wire);
      hash_combine(ret, hit.peakTick);
      return ret;
    }
  };
}

struct CRTHit {
  double x, y, z;
  /// tick and ns;
  double t;

  CRTHit(double x, double y, double z, double t) : x(x), y(y), z(z), t(t) {}

  friend std::ostream &operator<<(std::ostream &os, const CRTHit &hit) {
    os << "x: " << hit.x << ", y: " << hit.y << ", z: " << hit.z << ", t: " << hit.t;
    return os;
  }
};

namespace std {
  template<>
  struct hash<CRTHit> {
    std::size_t operator()(CRTHit const &hit) const noexcept {
      std::size_t ret = 0;
      hash_combine(ret, hit.x);
      hash_combine(ret, hit.y);
      hash_combine(ret, hit.z);
      hash_combine(ret, hit.t);
      return ret;
    }
  };
}

struct CRTStrip {
  int id;
  double x1, y1, z1;
  double x2, y2, z2;

  void draw() const {
    double xs[] {x1, x2},
        ys[] {y1, y2},
        zs[] {z1, z2};
    if constexpr (FLIP) {
      flip(ys, zs, 2);
    }
    double cx = (xs[0] + xs[1]) / 2.0,
        cy = (ys[0] + ys[1]) / 2.0,
        cz = (zs[0] + zs[1]) / 2.0;
    double hlx = std::abs(xs[0] - xs[1]) / 2.0,
        hly = std::abs(ys[0] - ys[1]) / 2.0,
        hlz = std::abs(zs[0] - zs[1]) / 2.0;
    auto *box = new TMarker3DBox(
        cx, cy, cz,
        hlx, hly, hlz,
        0.0, 0.0
    );
    box->Draw();
  }

  friend std::ostream &operator<<(std::ostream &os, const CRTStrip &strip) {
    os << "id: " << strip.id << ", x1: " << strip.x1 << ", y1: " << strip.y1 << ", z1: " << strip.z1 << ", x2: "
       << strip.x2 << ", y2: " << strip.y2 << ", z2: " << strip.z2;
    return os;
  }
};

struct CRTTrack {
  // top
  double x1, y1, z1;
  // middle
  double x2 {}, y2 {}, z2 {};
  // bottom
  double x3 {}, y3 {}, z3 {};

  CRTTrack(double x1, double y1, double z1) : x1(x1), y1(y1), z1(z1) {}

  CRTTrack(double x1,
           double y1,
           double z1,
           double x2,
           double y2,
           double z2
  ) : x1(x1), y1(y1), z1(z1), x2(x2), y2(y2), z2(z2) {}

  CRTTrack(double x1,
           double y1,
           double z1,
           double x2,
           double y2,
           double z2,
           double x3,
           double y3,
           double z3
  ) : x1(x1), y1(y1), z1(z1), x2(x2), y2(y2), z2(z2), x3(x3), y3(y3), z3(z3) {}

  CRTTrack(CRTTrack top, double x3, double y3, double z3)
      : x1(top.x1), y1(top.y1), z1(top.z1), x2(top.x2), y2(top.y2), z2(top.z2), x3(x3), y3(y3), z3(z3) {}

  void draw(Color_t color = kBlack) const {
    TPolyLine3D *line;
    if (x3 == 0 && y3 == 0 && z3 == 0) {
      double xs[] {x1, x2},
          ys[] {y1, y2},
          zs[] {z1, z2};
      if constexpr (FLIP) {
        flip(ys, zs, 2);
      }
      line = new TPolyLine3D(2, xs, ys, zs);
    } else if (x2 == 0 && y2 == 0 && z2 == 0) {
      double xs[] {x1, x3},
          ys[] {y1, y3},
          zs[] {z1, z3};
      if constexpr (FLIP) {
        flip(ys, zs, 2);
      }
      line = new TPolyLine3D(2, xs, ys, zs);
    } else {
      double xs[] {x1, x2, x3},
          ys[] {y1, y2, y3},
          zs[] {z1, z2, z3};
      if constexpr (FLIP) {
        flip(ys, zs, 3);
      }
      line = new TPolyLine3D(3, xs, ys, zs);
    }
    line->SetLineColor(color);
    line->Draw();
  }

  // alg from https://rosettacode.org/wiki/Find_the_intersection_of_a_line_with_a_plane
  [[nodiscard]] TVector3 projectY(double y3) const {
    TVector3 trackDir(x1 - x2, y1 - y2, z1 - z2),
        trackPt(x2, y2, z2),
        planeNorm(0.0, 1.0, 0.0),
        planePt(0.0, y3, 0.0);
    TVector3 diff = trackPt - planePt;
    double prod1 = diff.Dot(planeNorm);
    double prod2 = trackDir.Dot(planeNorm);
    double prod3 = prod1 / prod2;
    return trackPt - trackDir * prod3;
  }

  [[nodiscard]] bool sameStart(const CRTTrack &other) const {
    return x1 == other.x1 &&
           y1 == other.y1 &&
           z1 == other.z1 &&
           x2 == other.x2 &&
           y2 == other.y2 &&
           z2 == other.z2;
  }

  friend std::ostream &operator<<(std::ostream &os, const CRTTrack &track) {
    os << "x1: " << track.x1 << ", y1: " << track.y1 << ", z1: " << track.z1 << ", x2: " << track.x2 << ", y2: "
       << track.y2 << ", z2: " << track.z2 << ", x3: " << track.x3 << ", y3: " << track.y3 << ", z3: " << track.z3;
    return os;
  }
};

namespace std {
  template<>
  struct hash<CRTTrack> {
    std::size_t operator()(CRTTrack const &track) const noexcept {
      std::size_t ret = 0;
      hash_combine(ret, track.x1);
      hash_combine(ret, track.y1);
      hash_combine(ret, track.z1);
      hash_combine(ret, track.x2);
      hash_combine(ret, track.y2);
      hash_combine(ret, track.z2);
      hash_combine(ret, track.x3);
      hash_combine(ret, track.y3);
      hash_combine(ret, track.z3);
      return ret;
    }
  };
}

std::optional<TVector3> intersection(Wire plane2, Wire plane1, Wire plane0);

std::vector<Wire> parse_wires();

std::vector<CRTStrip> parse_strips();

void
tpcMatch(const std::vector<Wire> &wires, const std::vector<WireHit> &tcpPlane, std::vector<TVector3> &intersects,
         TPolyMarker3D *intersectMarks, int n);

void addPoint(TPolyMarker3D *pts, const TVector3 &vector);

TVector3 mean(const std::vector<TVector3> &pts);

double stddev(const std::vector<TVector3> &pts, const std::optional<TVector3> &omean);

std::string vecToStr(const TVector3 &vec);

void wires(int n = -1) {
  auto *c1 = new TCanvas("c1", "c1"); // to make root not print that this is created
  auto *file = TFile::Open("hitdumper_tree.root");
  auto *hitdumper = (TDirectoryFile *) file->Get("hitdumper");
  hitdumper->cd();
  auto *tree = (TTree *) hitdumper->Get("hitdumpertree");

  std::vector<Wire> wires = parse_wires();
  std::vector<CRTStrip> strips = parse_strips();

  // Wire Hits
  int channel[WHITS], cryo[WHITS], tpc[WHITS], plane[WHITS], wire[WHITS];
  double peakT[WHITS];
  int nwhits;
  tree->SetBranchAddress("hit_channel", &channel);
  tree->SetBranchAddress("hit_cryostat", &cryo);
  tree->SetBranchAddress("hit_tpc", &tpc);
  tree->SetBranchAddress("hit_plane", &plane);
  tree->SetBranchAddress("hit_wire", &wire);
  tree->SetBranchAddress("hit_peakT", &peakT);
  tree->SetBranchAddress("nhits", &nwhits);

  // CRT Hits
  double x[CHITS], y[CHITS], z[CHITS], t[CHITS];
  int nchits;
  tree->SetBranchAddress("chit_x", &x);
  tree->SetBranchAddress("chit_y", &y);
  tree->SetBranchAddress("chit_z", &z);
  tree->SetBranchAddress("chit_time", &t);
  tree->SetBranchAddress("nchits", &nchits);

  // draw hits
  auto *chitMarker = new TPolyMarker3D(CHITS);
  long low = (n == -1) ? 0 : n;
  long high = (n == -1) ? tree->GetEntries() : n + 1;
  for (long i = low; i < high; ++i) {
    tree->GetEntry(i);
    if (nchits > CHITS) {
      std::cout << "TOO FEW MAX CRT HITS, only " << CHITS << " of " << nchits << " hits are being used" << std::endl;
      nchits = CHITS;
    }
    if (nwhits > WHITS) {
      std::cout << "TOO FEW MAX WIRE HITS, only " << WHITS << " of " << nwhits << " hits are being used" << std::endl;
      nwhits = WHITS;
    }

    std::vector<CRTHit> chits;
    std::vector<WireHit> whits;
    for (int j = 0; j < nchits; ++j) {
      chits.emplace_back(x[j], y[j], z[j], t[j]);
      if constexpr (FLIP) {
        chitMarker->SetNextPoint(x[j], -z[j], y[j]);
      } else {
        chitMarker->SetNextPoint(x[j], y[j], z[j]);
      }
    }
    for (int j = 0; j < nwhits; ++j) {
      whits.emplace_back(channel[j], cryo[j], tpc[j], plane[j], wire[j], peakT[j]);
    }

    std::vector<CRTHit> top, mid, bot;
    for (CRTHit hit : chits) {
      if (near_top(hit.y)) {
        top.push_back(hit);
      } else if (near_mid(hit.y)) {
        mid.push_back(hit);
      } else {
        bot.push_back(hit);
      }
      /* else if (near_bot(hit.y)) {
        bot.push_back(hit);
      } else {
        std::cout << "Not on plane || to Y -> " << hit << std::endl;
      }*/
    }

//    std::vector<CRTTrack> tracks, maybeTracks;
//    hashset<CRTHit> usedHits;
//    if (!top.empty()) {
//      hashset<CRTHit> tmpUsed;
//      if (!mid.empty()) {
//        for (const auto &thit : top) {
//          // TODO: instead of just doing closest hit, should do (closest?) pair where the time and dist
//          //  between the points are the same (velocity)
//          std::sort(mid.begin(), mid.end(), [&](const CRTHit &a, const CRTHit &b) {
//            return std::hypot(a.x - thit.x, a.z - thit.z) < std::hypot(b.x - thit.x, b.z - thit.z);
//          });
//          double dist = std::hypot(thit.x - mid[0].x, thit.z - mid[0].z);
//          std::cout << "dist = " << dist << std::endl;
//
//          if (dist < 5) {
//            CRTTrack track(thit.x, thit.y, thit.z, mid[0].x, mid[0].y, mid[0].z);
//            tracks.push_back(track);
//            tmpUsed.insert({thit, mid[0]});
//          }
//        }
//      } else { // top !empty, mid empty
//        if (!bot.empty()) {
//          for (const auto &thit : top) {
//            // TODO: same as above
//            std::sort(bot.begin(), bot.end(), [&](const CRTHit &a, const CRTHit &b) {
//              return std::hypot(a.x - thit.x, a.z - thit.z) < std::hypot(b.x - thit.x, b.z - thit.z);
//            });
//            double dist = std::hypot(thit.x - bot[0].x, thit.z - bot[0].z);
//            std::cout << "dist = " << dist << std::endl;
//
//            if (dist < 30) {
//              CRTTrack track(thit.x, thit.y, thit.z, bot[0].x, bot[0].y, bot[0].z);
//              tracks.push_back(track);
//              tmpUsed.insert({thit, bot[0]});
//            }
//          }
//        }
//      }
//    } else { // top empty
//      if (!bot.empty()) {
//        for (const auto &thit : top) {
//          // TODO: same as above
//          std::sort(bot.begin(), bot.end(), [&](const CRTHit &a, const CRTHit &b) {
//            return std::hypot(a.x - thit.x, a.z - thit.z) < std::hypot(b.x - thit.x, b.z - thit.z);
//          });
//          double dist = std::hypot(thit.x - bot[0].x, thit.z - bot[0].z);
//          std::cout << "dist = " << dist << std::endl;
//
//          if (dist < 30) {
//            CRTTrack track(thit.x, thit.y, thit.z, bot[0].x, bot[0].y, bot[0].z);
//            tracks.push_back(track);
//            tmpUsed.insert({thit, bot[0]});
//          }
//        }
//      }
//    }

    std::vector<CRTTrack> tracks;
    if (!mid.empty()) {
      for (const auto &thit : top) {
        std::sort(mid.begin(), mid.end(), [&](const CRTHit &a, const CRTHit &b) {
          return std::hypot(a.x - thit.x, a.z - thit.z) < std::hypot(b.x - thit.x, b.z - thit.z);
        });
        double dist = std::hypot(thit.x - mid[0].x, thit.z - mid[0].z);
        std::cout << "dist = " << dist << std::endl;

        CRTTrack track(thit.x, thit.y, thit.z, mid[0].x, mid[0].y, mid[0].z);
        tracks.push_back(track);
      }
      auto *pts = new TPolyMarker3D();
      if (!bot.empty()) {
        for (int j = tracks.size(); j > 0; --j) {
          auto track = tracks[0];
          tracks.erase(tracks.begin());
          for (const auto &bhit : bot) {
            TVector3 proj = track.projectY(bhit.y);
            addPoint(pts, proj);
            TVector3 bvec(bhit.x, bhit.y, bhit.z);
            double dist = (bvec - proj).Mag();
            std::cout << "proj dist = " << dist << std::endl;
            if (dist < 150) { // TODO tune this a bit
              tracks.emplace_back(track, proj.x(), bhit.y, bhit.z);
            }
          }
        }
      }
      pts->SetMarkerStyle(kFullDotMedium);
      pts->SetMarkerColor(kGreen);
      pts->Draw();

//      for (const auto &track : tracks) {
//        track.draw();
//      }
    }/* else {
      auto *pts = new TPolyMarker3D();
      if (!bot.empty()) {
        for (const auto &thit : top) {
          std::sort(bot.begin(), bot.end(), [&](const CRTHit &a, const CRTHit &b) {
            return std::hypot(a.x - thit.x, a.z - thit.z) < std::hypot(b.x - thit.x, b.z - thit.z);
          });
          double dist = std::hypot(thit.x - bot[0].x, thit.z - bot[0].z);
          std::cout << "dist (b hit) = " << dist << std::endl;
          tracks.emplace_back(thit.x, thit.y, thit.z);

//          for (const auto &)
        }
      }
    }*/

    // Sort by time, with the vertical planes always first
    std::sort(whits.begin(), whits.end(), [](const WireHit &a, const WireHit &b) {
      if (a.peakTick < b.peakTick) {
        return true;
      } else if (a.peakTick == b.peakTick) {
        return a.plane > b.plane;
      } else return false;
    });

    std::vector<WireHit> tpc0, tpc1;
    for (WireHit whit : whits) {
      if (whit.tpc == 0) {
        tpc0.push_back(whit);
      } else {
        tpc1.push_back(whit);
      }
    }

    auto *intersectMarks = new TPolyMarker3D();
    std::vector<TVector3> intersects;
    tpcMatch(wires, tpc0, intersects, intersectMarks, n);
    tpcMatch(wires, tpc1, intersects, intersectMarks, n);
    std::cout << "num intersects = " << intersectMarks->Size() << std::endl;
    intersectMarks->SetMarkerStyle(kFullDotMedium);
    intersectMarks->SetMarkerColor(kRed);
    intersectMarks->Draw();

//    auto *trace = new TPolyMarker3D();
    std::vector<int> scores;
    for (const auto &track : tracks) {
      int score = 0;

      double slope = (track.z2 - track.z3) / (track.y2 - track.y3);
      int delta = 1;
      double planeX = (track.x1 > 0) ? 201.8 : -201.8;
      for (
          TVector3 pt {planeX, track.y2, track.z2};
          pt.Y() > -200; // -200 is the bottom of tcp plane
          pt = TVector3(pt.X(), pt.Y() - delta, pt.Z() - slope * delta)
          ) {
        if (pt.Y() > 200) continue; // 200 is top of tcp plane

        for (const auto &intersect : intersects) {
          double mag = (pt - intersect).Mag();
          int maxDist = 10;
          if (mag < maxDist) {
//            double y = intersect.Y(), z = intersect.Z();
//            flip(&y, &z, 1);
//            trace->SetNextPoint(intersect.X(), y, z);
            ++score;
          }
        }
      }

      std::cout << "score = " << score << std::endl;
      scores.push_back(score);
    }
//    trace->SetMarkerStyle(kFullDotMedium);
//    trace->SetMarkerColor(kGreen + 2);
//    std::cout << trace->Size() << std::endl;
//    trace->Draw();

    std::vector<int> bounds;
    bounds.push_back(0);
    for (int j = 1; j < tracks.size(); ++j) {
      if (!tracks[j - 1].sameStart(tracks[j])) {
        bounds.push_back(j);
      }
    }
    bounds.push_back(tracks.size());
    std::cout << "bounds = [";
    for (int bound : bounds) {
      std::cout << bound << ",";
    }
    std::cout << "]" << std::endl;
    for (int b = 0; b < bounds.size() - 1; ++b) {
      int maxIdx;
      int maxScore = std::numeric_limits<int>::min();
      for (int j = bounds[b]; j < bounds[b + 1]; ++j) {
        int score = scores[j];
        if (score >= maxScore) {
          maxScore = score;
          maxIdx = j;
        }
      }
      tracks[maxIdx].draw(kGreen + 4);
      std::cout << "maxIdx = " << maxIdx
                << ", track = " << tracks[maxIdx] << std::endl;
    }

    std::cout << std::endl;
  }
  chitMarker->SetMarkerStyle(kFullDotMedium);
  chitMarker->SetMarkerColor(kBlue);
  chitMarker->Draw();

//   draw wires
  for (int i = 0; i < wires.size(); i += 50) {
    wires[i].draw();
  }

  // draw strips
/*
  for (int i = 0; i < strips.size(); i += 3) {
    strips[i].draw();
  }
*/
}

void addPoint(TPolyMarker3D *pts, const TVector3 &vector) {
  if constexpr (FLIP) {
    pts->SetNextPoint(vector.x(), -vector.z(), vector.y());
  } else {
    pts->SetNextPoint(vector.x(), vector.y(), vector.z());
  }
}

void
tpcMatch(const std::vector<Wire> &wires, const std::vector<WireHit> &tcpPlane, std::vector<TVector3> &intersects,
         TPolyMarker3D *intersectMarks, int n) {
  std::vector<WireHit> plane2, plane1, plane0;
  for (WireHit whit : tcpPlane) {
    if (whit.plane == 2) {
      plane2.push_back(whit);
    } else if (whit.plane == 1) {
      plane1.push_back(whit);
    } else {
      plane0.push_back(whit);
    }
  }

  for (const WireHit &p2 : plane2) {
    double t2 = p2.peakTick;
    for (const WireHit &p1 : plane1) {
      double t1 = p1.peakTick;
      double diff = t2 - t1;
      if (diff < 3 || diff > 4) continue;
      for (const WireHit &p0 : plane0) {
        double t0 = p0.peakTick;
        double diff = t1 - t0;
        if (diff < 3 || diff > 4) continue;
        auto intersect = intersection(wires[p2.channel], wires[p1.channel], wires[p0.channel]);
        if (intersect) {
          intersects.push_back(*intersect);
          addPoint(intersectMarks, *intersect);
        }
      }
    }
  }
//  TVector3 center = mean(intersects);
//  double dist = stddev(intersects, std::optional(center));
//
//  std::cout << "center = " << vecToStr(center)
//            << ", dist = " << dist
//            << std::endl;
//
//  double pt[] {center.X(), center.Y(), center.Z()};
//  flip(&pt[1], &pt[2], 1);
//  auto pts = new TPolyMarker3D(1, pt);
//  pts->SetMarkerColor(kGreen);
//  pts->SetMarkerStyle(kFullDotLarge);
//  pts->Draw();
//  double xs[] {center.X(), center.X(), center.X(), center.X(), center.X()},
//      ys[] {center.Y() + dist, center.Y(), center.Y() - dist, center.Y(), center.Y() + dist},
//      zs[] {center.Z(), center.Z() + dist, center.Z(), center.Z() - dist, center.Z()};
//  flip(ys, zs, 5);
//  auto bound = new TPolyLine3D(5, xs, ys, zs);
//  bound->SetLineColor(kGreen);
//  bound->Draw();
}

std::optional<TVector3> intersection(Wire plane2, Wire plane1, Wire plane0) {
  if (plane2.plane != 2 || plane1.plane != 1 || plane0.plane != 0) {
    return std::optional<TVector3>();
  }

  double x = plane0.x1; // entire tpc-plane is same x

  // slope: always 60 degrees, have to pick + or -
  double m1 = 1.0 / sqrt(3.0);
  auto[p0y1, p0y2] = (plane0.z1 < plane0.z2)
                     ? std::make_pair(plane0.y1, plane0.y2)
                     : std::make_pair(plane0.y2, plane0.y1);
  if (p0y2 > p0y1) {
    m1 *= -1;
  }

  // intersection of wire 1 and wire 0
  double z = (m1 * (plane1.z1 + plane0.z1) + plane0.y1 - plane1.y1) / (2 * m1);

//  if (std::abs(z - plane2.z1) > 1) {
//    return std::optional<TVector3>();
//  }

//  plane2.draw();
//  plane1.draw();
//  plane0.draw();

  double y = m1 * (z - plane1.z1) + plane1.y1;

  // how different are the Z values??
//  double p[] {x, y, plane2.z1};
//  flip(&p[1], &p[2], 1);
//  auto *pt = new TPolyMarker3D(1, p);
//  pt->SetMarkerColor(kGreen);
//  pt->SetMarkerStyle(kFullDotMedium);
//  pt->Draw();

  return std::optional<TVector3>(TVector3(x, y, plane2.z1));
//  return std::optional<TVector3>(TVector3(x, y, z));
}

std::vector<Wire> parse_wires() {
  std::vector<Wire> wires;

  std::ifstream file("../WireDumpSBND.txt");
//    std::ifstream file(childDir.append("/..").c_str());
  std::string line;
  while (std::getline(file, line)) {
    Wire wire {};
    std::stringstream ss(line);
    for (int i = 0; i < 11; ++i) {
      switch (i) {
        case 0:
          ss >> wire.channel;
          break;
        case 1:
          ss >> wire.cryo;
          break;
        case 2:
          ss >> wire.tpc;
          break;
        case 3:
          ss >> wire.plane;
          break;
        case 4:
          ss >> wire.wire;
          break;
        case 5:
          ss >> wire.x1;
          break;
        case 6:
          ss >> wire.y1;
          break;
        case 7:
          ss >> wire.z1;
          break;
        case 8:
          ss >> wire.x2;
          break;
        case 9:
          ss >> wire.y2;
          break;
        case 10:
          ss >> wire.z2;
          break;
        default:
          break;
      }
    }
    wires.push_back(wire);
  }
  return wires;
}

std::vector<CRTStrip> parse_strips() {
  std::vector<CRTStrip> strips;

  std::ifstream file("../StripDumpSBND.txt");
  std::string line;
  while (std::getline(file, line)) {
    CRTStrip strip {};
    std::stringstream ss(line);
    for (int i = 0; i < 7; ++i) {
      switch (i) {
        case 0:
          ss >> strip.id;
          break;
        case 1:
          ss >> strip.x1;
          break;
        case 2:
          ss >> strip.y1;
          break;
        case 3:
          ss >> strip.z1;
          break;
        case 4:
          ss >> strip.x2;
          break;
        case 5:
          ss >> strip.y2;
          break;
        case 6:
          ss >> strip.z2;
          break;
        default:
          break;
      }
    }
    strips.push_back(strip);
  }
  return strips;
}

std::string vecToStr(const TVector3 &vec) {
  std::ostringstream s;
  s << "(" << vec.x() << ", " << vec.y() << ", " << vec.z() << ")";
  return s.str();
}

TVector3 mean(const std::vector<TVector3> &pts) {
  auto sum = pts[0];
  for (int i = 1; i < pts.size(); ++i) {
    sum += pts[i];
  }
  return sum * (1.0 / pts.size());
}

double stddev(const std::vector<TVector3> &pts, const std::optional<TVector3> &omean = std::nullopt) {
  TVector3 sum;
  auto avg = omean ? *omean : mean(pts);
  for (const auto &pt : pts) {
    const TVector3 &diff = pt - avg;
//    std::cout << "mag = " << diff.Mag() << std::endl;
    sum += TVector3(
        std::pow(diff.X(), 2),
        std::pow(diff.Y(), 2),
        std::pow(diff.Z(), 2)
    );
  }
  sum *= 1.0 / pts.size();
  return sum.Mag();
//  double x = std::sqrt(sum.X());
//  double y = std::sqrt(sum.Y());
//  double z = std::sqrt(sum.Z());
//  let stdev =
//  return TVector3()
}

#pragma clang diagnostic pop
#pragma clang diagnostic pop