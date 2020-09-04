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
#include <TH2D.h>

#pragma clang diagnostic push
#pragma ide diagnostic ignored "OCSimplifyInspection"

template<class T>
using hashset = std::unordered_set<T>;

template<class K, class V>
using hashmap = std::unordered_map<K, V>;

constexpr bool DRAW = true;
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

/// rotate coordinates so that display in ROOT has Y as vertical
void flip(double *ys, double *zs, int len) {
  for (int i = 0; i < len; ++i) {
    double temp = ys[i];
    ys[i] = -zs[i];
    zs[i] = temp;
  }
}

bool epsEqual(double n1, double n2, double eps = 0.01) {
  return std::abs(n1 - n2) < eps;
}

struct Wire {
  int channel;
  int cryo, tpc, plane, wire;
  double x1, y1, z1;
  double x2, y2, z2;

  void draw() const {
    if constexpr (!DRAW) return;
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
       << wire.plane << ", wire: " << wire.wire << ", xt: " << wire.x1 << ", yt: " << wire.y1 << ", zt: "
       << wire.z1 << ", xm: " << wire.x2 << ", ym: " << wire.y2 << ", zm: " << wire.z2;
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

/// from https://stackoverflow.com/questions/2590677/how-do-i-combine-hash-values-in-c0x#2595226
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
    if constexpr (!DRAW) return;
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
    os << "id: " << strip.id << ", xt: " << strip.x1 << ", yt: " << strip.y1 << ", zt: " << strip.z1 << ", xm: "
       << strip.x2 << ", ym: " << strip.y2 << ", zm: " << strip.z2;
    return os;
  }
};

/// reconstructed CRTTrack, will always have at least two points but handles any pair of two,
/// including uncertainty in xb due to the single plane of strips on bottom
struct CRTTrack {
  // top
  double xt, yt, zt;
  // middle
  double xm, ym, zm;
  // bottom
  double xb, yb, zb;

  CRTTrack(double xt, double yt, double zt,
           double xm, double ym, double zm,
           double xb, double yb, double zb
  ) : xt(xt), yt(yt), zt(zt), xm(xm), ym(ym), zm(zm), xb(xb), yb(yb), zb(zb) {}

  CRTTrack(CRTTrack top, double xb, double yb, double zb)
      : xt(top.xt), yt(top.yt), zt(top.zt), xm(top.xm), ym(top.ym), zm(top.zm), xb(xb), yb(yb), zb(zb) {}

  [[nodiscard]] std::optional<TVector3> top() const {
    if (xt == 0 && yt == 0 && zt == 0) {
      return std::nullopt;
    } else {
      return TVector3(xt, yt, zt);
    }
  }

  [[nodiscard]] std::optional<TVector3> mid() const {
    if (xm == 0 && ym == 0 && zm == 0) {
      return std::nullopt;
    } else {
      return TVector3(xm, ym, zm);
    }
  }

  [[nodiscard]] std::optional<TVector3> bot() const {
    if (xb == 0 && yb == 0 && zb == 0) {
      return std::nullopt;
    } else {
      return TVector3(xb, yb, zb);
    }
  }

  [[nodiscard]] int numPts() const {
    int n = 0;
    if (top()) ++n;
    if (mid()) ++n;
    if (bot()) ++n;
    return n;
  }

  static CRTTrack topMid(double xt, double yt, double zt,
                         double xm, double ym, double zm) {
    return CRTTrack(xt, yt, zt, xm, ym, zm, 0, 0, 0);
  }

  static CRTTrack topBot(double xt, double yt, double zt,
                         double xb, double yb, double zb) {
    return CRTTrack(xt, yt, zt, 0, 0, 0, xb, yb, zb);
  }

  static CRTTrack midBot(double xm, double ym, double zm,
                         double xb, double yb, double zb) {
    return CRTTrack(0, 0, 0, xm, ym, zm, xb, yb, zb);
  }

  void draw(Color_t color = kBlack) const {
    if constexpr (!DRAW) return;
    TPolyLine3D *line;
    if (xb == 0 && yb == 0 && zb == 0) {
      double xs[] {xt, xm},
          ys[] {yt, ym},
          zs[] {zt, zm};
      if constexpr (FLIP) {
        flip(ys, zs, 2);
      }
      line = new TPolyLine3D(2, xs, ys, zs);
    } else if (xm == 0 && ym == 0 && zm == 0) {
      double xs[] {xt, xb},
          ys[] {yt, yb},
          zs[] {zt, zb};
      if constexpr (FLIP) {
        flip(ys, zs, 2);
      }
      line = new TPolyLine3D(2, xs, ys, zs);
    } else if (xt == 0 && yt == 0 && zt == 0) {
      double xs[] {xm, xb},
          ys[] {ym, yb},
          zs[] {zm, zb};
      if constexpr (FLIP) {
        flip(ys, zs, 2);
      }
      line = new TPolyLine3D(2, xs, ys, zs);
    } else {
      double xs[] {xt, xm, xb},
          ys[] {yt, ym, yb},
          zs[] {zt, zm, zb};
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
    TVector3 trackDir(xt - xm, yt - ym, zt - zm),
        trackPt(xm, ym, zm),
        planeNorm(0.0, 1.0, 0.0),
        planePt(0.0, y3, 0.0);
    TVector3 diff = trackPt - planePt;
    double prod1 = diff.Dot(planeNorm);
    double prod2 = trackDir.Dot(planeNorm);
    double prod3 = prod1 / prod2;
    return trackPt - trackDir * prod3;
  }

  [[nodiscard]] bool topEq(const CRTTrack &other) const {
    return epsEqual(xt, other.xt) &&
        epsEqual(yt, other.yt) &&
        epsEqual(zt, other.zt);
  }

  [[nodiscard]] bool midEq(const CRTTrack &other) const {
    return epsEqual(xm, other.xm) &&
        epsEqual(ym, other.ym) &&
        epsEqual(zm, other.zm);
  }

  [[nodiscard]] bool botEq(const CRTTrack &other) const {
    // ignore x because the hit doesn't tell us x
    return epsEqual(yb, other.yb) &&
        epsEqual(zb, other.zb);
  }

  [[nodiscard]] TVector3 topmostPt() const {
    return (xt == 0 && yt == 0 && zt == 0)
           ? TVector3(xm, ym, zm)
           : TVector3(xt, yt, zt);
  }

  friend std::ostream &operator<<(std::ostream &os, const CRTTrack &track) {
    os << "xt: " << track.xt << ", yt: " << track.yt << ", zt: " << track.zt << ", xm: " << track.xm << ", ym: "
        << track.ym << ", zm: " << track.zm << ", xb: " << track.xb << ", yb: " << track.yb << ", zb: " << track.zb;
    return os;
  }

  [[nodiscard]] bool containsPt(const CRTHit &hit) const {
    return (epsEqual(hit.x, xt) && epsEqual(hit.y, yt) && epsEqual(hit.z, zt)) ||
        (epsEqual(hit.x, xm) && epsEqual(hit.y, ym) && epsEqual(hit.z, zm)) ||
        (epsEqual(hit.x, xb) && epsEqual(hit.y, yb) && epsEqual(hit.z, zb));
  }

  bool operator==(const CRTTrack &rhs) const {
    return xt == rhs.xt &&
        yt == rhs.yt &&
        zt == rhs.zt &&
        xm == rhs.xm &&
        ym == rhs.ym &&
        zm == rhs.zm &&
        xb == rhs.xb &&
        yb == rhs.yb &&
        zb == rhs.zb;
  }

  bool operator!=(const CRTTrack &rhs) const {
    return !(rhs == *this);
  }
};

namespace std {
  template<>
  struct hash<CRTTrack> {
    std::size_t operator()(const CRTTrack &track) const noexcept {
      std::size_t ret = 0;
      hash_combine(ret, track.xt);
      hash_combine(ret, track.yt);
      hash_combine(ret, track.zt);
      hash_combine(ret, track.xm);
      hash_combine(ret, track.ym);
      hash_combine(ret, track.zm);
      hash_combine(ret, track.xb);
      hash_combine(ret, track.yb);
      hash_combine(ret, track.zb);
      return ret;
    }
  };
}

namespace std {
  template<>
  struct hash<TVector3> {
    std::size_t operator()(const TVector3 &vec) const noexcept {
      std::size_t ret = 0;
      hash_combine(ret, vec.X());
      hash_combine(ret, vec.Y());
      hash_combine(ret, vec.Z());
      return ret;
    }
  };
}

/// coordinate where the three wires intersect
std::optional<TVector3> intersection(Wire plane2, Wire plane1, Wire plane0);

/// load wire geometry
std::vector<Wire> parse_wires();

/// load crt strip geometry
std::vector<CRTStrip> parse_strips();

/// finds all intersections points of WireHits on one TPC wire plane
void wireHitIntersections(const std::vector<Wire> &wires, const std::vector<WireHit> &tpcPlane,
                          std::vector<TVector3> &intersects,
                          TPolyMarker3D *intersectMarks);

void addPoint(TPolyMarker3D *pts, const TVector3 &vector);

double mean(const std::vector<double> &data);

TVector3 mean(const std::vector<TVector3> &pts);

double stddev(const std::vector<double> &data, const std::optional<double> &omean);

double stddev(const std::vector<TVector3> &pts, const std::optional<TVector3> &omean);

std::string vecToStr(const TVector3 &vec);

/// main method. n is the event to look at, or -1 to look at all events
/// makes some assumptions about directory structure: the dir this is run from has a `hitdumper_tree.root`,
/// and the dir above this has `WireDumpSBND.txt` and `StripDumpSBND.txt`
void wires(int n) {
  // initialize ROOT stuff
  auto *c1 = new TCanvas("c1", "c1"); // to make root not print that this is created
  auto *file = TFile::Open("hitdumper_tree.root");
  auto *hitdumper = (TDirectoryFile *) file->Get("hitdumper");
  hitdumper->cd();
  auto *tree = (TTree *) hitdumper->Get("hitdumpertree");

  // load geometry
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

  int crtHitTotal = 0;
  int crtHitMatches = 0;

  // draw crt hits
  auto *chitMarker = new TPolyMarker3D(CHITS);
  long low = (n == -1) ? 0 : n;
  long high = (n == -1) ? tree->GetEntries() : n + 1;
  for (long i = low; i < high; ++i) {
    tree->GetEntry(i);
    if (nchits > CHITS) {
      std::cout << "TOO FEW MAX CRT HITS, only " << CHITS << " of " << nchits << " hits are being used" << std::endl;
      nchits = CHITS;
    }
    crtHitTotal += nchits;
    if (nwhits > WHITS) {
      std::cout << "TOO FEW MAX WIRE HITS, only " << WHITS << " of " << nwhits << " hits are being used" << std::endl;
      nwhits = WHITS;
    }

    // read CRT/Wire hits from ROOT tree
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

    // top and mid are the planes above detector, bot is below
    std::vector<CRTHit> top, mid, bot;
    for (CRTHit hit : chits) {
      if (near_top(hit.y)) {
        top.push_back(hit);
      } else if (near_mid(hit.y)) {
        mid.push_back(hit);
      } else {
        bot.push_back(hit);
      }
    }

    // reconstruct CRTTracks (neither of the two in the ROOT tree worked for this, one because it was always empty)
    std::vector<CRTTrack> tracks;
    for (const auto &thit : top) {
      for (const auto &mhit : mid) {
        tracks.push_back(CRTTrack::topMid(thit.x, thit.y, thit.z, mhit.x, mhit.y, mhit.z));
      }
    }
    // the bottom plane doesn't have x coord because of it only has one plane of strips, so the actual
    // x coordinate is projected down from the top planes. Replaces all top->mid tracks with (potentially multiple)
    // top->mid->bot tracks.
    auto *projPts = new TPolyMarker3D();
    if (!bot.empty()) {
      for (int j = tracks.size(); j > 0; --j) {
        auto track = tracks[0];
        tracks.erase(tracks.begin());
        for (const auto &bhit : bot) {
          TVector3 proj = track.projectY(bhit.y);
          TVector3 bvec(bhit.x, bhit.y, bhit.z);
          double dist = (bvec - proj).Mag();
          std::cout << "proj dist = " << dist << std::endl;
          bool sameSide = std::signbit(proj.x()) == std::signbit(bhit.x);
          if (sameSide && dist < 600) { // TODO tune this a bit (maybe instead check that proj is inside CRT?)
            addPoint(projPts, TVector3(proj.x(), bhit.y, bhit.z));
            tracks.emplace_back(track, proj.x(), bhit.y, bhit.z);
          }
        }
      }
    }
    projPts->SetMarkerStyle(kFullDotMedium);
    projPts->SetMarkerColor(kGreen);
    if constexpr (DRAW) {
      projPts->Draw();
    }

    unsigned long firstIncomplete = tracks.size();
    // every permutation of incomplete top->bot and mid->bot tracks
    for (const auto &bhit : bot) {
      for (const auto &thit : top) {
        tracks.push_back(CRTTrack::topBot(thit.x, thit.y, thit.z, bhit.x, bhit.y, bhit.z));
      }
      for (const auto &mhit : mid) {
        tracks.push_back(CRTTrack::midBot(mhit.x, mhit.y, mhit.z, bhit.x, bhit.y, bhit.z));
      }
    }

    // remove duplicate top->bot or mid->bot tracks that are already part of a top->mid->bot track
    if (firstIncomplete != 0) {
      for (unsigned long j = tracks.size() - 1; j >= firstIncomplete; --j) {
        auto incomp = tracks[j];
        for (unsigned long k = 0; k < firstIncomplete; ++k) {
          auto comp = tracks[k];
          if (incomp.botEq(comp) && (incomp.topEq(comp) || incomp.midEq(comp))) {
            tracks.erase(tracks.begin() + j);
            break;
          }
        }
      }
    }
    std::cout << "#tracks = " << tracks.size() << std::endl;

    // Sort by time then by plane (vertical plane first)
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
    wireHitIntersections(wires, tpc0, intersects, intersectMarks);
    wireHitIntersections(wires, tpc1, intersects, intersectMarks);
    std::cout << "num intersects = " << intersectMarks->Size() << std::endl;
    intersectMarks->SetMarkerStyle(kFullDotMedium);
    intersectMarks->SetMarkerColor(kRed);
    if constexpr (DRAW) {
      intersectMarks->Draw();
    }

    // score each track: walk down track, project it onto the wire plane, and count nearby wire hits
    hashmap<CRTTrack, double> scores;
    for (const auto &track : tracks) {
      double score = 0;
      int tot = 0;

      hashset<TVector3> usedIntersects;
      auto topmost = track.topmostPt();
      double slopeZ = (topmost.z() - track.zb) / (topmost.y() - track.yb);
      double slopeX = (topmost.x() - track.xb) / (topmost.y() - track.yb);
      int delta = 1;
      for (
          TVector3 pt = track.topmostPt();
          pt.Y() > -200; // -200 is the bottom of tcp plane
          pt -= TVector3(slopeX * delta, delta, slopeZ * delta)
          ) {
        if (pt.Y() > 200) continue; // 200 is top of tcp plane
        for (const auto &intersect : intersects) {
          if (std::signbit(intersect.X()) != std::signbit(pt.X())) {
            continue; // ensure intersect and track are at same wire plane
          }
          ++tot;
          auto projected = pt;
          projected.SetX(intersect.x());
          double mag = (projected - intersect).Mag();
          int maxDist = 5;
          if (mag < maxDist) {
            ++score;
            usedIntersects.insert(intersect);
          }
        }
      }

      if (tot == 0) {
        score = 0;
      } else {
        score /= tot;
      }
      std::cout << "score = " << score << std::endl;
      scores.emplace(track, score);
    }

    hashmap<CRTTrack, std::vector<TVector3>> trackHits; // unused
    // each point to a list of each track that includes that point
    hashmap<TVector3, std::vector<CRTTrack>> hitTracks;
    for (const CRTTrack track : tracks) {
      trackHits.emplace(track, std::vector<TVector3>());
      std::optional<TVector3> top = track.top();
      if (top) {
        trackHits[track].push_back(*top);
        if (hitTracks.find(*top) == hitTracks.end()) { // not in map
          hitTracks.insert({*top, std::vector<CRTTrack>()});
        }
        hitTracks[*top].push_back(track);
      }
      std::optional<TVector3> mid = track.mid();
      if (mid) {
        trackHits[track].push_back(*mid);
        if (hitTracks.find(*mid) == hitTracks.end()) { // not in map
          hitTracks.insert({*mid, std::vector<CRTTrack>()});
        }
        hitTracks[*mid].push_back(track);
      }
      std::optional<TVector3> bot = track.bot();
      if (bot) {
        trackHits[track].push_back(*bot);
        if (hitTracks.find(*bot) == hitTracks.end()) { // not in map
          hitTracks.insert({*bot, std::vector<CRTTrack>()});
        }
        hitTracks[*bot].push_back(track);
      }
    }

    // only the tracks that match a wire intersects "track"
    std::vector<CRTTrack> matches;
    for (const auto &hitTrack : hitTracks) {
      const auto&[hit, htracks] = hitTrack;
//      std::cout << "vecToStr(hit) = " << vecToStr(hit) << std::endl;
      int maxIdx;
      double maxScore = -1;
      // get best match for TVector3 track point
      for (int j = 0; j < htracks.size(); ++j) {
        double score = scores[htracks[j]];
//        std::cout << "htracks[j]: score = " << score << ",  = " << htracks[j] << std::endl;
        if (score >= maxScore) {
          maxScore = score;
          maxIdx = j;
        }
      }
      if (maxScore > 0.00001) {
        auto newTrack = htracks[maxIdx];
        std::vector<int> oldTrackIdxs;
        for (int j = 0; j < matches.size(); ++j) {
          auto track = matches[j];
          if (track.topEq(newTrack) || track.midEq(newTrack) || track.botEq(newTrack)) {
            oldTrackIdxs.push_back(j);
          }
        }
        if (!oldTrackIdxs.empty()) { // if new track shares point with any track already in matches, only keep the highest scored one
//          std::cout << "NOT NEW: " << oldTrackIdxs.size() << " overlaps!" << std::endl;
          bool newBest = true;
          for (unsigned long j = oldTrackIdxs.size() - 1; j > 0; --j) {
            int oldTrackIdx = oldTrackIdxs[j];
            auto oldTrack = matches[oldTrackIdx];
            double oldScore = scores[oldTrack];
            double newScore = scores[newTrack];
            if (newScore < oldScore) {
              newBest = false;
            } else {
//              std::cout << "REMOVING: " << oldTrackIdx << std::endl;
              matches.erase(matches.begin() + oldTrackIdx);
            }
          }
          if (newBest) { // erase oldTrack, put in newTrack
//            std::cout << "UPDATING " << std::endl;
            matches.push_back(newTrack);
          }
        } else {
//          std::cout << "NEW" << std::endl;
          matches.push_back(newTrack);
        }
      }
    }

    std::cout << "matches.size() = " << matches.size() << std::endl;
    for (const auto &match : matches) {
      match.draw();
    }

    int nmatch = 0;
    for (const auto &chit : chits) {
      for (const auto &match : matches) {
        if (match.containsPt(chit)) {
          ++nmatch;
        }
      }
    }
    crtHitMatches += nmatch;

    std::cout << std::endl;
  } // end of events for loop

  chitMarker->SetMarkerStyle(kFullDotMedium);
  chitMarker->SetMarkerColor(kBlue);
  if constexpr (DRAW) {
    chitMarker->Draw();
  }

  std::cout << "crtHitTotal = " << crtHitTotal << std::endl;
  std::cout << "crtHitMatches = " << crtHitMatches << std::endl;
  std::cout << "real% = " << ((double) crtHitMatches) / (crtHitTotal) * 100.0 << std::endl;

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

void wireHitIntersections(const std::vector<Wire> &wires, const std::vector<WireHit> &tpcPlane,
                          std::vector<TVector3> &intersects,
                          TPolyMarker3D *intersectMarks) {
  std::vector<WireHit> plane2, plane1, plane0;
  for (WireHit whit : tpcPlane) {
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

  double y = m1 * (z - plane1.z1) + plane1.y1;

  return std::optional<TVector3>(TVector3(x, y, plane2.z1));
//  return std::optional<TVector3>(TVector3(x, y, z));
}

std::vector<Wire> parse_wires() {
  std::vector<Wire> wires;

  std::ifstream file("../WireDumpSBND.txt");
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

double mean(const std::vector<double> &data) {
  auto sum = data[0];
  for (int i = 1; i < data.size(); ++i) {
    sum += data[i];
  }
  return sum * (1.0 / data.size());
}

double stddev(const std::vector<TVector3> &pts, const std::optional<TVector3> &omean = std::nullopt) {
  TVector3 sum;
  auto avg = omean ? *omean : mean(pts);
  for (const auto &pt : pts) {
    const TVector3 &diff = pt - avg;
    sum += TVector3(
        std::pow(diff.X(), 2),
        std::pow(diff.Y(), 2),
        std::pow(diff.Z(), 2)
    );
  }
  sum *= 1.0 / pts.size();
  return sum.Mag();
}

double stddev(const std::vector<double> &data, const std::optional<double> &omean) {
  double sum;
  auto avg = omean ? *omean : mean(data);
  for (const auto &datum : data) {
    sum += std::pow(datum - avg, 2);
  }
  sum /= data.size();
  return std::sqrt(sum);
}

#pragma clang diagnostic pop
#pragma clang diagnostic pop