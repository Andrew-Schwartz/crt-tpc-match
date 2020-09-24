#pragma clang diagnostic push
#pragma ide diagnostic ignored "Simplify" // IDE was angry that OUTPUT is const and wanted to simplify ifs with it

#include <iostream>
#include <TPolyLine3D.h>
#include <TMarker3DBox.h>
#include "utils.h"

/// TPC wire geometry with identifying information
struct Wire {
  int channel;
  int cryo, tpc, plane, wire;
  double x1, y1, z1;
  double x2, y2, z2;

  void draw() const {
    if constexpr (OUTPUT != Output::Draw) return;
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

/// Represents a reconstructed hit on a TPC wire at some tick time
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

  friend std::ostream &operator<<(std::ostream &os, const WireHit &hit) {
    os << "channel: " << hit.channel << ", cryo: " << hit.cryo << ", tpc: " << hit.tpc << ", plane: " << hit.plane
       << ", wire: " << hit.wire << ", peakTick: " << hit.peakTick;
    return os;
  }
};

//todo remove this?
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

/// Represents a CRT hit at some time (ns or tick)
struct CRTHit {
  double x, y, z;
  /// 1 tick == 1 ns
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

// CRT strip geometry with identifying id
struct CRTStrip {
  int id;
  double x1, y1, z1;
  double x2, y2, z2;

  void draw() const {
    if constexpr (OUTPUT != Output::Draw) return;
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
    if constexpr (OUTPUT != Output::Draw) return;
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

  /// Projects a CRTTrack with top and mid points onto the x-z plane at y = y3
  /// Alg from https://rosettacode.org/wiki/Find_the_intersection_of_a_line_with_a_plane
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

// todo is this needed?
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

// todo is this needed?
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

/// load wire geometry
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

/// load crt strip geometry
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

#pragma clang diagnostic pop