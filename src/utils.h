#ifndef SCRATCH_UTILS_H
#define SCRATCH_UTILS_H

template<class T>
using hashset = std::unordered_set<T>;

template<class K, class V>
using hashmap = std::unordered_map<K, V>;

enum class Output {
  Draw, Histogram
};

// config variables
constexpr Output OUTPUT = Output::Histogram;
constexpr bool FLIP = true;
constexpr int CHITS = 100;
constexpr int WHITS = 30000;
constexpr size_t NHeightBins = 50;

// detector geometry
constexpr double ytop = 866.42;
constexpr double ymid = 662.42;
constexpr double ybot = -358.9;

/// rotate coordinates so that display in ROOT has Y as vertical
void flip(double *ys, double *zs, int len) {
  for (int i = 0; i < len; ++i) {
    double temp = ys[i];
    ys[i] = -zs[i];
    zs[i] = temp;
  }
}

/// from https://stackoverflow.com/questions/2590677/how-do-i-combine-hash-values-in-c0x#2595226
template<class T>
inline void hash_combine(std::size_t &s, const T &v) {
  std::hash<T> h;
  s ^= h(v) + 0x9e3779b9 + (s << 6) + (s >> 2);
}

bool epsEqual(double n1, double n2, double eps = 0.01) {
  return std::abs(n1 - n2) < eps;
}

/// mean of a list of data
template<class T>
T mean(const std::vector<T> &data) {
  auto sum = data[0];
  for (size_t i = 0; i < data.size(); ++i) {
    sum += data[i];
  }
  // multiplying by reciprocal instead of dividing to support TVector3
  return sum * (1.0 / data.size());
}

double stddev(const std::vector<TVector3> &pts, const std::optional<TVector3> &omean = std::nullopt) {
  TVector3 sum;
  auto avg = omean ? *omean : mean(pts);
  for (const auto &pt : pts) {
    TVector3 diff = pt - avg;
    sum += TVector3(
        std::pow(diff.X(), 2),
        std::pow(diff.Y(), 2),
        std::pow(diff.Z(), 2)
    );
  }
  sum *= 1.0 / pts.size();
  return sum.Mag();
}

double stddev(const std::vector<double> &data, const std::optional<double> &omean = std::nullopt) {
  double sum = 0;
  auto avg = omean ? *omean : mean(data);
  for (const auto &datum : data) {
    sum += std::pow(datum - avg, 2);
  }
  sum /= data.size();
  return std::sqrt(sum);
}

/// sometimes you just want to print out a TVector3
std::string vecToStr(const TVector3 &vec) {
  std::ostringstream s;
  s << "(" << vec.x() << ", " << vec.y() << ", " << vec.z() << ")";
  return s.str();
}

/// convenience method to add a point to a TPolyMarker3d and flip axis if needed
void addPoint(TPolyMarker3D *to, const TVector3 &from) {
  if constexpr (FLIP) {
    to->SetNextPoint(from.x(), -from.z(), from.y());
  } else {
    to->SetNextPoint(from.x(), from.y(), from.z());
  }
}

#endif //SCRATCH_UTILS_H