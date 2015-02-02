#include "stdafx.h"
#include "vertex.hpp"
#include "mesh.hpp"
//#include <omp.h>
#include <fstream>
#include <map>
#include <set>
#include <vector>
#include <string>
#include <numeric>
#include <cmath>
#include <fstream>
#include <limits>
#include <algorithm>
#include <d3dx10.h>
#include <boost/tuple/tuple.hpp>
#include <boost/tuple/tuple_comparison.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/scoped_array.hpp>
#include <boost/format.hpp>
#include <boost/array.hpp>
#include <boost/utility.hpp>
#include <boost/foreach.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/triangular.hpp>
#include <boost/numeric/ublas/lu.hpp>

using namespace std;
using namespace boost;
using namespace boost::algorithm;
namespace ublas = boost::numeric::ublas;
#define foreach BOOST_FOREACH

Mesh::Mesh()
     : adjvert_()
{
}

inline float concentrate(float v, float center, float mag) {
  return (v - center) * mag + center;
}

#ifdef max
#undef max
#endif
#ifdef V
#undef V
#endif
template<typename In, typename Out>
void safe_inverse(const In& A_, Out& B) {
  typedef typename In::value_type value_type;
  ublas::matrix<value_type> OrigA(A_);
  ublas::matrix<value_type> A(A_);
  ublas::matrix<value_type> U(ublas::identity_matrix<value_type>(A_.size1()));

  value_type threshold = 0.0;
  for(int t = 0; t < 100; ++t) {
    int maxi = 0, maxj = 0;
    value_type maxval = 0;
    for(int i = 0; i < A.size1(); ++i) {
      for(int j = i+1; j < A.size2(); ++j) {
        if(fabs(A(i, j)) > maxval) {
          maxval = fabs(A(i, j));
          maxi = i;
          maxj = j;
        }
      }
    }

    if(t == 0) {
      threshold = maxval * 0.00001;
    }
    if(maxval <= threshold) {
      break;
    }

    ublas::matrix<value_type> V(ublas::identity_matrix<value_type>(A_.size1()));
    float theta;
    if(A(maxi, maxi) == A(maxj, maxj)) {
      theta = 3.141592 / 4;
    }
    else {
      theta = atan(-2*A(maxi, maxj) / (A(maxi, maxi) - A(maxj, maxj))) / 2; 
    }
    value_type c = cos(theta), s = sin(theta);
    V(maxi, maxi) = V(maxj, maxj) = c;
    V(maxi, maxj) =  s;
    V(maxj, maxi) = -s;

    U = ublas::prod(U, V);
    ublas::matrix<value_type> tVA(ublas::prod(ublas::trans(V), A));
    A = ublas::prod(tVA, V);
  }

  {
    value_type maxval = -numeric_limits<value_type>::max();
    int maxi;
    for(int i = 0; i < A.size1(); ++i) {
      if(A(i, i) > maxval) {
        maxval = A(i, i);
        maxi = i;
      }
    }
    for(int i = 0; i < A.size1(); ++i) {
      if(A(i, i) >= maxval * 0.0001) {
        A(i, i) = 1 / A(i, i);
      }
    }
    ublas::matrix<value_type> tUA(ublas::prod(A, ublas::trans(U)));
    B = ublas::prod(U, tUA);
  }
}

// sorted tuple
template<typename T, typename U>
tuple<int, int> make_stuple(const T& i, const U& j) {
  if(i > j) {
    return make_tuple(j, i);
  }
  else {
    return make_tuple(i, j);
  }
}

inline float sign(float v) {
  if(v < 0.f) { return -1.f; }
  else if(v > 0.f) { return 1.f; }
  else { return 0.f; }
}

inline float compute_dir_curv(float* tensor, const D3DXVECTOR2& v) {
  return tensor[0]*v.x*v.x + 2*tensor[1]*v.x*v.y + tensor[2]*v.y*v.y;
}
inline float compute_dir_curv_deriv(float* tensor, const D3DXVECTOR2& v) {
  return tensor[0]*v.x*v.x*v.x +
    3*tensor[1]*v.x*v.x*v.y + 3*tensor[2]*v.x*v.y*v.y +
    tensor[3]*v.y*v.y*v.y;
}

// parse obj file
void Mesh::load_obj(istream& ifs) {
  vert_.clear();
  face_.clear();

  string line;
  vector<string> tokens;
  D3DXVECTOR3 maxpos(-FLT_MAX, -FLT_MAX, -FLT_MAX), minpos(FLT_MAX, FLT_MAX, FLT_MAX);
  while(ifs) {
    getline(ifs, line);
    static int count = 0;

    // empty line, or comment line.
    if(line.empty()) { continue; }
    if(line[0] == '#') { continue; }

    // from boost::string_algo, splits input line by whitespace.
    split(tokens, line, boost::algorithm::is_space(), token_compress_on);

    if(tokens[0] == "g") {
      // ignore it!
    }
    else if(tokens[0] == "v") {
      vert_type v;
      D3DXVECTOR3 p(lexical_cast<float>(tokens[1]),
                    lexical_cast<float>(tokens[2]),
                    lexical_cast<float>(tokens[3]));
      v.p = p;
      maxpos = D3DXVECTOR3(max(maxpos.x, p.x), max(maxpos.y, p.y), max(maxpos.z, p.z));
      minpos = D3DXVECTOR3(min(minpos.x, p.x), min(minpos.y, p.y), min(minpos.z, p.z));
      vert_.push_back(v);
    }
    else if(tokens[0] == "f" && tokens.size() >= 3+1) {
      vector<int> vs;

      // parse verts
      for(int i = 1; i < tokens.size(); ++i) {
        vector<string> vinds;
        split(vinds, tokens[i], is_any_of("/"), token_compress_off);
        vs.push_back(lexical_cast<int>(vinds[0]) - 1);
      }

      // calc face's norm
      D3DXVECTOR3 norm;
      {
        const D3DXVECTOR3& p0 = vert_[vs[0]].p;
        const D3DXVECTOR3& p1 = vert_[vs[1]].p;
        const D3DXVECTOR3& p2 = vert_[vs[2]].p;
        const D3DXVECTOR3 v1 = p1 - p0;
        const D3DXVECTOR3 v2 = p2 - p1;
        D3DXVec3Cross(&norm, &v1, &v2);
        D3DXVec3Normalize(&norm, &norm);
      }

      // make triangles
      for(int i = 0; i <= vs.size() - 3; ++i) {
        face_type f;
        {
          f.vpos[0]  = vs[0];
        }
        for(int j = 1; j < 3; ++j) {
          f.vpos[j]  = vs[j + i];
        }
        f.norm = norm;
        face_.push_back(f);
      }
    }
  }

  D3DXVECTOR3 center = (maxpos + minpos) / 2;
  maxpos -= center;
  float scale = 1.f / max(maxpos.x, max(maxpos.y, maxpos.z));
  // scaling
  foreach(vert_type& v, vert_) {
    v.p -= center;
    v.p *= scale;
  }
}

void Mesh::load_ifs(istream& ifs) {
  struct {
    void operator()(istream& is, const char* expected = 0) const {
      int len;
      is.read((char*)&len, 4);
      vector<char> buf(len);
      is.read(&buf[0], len);
      if(expected && strcmp(expected, &buf[0]) != 0) {
        throw runtime_error(string("ifs error: expected ") + expected + ", but " + &buf[0]);
      }
    }
  } eat_string;
  struct {
    float operator()(istream& is) const {
      float v; is.read((char*)&v, 4); return v;
    }
  } read_float;
  struct {
    int operator()(istream& is) const {
      int v; is.read((char*)&v, 4); return v;
    }
  } read_int;

  eat_string(ifs, "IFS");
  read_float(ifs);
  eat_string(ifs);
  eat_string(ifs, "VERTICES");
  int vnum = read_int(ifs);
  vert_.reserve(vnum);
  for(int i = 0; i < vnum; ++i) {
    vert_type v;
    v.p.x = read_float(ifs);
    v.p.y = read_float(ifs);
    v.p.z = read_float(ifs);
    vert_.push_back(v);
  }
  eat_string(ifs, "TRIANGLES");
  int fnum = read_int(ifs);
  face_.reserve(fnum);
  for(int i = 0; i < fnum; ++i) {
    face_type f;
    f.vpos[0] = read_int(ifs);
    f.vpos[1] = read_int(ifs);
    f.vpos[2] = read_int(ifs);
    D3DXVECTOR3 norm;
    const D3DXVECTOR3& p0 = vert_[f.vpos[0]].p;
    const D3DXVECTOR3& p1 = vert_[f.vpos[1]].p;
    const D3DXVECTOR3& p2 = vert_[f.vpos[2]].p;
    const D3DXVECTOR3 v1 = p1 - p0;
    const D3DXVECTOR3 v2 = p2 - p1;
    D3DXVec3Cross(&norm, &v1, &v2);
    D3DXVec3Normalize(&f.norm, &norm);
    face_.push_back(f);
  }
}

void Mesh::make_adjinfo() {
  adjface_.resize(vert_.size());
  adjvert_elem_type defaultval;
  adjvert_.resize(vert_.size(), defaultval);
  for(int i = 0; i < face_.size(); ++i) {
    face_type& f = face_[i];
    for(int j = 0; j < 3; ++j) {
      // make adjacent info
      int vind = f.vpos[j];
      vert_type& v = vert_[vind];
      adjface_[vind].push_back(make_tuple(i, j));
      adjvert_[vind].push_back(f.vpos[(j+1)%3]);
      adjvert_[vind].push_back(f.vpos[(j+2)%3]);
    }
  }

  // uniquen adjverts
  foreach(adjvert_elem_type& adjs, adjvert_) {
    sort(adjs.begin(), adjs.end());
    adjs.erase(unique(adjs.begin(), adjs.end()), adjs.end());
  }

}

void Mesh::calc_normal() {
/*
 * Please refer the following paper for computation of normal and curvature.
 * [Meyer '02]
 * Discrete Differential-Geometry Operators for Triangulated 2-Manifolds
 * Mark Meyer, Mathieu Desbrun, Peter Schroder, and Alan H. Barr
 *
 * with, fix of following paper.
 * [Goldfeather '04]
 * A Novel Cubic-Order Algorithm for Approximating Principal Direction Vectors
 * Jack Goldfeather and Victoria Interrante
 */
  // initialize
  vector<D3DXVECTOR3> k(vert_.size(), D3DXVECTOR3(0.f, 0.f, 0.f));
  vector<float> total_angle(vert_.size(), 0.f);

  // for each faces
  foreach(face_type& f, face_) {
    vert_type& v0 = vert_[f.vpos[0]];
    vert_type& v1 = vert_[f.vpos[1]];
    vert_type& v2 = vert_[f.vpos[2]];
    const D3DXVECTOR3 e01 = v1.p - v0.p;
    const D3DXVECTOR3 e12 = v2.p - v1.p;
    const D3DXVECTOR3 e20 = v0.p - v2.p;
    const float e01_len = D3DXVec3Length(&e01);
    const float e12_len = D3DXVec3Length(&e12);
    const float e20_len = D3DXVec3Length(&e20);
    D3DXVECTOR3 cross0, cross1, cross2;
    D3DXVec3Cross(&cross0, &e01, &e20);
    D3DXVec3Cross(&cross1, &e12, &e01);
    D3DXVec3Cross(&cross2, &e20, &e12);
    const float cotan0 = -D3DXVec3Dot(&e01, &e20) / D3DXVec3Length(&cross0);
    const float cotan1 = -D3DXVec3Dot(&e12, &e01) / D3DXVec3Length(&cross1);
    const float cotan2 = -D3DXVec3Dot(&e20, &e12) / D3DXVec3Length(&cross2);
    const float theta0 = acos(-D3DXVec3Dot(&e01, &e20) / e01_len / e20_len);
    const float theta1 = acos(-D3DXVec3Dot(&e12, &e01) / e12_len / e01_len);
    const float theta2 = acos(-D3DXVec3Dot(&e20, &e12) / e20_len / e12_len);

    // calc normal, see ref.paper expression 8.
    k[f.vpos[0]] += e01 * cotan2 - e20 * cotan1;
    k[f.vpos[1]] += e12 * cotan0 - e01 * cotan2;
    k[f.vpos[2]] += e20 * cotan1 - e12 * cotan0;
    total_angle[f.vpos[0]] += theta0;
    total_angle[f.vpos[1]] += theta1;
    total_angle[f.vpos[2]] += theta2;
  }

  // for each vertices
  for(int i = 0; i < vert_.size(); ++i) {
    vert_type& v = vert_[i];

    if(fabs(total_angle[i] - 2*3.141592f) < 0.1f) {
      // in almost flat face, calclation of normal would be incorrect because of float error.
      // we use the average of adjacent faces' normals.
      v.n = D3DXVECTOR3(0.f, 0.f, 0.f);
      typedef tuple<int, int> adjtype;
      foreach(adjtype& adjf, adjface_[i]) {
        v.n += face_[adjf.get<0>()].norm;
      }
      D3DXVec3Normalize(&v.n, &v.n);
    }
    else {
      D3DXVec3Normalize(&v.n, &k[i]);
    }

    // fitting the direction of normal to the adj face's normal
    if(!adjface_[i].empty()) {
      float signtotal = 0;
      typedef tuple<int, int> adjf_type;
      foreach(adjf_type& t, adjface_[i]) {
        signtotal += sign(D3DXVec3Dot(&face_[t.get<0>()].norm, &v.n));
      }
      if(signtotal < 0) {
        v.n = -v.n;
      }
    }
  }
}

// normal, tangent, binormal
void Mesh::make_ntb() {
  foreach(vert_type& v, vert_) {
    // make tangent
    if(v.n.y == 0.f && v.n.z == 0.f) {
      v.t.y = (float)_copysign(1.0, v.n.x);
      v.t.x = v.t.z = 0.f;
    }
    else {
      v.t.x = 0.f;
      v.t.y =  v.n.z;
      v.t.z = -v.n.y;
      D3DXVec3Normalize(&v.t, &v.t);
    }

    // make binormal
    D3DXVec3Cross(&v.b, &v.n, &v.t);
  }
}

void Mesh::calc_curv_dir() {
  // implementing ref.paper 5.3, with fix of [Goldfeather '04].
  for(int i = 0; i < vert_.size(); ++i) {
    // using 2-adj-verts, for smooth result
    set<int> adjverts;
    foreach(int v0, adjvert_[i]) {
      adjverts.insert(v0);
      foreach(int v1, adjvert_[v0]) {
        adjverts.insert(v1);
      }
    }
    adjverts.erase(i);

    vert_type& vi = vert_[i];

    int adjnum = adjverts.size();
    ublas::matrix<float> D(adjnum, 3);
    ublas::matrix_column<ublas::matrix<float> > uu(D, 0), uv(D, 1), vv(D, 2);
    ublas::vector<float> kappa(adjnum);

    int j = 0;
    foreach(int vind, adjverts) {
      vert_type& vj = vert_[vind];
      D3DXVECTOR3 dij3 = vi.p - vj.p;
      D3DXVECTOR2 dij(D3DXVec3Dot(&dij3, &vi.t), D3DXVec3Dot(&dij3, &vi.b));
      D3DXVec2Normalize(&dij, &dij);
      uu(j) = dij.x * dij.x;
      uv(j) = dij.x * dij.y * 2;
      vv(j) = dij.y * dij.y;
      kappa(j) = 2 * D3DXVec3Dot(&dij3, &vi.n) / D3DXVec3LengthSq(&dij3);
      ++j;
    }

    ublas::matrix<float> tDD_inv(3, 3);
    safe_inverse(ublas::prod(ublas::trans(D), D), tDD_inv);
    ublas::vector<float> abc = ublas::prod(ublas::prod(tDD_inv, ublas::trans(D)), kappa);

    // now, we got B = ((a, b), (b, c)).
    float a = abc(0), b = abc(1), c = abc(2);
    if(!(a+c <= 0) && !(a+c >= 0)) {
      cerr << "failed to compute curvature" << endl;
      //DebugBreak();
      a = b = c = 0;
    }
    vi.curv_tensor[0] = a;
    vi.curv_tensor[1] = b;
    vi.curv_tensor[2] = c;

    // computing eigenvector/values of B.
    const float curv_mean = (a + c) / 2;
    const float curv_gauss = a*c - b*b;
    const float discriminant = sqrt(fabs(curv_mean * curv_mean - curv_gauss));
    const float curv_min = curv_mean - discriminant;
    const float curv_max = curv_mean + discriminant;

    if(fabs(b) < 0.00001f) {
      vi.curv_minmax_dir[0] = D3DXVECTOR2(0.f, 1.f);
      vi.curv_minmax_dir[1] = D3DXVECTOR2(1.f, 0.f);
      if(c > a) { swap(vi.curv_minmax_dir[0], vi.curv_minmax_dir[1]); }
    }
    else {
      const float d = sqrt((a-c)*(a-c) + 4*b*b);
      D3DXVECTOR2 mindir = D3DXVECTOR2((a - c + d) / 2 / b, 1.f);
      D3DXVECTOR2 maxdir = D3DXVECTOR2((a - c - d) / 2 / b, 1.f);
      D3DXVec2Normalize(&mindir, &mindir);
      D3DXVec2Normalize(&maxdir, &maxdir);
      const float curv0 = compute_dir_curv(vi.curv_tensor, mindir);
      const float curv1 = compute_dir_curv(vi.curv_tensor, maxdir);
      if(curv0 > curv1) {
        swap(mindir, maxdir);
      }

      vi.curv_minmax_dir[0] = mindir;
      vi.curv_minmax_dir[1] = maxdir;
    }
  }
}

void Mesh::smooth_curv() {
  vector<array<float, 3> > oldcurvs(vert_.size());
  for(int i = 0; i < vert_.size(); ++i) {
    copy(vert_[i].curv_tensor, vert_[i].curv_tensor + 3, oldcurvs[i].begin());
  }

  for(int i = 0; i < vert_.size(); ++i) {
    set<int> adjverts;
    foreach(int v0, adjvert_[i]) {
      adjverts.insert(v0);
    }
    adjverts.erase(i);

    vert_type& v = vert_[i];

    foreach(int j, adjverts) {
      v.curv_tensor[0] += oldcurvs[j][0];
      v.curv_tensor[1] += oldcurvs[j][1];
      v.curv_tensor[2] += oldcurvs[j][2];
    }
    int denom = 1 + adjverts.size();
    v.curv_tensor[0] /= denom;
    v.curv_tensor[1] /= denom;
    v.curv_tensor[2] /= denom;
  }
}


Mesh themesh;
struct octtree_node;

struct octtree {
  vector<int> vert;
  auto_ptr<octtree_node> root;

  template<typename Functor>
  void search(Mesh& m, Functor f, const D3DXVECTOR3& center, float distance); // f(vert_type& v, float distance)
};

struct octtree_node {
  vector<int>::iterator from, to;
  D3DXVECTOR3 v0, v1;
  D3DXVECTOR3 center;
  float diagonal;
  array<auto_ptr<octtree_node>, 8> children; // arary of auto_ptr. IS THIS OK?

  template<typename Functor>
  void search(Mesh& m, Functor f, const D3DXVECTOR3& center, float distance);
};

template<float _D3DVECTOR::* memptr>
struct veccoord_cmp : unary_function<bool, int> {
  veccoord_cmp(Mesh& mesh, float threshold) : m(mesh), t(threshold) {}
  Mesh& m;
  float t;
  bool operator()(int i) const {
    const D3DXVECTOR3& v = m.vert_[i].p;
    return v.*memptr < t;
  }
};

template<typename Iterator>
auto_ptr<octtree_node> create_octtree_node(
  Mesh& m, Iterator from, Iterator to,
  const D3DXVECTOR3& v0, const D3DXVECTOR3& v1, int depth)
{
  // typeof(*from) == typeof(int), index of Mesh::vert_.
  auto_ptr<octtree_node> n(new octtree_node);
  n->from = from;
  n->to   = to;
  n->v0   = v0;
  n->v1   = v1;
  n->center   = (v0 + v1) / 2;
  D3DXVECTOR3 diag(v1 - v0);
  n->diagonal = D3DXVec3Length(&diag);
  if(/*depth != 0 && */to - from >= 16) {
    float vz[3];
    Iterator zit[3];
    vz[0] = v0.z; vz[2] = v1.z;
    vz[1] = (v0.z + v1.z) / 2;
    zit[0] = from; zit[2] = to;
    zit[1] = partition(from, to, veccoord_cmp<&D3DXVECTOR3::z>(m, vz[1]));
    for(int z = 0; z < 2; ++z) {
      float vy[3];
      Iterator yit[3];
      vy[0] = v0.y; vy[2] = v1.y;
      vy[1] = (v0.y + v1.y) / 2;
      yit[0] = zit[z]; yit[2] = zit[z+1];
      yit[1] = partition(zit[z], zit[z+1], veccoord_cmp<&D3DXVECTOR3::y>(m, vy[1]));
      for(int y = 0; y < 2; ++y) {
        float vx[3];
        Iterator xit[3];
        vx[0] = v0.x; vx[2] = v1.x;
        vx[1] = (v0.x + v1.x) / 2;
        xit[0] = yit[y]; xit[2] = yit[y+1];
        xit[1] = partition(yit[y], yit[y+1], veccoord_cmp<&D3DXVECTOR3::x>(m, vx[1]));
        for(int x = 0; x < 2; ++x) {
          D3DXVECTOR3 newv0(vx[x  ], vy[y  ], vz[z  ]);
          D3DXVECTOR3 newv1(vx[x+1], vy[y+1], vz[z+1]);
          n->children[z*4+y*2+x] = create_octtree_node(m, xit[x], xit[x+1], newv0, newv1, depth-1);
        }
      }
    }
  }
  return n;
}

auto_ptr<octtree> create_octtree(Mesh& m, const D3DXVECTOR3& v0, const D3DXVECTOR3& v1) {
  auto_ptr<octtree> o(new octtree);
  o->vert.resize(m.vert_.size());
  for(int i = 0; i < o->vert.size(); ++i) {
    o->vert[i] = i;
  }
  int depth = (int)(log((double)m.vert_.size()) - log(8.0));
  depth += 0;
  o->root = create_octtree_node(m, o->vert.begin(), o->vert.end(), v0, v1, depth);
  return o;
}

template<typename Functor>
void octtree_node::search(Mesh& m, Functor f, const D3DXVECTOR3& center, float distance) {
  D3DXVECTOR3 d3(center - this->center);
  const float d = D3DXVec3Length(&d3);
  if(d + diagonal < distance) {
    // all included
    foreach(int i, make_iterator_range(from, to)) {
      Mesh::vert_type& v = m.vert_[i];
      D3DXVECTOR3 vtoc(center - v.p);
      f(v, D3DXVec3Length(&vtoc));
    }
  }
  else if(d - diagonal > distance) {
    // not included
  }
  else {
    // may included
    if(children[0].get()) {
      // if I have children, give over.
      for(int i = 0; i < 8; ++i) {
        children[i]->search(m, f, center, distance);
      }
    }
    else {
      // if not, calc on each verts.
      foreach(int i, make_iterator_range(from, to)) {
        Mesh::vert_type& v = m.vert_[i];
        D3DXVECTOR3 vtoc(center - v.p);
        float d = D3DXVec3Length(&vtoc);
        if(d < distance) { f(v, d); }
      }
    }
  }
}

template<typename Functor>
void octtree::search(Mesh& m, Functor f, const D3DXVECTOR3& center, float distance) {
  root->search(m, f, center, distance);
}

struct sal_func {
  sal_func(float epsilon, double* sigma, double* sigma_denom)
       : eps(epsilon), sig(sigma), sigd(sigma_denom)
    {
      for(int i = 0; i < 13; ++i) {
        exp_rcp2_s[i] = exp(1.f / (i*i));
      }
    }
  double *sig, *sigd;
  float eps;
  float exp_rcp2_s[13];
  void operator()(Mesh::vert_type& vj, float dist) const {
    if(dist == 0.f) { return; }
    float distance2 = dist*dist;
    int distancei = (int)(fabs(dist / eps));
    float e = exp(-distance2 / (eps*eps) / 2);

    for(int s = distancei; s < 13; ++s) {
      float e2 = e * exp_rcp2_s[s];
      float curv_mean = (vj.curv_tensor[0] + vj.curv_tensor[2]) / 2;
      sig[s]  += e2 * curv_mean;
      sigd[s] += e2;
    }
  }
};


void Mesh::calc_saliency() {
  // make bounding-box
  float epsilon;
  D3DXVECTOR3 bound_min( FLT_MAX,  FLT_MAX,  FLT_MAX);
  D3DXVECTOR3 bound_max(-FLT_MAX, -FLT_MAX, -FLT_MAX);
  {
    for(int i = 0; i < vert_.size(); ++i) {
      vert_type& v = vert_[i];
      D3DXVec3Minimize(&bound_min, &bound_min, &v.p);
      D3DXVec3Maximize(&bound_max, &bound_max, &v.p);
    }
    D3DXVECTOR3 diagonal(bound_max-bound_min);
    epsilon = D3DXVec3Length(&diagonal) * 0.003;
  }
  auto_ptr<octtree> ot(create_octtree(*this, bound_min, bound_max));

  array<float, 5> zeros = { };
  vector<array<float, 5> > saliencies(vert_.size(), zeros);

  // compute saliencies on each scale level.
//#pragma omp parallel for
  for(int i = 0; i < vert_.size(); ++i) {
    // print progress
    if(i % (vert_.size() / 10) == 0) {
      cout << "saliency " << (i / (vert_.size() / 10)) << "0%" << endl;
    }

    vert_type& vi = vert_[i];
    array<double, 13> sigma = {};
    array<double, 13> sigma_denom = {};
    sal_func f(epsilon, sigma.c_array(), sigma_denom.c_array());
    ot->search(*this, f, vi.p, epsilon * 12);

    for(int s = 0; s < 13; ++s) {
      if(sigma_denom[s] == 0.0) { sigma_denom[s] += 1.0; }
      sigma[s] /= sigma_denom[s];
    }

    for(int e = 2; e <= 6; ++e) {
      saliencies[i][e-2] = fabs(sigma[e] - sigma[e*2]);
    }
  }

  // normalize saliencies
  vector<float> max_saliencies(5, 0.0);
  for(int i = 0; i < vert_.size(); ++i) {
    for(int e = 0; e < 5; ++e) {
      max_saliencies[e] = max(max_saliencies[e], saliencies[i][e]);
    }
  }
  for(int i = 0; i < vert_.size(); ++i) {
    for(int e = 0; e < 5; ++e) {
      saliencies[i][e] /= max_saliencies[e];
    }
  }

  // find local maximams
  array<float, 5> local_max = {};
  array<float, 5> local_max_denom = {};
  for(int i = 0; i < vert_.size(); ++i) {
    float maxs[] = { -FLT_MAX, -FLT_MAX, -FLT_MAX, -FLT_MAX, -FLT_MAX };
    for(int j = 0; j < adjvert_[i].size(); ++j) {
      for(int e = 0; e < 5; ++e) {
        maxs[e] = max(maxs[e], saliencies[j][e]);
      }
    }

    for(int e = 0; e < 5; ++e) {
      float sal = saliencies[i][e];
      if(maxs[e] == sal && sal != 1.0) {
        local_max[e] += sal;
        local_max_denom[e] += 1.f;
      }
    }
  }
  for(int e = 0; e < 5; ++e) {
    if(local_max_denom[e] == 0.0) { local_max_denom[e] = 1.0; }
    local_max[e] /= local_max_denom[e];
  }

  // compute final saliency
  float scale_coeff[5];
  for(int e = 0; e < 5; ++e) {
    scale_coeff[e] = pow(1.0 - local_max[e], 2);
  }
  vector<float> total_saliencies(vert_.size(), 0.f);
  for(int i = 0; i < vert_.size(); ++i) {
    for(int e = 0; e < 5; ++e) {
      total_saliencies[i] += saliencies[i][e] * scale_coeff[e];
    }
  }

  foreach(float& sal, total_saliencies) {
    if(!(sal <= 0.f || 0.f <= sal)) {
      // v.saliency is NaN
      sal = 1.f;
    }
  }

  // smooth
  for(int i = 0; i < vert_.size(); ++i) {
    vert_type& v = vert_[i];
    v.saliency = total_saliencies[i];
    for(int j = 0; j < adjvert_[i].size(); ++j) {
      v.saliency += total_saliencies[adjvert_[i][j]];
    }
    v.saliency /= adjvert_[i].size() + 1;
    //printf("%f ",v.saliency);
  }
}

