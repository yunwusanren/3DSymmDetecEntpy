#pragma once

#include "vertex.hpp"
#include <d3dx10.h>
#include <vector>
#include <istream>
#include <boost/tuple/tuple.hpp>

void loadMesh(ID3D10Device* pdevice, ID3D10Buffer** ppvbuf, ID3D10Buffer** ppibuf,
              int* numvertex, int* numprimitive, const std::string& objfile);

class Mesh {
public:
  Mesh();

  void load_obj(std::istream& ifs);
  void load_ifs(std::istream& ifs);
  void make_adjinfo();
  void calc_normal();
  void make_ntb();
  void calc_curv_dir();
  void smooth_curv();
  void calc_saliency();

  HRESULT make_D3Dvert(ID3D10Device* pdevice,
                       ID3D10Buffer** ppverts,
                       ID3D10Buffer** ppindices,
                       int* pnumvertex,
                       int* pnumprimitive);

  typedef d3dvertex vert_type;
  std::vector<vert_type> vert_;

  typedef std::vector<int> adjvert_elem_type;
  typedef std::vector<adjvert_elem_type> adjvert_type;
  adjvert_type adjvert_;
  std::vector<std::vector<boost::tuple<int, int> > > adjface_;

  struct face_type {
    int vpos[3];
    bool hasadj;
    D3DXVECTOR3 norm;
  };
  std::vector<face_type> face_;
};
