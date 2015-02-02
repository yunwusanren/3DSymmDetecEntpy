
#pragma once

#include <windows.h>
#include <d3d10.h>
#include <d3dx10.h>

struct d3dvertex {
  D3DXVECTOR3 p;

  D3DXVECTOR3 n, t, b; // smooth surface,s normal, and tangent plane.

  D3DXVECTOR2 curv_minmax_dir[2];

  /*
   * |c[0] c[1]|
   * |c[1] c[2]|, this matrix defines the curvature of the point.
   */
  float curv_tensor[3];

  float saliency;
};
