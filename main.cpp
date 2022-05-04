#include "XNumeric.hpp"
#include "XRandom.hpp"
#include "XDiffEq.hpp"
#include <vector>
#include <cmath>
#include <ctime>
#include <fstream>

int N;
int P;
int d = 2;
XNum m = 1;
XNum kB = 1;
XNum T;
XNum beta;
XNum hBar = 1.0;
XNum omgP;
XNum L = 30.0;
XNum omg0 = 1.0 / hBar;
//XNum lamda = 0.0;
XNum s = 3.0;
XNum g = 0.0;
int cor_num = 60;
XNum incre = 5.0 / cor_num;
XNum vi = 0.0;

XNum minimum_image(XNum d) {
  if (std::abs(d) > L / 2)
    return L - std::abs(d);
  return d;
}

XNum minimum_image2(XNum d) {
  if (std::abs(d) > L / 2) {
    if (d < 0)
      return L - std::abs(d);
    else
      return -(L - std::abs(d));
  }
  return d;
}

XNum distance(int i, int j, const std::vector<XNum> &Rs) {
  XNum a = minimum_image(Rs[d*i]-Rs[d*j]);
  XNum b = minimum_image(Rs[d*i+1]-Rs[d*j+1]);
  return a*a+b*b;
}

std::vector<XNum> Aj(const std::vector<XNum> &Rs) {
  int j, k, l;
  std::vector<XNum> res;
  for (j = 0; j < N; ++j)
    for (l = 0; l < P; ++l) {
      XNum sum = 0, sum2 = 0;
      /*for (k = 0; k < N; ++k) {
        if (k == j) continue;
        int index = j*P+l;
        int index2 = k*P+l;
        XNum dis = distance(index, index2, Rs);
        sum += vi*minimum_image2(Rs[d*index+1]-Rs[d*index2+1])/dis;
        sum2 += -vi*minimum_image2(Rs[d*index]-Rs[d*index2])/dis;
      }*/
      int in = j*P+l;
      XNum a = Rs[d*in]-L/2;
      XNum b = Rs[d*in+1]-L/2;
      XNum dis = minimum_image(a)*minimum_image(a)+minimum_image(b)*minimum_image(b);
      sum += vi*minimum_image2(b)/dis;
      sum2 += -vi*minimum_image2(a)/dis;
      /*sum += -vi*minimum_image2(b);
      sum2 += vi*minimum_image2(a);*/
      res.push_back(sum);
      res.push_back(sum2);
    }
  return res;
}

/*XNum ANk(int N2, int k, const std::vector<XNum> &Rs, const std::vector<XNum> &a) {
  XNum res = 0;
  int l, j;
  for (l = N2 - k + 1; l <= N2; ++l)
    for (j = 1; j <= P; ++j) {
      int index = (l - 1) * P + j - 1;
      int index2 = (l - 1) * P + j;
      if (j == P) {
        if (l == N2)
          index2 = (N2 - k) * P;
        else
          index2 = l * P;
      }
      XNum a1 = minimum_image2(Rs[d*index2]-L/2);
      XNum a2 = minimum_image2(Rs[d*index2+1]-L/2);
      XNum b1 = minimum_image2(Rs[d*index]-L/2);
      XNum b2 = minimum_image2(Rs[d*index+1]-L/2);
      XNum angle = atan2(b1*a2-b2*a1, a1*b1+a2*b2);
      res += angle;
      //res += minimum_image2(Rs[d*index2]-Rs[d*index])*(a[2*index]+a[2*index2])/2;
      //res += minimum_image2(Rs[d*index2+1]-Rs[d*index+1])*(a[2*index+1]+a[2*index2+1])/2;
    }
  for (j = 0; j < P; ++j) {
    int index = j, index2 = (j+1)%P;
    XNum a1 = minimum_image2(Rs[d*index2]-L/2);
    XNum a2 = minimum_image2(Rs[d*index2+1]-L/2);
    XNum b1 = minimum_image2(Rs[d*index]-L/2);
    XNum b2 = minimum_image2(Rs[d*index+1]-L/2);
    XNum angle = atan2(b1*a2-b2*a1, a1*b1+a2*b2);
    res += angle;
  }
  res /= 2 * M_PI;
  return res;
}*/

/*XNum ANk(int N2, int k, const std::vector<XNum> &Rs) {
  XNum res = 0;
  int l, j, l2;
  for (l = N2 - k + 1; l <= N2; ++l)
    for (j = 1; j <= P; ++j) {
      for (l2 = N2 - k + 1; l2 <= N2; ++l2) {
        if (l == l2) continue;
        int index = (l - 1) * P + j - 1;
        int index2 = (l - 1) * P + j;
        if (j == P) {
          if (l == N2)
            index2 = (N2 - k) * P;
          else
            index2 = l * P;
          //index2 = (l - 1) * P;
        }
        XNum r1x, r1y, r2x, r2y;
        r1x = Rs[d*index];
        r1y = Rs[d*index+1];
        r2x = Rs[d*index2];
        r2y = Rs[d*index2+1];
        int index3 = (l2 - 1) * P + j - 1;
        int index4 = (l2 - 1) * P + j;
        if (j == P) {
          if (l2 == N2)
            index4 = (N2 - k) * P;
          else
            index4 = l2 * P;
          index4 = (l2 - 1) * P;
        }
        XNum r3x, r3y, r4x, r4y;
        r3x = Rs[d*index3];
        r3y = Rs[d*index3+1];
        r4x = Rs[d*index4];
        r4y = Rs[d*index4+1];
        XNum a1 = minimum_image2(r2x-r4x);
        XNum a2 = minimum_image2(r2y-r4y);
        XNum b1 = minimum_image2(r1x-r3x);
        XNum b2 = minimum_image2(r1y-r3y);
        XNum angle = atan2(b1*a2-b2*a1, a1*b1+a2*b2);
        res += angle;
      }
    }
  res /= 2 * M_PI;
  return res;
}*/

XNum ANk(int N2, int k, const std::vector<XNum> &Rs) {
  XNum res = 0;
  int l, j, l2;
  for (l = N2 - k + 1; l <= N2; ++l)
    for (j = 1; j <= P; ++j) {
      int l3 = l;
      /*if (j == P) {
        if (l == N2)
          l3 = N2 - k + 1;
        else
          l3 = l + 1;
      }*/
      for (l2 = N2 - k + 1; l2 <= N2; ++l2) {
        if (l2 == l3) continue;
        int index = (l - 1) * P + j - 1;
        int index2 = (l - 1) * P + j;
        if (j == P) {
          if (l == N2)
            index2 = (N2 - k) * P;
          else
            index2 = l * P;
          index2 = (l - 1) * P;
        }
        XNum r1x, r1y, r2x, r2y;
        r1x = Rs[d*index];
        r1y = Rs[d*index+1];
        r2x = Rs[d*index2];
        r2y = Rs[d*index2+1];
        int index3 = (l2 - 1) * P + j - 1;
        int index4 = (l2 - 1) * P + j;
        if (j == P) {
          if (l2 == N2)
            index4 = (N2 - k) * P;
          else
            index4 = l2 * P;
          index4 = (l2 - 1) * P;
        }
        XNum r3x, r3y, r4x, r4y;
        r3x = Rs[d*index3];
        r3y = Rs[d*index3+1];
        r4x = Rs[d*index4];
        r4y = Rs[d*index4+1];
        XNum a1 = minimum_image2(r2x-r4x);
        XNum a2 = minimum_image2(r2y-r4y);
        XNum b1 = minimum_image2(r1x-r3x);
        XNum b2 = minimum_image2(r1y-r3y);
        XNum angle = atan2(b1*a2-b2*a1, a1*b1+a2*b2);
        res += angle;
      }
    }
  res /= 2 * M_PI;
  return res/2;
}

std::pair<XNum, XNum> APhase(const std::vector<XNum> &Rs) {
  std::vector<XNum> RP, IP;
  //std::vector<XNum> a = Aj(Rs);
  /*RP.push_back(1);
  IP.push_back(0);
  int N2, k;
  for (N2 = 1; N2 <= N; ++N2) {
    XNum tmpr = 0, tmpi = 0;
    for (k = 1; k <= N2; ++k) {
      XNum tmp = ANk(N2, k, Rs);
      //tmpr += RP[N2-k]*cos(tmp)-IP[N2-k]*sin(tmp);
      //tmpi += RP[N2-k]*sin(tmp)+IP[N2-k]*cos(tmp);
      tmpr += RP[N2-k]*cos(vi*M_PI*tmp)-IP[N2-k]*sin(vi*M_PI*tmp);
      tmpi += RP[N2-k]*sin(vi*M_PI*tmp)+IP[N2-k]*cos(vi*M_PI*tmp);
    }
    RP.push_back(tmpr/N2);
    IP.push_back(tmpi/N2);
  }
  std::pair<XNum, XNum> res(RP[N], IP[N]);*/
  std::pair<XNum, XNum> res(cos(vi*M_PI*ANk(N,N,Rs)), sin(vi*M_PI*ANk(N,N,Rs)));
  return res;
}

std::vector<XNum> AjCache;

XNum ENk(int N2, int k, const std::vector<XNum> &Rs) {
  XNum res = 0;
  int l, j;
  for (l = N2 - k + 1; l <= N2; ++l)
    for (j = 1; j <= P; ++j) {
      int index = (l - 1) * P + j - 1;
      int index2 = (l - 1) * P + j;
      if (j == P) {
        if (l == N2)
          index2 = (N2 - k) * P;
        else
          index2 = l * P;
      }
      res += distance(index2, index, Rs);
    }
  return 0.5*m*omgP*omgP*res;
  /*XNum res2 = 0;
  for (l = N2 - k + 1; l <= N2; ++l)
    for (j = 1; j <= P; ++j) {
      int index = (l - 1) * P + j - 1;
      int index2 = (l - 1) * P + j;
      if (j == P) {
        if (l == N2)
          index2 = (N2 - k) * P;
        else
          index2 = l * P;
      }
      XNum a = AjCache[2*index2]-AjCache[2*index];
      XNum b = AjCache[2*index2+1]-AjCache[2*index+1];
      res2 += a*a+b*b;
    }
  return 0.5*m*omgP*omgP*res+res2/(8*m*P);*/
}

XNum ENk2(int N2, int k, const std::vector<XNum> &Rs) {
  XNum res = 0;
  int l, j;
  for (l = N2 - k + 1; l <= N2; ++l)
    for (j = 1; j <= P; ++j) {
      int index = (l - 1) * P + j - 1;
      int index2 = (l - 1) * P + j;
      if (j == P) {
        if (l == N2)
          index2 = (N2 - k) * P;
        else
          index2 = l * P;
      }
      res += distance(index2, index, Rs);
    }
  return 0.5*m*omgP*omgP*res;
}

std::vector<XVecD> ENkCache;
std::vector<XVecD> ENkCache2;
std::vector<XNum> VBCache;

std::vector<XNum> expVB(int N, const std::vector<XNum> &Rs) {
  std::vector<XNum> res;
  res.push_back(0);
  int N2, k;
  for (N2 = 1; N2 <= N; ++N2) {
    XNum sum = 0;
    XNum tmp = 0.5*(ENkCache[N2-1][N2-1]+res[N2-1]);
    for (k = 1; k <= N2; ++k)
      sum += std::exp(-beta*(ENkCache[N2-1][k-1]+res[N2-k]-tmp));
    res.push_back(tmp-(1/beta)*std::log(sum / N2));
  }
  return res;
}

std::vector<XNum> gradient(int N2, int k, int l, int j, const std::vector<XNum> &Rs) {
  std::vector<XNum> res;
  res.push_back(0);
  res.push_back(0);
  if (l >= N2 - k + 1 && l <= N2) {
    int index = (l - 1) * P + j - 1;
    int index2 = (l - 1) * P + j;
    if (j == P) {
      if (l == N2)
        index2 = (N2 - k) * P;
      else
        index2 = l * P;
    }
    int index3 = (l - 1) * P + j - 2;
    if (j == 1) {
      if (l == N2 - k + 1)
        index3 = (N2 - 1) * P + P - 1;
      else
        index3 = (l - 2) * P + P - 1;
    }
    XNum a = Rs[d*index] - Rs[d*index2];
    XNum b = Rs[d*index] - Rs[d*index3];
    res[0] = m*omgP*omgP*(minimum_image2(a)+minimum_image2(b));
    a = Rs[d*index+1] - Rs[d*index2+1];
    b = Rs[d*index+1] - Rs[d*index3+1];
    res[1] = m*omgP*omgP*(minimum_image2(a)+minimum_image2(b));
    /*a = minimum_image2(Rs[d*index]-L/2);
    b = minimum_image2(Rs[d*index+1]-L/2);
    XNum daxx = -2*vi*a*b/std::pow(a*a+b*b,2);
    XNum dayx = 2*vi*a*a/std::pow(a*a+b*b,2)-vi/(a*a+b*b);
    //XNum daxx = 0;
    //XNum dayx = vi;
    res[0] += (AjCache[2*index]-AjCache[2*index2])*daxx/(4*m*P);
    res[0] += (AjCache[2*index+1]-AjCache[2*index2+1])*dayx/(4*m*P);
    res[0] += (AjCache[2*index]-AjCache[2*index3])*daxx/(4*m*P);
    res[0] += (AjCache[2*index+1]-AjCache[2*index3+1])*dayx/(4*m*P);
    XNum daxy = -2*vi*b*b/std::pow(a*a+b*b,2)+vi/(a*a+b*b);
    XNum dayy = 2*vi*a*b/std::pow(a*a+b*b,2);
    //XNum daxy = -vi;
    //XNum dayy = 0;
    res[1] += (AjCache[2*index]-AjCache[2*index2])*daxy/(4*m*P);
    res[1] += (AjCache[2*index+1]-AjCache[2*index2+1])*dayy/(4*m*P);
    res[1] += (AjCache[2*index]-AjCache[2*index3])*daxy/(4*m*P);
    res[1] += (AjCache[2*index+1]-AjCache[2*index3+1])*dayy/(4*m*P);*/
  }
  return res;
}

std::vector<XNum> XForceVBCache;

void forceVB(int N, const std::vector<XNum> &Rs) {
  int N2, k, l, j;
  XForceVBCache.clear();
  std::vector<XNum> tmp;
  for (l = 1; l <= N; ++l)
    for (j = 1; j <= P; ++j) {
      XVecD res;
      res.push_back(0);
      res.push_back(0);
      for (N2 = 1; N2 <= N; ++N2) {
        int i;
        XNum sum2 = 0;
        XNum tmp2 = 0.5*(ENkCache[N2-1][N2-1]+VBCache[N2-1]);
        for (k = 1; k <= N2; ++k)
          sum2 += exp(-beta*(ENkCache[N2-1][k-1]+VBCache[N2-k]-tmp2));
        for (i = 0; i < d; ++i) {
          XNum sum = 0;
          for (k = 1; k <= N2; ++k) {
            tmp = gradient(N2, k, l, j, Rs);
            sum += (tmp[i] + res[d*(N2-k)+i])*exp(-beta*(ENkCache[N2-1][k-1]+VBCache[N2-k]-tmp2));
          }
          res.push_back(sum / sum2);
        }
      }
      XForceVBCache.push_back(res[d*N]);
      XForceVBCache.push_back(res[d*N+1]);
    }
}

XVecD force(XNum t, const XVecD &Rs) {
  XVecD res;
  std::vector<XNum> tmp;
  int l, j;
  AjCache = Aj(Rs);
  if (ENkCache.empty()) {
    for (l = 1; l <= N; ++l) {
      XVecD tmp;
      for (j = 1; j <= l; ++j)
        tmp.push_back(0);
      ENkCache.push_back(tmp);
    }
  }
  for (l = 1; l <= N; ++l)
    for (j = 1; j <= l; ++j)
      ENkCache[l-1][j-1] = ENk(l, j, Rs);
  if (ENkCache2.empty()) {
    for (l = 1; l <= N; ++l) {
      XVecD tmp;
      for (j = 1; j <= l; ++j)
        tmp.push_back(0);
      ENkCache2.push_back(tmp);
    }
  }
  for (l = 1; l <= N; ++l)
    for (j = 1; j <= l; ++j)
      ENkCache2[l-1][j-1] = ENk2(l, j, Rs);
  VBCache = expVB(N, Rs);
  forceVB(N, Rs);
  for (l = 1; l <= N; ++l)
    for (j = 1; j <= P; ++j) {
      int index = (l - 1) * P + j - 1;
      XNum a = Rs[d*index] - L/2;
      res.push_back(-XForceVBCache[d*index]/m-minimum_image2(a)*omg0*omg0/P);
      a = Rs[d*index+1] - L/2;
      res.push_back(-XForceVBCache[d*index+1]/m-minimum_image2(a)*omg0*omg0/P);
    }
  int k;
  for (j = 1; j <= P; ++j)
    for (l = 1; l <= N; ++l) {
      int index = (l - 1) * P + j - 1;
      for (k = 1; k <= N; ++k) {
        if (l == k) continue;
        int index2 = (k - 1) * P + j - 1;
        /*XNum inter = lamda / std::pow(distance(index, index2, Rs), 1.5);
        XNum a = Rs[d*index]-Rs[d*index2];
        res[d*index] += minimum_image2(a)*inter/P;
        a = Rs[d*index+1]-Rs[d*index2+1];
        res[d*index+1] += minimum_image2(a)*inter/P;*/
        XNum inter = (g/(M_PI*s*s))*std::exp(-distance(index, index2, Rs)/(s*s));
        XNum a = Rs[d*index]-Rs[d*index2];
        res[d*index] += (2*minimum_image2(a)/(s*s))*inter/P;
        a = Rs[d*index+1]-Rs[d*index2+1];
        res[d*index+1] += (2*minimum_image2(a)/(s*s))*inter/P;
      }
      /*XNum a = minimum_image2(Rs[d*index]-L/2);
      XNum b = minimum_image2(Rs[d*index+1]-L/2);
      res[d*index] += 2*vi*vi*a*b*b/std::pow(a*a+b*b,3)/P;
      res[d*index+1] += 2*vi*vi*b*a*a/std::pow(a*a+b*b,3)/P;*/
      XNum a = minimum_image2(Rs[d*index]-L/2);
      XNum b = minimum_image2(Rs[d*index+1]-L/2);
      //res[d*index] += 2*0.05*a/std::pow(a*a+b*b,2)/P;
      //res[d*index+1] += 2*0.05*b/std::pow(a*a+b*b,2)/P;
      res[d*index] += 0.05*a/std::pow(a*a+b*b,1.5)/P;
      res[d*index+1] += 0.05*b/std::pow(a*a+b*b,1.5)/P;
      /*XNum a = minimum_image2(Rs[d*index]-L/2);
      XNum b = minimum_image2(Rs[d*index+1]-L/2);
      XNum g = ::g / 5, s = ::s / 5;
      XNum inter = (g/(M_PI*s*s))*std::exp(-(a*a+b*b)/(s*s));
      res[d*index] += (2*a/(s*s))*inter/P;
      res[d*index+1] += (2*b/(s*s))*inter/P;*/
    }
  return res;
}

XNum GauEnergy(const std::vector<XNum> &Rs) {
  XNum res = 0;
  int j, k, l;
  for (j = 1; j <= P; ++j)
    for (l = 1; l <= N; ++l) {
      int index = (l - 1) * P + j - 1;
      for (k = 1; k <= N; ++k) {
        if (l == k) continue;
        int index2 = (k - 1) * P + j - 1;
        //XNum inter = lamda / std::pow(distance(index, index2, Rs), 0.5);
        XNum inter = (g/(M_PI*s*s))*std::exp(-distance(index, index2, Rs)/(s*s));
        res += 0.5*inter;
      }
      /*XNum a = minimum_image2(Rs[d*index]-L/2);
      XNum b = minimum_image2(Rs[d*index+1]-L/2);
      res += vi*vi*b*b/(2*std::pow(a*a+b*b,2));
      res += vi*vi*a*a/(2*std::pow(a*a+b*b,2));*/
      XNum a = minimum_image2(Rs[d*index]-L/2);
      XNum b = minimum_image2(Rs[d*index+1]-L/2);
      //res += 0.05/(a*a+b*b);
      res += 0.05/std::sqrt(a*a+b*b);
      //XNum g = ::g / 5, s = ::s / 5;
      //res += (g/(M_PI*s*s))*std::exp(-(a*a+b*b)/(s*s));
    }
  return res/P;
}

XNum energy(const std::vector<XNum> &Rs) {
  std::vector<XNum> res;
  res.push_back(0);
  int N2, k;
  for (N2 = 1; N2 <= N; ++N2) {
    XNum tmp2 = 0.5*(ENkCache[N2-1][N2-1]+VBCache[N2-1]);
    XNum sum2 = 0;
    for (k = 1; k <= N2; ++k)
      sum2 += exp(-beta*(ENkCache[N2-1][k-1]+VBCache[N2-k]-tmp2));
    XNum sum = 0;
    for (k = 1; k <= N2; ++k)
      sum += (res[N2-k]+ENkCache[N2-1][k-1]-2*ENkCache2[N2-1][k-1])*exp(-beta*(ENkCache[N2-1][k-1]+VBCache[N2-k]-tmp2));
    res.push_back(sum / sum2);
  }
  return res[N];
}

void period_boundary(std::vector<XNum> &Rs) {
  int i, j, k;
  for (i = 0; i < N; ++i)
    for (j = 0; j < P; ++j) {
      int index = i * P + j;
      for (k = 0; k < d; ++k) {
        if (Rs[d*index+k] < 0) {
          int n = (int)(std::abs(Rs[d*index+k]) / L);
          Rs[d*index+k] += (n + 1) * L;
        }
        else {
          int n = (int)(std::abs(Rs[d*index+k]) / L);
          Rs[d*index+k] -= n * L;
        }
      }
    }
}

void initial(std::vector<XNum> &Rs, std::vector<XNum> &Vs) {
  int i, j;
  XNum v1 = 0, v2 = 0;
  for (i = 0; i < N; ++i)
    for (j = 0; j < P; ++j) {
      XNum tmp = XRandGauss() * std::sqrt(1 / (m * beta));
      v1 += tmp;
      Vs.push_back(tmp);
      tmp = XRandGauss() * std::sqrt(1 / (m * beta));
      v2 += tmp;
      Vs.push_back(tmp);
      Rs.push_back(L / 2 + 1.0 * (XRandFloat()-0.5));
      Rs.push_back(L / 2 + 1.0 * (XRandFloat()-0.5));
    }
  v1 /= N * P;
  v2 /= N * P;
  for (i = 0; i < N; ++i)
    for (j = 0; j < P; ++j) {
      int index = i * P + j;
      Vs[d*index] -= v1;
      Vs[d*index+1] -= v2;
    }
}

void test() {
  std::vector<XNum> Rs, Vs;
  int i, j;
  initial(Rs, Vs);
  for (i = 0; i < N; ++i)
    for (j = 0; j < P; ++j) {
      int index = i * P + j;
      std::cout << Rs[index * d] << " " << Rs[index * d + 1] << std::endl;
    }
  XVecD f = force(0, Rs);
  for (i = 0; i < f.size(); ++i)
    std::cout << f[i] << " ";
  std::cout << std::endl;
}

void logR(const std::vector<XNum> &Rs, std::ostream &out) {
  out << "{";
  for (int i = 0; i < N * P; ++i) {
    out << "{" << Rs[d*i] << "," << Rs[d*i+1] << "}";
    if (i != N * P - 1)
      out << ",";
  }
  out << "}";
}

void velocity_rescale(std::vector<XNum> &Vs) {
  XNum sum = 0;
  int i;
  for (i = 0; i < N * P; ++i)
    sum += Vs[d*i]*Vs[d*i]+Vs[d*i+1]*Vs[d*i+1];
  XNum lam = std::sqrt((N*P)*d*(1/beta)/(m*sum));
  for (i = 0; i < N * P; ++i) {
    Vs[d*i] *= lam;
    Vs[d*i+1] *= lam;
  }
}

XNum temperature(const std::vector<XNum> &Vs) {
  XNum sum = 0;
  int i;
  for (i = 0; i < N * P; ++i)
    sum += Vs[d*i]*Vs[d*i]+Vs[d*i+1]*Vs[d*i+1];
  return m*sum/(d*(N*P));
}

XNum trapEnergy(const std::vector<XNum> &Rs) {
  XNum res = 0;
  int i;
  for (i = 0; i < N * P; ++i) {
    XNum a = Rs[d*i] - L/2;
    res += 0.5 * m * omg0 * omg0 * minimum_image(a) * minimum_image(a);
    a = Rs[d*i+1] - L/2;
    res += 0.5 * m * omg0 * omg0 * minimum_image(a) * minimum_image(a);
  }
  return res/P;
}

void XMVerletNHChain(XNum t, XNum h, XVecD &x, XVecD &v, std::vector<XVecD> &theta, std::vector<XVecD> &vtheta, XNum beta, XNum m, std::vector<XVecD> &Q, XForceFunc *func) {
  int i, j;
  int N2 = N;
  int N = d;
  int M = vtheta[0].size();
  if (XForceCache.empty())
    XForceCache = func(t, x);
  for (j = 0; j < N2*P; ++j) {
    XVecD v2;
    for (i = 0; i < d; ++i)
      v2.push_back(v[d*j+i]);
    XVecD f2 = NHForce(beta, m, d, Q[j], v2, vtheta[j]);
    for (i = 0; i < N; ++i)
      v[d*j+i] = v[d*j+i]*std::exp(-0.5*h*vtheta[j][0])+0.5*h*XForceCache[d*j+i]*std::exp(-0.25*h*vtheta[j][0]);
    int M2 = M / 2;
    for (i = 1; i <= M2; ++i)
      theta[j][2*i-2] = theta[j][2*i-2]+h*vtheta[j][2*i-2]/2;
    for (i = 1; i <= M2; ++i)
      vtheta[j][2*i-1] = vtheta[j][2*i-1]*std::exp(-0.5*h*((i==M2)?0:vtheta[j][2*i]))+0.5*h*f2[2*i-1]*std::exp(-0.25*h*((i==M2)?0:vtheta[j][2*i]));
    for (i = 0; i < N; ++i)
      x[d*j+i] = x[d*j+i] + h*v[d*j+i];
    for (i = 1; i <= M2; ++i)
      theta[j][2*i-1] = theta[j][2*i-1]+h*vtheta[j][2*i-1];
    for (i = 0; i < d; ++i)
      v2[i] = v[d*j+i];
    XVecD f3 = NHForce(beta, m, d, Q[j], v2, vtheta[j]);
    for (i = 1; i <= M2; ++i)
      vtheta[j][2*i-2] = vtheta[j][2*i-2]*std::exp(-h*vtheta[j][2*i-1])+h*f3[2*i-2]*std::exp(-0.5*h*vtheta[j][2*i-1]);
  }
  XForceCache = func(t+h, x);
  for (j = 0; j < N2*P; ++j) {
    int M2 = M / 2;
    for (i = 0; i < N; ++i)
      v[d*j+i] = v[d*j+i]*std::exp(-0.5*h*vtheta[j][0])+0.5*h*XForceCache[d*j+i]*std::exp(-0.25*h*vtheta[j][0]);
    for (i = 1; i <= M2; ++i)
      theta[j][2*i-2] = theta[j][2*i-2]+h*vtheta[j][2*i-2]/2;
    XVecD v2;
    for (i = 0; i < d; ++i)
      v2.push_back(v[d*j+i]);
    XVecD f3 = NHForce(beta, m, d, Q[j], v2, vtheta[j]);
    for (i = 1; i <= M2; ++i)
      vtheta[j][2*i-1] = vtheta[j][2*i-1]*std::exp(-0.5*h*((i==M2)?0:vtheta[j][2*i]))+0.5*h*f3[2*i-1]*std::exp(-0.25*h*((i==M2)?0:vtheta[j][2*i]));
  }
}

void XMMVerletNHChain(XNum t, XNum h, XVecD &x, XVecD &v, std::vector<XVecD> &theta, std::vector<XVecD> &vtheta, XNum beta, XNum m, std::vector<XVecD> &Q, XForceFunc *func) {
  int i, j;
  int N2 = N;
  int d2 = d;
  int d = 1;
  int N = d;
  int M = vtheta[0].size();
  if (XForceCache.empty())
    XForceCache = func(t, x);
  for (j = 0; j < N2*P*d2; ++j) {
    XVecD v2;
    for (i = 0; i < d; ++i)
      v2.push_back(v[d*j+i]);
    XVecD f2 = NHForce(beta, m, d, Q[j], v2, vtheta[j]);
    for (i = 0; i < N; ++i)
      v[d*j+i] = v[d*j+i]*std::exp(-0.5*h*vtheta[j][0])+0.5*h*XForceCache[d*j+i]*std::exp(-0.25*h*vtheta[j][0]);
    int M2 = M / 2;
    for (i = 1; i <= M2; ++i)
      theta[j][2*i-2] = theta[j][2*i-2]+h*vtheta[j][2*i-2]/2;
    for (i = 1; i <= M2; ++i)
      vtheta[j][2*i-1] = vtheta[j][2*i-1]*std::exp(-0.5*h*((i==M2)?0:vtheta[j][2*i]))+0.5*h*f2[2*i-1]*std::exp(-0.25*h*((i==M2)?0:vtheta[j][2*i]));
    for (i = 0; i < N; ++i)
      x[d*j+i] = x[d*j+i] + h*v[d*j+i];
    for (i = 1; i <= M2; ++i)
      theta[j][2*i-1] = theta[j][2*i-1]+h*vtheta[j][2*i-1];
    for (i = 0; i < d; ++i)
      v2[i] = v[d*j+i];
    XVecD f3 = NHForce(beta, m, d, Q[j], v2, vtheta[j]);
    for (i = 1; i <= M2; ++i)
      vtheta[j][2*i-2] = vtheta[j][2*i-2]*std::exp(-h*vtheta[j][2*i-1])+h*f3[2*i-2]*std::exp(-0.5*h*vtheta[j][2*i-1]);
  }
  XForceCache = func(t+h, x);
  for (j = 0; j < N2*P*d2; ++j) {
    int M2 = M / 2;
    for (i = 0; i < N; ++i)
      v[d*j+i] = v[d*j+i]*std::exp(-0.5*h*vtheta[j][0])+0.5*h*XForceCache[d*j+i]*std::exp(-0.25*h*vtheta[j][0]);
    for (i = 1; i <= M2; ++i)
      theta[j][2*i-2] = theta[j][2*i-2]+h*vtheta[j][2*i-2]/2;
    XVecD v2;
    for (i = 0; i < d; ++i)
      v2.push_back(v[d*j+i]);
    XVecD f3 = NHForce(beta, m, d, Q[j], v2, vtheta[j]);
    for (i = 1; i <= M2; ++i)
      vtheta[j][2*i-1] = vtheta[j][2*i-1]*std::exp(-0.5*h*((i==M2)?0:vtheta[j][2*i]))+0.5*h*f3[2*i-1]*std::exp(-0.25*h*((i==M2)?0:vtheta[j][2*i]));
  }
}

void pair_correlation3(const std::vector<XNum> &Rs, std::vector<XNum> &corrR, std::vector<XNum> &corrI, XNum phaR, XNum phaI) {
  int n = corrR.size();
  int i, j, k;
  for (k = 0; k < P; ++k)
    for (i = 0; i < N; ++i) {
      int in = i*P+k;
      XNum a = Rs[d*in]-L/2;
      XNum b = Rs[d*in+1]-L/2;
      XNum dis = std::sqrt(minimum_image(a)*minimum_image(a)+minimum_image(b)*minimum_image(b));
      //XNum dis = std::sqrt(minimum_image(b)*minimum_image(b));
      int index = int(dis / incre);
      if (index < 0)
        index = 0;
      if (index >= n)
        index = n - 1;
      corrR[index] += phaR;
      corrI[index] += phaI;
    }
  /*int n = corrR.size();
  int i, j, k;
  for (k = 0; k < P; ++k)
    for (i = 0; i < N; ++i)
      for (j = i + 1; j < N; ++j) {
        XNum dis = std::sqrt(distance(i*P+k, j*P+k, Rs));
        int index = int(dis / incre);
        if (index < 0)
          index = 0;
        if (index >= n)
          index = n - 1;
        corrR[index] += phaR;
        corrI[index] += phaI;
      }*/
}

std::ofstream out("data.txt");

void phase_dis(XNum p, std::vector<XNum> &corr) {
  int n = corr.size();
  /*if (p > 0)
    p = p-2*M_PI*int(p/(2*M_PI));
  else
    p = p+2*M_PI*(int(-p/(2*M_PI))+1);
  XNum incre = 2*M_PI/n;*/
  p += 4.5;
  XNum incre = 9.0/n;
  int index = int(p / incre);
  if (index < 0)
    index = 0;
  if (index >= n)
    index = n - 1;
  corr[index] += 1;
}

XNum simulation() {
  XForceCache.clear();
  ENkCache.clear();
  ENkCache2.clear();
  VBCache.clear();
  AjCache.clear();
  std::vector<XNum> Rs, Vs;
  std::vector<XNum> corrR;
  int i, j;
  initial(Rs, Vs);
  logR(Rs, std::cout);
  std::cout << std::endl;
  std::cout << ANk(N, N, Rs) << std::endl;
  velocity_rescale(Vs);
  XNum t = 0, h = 0.001;
  int skip = 100000;
  for (i = 0; i < cor_num; ++i)
    corrR.push_back(0);
  std::vector<XNum> pd = corrR;
  std::vector<XNum> corrI = corrR;
  std::vector<XVecD> theta, vtheta, Q;
  for (i = 0; i < d*N*P; ++i) {
    XVecD tmp;
    for (j = 0; j < 4; ++j) {
      tmp.push_back(0);
      tmp.push_back(0);
    }
    theta.push_back(tmp);
    tmp.clear();
    for (j = 0; j < 4; ++j) {
      tmp.push_back(1);
      tmp.push_back(1);
    }
    vtheta.push_back(tmp);
    tmp.clear();
    for (j = 0; j < 4; ++j) {
      tmp.push_back(1+0.5*j);
      tmp.push_back(1+0.5*j);
    }
    Q.push_back(tmp);
  }
  for (i = 0; i < skip; ++i) {
    XMMVerletNHChain(t, h, Rs, Vs, theta, vtheta, beta, m, Q, force);
    period_boundary(Rs);
    t += h;
  }
  logR(Rs, std::cout);
  logR(Rs, out);
  std::cout << std::endl;
  out << std::endl;
  int count = 0;
  int steps = 5000000;
  XNum e1 = 0, e2 = 0, s1 = 0, s2 = 0;
  XNum temp = 0;
  XVecD es1, es2, ss1, ss2;
  for (i = 0; i <= 10; ++i) {
    es1.push_back(0);
    es2.push_back(0);
    ss1.push_back(0);
    ss2.push_back(0);
  }
  for (i = 0; i < steps; ++i) {
    XMMVerletNHChain(t, h, Rs, Vs, theta, vtheta, beta, m, Q, force);
    period_boundary(Rs);
    if (i % 10 == 0) {
      temp += temperature(Vs);
      XNum te = P * d * N / (2 * beta);
      te += energy(Rs);
      te += trapEnergy(Rs);
      te += GauEnergy(Rs);
      std::pair<XNum, XNum> pha = APhase(Rs);
      e1 += pha.first*te;
      e2 += pha.second*te;
      s1 += pha.first;
      s2 += pha.second;
      XNum tvi = vi;
      for (j = 0; j <= 10; ++j) {
        vi = 0.0+0.1*j;
        std::pair<XNum, XNum> pha = APhase(Rs);
        es1[j] += pha.first*te;
        es2[j] += pha.second*te;
        ss1[j] += pha.first;
        ss2[j] += pha.second;
      }
      vi = tvi;
      //phase_dis(std::acos(pha.first), pd);
      phase_dis(ANk(N, N, Rs), pd);
      pair_correlation3(Rs, corrR, corrI, pha.first, pha.second);
      ++count;
    }
    if (i % 100000 == 0) {
      std::cout << "i=" << i << std::endl;
      out << "i=" << i << std::endl;
    }
    t += h;
  }
  temp /= count;
  std::cout << temp << std::endl;
  out << temp << std::endl;
  e1 /= count;
  e2 /= count;
  s1 /= count;
  s2 /= count;
  std::cout << e1 << "+" << e2 << "i" << std::endl;
  out << e1 << "+" << e2 << "i" << std::endl;
  std::cout << s1 << "+" << s2 << "i" << std::endl;
  out << s1 << "+" << s2 << "i" << std::endl;
  for (j = 0; j <= 10; ++j) {
    es1[j] /= count;
    es2[j] /= count;
    ss1[j] /= count;
    ss2[j] /= count;
    std::cout << "vi=" << 0.0+0.1*j << std::endl;
    std::cout << es1[j] << "+" << es2[j] << "i" << std::endl;
    std::cout << ss1[j] << "+" << ss2[j] << "i" << std::endl;
    out << "vi=" << 0.0+0.1*j << std::endl;
    out << es1[j] << "+" << es2[j] << "i" << std::endl;
    out << ss1[j] << "+" << ss2[j] << "i" << std::endl;
  }
  XNum norm = 0;
  for (i = 0; i < cor_num; ++i)
    pd[i] /= count;
  for (i = 0; i < cor_num; ++i)
    out << "{" << -4.5+i*9.0/cor_num << "," << pd[i] << "},";
    //out << "{" << i*2*M_PI/cor_num << "," << pd[i] << "},";
  out << std::endl;
  for (i = 0; i < cor_num; ++i)
    std::cout << "{" << -4.5+i*9.0/cor_num << "," << pd[i] << "},";
    //std::cout << "{" << i*2*M_PI/cor_num << "," << pd[i] << "},";
  std::cout << std::endl;
  XVecD denR;
  for (i = 0; i < cor_num; ++i) {
    XNum r = (corrR[i]*s1+corrI[i]*s2)/(s1*s1+s2*s2);
    denR.push_back(r);
  }
  norm = 0;
  for (i = 0; i < cor_num; ++i)
    norm += denR[i]*incre;
  for (i = 0; i < cor_num; ++i)
    denR[i] /= norm;
  for (i = 0; i < cor_num; ++i)
    denR[i] /= 2*M_PI*(i+0.5)*incre;
  for (i = 0; i < cor_num; ++i)
    std::cout << "{" << i*incre << "," << denR[i] << "},";
  std::cout << std::endl;
  for (i = 0; i < cor_num; ++i)
    out << "{" << i*incre << "," << denR[i] << "},";
  out << std::endl;
  logR(Rs, std::cout);
  logR(Rs, out);
  std::cout << std::endl;
  out << std::endl;
  return e1;
}

int main() {
  clock_t t;
  t = clock();
  N = 2;
  T = 1.0 / 1.0;
  beta = 1 / (kB * T);
  XSetRandSeed(7);
  for (P = 24; P <= 24; P += 3) {
    omgP = std::sqrt(P) / (beta * hBar);
    XNum e = simulation();
    std::cout << "{" << P << ", " << e << "}, ";
    out << "{" << P << ", " << e << "}, ";
  }
  std::cout << std::endl;
  out << std::endl;
  t = clock() - t;
  std::cout << (int)(((double)1000 * t) / CLOCKS_PER_SEC) << std::endl;
  out << (int)(((double)1000 * t) / CLOCKS_PER_SEC) << std::endl;
  return 0;
}
