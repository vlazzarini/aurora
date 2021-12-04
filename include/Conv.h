// Conv.h
// Partitioned convolution
//
// (c) V Lazzarini, 2021
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
// 1. Redistributions of source code must retain the above copyright notice,
// this list of conditions and the following disclaimer.
// 2. Redistributions in binary form must reproduce the above copyright notice,
// this list of conditions and the following disclaimer in the documentation
// and/or other materials provided with the distribution.
// 3. Neither the name of the copyright holder nor the names of its contributors
// may be used to endorse or promote products derived from this software without
// specific prior written permission.
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
// LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
// CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
// SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
// INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
// CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
// ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
// POSSIBILITY OF SUCH DAMAGE

#include "FFT.h"
#include "SndBase.h"
#include <complex>
#include <vector>

namespace Aurora {
const std::size_t def_psize = 2048;

/** IR class \n
    Impulse response table
*/
template <typename S> class IR {
  std::vector<std::vector<std::complex<S>>> parts;

  void create(const std::vector<S> &s, std::size_t psize) {
    std::vector<S> buffer(psize);
    std::size_t n = 0, k = 0, cnt = s.size();
    FFT<S> fft(2 * psize, !packed, 1);

    while (cnt) {
      cnt = cnt > psize ? psize : cnt;
      std::copy(s.begin() + n, s.begin() + n + cnt, buffer.begin());
      fft.transform(buffer);
      auto v = fft.vector();
      std::copy(v.begin(), v.end(), parts[k++].begin());
      std::fill(buffer.begin(), buffer.end(), 0);
      n += cnt;
      cnt = s.size() - n;
    }
  }

public:
  /** Constructor \n
   s: impulse response \n
   psize: partition size
  */
  IR(const std::vector<S> &s, std::size_t psize = def_psize)
      : parts(ceil((float)s.size() / psize),
              std::vector<std::complex<S>>(psize + 1))

  {
    create(s, psize);
  }

  /** IR spectrum access \n
      returns a vector of impulse response partitions
  */
  const std::vector<std::vector<std::complex<S>>> &spectrum() const {
    return parts;
  }

  /** IR partition size \n
      returns the size of each partition
  */
  const std::size_t psize() const { return parts[0].size() - 1; }

  /** IR length \n
      returns the number of partitions
  */
  const std::size_t nparts() const { return parts.size(); }

  /** Reset the IR table \n
     s: impulse response \n
     psize: partition size
  */
  void reset(const std::vector<S> &s, std::size_t psize = def_psize) {
    parts.clear();
    parts.resize(ceil((float)s.size() / psize));
    for (auto &part : parts)
      part(std::vector<std::complex<S>>(psize + 1));
    create(s, psize);
  }
};

/** Overlap-save function for Conv \n
    in: input  \n
    bufin: convolution input buffer \n
    bufout: convolution output buffer \n
    unused: not used \n
    p: read/write pos \n
    psz: partition size \n
    returns convolution sample
*/
template <typename S>
S ols(S in, S *bufin, const S *bufout, S *unused, std::size_t p,
      std::size_t psz) {
  (void)unused;
  auto s = bufout[p + psz];
  bufin[p] = bufin[p + psz];
  bufin[p + psz] = in;
  return s;
}

/** Overlap-add function for Conv \n
    in: input \n
    bufin: convolution input buffer \n
    bufout: convolution output buffer \n
    olabuf: overlap-add buffer \n
    p: read/write pos \n
    psz: partition size \n
    returns convolution sample
*/
template <typename S>
S ola(S in, S *bufin, const S *bufout, S *olabuf, std::size_t cnt,
      std::size_t psz) {
  auto s = bufout[cnt] + olabuf[cnt];
  bufin[cnt] = in;
  olabuf[cnt] = bufout[cnt + psz];
  return s;
}

/** Conv class \n
    Partitioned convolution
*/
template <typename S> class Conv : public SndBase<S> {
  using SndBase<S>::process;
  const IR<S> *ir;
  std::vector<std::vector<std::complex<S>>> del;
  std::vector<std::complex<S>> mix;
  std::vector<S> inbuf;
  std::vector<S> olabuf;
  std::size_t p, sn;
  FFT<S> fft;
  std::function<S(S, S *, const S *, S *, std::size_t, std::size_t)> fun;

  void convol(const std::vector<S> &in) {
    fft.transform(in);
    auto &v = fft.vector();
    std::copy(v.begin(), v.end(), del[p].begin());
    std::fill(mix.begin(), mix.end(), 0);
    p = p == del.size() - 1 ? 0 : p + 1;
    auto dl = del.begin() + p;
    auto &part = ir->spectrum();
    for (auto prt = part.rbegin(); prt != part.rend(); prt++, dl++) {
      if (dl == del.end())
        dl = del.begin();
      auto dsamp = dl->begin();
      auto psamp = prt->begin();
      for (auto &mx : mix)
        mx += (*dsamp++ * *psamp++);
    }
    fft.transform(mix);
  }

public:
  /** Constructor \n
      IR: impulse response table\n
      f: convolution function \n
      vsize: vector size
  */
  Conv(const IR<S> *imp,
       const std::function<S(S, S *, const S *, S *, std::size_t, std::size_t)>
           f = ols<S>,
       std::size_t vsize = def_vsize)
      : SndBase<S>(vsize), ir(imp),
        del(imp->nparts(), std::vector<std::complex<S>>(imp->psize() + 1)),
        mix(imp->psize() + 1), inbuf(2 * imp->psize()), olabuf(imp->psize()),
        p(0), sn(0), fft(imp->psize() * 2, !packed, 1), fun(f){};

  /** Convolution \n
    in: input \n
    scal: output ampltude scaling
  */
  const std::vector<S> &operator()(const std::vector<S> &in, S scal) {
    std::size_t n = 0;
    S *bufin = inbuf.data();
    S *obuff = olabuf.data();
    std::size_t sz = ir->psize();
    return process(
        [&]() {
          auto s = fun(in[n++], bufin, fft.data(), obuff, sn, sz);
          if (++sn == sz) {
            convol(inbuf);
            sn = 0;
          }
          return s * scal;
        },
        in.size());
  }
};
} // namespace Aurora
