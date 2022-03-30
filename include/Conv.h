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
const bool ola = 0;
const bool ols = 1;
/** IR class \n
    Impulse response table
*/
template <typename S> class IR {
  std::vector<std::vector<std::complex<S>>> parts;

  void create(const std::vector<S> &s, std::size_t psize) {
    std::vector<S> buffer(psize);
    std::size_t n = 0, k = 0, cnt = s.size();
    FFT<S> fft(2 * psize, !packed, inverse);

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
              std::move(std::vector<std::complex<S>>(psize + 1)))

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
      part = std::move(std::vector<std::complex<S>>(psize + 1));
    create(s, psize);
  }
};

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
  bool meth;

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

  S oladd(S in, S *bufin, const S *bufout, S *olabuf, std::size_t cnt,
          std::size_t psz) {
    auto s = bufout[cnt] + olabuf[cnt];
    bufin[cnt] = in;
    olabuf[cnt] = bufout[cnt + psz];
    return s;
  }

  S olsave(S in, S *bufin, const S *bufout, S *unused, std::size_t p,
           std::size_t psz) {
    (void)unused;
    auto s = bufout[p + psz];
    bufin[p] = bufin[p + psz];
    bufin[p + psz] = in;
    return s;
  }

public:
  /** Constructor \n
      IR: impulse response table\n
      algo: convolution algorithm (ols = 0, ola = 1) \n
      vsize: vector size
  */
  Conv(const IR<S> *imp, bool algo = ols, std::size_t vsize = def_vsize)
      : SndBase<S>(vsize), ir(imp),
        del(imp->nparts(), std::vector<std::complex<S>>(imp->psize() + 1)),
        mix(imp->psize() + 1), inbuf(2 * imp->psize()), olabuf(imp->psize()),
        p(0), sn(0), fft(imp->psize() * 2, !packed, inverse), meth(algo){};

  /** Convolution \n
    in: input \n
    scal: output ampltude scaling
  */
  const std::vector<S> &operator()(const std::vector<S> &in, S scal) {
    std::size_t n = 0;
    S *bufin = inbuf.data();
    S *obuff = olabuf.data();
    std::size_t sz = ir->psize();
    return meth ? process(
                      [&]() {
                        auto s =
                            oladd(in[n++], bufin, fft.data(), obuff, sn, sz);
                        if (++sn == sz) {
                          convol(inbuf);
                          sn = 0;
                        }
                        return s * scal;
                      },
                      in.size())
                : process(
                      [&]() {
                        auto s =
                            olsave(in[n++], bufin, fft.data(), obuff, sn, sz);
                        if (++sn == sz) {
                          convol(inbuf);
                          sn = 0;
                        }
                        return s * scal;
                      },
                      in.size());
  }
};



/** TVConv class \n
    Time-Varying Partitioned convolution
*/
template <typename S> class TVConv : public SndBase<S> {
  using SndBase<S>::process;
  std::vector<std::vector<std::complex<S>>> del1, del2;
  std::vector<std::complex<S>> mix;
  std::vector<S> inbuf1, inbuf2;
  std::vector<S> olabuf;
  std::size_t p, sn;
  FFT<S> fft;
  std::size_t sz;

  void convol(const std::vector<S> &in1, const std::vector<S> &in2) {
    fft.transform(in1);
    auto &v1 = fft.vector();
    std::copy(v1.begin(), v1.end(), del1[p].begin());
    fft.transform(in2);
    auto &v2 = fft.vector();
    std::copy(v2.begin(), v2.end(), del2[p].begin());
    std::fill(mix.begin(), mix.end(), 0);
    p = p == del1.size() - 1 ? 0 : p + 1;
    auto dl = del1.begin() + p;
    auto &part = del2;
    for (auto prt = part.rbegin(); prt != part.rend(); prt++, dl++) {
      if (dl == del1.end())
        dl = del1.begin();
      auto dsamp = dl->begin();
      auto psamp = prt->begin();
      for (auto &mx : mix)
        mx += (*dsamp++ * *psamp++);
    }
    fft.transform(mix);
  }

  S oladd(S in1, S in2, S *bufin1,  S *bufin2,
	  const S *bufout, S *olabuf, std::size_t cnt,
          std::size_t psz) {
    auto s = bufout[cnt] + olabuf[cnt];
    bufin1[cnt] = in1;
    bufin2[cnt] = in2;
    olabuf[cnt] = bufout[cnt + psz];
    return s;
  }

public:
  /** Constructor \n
      len: filter len \n
      psize: partition size \n
      vsize: vector size
  */
 TVConv(std::size_t len, std::size_t psize = def_psize, bool algo = ols, std::size_t vsize = def_vsize)
      : SndBase<S>(vsize),
        del1(len/psize, std::vector<std::complex<S>>(psize + 1)),
        del2(len/psize, std::vector<std::complex<S>>(psize + 1)),
        mix(psize + 1), inbuf1(2 * psize), inbuf2(psize),
        olabuf(psize),
        p(0), sn(0), fft(psize * 2, !packed, inverse), sz(psize){};

  /** Convolution \n
    in: input \n
    scal: output ampltude scaling
  */
  const std::vector<S> &operator()(const std::vector<S> &in1, const std::vector<S> &in2, S scal) {
    std::size_t n = 0;
    S *bufin1 = inbuf1.data();
    S *bufin2 = inbuf2.data();
    S *obuff = olabuf.data();
    return process( [&]() {
                        auto s =
			  oladd(in1[n], in2[n], bufin1, bufin2, fft.data(), obuff, sn, sz);
                        if (++sn == sz) {
                          convol(inbuf1,inbuf2);
                          sn = 0;
                        }
			n++;
                        return s *scal;
                      },
                      in1.size());
  }

  void reset(std::size_t len, std::size_t psize = def_psize) {
    del1 = std::vector<std::vector<std::complex<S>>>(len/psize, std::vector<std::complex<S>>(psize + 1));
    del2 = std::vector<std::vector<std::complex<S>>>(len/psize, std::vector<std::complex<S>>(psize + 1));
    mix = std::vector<std::complex<S>>(psize + 1);
    inbuf1 = std::vector<S>(2 * psize);
    inbuf2 = std::vector<S>(psize);
    olabuf = std::vector<S>(psize);
    p = 0;
    sn = 0;
    fft = FFT<S>(psize * 2, !packed, inverse);
    sz = psize;     
  }
           
};

 

 
} // namespace Aurora
