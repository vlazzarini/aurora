// BlOsc.h
// Wavetable sets
// Bandlimited Wavetable Oscillator
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

#ifndef _AURORA_BLOSC_
#define _AURORA_BLOSC_

#include "FFT.h"
#include "Osc.h"

namespace Aurora {
enum { SAW = 0, SQUARE, TRIANGLE, PULSE };
const double def_base = 16.;

/** TableSet class \n
    Creates a set of tables for BlOsc. \n
    S: sample type
*/
template <typename S> class TableSet {
  std::size_t tlen;
  std::vector<std::vector<S>> waves;
  S base;

  void norm(std::vector<S> &wave) {
    S max = 0.;
    for (auto s : wave)
      if (s > max)
        max = s;
    for (auto &s : wave)
      s *= 1. / max;
  }

  void fourier(const std::vector<S> &src, S fs, int32_t type = -1) {
    std::size_t nh;
    uint32_t k = 0;
    FFT<S> fft(tlen);
    std::vector<std::complex<S>> blsp(fft.size() / 2);
    if (type < 0) {
      auto sp = fft.transform(src);
      std::copy(sp, sp + blsp.size(), blsp.begin());
    } else {
      std::size_t n = 0;
      for (auto &s : blsp) {
        if (type > 1)
          s.real(src[n++]);
        else
          s.imag(-src[n++]);
      }
    }
    for (auto &wave : waves) {
      double fr = base * std::pow(2, k++);
      if (fr > fs * .375)
        nh = 2;
      else
        nh = .375 * fs / fr + 1;
      std::fill(blsp.begin() + nh, blsp.end(), std::complex<S>(0, 0));
      auto wv = fft.transform(blsp);
      std::copy(wv, wv + tlen, wave.begin());
      norm(wave);
    }
  }

  void create(S fs, uint32_t type) {
    std::vector<S> src(tlen / 2);
    if (type >= 0) {
      std::size_t n = 0;
      for (auto &s : src) {
        switch (type) {
        case SAW:
          if (n)
            s = (S)1. / n;
          break;
        case SQUARE:
          if (n % 2)
            s = (S)1. / n;
          break;
        case TRIANGLE:
          if (n % 2)
            s = (S)1. / (n * n);
          break;
        default:
          s = (S)1.;
        }
        n++;
      }
    }
    fourier(src, fs, type);
  }

  void resize(std::size_t len, S fs) {
    tlen = len;
    waves.clear();
    waves.resize((int)std::log2(fs / base));
    for (auto &w : waves) {
      w = std::vector<S>(tlen);
    }
  }

public:
  /** Constructor \n
      type: wave type (SAW, SQUARE, TRIANGLE, PULSE) \n
      fs: sampling rate for which these will be built \n
      len: table length
  */
  TableSet(uint32_t type, S fs = def_sr, std::size_t len = def_ftlen)
      : tlen(len), waves((int)std::log2(fs / def_base), std::vector<S>(len)),
        base((S)def_base) {
    create(fs, type);
  }

  /** Constructor \n
      src: source wave table \n
      b:  source table base frequency \n
      fs: source table sampling rate
  */
  TableSet(const std::vector<S> &src, S b = def_base, S fs = def_sr)
      : tlen(src.size()), waves((int)std::log2(fs / b), std::vector<S>(tlen)),
        base(1 / (b * tlen / fs)) {
    fourier(src, fs);
  }

  /** table selection \n
      f: fundamental frequency used for playback
  */
  const std::vector<S> &func(S f) const {
    std::size_t num = f > base ? (int32_t)round(std::log2(f / base)) : 0;
    return num < waves.size() ? waves[num] : waves.back();
  }

  /** reset the table set \n
      type: wave type (SAW, SQUARE, TRIANGLE, PULSE) \n
      fs: sampling rate for which these will be built
      tlen: table size
   */
  void reset(uint32_t type, S fs, std::size_t len = def_ftlen) {
    resize(len, fs);
    create(fs, type);
  }
  /** reset the table set \n
     src: source wave table \n
     b:  source table base frequency \n
     fs: source table sampling rate \n
     tlen: table size
  */
  void reset(const std::vector<S> &src, S b, S fs) {
    base(1 / (b * src.size() / fs));
    resize(src.size(), fs);
    fourier(src, fs);
  }

  void guardpoint() {
    for (auto &w : waves)
      w.push_back(w[0]);
  }   

};

/** BlOsc class \n
    Bandlimited wavetable oscillator. \n
    S: sample type \n
    FN: oscillator function
*/
template <typename S, S (*FN)(double, const std::vector<S> *) = lookupi<S>>
class BlOsc : public Osc<S, FN> {
  using Osc<S, FN>::ph;
  using Osc<S, FN>::ts;
  using Osc<S, FN>::tab;
  const TableSet<S> *tset;
  S ff;

  S synth(S a, S f, double &phs, const std::vector<S> *t, S pm = 0) override {
    t = ff != f || !t ? &tset->func(f) : t;
    ff = f;
    return Osc<S, FN>::synth(a, f, phs, t, pm);
  }

  void reset_obj(S fs) override {
    ff = 0;
    Osc<S, FN>::reset_obj(fs);
  }

public:
  /** Constructor \n
      t: wavetable set \n
      fs: sampling rate \n
      vsize: vector size
  */
  BlOsc(const TableSet<S> *t, S fs = (S)def_sr, std::size_t vsize = def_vsize)
      : Osc<S, FN>(nullptr, fs, vsize), tset(t), ff(-fs){};

  /** Change the wavetable set
      t: wavetable set
  */
  void waveset(const TableSet<S> *t) { tset = t; }
};
} // namespace Aurora

#endif // _AURORA_BLOSC_
