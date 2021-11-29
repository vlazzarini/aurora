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
#include <iostream>

namespace Aurora {
enum { SAW = 0, SQUARE, TRIANGLE, PULSE };
const int def_octs = 14;
const double def_base = 1.f;
const int def_ftlen = 16384;

/** TableSet class \n
    Creates a set of tables for BlOsc.
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
    uint32_t k = 1;
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
      if (fr > fs / 2)
        nh = 2;
      else
        nh = fs / (2 * fr) + 1;
      std::fill(blsp.begin() + nh, blsp.end(), std::complex<S>(0, 0));
      auto wv = fft.transform(blsp);
      std::copy(wv, wv + tlen, wave.begin());
      norm(wave);
    }
  }

  const std::vector<S> &select(S f) const {
    int32_t num = f > base ? (int32_t)std::log2(f / base) : 0;
    return num < waves.size() ? waves[num] : waves.back();
  }


public:
  /** Constructor \n
      type: wave type (SAW, SQUARE, TRIANGLE, PULSE) \n
      fs: sampling rate for which these will be built \n
      len: table length
  */
  TableSet(uint32_t type, S fs = def_sr, std::size_t len = def_ftlen)
      : tlen(len), waves(def_octs, std::vector<S>(len)), base((S)def_base) {

    std::vector<S> src(tlen / 2);
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
    fourier(src, fs, type);
  }

  /** Constructor \n
      src: source wave table \n
      b: base frequency for wavetables \n
      octs: number of octaves to generate \n
      fs: sampling rate for which these will be built
  */
  TableSet(const std::vector<S> &src, S b = def_base, S fs = def_sr)
      : tlen(src.size()), waves((int)std::log2(fs / b), std::vector<S>(tlen)),
        base(1 / (b * tlen / fs)) {
    fourier(src, fs);
  }

  /** Function selection \n
      f: fundamental frequency used for playback
  */
  std::function<S(S)> func(S f) const {
    return lookupi_gen(select(f));
  }

};

/** BlOsc class \n
    Bandlimited wavetable oscillator.
*/
template <typename S> class BlOsc : public Osc<S> {
  using Osc<S>::ph;
  using Osc<S>::ts;
  using Osc<S>::fun;
  const TableSet<S> *tset;
  S ff;

  virtual S synth(S a, S f, double &phs) {
    fun = ff != f ? tset->func(f) : fun;
    return Osc<S>::synth(a, f, phs);
  }

public:
  /** Constructor \n
      t: wavetable set \n
      fs: sampling rate \n
      vsize: vector size
  */
  BlOsc(const TableSet<S> *t, S fs = (S)def_sr, std::size_t vsize = def_vsize)
    : Osc<S>(nullptr, fs, vsize), tset(t), ff(0) {};

  /** Change the wavetable set
      t: wavetable set
   */
  void waveset(const TableSet<S> *t) { tset = t; }
};
} // namespace Aurora

#endif // _AURORA_BLOSC_
