// fft.h
// Fast Fourier transform functions
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

#ifndef _AURORA_FFT_
#define _AURORA_FFT_

#include <cmath>
#include <complex>

namespace Aurora {

static inline uint32_t np2(uint32_t n) {
  uint32_t v = 2;
  while (v < n)
    v <<= 1;
  return v;
}

/** FFT class  \n
    Radix-2 fast Fourier transform
*/
template <typename S> class FFT {

  static constexpr bool forward = true;
  static constexpr bool inverse = false;
  std::vector<std::complex<S>> c;
  bool pckd;

  void reorder(std::vector<std::complex<S>> &s) {
    uint32_t N = s.size();
    uint32_t j = 0, m;
    for (uint32_t i = 0; i < N; i++) {
      if (j > i) {
        std::swap(s[i], s[j]);
      }
      m = N / 2;
      while (m >= 2 && j >= m) {
        j -= m;
        m /= 2;
      }
      j += m;
    }
  }

public:
  /** Constructor \n
      N: transform size \n
      packed: complex data format in packed form (true) or not
  */
  FFT(std::size_t N, bool packed = true)
      : c(packed ? np2(N) / 2 : np2(N) / 2 + 1), pckd(packed){};

  std::size_t size() { return 2 * c.size(); }

  /** In-place complex-to-complex FFT \n
      data: data to be transformed  \n
      dir: true for forward operation, false for inverse \n
      pckd: true for packed data (first point is DC,Nyq)
      Size of data is expected to be a power-of-two.
  */
  void transform(std::vector<std::complex<S>> &s, bool dir) {
    uint32_t N = s.size();
    std::complex<S> wp, w, even, odd;
    S o;
    uint32_t i;
    reorder(s);
    for (uint32_t n = 1; n < N; n *= 2) {
      o = dir == forward ? -M_PI / n : M_PI / n;
      wp.real(std::cos(o)), wp.imag(std::sin(o));
      w = 1.;
      for (uint32_t m = 0; m < n; m++) {
        for (uint32_t k = m; k < N; k += n * 2) {
          i = k + n;
          even = s[k];
          odd = w * s[i];
          s[k] = even + odd;
          s[i] = even - odd;
        }
        w *= wp;
      }
    }
    if (dir == forward)
      for (uint32_t n = 0; n < N; n++)
        s[n] /= N;
  }

  /** Real-to-complex forward FFT \n
      r: input vector (real)
      returns pointer to complex data containing the transform
  */
  const std::complex<S> *transform(const std::vector<S> &r) {
    using namespace std::complex_literals;
    uint32_t N = c.size() - (pckd ? 0 : 1);
    std::complex<S> wp, w = 1., even, odd;
    S o, zro, nyq;
    S *s = reinterpret_cast<S *>(c.data());
    std::fill(c.begin(), c.end(), std::complex<S>(0, 0));
    std::copy(r.begin(), r.end(), s);
    if (!pckd)
      c.resize(N);
    transform(c, forward);
    zro = c[0].real() + c[0].imag();
    nyq = c[0].real() - c[0].imag();
    c[0].real(zro * .5), c[0].imag(nyq * .5);
    o = -M_PI / N;
    wp.real(std::cos(o)), wp.imag(std::sin(o));
    w *= wp;
    for (uint32_t i = 1, j = 0; i < N / 2; i++) {
      j = N - i;
      even = S(.5) * (c[i] + conj(c[j]));
      odd = std::complex<S>(.5i) * (conj(c[j]) - c[i]);
      c[i] = even + w * odd;
      c[j] = conj(even - w * odd);
      w *= wp;
    }
    if (!pckd) {
      c.resize(N + 1);
      c[N].real(c[0].imag());
      c[0].imag(0.);
      c[N].imag(0.);
    }
    return c.data();
  }

  /** Complex-to-real inverse FFT \n
      sp - input vector (complex) \n
      returns pointer to real data containing the transform
  */
  const S *transform(const std::vector<std::complex<S>> &sp, bool pckd = true) {
    using namespace std::complex_literals;
    uint32_t N = c.size() - (pckd ? 0 : 1);
    std::complex<S> wp, w = 1., even, odd;
    S o, zro, nyq;
    S *s = reinterpret_cast<S *>(c.data());
    std::copy(sp.begin(), sp.end(), c.begin());
    if (pckd)
      zro = c[0].real(), nyq = c[0].imag();
    else
      zro = c[0].real(), nyq = c[N].real();
    c[0].real(zro + nyq), c[0].imag(zro - nyq);
    o = M_PI / N;
    wp.real(std::cos(o)), wp.imag(std::sin(o));
    w *= wp;
    int j;
    for (uint32_t i = 1; i < N / 2; i++) {
      j = N - i;
      even = S(.5) * (c[i] + conj(c[j]));
      odd = std::complex<S>(.5i) * (c[i] - conj(c[j]));
      c[i] = even + w * odd;
      c[j] = conj(even - w * odd);
      w *= wp;
    }
    if (!pckd)
      c.resize(N);
    transform(c, inverse);
    return s;
  }
};
} // namespace Aurora

#endif // _AURORA_FFT_
